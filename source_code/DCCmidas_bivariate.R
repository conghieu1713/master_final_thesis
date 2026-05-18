# ============================================================================== #
# ------------                                                       ----------- #
# ------------    ƯỚC LƯỢNG MÔ HÌNH DCC-MIDAS VÀ DCC-MIDAS-X         ----------- #
# ------------            (CHO CẶP TÀI SẢN TÙY CHỌN)                 ----------- #
# ------------             Sử dụng thư viện dccmidas                 ----------- #
# ============================================================================== #

# Tác giả: Công Hiếu, NGUYỄN (524102110660)
# Học viên Cao học, Đại học Kinh tế TP. Hồ Chí Minh


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 1. LÀM SẠCH MÔI TRƯỜNG ---------------------------------------------------- #
# ------------------------------------------------------------------------------ #
  rm(list = ls(all.names = TRUE))
  graphics.off()
  library(xts)     # Khai báo một thư viện để lệnh gỡ bỏ thư viện không gặp lỗi
  invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), 
                 detach, character.only = TRUE, unload = TRUE, force = TRUE))


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 2. THIẾT LẬP CÁC THAM SỐ CHO MÔ HÌNH -------------------------------------- #
# ------------------------------------------------------------------------------ #
  # Đường dẫn thư mục làm việc và tệp dữ liệu
    dir <- "D:/Post-graduate programme/Final Thesis/Model"   
    setwd(dir)
    file_path <- "official.xlsx"

  # CHỈ ĐỊNH TÊN SHEET DỮ LIỆU
    rate_sheet_spec = "rate"
    shocks_sheet_spec = "shocks_lg"
    midas_sheet_spec = "rv"

  # CHỈ ĐỊNH 2 TÀI SẢN BẠN MUỐN PHÂN TÍCH
    asset1 = "VNIndex"
    asset2 = "GC_F"
    
  # CHỈ ĐỊNH YẾU TỐ MIDAS (tần suất thấp)
  # Chỉ định tên cột tương ứng với asset1 và asset2
    midas1_col_spec = "RV_VNIndex"
    midas2_col_spec = "RV_GC_F"

  # CHỈ ĐỊNH BIẾN VĨ MÔ (X) CHO MÔ HÌNH DCC-MIDAS-X
    shock_name = "GPRT"

  # CÁC THAM SỐ KHÁC
    K <- 36; N_c <- 21; K_c <- 756; K_x <- 36
    do_arma = FALSE
    univ_model = "GM_noskew"
    distribution = "std"
    lag_fun = "Beta"


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 3. KHAI BÁO THƯ VIỆN ------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(dccmidas)
  library(rumidas)
  library(ggplot2)
  library(patchwork)
  library(zoo)
  library(xts)
  library(FinTS)
  library(forecast)
  library(Rsolnp)

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 4. NHẬP DỮ LIỆU ----------------------------------------------------------- #
# ------------------------------------------------------------------------------ #
  rate_data <- read_excel(file_path, sheet = rate_sheet_spec, guess_max = 100000) %>% mutate(Date = as.Date(Date))
  shocks_data <- read_excel(file_path, sheet = shocks_sheet_spec, guess_max = 100000) %>% mutate(Date = as.Date(Date))
  midas_data <- read_excel(file_path, sheet = midas_sheet_spec, guess_max = 100000) %>% mutate(Date = as.Date(Date))


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 5. HÀM LẤY DỮ LIỆU SONG BIẾN ---------------------------------------------- #
# ------------------------------------------------------------------------------ #
get_bivariate_data <- function(asset1, asset2, df_rate, df_midas, midas1_col, midas2_col, K, do_arma = FALSE) {
  # ============================================
  # BỘ ĐỆM ĐỘNG (DYNAMIC BUFFER)
  # Tránh lỗi Hessian (At least one parameter must not be fixed)
  # Loại bỏ các giá trị NA
  # ============================================
    valid_midas_data <- df_midas %>% 
      select(Date, all_of(midas1_col), all_of(midas2_col)) %>% 
      drop_na()
  # Ngày bắt đầu thực tế
    start_mv <- min(valid_midas_data$Date)
  # Bộ đệm cho K
    start_rate_limit <- start_mv %m+% months(K+1)
  # =============================================
  # DỮ LIỆU TẦN SỐ CAO (HIGH-FREQ)
  # =============================================
    rate_pair <- df_rate %>% 
      select(Date, all_of(asset1), all_of(asset2)) %>% 
      filter(Date >= start_rate_limit) %>% 
      drop_na()
  
  # Giữ lại tên cột: dùng drop = FALSE để tránh lỗi "replacement has length zero"
    xts_asset1 <- xts(rate_pair[, asset1, drop = FALSE], order.by = rate_pair$Date)
    xts_asset2 <- xts(rate_pair[, asset2, drop = FALSE], order.by = rate_pair$Date)
  
  # Danh sách r_t: Tránh lỗi "the condition has length > 1" của dcc_fit
  # Thêm nhiễu trắng (sd = 1e-6) cho SJC nếu có
    if (asset1 == "SJC") {
      set.seed(123) # Cố định seed để kết quả không đổi
      noise <- rnorm(length(xts_asset1), mean = 0, sd = 1e-6)
      xts_asset1 <- xts_asset1 + noise
    }
    if (asset2 == "SJC") {
      set.seed(456) # Dùng seed khác để tránh trùng lặp nếu cả 2 là SJC
      noise <- rnorm(length(xts_asset2), mean = 0, sd = 1e-6)
      xts_asset2 <- xts_asset2 + noise
    }
  
  # Loại bỏ trung bình (demean) trước khi lọc ARMA hoặc chạy GARCH
    mean_asset1 <- mean(xts_asset1, na.rm = TRUE)
    mean_asset2 <- mean(xts_asset2, na.rm = TRUE)
    xts_asset1 <- xts_asset1 - mean_asset1
    xts_asset2 <- xts_asset2 - mean_asset2
    
  # LỌC TRƯỚC (PRE-FILTERING) bằng ARMA(1,1)
    if (do_arma) {
      # Áp dụng bộ lọc ARMA(1,1) cho chuỗi lợi nhuận
      arma_asset1 <- Arima(xts_asset1, order = c(1,0,1), include.mean = TRUE)
      arma_asset2 <- Arima(xts_asset2, order = c(1,0,1), include.mean = TRUE)
      xts_asset1 <- xts(as.numeric(residuals(arma_asset1)), order.by = index(xts_asset1))
      colnames(xts_asset1) <- asset1
      xts_asset2 <- xts(as.numeric(residuals(arma_asset2)), order.by = index(xts_asset2))
      colnames(xts_asset2) <- asset2
    }
    
    r_t_list <- list(xts_asset1, xts_asset2)  
    
  # ===========================================
  # DỮ LIỆU TẦN SỐ THẤP (LOW-FREQ)
  # ===========================================
    mv_xts <- xts(valid_midas_data[, -1], order.by = valid_midas_data$Date)
  # Chuyển đổi MV thành ma trận
    mv_asset1 <- mv_into_mat(r_t_list[[1]], mv_xts[, midas1_col], K = K, type = "monthly")
    mv_asset2 <- mv_into_mat(r_t_list[[2]], mv_xts[, midas2_col], K = K, type = "monthly")
  
    Z_list <- list(mv_asset1, mv_asset2)
    
  # ===========================================
  # TRẢ VỀ KẾT QUẢ
  # ===========================================
    return(list(
      r_t = r_t_list,
      MV = Z_list,
      dates = index(xts_asset1),
      asset1 = asset1,
      asset2 = asset2,
      K_used = K
    ))
}

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 6. HÀM ƯỚC LƯỢNG DCC-MIDAS ------------------------------------------------ #
# ------------------------------------------------------------------------------ #
run_dcc_estimation <- function(asset1, asset2, rate_data, midas_data, midas1_col, midas2_col,
                               K = 36, N_c = 21, K_c = 504,
                               do_arma = FALSE,
                               univ_model = "GM_noskew",
                               distribution = "norm",
                               lag_fun = "Beta") {
  # ====================================================
  # CHUẨN BỊ DỮ LIỆU, sử dụng hàm get_bivariate_data
  # ====================================================
    prepared_data <- get_bivariate_data(asset1, asset2, rate_data, midas_data, midas1_col, midas2_col, K = K, do_arma = do_arma)
  
  if (is.null(prepared_data)) {
    cat(">> [LỖI] Không thể xử lý dữ liệu cho cặp", asset1, "-", asset2, "\n")
    return(NULL)
  }
  
  # ====================================================
  # ƯỚC LƯỢNG MÔ HÌNH
  # ====================================================
  cat(sprintf("\n[ƯỚC LƯỢNG] DCC-MIDAS cho %s và %s\n", asset1, asset2))
  cat(sprintf("K=%d, N_c=%d, Kc=%d, Model=%s, distribution=%s, lag_fun=%s\n",
              K, N_c, K_c, univ_model, distribution, lag_fun))
  
  fit <- tryCatch({
    dcc_fit(
      r_t = prepared_data$r_t, 
      MV = prepared_data$MV, 
      univ_model = univ_model, 
      corr_model = "DCCMIDAS", 
      K = prepared_data$K_used,
      N_c = N_c, 
      K_c = K_c, 
      lag_fun = "Beta",
      distribution = "norm"
    )
  }, error = function(e) {
    cat(">> [LỖI ƯỚC LƯỢNG]:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(fit)) {
    cat(">> [ƯỚC LƯỢNG] Thành công!\n")
  }
  
  # ====================================================
  # TRẢ VỀ KẾT QUẢ
  # ====================================================
  return(list(
    fit = fit,               
    data = prepared_data
  ))
}

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 7. CÁC HÀM TIỆN ÍCH (VẼ BIỂU ĐỒ, TRÍCH XUẤT) ------------------------------ #
# ------------------------------------------------------------------------------ #

# -- 7.1. HÀM VẼ BIỂU ĐỒ KẾT QUẢ DCC ------------------------------------------- #
plot_dcc_results <- function(estimation_result, filename = NULL) {
  # ====================================================
  # TRÍCH XUẤT DỮ LIỆU TỪ KẾT QUẢ ƯỚC LƯỢNG
  # ====================================================
    fit <- estimation_result$fit
    data_info <- estimation_result$data
    asset1_name <- data_info$asset1
    asset2_name <- data_info$asset2
    plot_dates <- data_info$dates
  
  # ====================================================
  # TÍNH TOÁN BIẾN ĐỘNG VÀ TƯƠNG QUAN
  # ====================================================
    vol_asset1 <- sqrt(fit$H_t[1, 1, ])
    vol_asset2 <- sqrt(fit$H_t[2, 2, ])
    rho_sr <- fit$R_t[1, 2, ]
    rho_lr <- fit$R_t_bar[1, 2, ]
  # ====================================================
  # TẠO DATA FRAME ĐỂ VẼ BIỂU ĐỒ
  # ====================================================
    plot_df <- data.frame(
      Date = plot_dates,
      Vol_Asset1 = as.numeric(vol_asset1), # Biến động có điều kiện của tài sản 1
      Vol_Asset2 = as.numeric(vol_asset2), # Biến động có điều kiện của tài sản 2
      SR_Corr = as.numeric(rho_sr),      # Tương quan ngắn hạn
      LR_Corr = as.numeric(rho_lr)       # Tương quan dài hạn
    )
  
  # ====================================================
  # LÀM SẠCH DỮ LIỆU - LOẠI BỎ GIAI ĐOẠN "BURN-IN"
  # ====================================================
  # Tự động xác định điểm bắt đầu (loại bỏ các dòng burn-in)
    burn_in_index <- which(plot_df$SR_Corr != 0 & !is.na(plot_df$SR_Corr))[1]
    if(is.na(burn_in_index)) burn_in_index <- 1
    plot_df_clean <- plot_df %>% slice(burn_in_index:nrow(plot_df))
  
  # ====================================================
  # THEME SETTING UP
  # ====================================================
    theme1 <- theme_bw() + 
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        # Cài đặt lưới
          panel.grid.major = element_line(linetype = "dotted", color = "darkgray"),
          panel.grid.minor.x = element_line(linetype = "dotted", color = "darkgray"),
          panel.grid.minor.y = element_blank(),
        # Ẩn tên trục
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
        # Xoay và căn chỉnh văn bản trên trục
          axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(r = 5)),
          axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1, margin = margin(t = 5)),
        # Đưa dấu tick vào trong
          axis.ticks.length = unit(-0.15, "cm"),
        # Viền đen sắc nét
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
      )
    
  # ====================================================
  # VẼ TỪNG BIỂU ĐỒ CON (sử dụng ggplot2)
  # ====================================================
  # Biểu đồ A: Biến động tài sản 1
    p1 <- ggplot(plot_df_clean, aes(x = Date, y = Vol_Asset1)) +
      geom_line(color = "black", linewidth = 0.75) +
      scale_x_date(date_labels = "%m/%Y", date_breaks = "2 years", date_minor_breaks = "1 year") +
      labs(title = asset1_name, y = "Conditional Volatility", x = "") +
      theme1 + ylim(0,8)
  # Panel B: Tương quan động
    p2 <- ggplot(plot_df_clean, aes(x = Date)) +
      geom_line(aes(y = SR_Corr), color = "black", alpha = 0.75, linewidth = 0.6) +
      geom_line(aes(y = LR_Corr), color = "red", linewidth = 1) +
      scale_x_date(date_labels = "%m/%Y", date_breaks = "2 years", date_minor_breaks = "1 year") +
      labs(title = paste(asset1_name, "-", asset2_name), y = "Dynamic Correlation", x = "") +
      theme1 + ylim(-0.3,0.3)
  # Panel C: Biến động tài sản 2
    p3 <- ggplot(plot_df_clean, aes(x = Date, y = Vol_Asset2)) +
      geom_line(color = "black", linewidth = 0.75) +
      scale_x_date(date_labels = "%m/%Y", date_breaks = "2 years", date_minor_breaks = "1 year") +
      labs(title = asset2_name, y = "Conditional Volatility", x = "") +
      theme1 + ylim(0,8)
  
  # ====================================================
  # KẾT HỢP VÀ XUẤT BIỂU ĐỒ
  # ====================================================
    des <- c(
      area(1,1),
      area(1,2),
      area(2,2)
    )
    final_plot <- p1 + p2 + p3 + 
      plot_layout(design = des) +
      plot_annotation(
        theme = theme(plot.caption = element_text(hjust = 0.5))
      )
    
  # In ra màn hình RStudio
    print(final_plot)
    
  # Nếu có tên tệp, lưu lại dưới dạng PNG chất lượng cao
    if (!is.null(filename)) {
      dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
      ggsave(filename, plot = final_plot, width = 12, height = 8, dpi = 300)
      cat(sprintf(">> Đã lưu biểu đồ tại: %s\n", filename))
    }
  
  return(final_plot)
}

# -- 7.2. HÀM TRÍCH XUẤT PHẦN DƯ CHUẨN HÓA -------------------------------------- #
  extract_std_residuals <- function(estimation_result) {
    # Trích xuất các thành phần từ kết quả ước lượng
      fit <- estimation_result$fit
      data_info <- estimation_result$data
    # Tên tài sản để đặt tên cho cột kết quả
      asset1_name <- data_info$asset1
      asset2_name <- data_info$asset2
    # LẤY LỢI NHUẬN - r_t
      r_asset1 <- as.numeric(data_info$r_t[[1]])
      r_asset2 <- as.numeric(data_info$r_t[[2]])
    # BIẾN ĐỘNG CÓ ĐIỀU KIỆN - sqrt(h_t)
      vol_asset1 <- as.numeric(sqrt(fit$H_t[1, 1, ]))
      vol_asset2 <- as.numeric(sqrt(fit$H_t[2, 2, ]))
    # TÍNH TOÁN PHẦN DƯ CHUẨN HÓA
      std_resid_asset1 <- r_asset1 / vol_asset1
      std_resid_asset2 <- r_asset2 / vol_asset2
    # TRẢ VỀ KẾT QUẢ
      std_resid_df <- data.frame(
        Date = data_info$dates,
        res1 = std_resid_asset1,
        res2 = std_resid_asset2
      )
      colnames(std_resid_df)[2] <- asset1_name
      colnames(std_resid_df)[3] <- asset2_name
      
      cat(sprintf("\n>> Đã trích xuất %d phần dư chuẩn hóa cho %s và %s.", 
                  nrow(std_resid_df), asset1_name, asset2_name), "\n")
    
    return(std_resid_df)
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 8. CÁC HÀM CHO MÔ HÌNH DCC-MIDAS-X ---------------------------------------- #
# ------------------------------------------------------------------------------ #

# -- 8.1. HÀM CHUẨN BỊ DỮ LIỆU CHO DCC-MIDAS-X --------------------------------- #
prepare_dcc_x_data <- function(estimation_result, macro_data, shock_name, K_x) {
  
  # Kiểm tra các cột cần thiết
    if (!"Date" %in% colnames(macro_data)) {
      stop("Lỗi: Dữ liệu vĩ mô phải có cột 'Date'.")
    }
    if (!shock_name %in% colnames(macro_data)) {
      stop("Lỗi: Không tìm thấy biến '", shock_name, "' trong dữ liệu vĩ mô.")
    }
  
  # TRÍCH XUẤT PHẦN DƯ CHUẨN HÓA (gọi lại hàm extract_std_residuals)
    std_resid_df <- extract_std_residuals(estimation_result)
    res_matrix <- as.matrix(std_resid_df[, 2:3]) # Lấy cột 2 tài sản
    rownames(res_matrix) <- as.character(std_resid_df$Date)
  
  # TẠO MA TRẬN MIDAS CHO BIẾN VĨ MÔ
    shock_xts <- xts(as.numeric(macro_data[[shock_name]]), order.by = macro_data$Date)
    X_mat <- mv_into_mat(
      x = estimation_result$data$r_t[[1]], 
      mv = shock_xts,  
      K = K_x, 
      type = "monthly"
    )
    
  # ĐỒNG BỘ DỮ LIỆU
  # Đảm bảo trục thời gian của phần dư hàng ngày và biến vĩ mô hàng tháng khớp nhau
    X_mat_temp <- t(X_mat)
    daily_dates <- as.Date(index(estimation_result$data$r_t[[1]]))
    row.names(X_mat_temp) <- as.character(daily_dates)
    common_dates <- intersect(row.names(res_matrix), row.names(X_mat_temp))
    res_matrix_final <- res_matrix[common_dates, ]
    X_mat_final <- X_mat_temp[common_dates, ]
    cat(sprintf("\n>> [HOÀN TẤT] Đã chuẩn bị %d quan sát cho biến %s.\n", 
                nrow(res_matrix_final), shock_name), "\n")
  
  # Trả về một danh sách chứa hai ma trận đã được làm sạch để đưa vào solnp
  return(list(
    res_matrix = res_matrix_final,
    X_mat = X_mat_final
  ))
}

# -- 8.2. HÀM TÍNH TRỌNG SỐ BETA MIDAS ----------------------------------------- #
# Hàm trọng số Beta (với tham số w1 = 1 để đảm bảo trọng số giảm dần)
beta_weights <- function(K, w) {
  k <- 1:K
  weights <- (1 - k/(K+1))^(w - 1)
  weights <- weights / sum(weights)
  return(weights)
}

# -- 8.3. HÀM LOG-LIKELIHOOD CHO DCC-MIDAS-X ------------------------------------ #
dcc_midas_x_loglik <- function(params, res_data, X_data) {
  a <- params[1]
  b <- params[2]
  m <- params[3]
  theta <- params[4]
  w <- params[5]
  
  # --- BẢO VỆ 1: Ràng buộc tham số ---
  # Nếu numDeriv cố gắng di chuyển tham số ra ngoài vùng an toàn, trả về một hình phạt lớn
    if (a <= 0 || b <= 0 || (a + b) >= 0.999 || w <= 1.001) {
      return(1e6) 
    }
    T_obs <- nrow(res_data)
    K_lags <- ncol(X_data)
    weights <- beta_weights(K_lags, w)
    Z_t <- m + theta * (X_data %*% weights)
    rho_bar_t <- (exp(2 * Z_t) - 1) / (exp(2 * Z_t) + 1)
    Q_prev <- cov(res_data)
    LL <- 0
  
    for (t in 2:T_obs) {
      R_bar_t <- matrix(c(1, rho_bar_t[t], 
                          rho_bar_t[t], 1), nrow = 2, ncol = 2)
      
      eps_prev <- matrix(res_data[t-1, ], ncol = 1)
      Q_t <- (1 - a - b) * R_bar_t + a * (eps_prev %*% t(eps_prev)) + b * Q_prev
    
    # --- BẢO VỆ 2: Ngăn chặn ma trận Q_t suy biến ---
      if (any(is.na(Q_t)) || Q_t[1,1] <= 0 || Q_t[2,2] <= 0) return(1e6)
      D_t <- diag(1 / sqrt(diag(Q_t)))
      R_t <- D_t %*% Q_t %*% D_t
      eps_t <- matrix(res_data[t, ], ncol = 1)
      det_R <- det(R_t)
    
    # --- BẢO VỆ 3: Định thức ma trận R_t ---
    # Thay vì gán 1e-8, chúng ta phạt trực tiếp để thay đổi hướng của Hessian
      if (is.na(det_R) || det_R <= 0) return(1e6 + sum(params^2)*1000) 
      inv_R <- solve(R_t)
      LL_t <- -0.5 * (log(det_R) + t(eps_t) %*% inv_R %*% eps_t - t(eps_t) %*% eps_t)
      LL <- LL + LL_t
      Q_prev <- Q_t
  }
  return(-as.numeric(LL)) 
}

# -- 8.4. HÀM TỔNG HỢP KẾT QUẢ DCC-MIDAS-X ------------------------------------- #
  report_dcc_midas_results <- function(opt_res, res_matrix, x_matrix, model_params) {
    # 1. Trích xuất các tham số cơ bản
      params <- opt_res$pars
      n_param <- length(params)
      n_obs <- nrow(res_matrix)
    
    # --- Tên tham số tương ứng
      param_names <- c("a (alpha)", "b (beta)", "m (mu)", "theta", "w (omega)")
    
    # 2. Tính ma trận Hessian và Sai số chuẩn (Std. Error)
    # --- Tính ma trận Hessian
      h <- numDeriv::hessian(func = dcc_midas_x_loglik, x = params, 
                             res_data = res_matrix, X_data = x_matrix)
    
    # --- Tính ma trận hiệp phương sai bằng cách nghịch đảo Hessian
      cov_mat <- tryCatch({
        solve(h)
      }, error = function(e) {
        cat("\n[CẢNH BÁO] Ma trận Hessian suy biến, hàm solve() thất bại!\n")
        cat("Chi tiết:", conditionMessage(e), "\n")
        # Trả về ma trận NA nếu thất bại để không làm dừng mã
        return(matrix(NA, nrow = length(params), ncol = length(params)))
      })
      
      std_err <- sqrt(abs(diag(cov_mat)))
    
    # 3. Tính toán t-stat và P-value
      t_stat <- params / std_err
      p_val <- 2 * (1 - pnorm(abs(t_stat)))
    
    # 4. Tính toán các chỉ số lựa chọn mô hình
      ll <- -tail(opt_res$values, 1)
      aic <- 2 * n_param - 2 * ll
      bic <- n_param * log(n_obs) - 2 * ll
    
    # 5. Tạo bảng kết quả
      results_df <- data.frame(
        Param = param_names,
        Estimate = round(params, 6),
        `Std. Error` = round(std_err, 6),
        `t value` = round(t_stat, 4),
        `Pr(>|t|)` = format.pval(p_val, digits = 4, eps = 0.0001),
        `Sig.` = ifelse(p_val < 0.01, "***", ifelse(p_val < 0.05, "**", ifelse(p_val < 0.1, "*", "")))
      )
    
    # 6. In bảng kết quả ra Console
      cat("\n======================================================================\n")
      cat("      KẾT QUẢ ƯỚC LƯỢNG: MÔ HÌNH DCC-MIDAS-X\n")
      cat("----------------------------------------------------------------------\n")
      cat(sprintf(" Cặp tài sản: %s - %s\n", model_params$asset1, model_params$asset2))
      cat(sprintf(" Biến vĩ mô (X): %s\n", model_params$shock_name))
      cat("----------------------------------------------------------------------\n")
      cat(sprintf(" Tham số: K=%d, N_c=%d, K_c=%d, K_x=%d\n", 
                  model_params$K, model_params$N_c, model_params$K_c, model_params$K_x))
      cat(sprintf(" Cấu hình: univ_model=%s, distribution=%s, lag_fun=%s\n",
                  model_params$univ_model, model_params$distribution, model_params$lag_fun))
      cat("======================================================================\n")
      print(format(results_df, justify = "left"), row.names = FALSE)
      cat("----------------------------------------------------------------------\n")
      cat(sprintf(" Log-Likelihood: %0.4f\n", ll))
      cat(sprintf(" AIC           : %0.4f\n", aic))
      cat(sprintf(" BIC           : %0.4f\n", bic))
      cat(sprintf(" Số quan sát   : %d\n", n_obs))
      cat("======================================================================\n")
      cat(" Ghi chú: *** p<0.01, ** p<0.05, * p<0.1\n")
    
    # Trả về một đối tượng để sử dụng nếu cần
    return(invisible(list(estimates = results_df, stats = c(LL=ll, AIC=aic, BIC=bic))))
  }
  
# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 9. THỰC THI MÔ HÌNH VÀ PHÂN TÍCH ------------------------------------------ #
# ------------------------------------------------------------------------------ #

# -- 9.1. Ước lượng mô hình DCC-MIDAS ------------------------------------------- #
    res <- run_dcc_estimation(
      asset1 = asset1,
      asset2 = asset2,
      rate_data = rate_data,
      midas_data = midas_data,
      midas1_col = midas1_col_spec,
      midas2_col = midas2_col_spec,
      do_arma = do_arma,
      K = K, N_c = N_c, K_c = K_c,
      univ_model = univ_model,
      distribution = distribution,
      lag_fun = lag_fun
    )
  
  # In kết quả tóm tắt và vẽ biểu đồ
    cat("\n======================================================================\n")
    cat("      KẾT QUẢ ƯỚC LƯỢNG: MÔ HÌNH DCC-MIDAS\n")
    cat("----------------------------------------------------------------------\n")
    cat(sprintf(" Cặp tài sản: %s - %s\n", asset1, asset2))
    cat(sprintf(" Yếu tố MIDAS: '%s' & '%s' (từ sheet '%s')\n", midas1_col_spec, midas2_col_spec, midas_sheet_spec))
    cat("----------------------------------------------------------------------\n")
    cat(sprintf(" Tham số: K=%d, N_c=%d, K_c=%d\n", K, N_c, K_c))
    cat(sprintf(" Cấu hình: univ_model=%s, distribution=%s, lag_fun=%s\n",
                univ_model, distribution, lag_fun))
    cat("======================================================================\n")
    print(summary.dccmidas(res$fit))
    plot_dcc_results(estimation_result = res,
                     filename = file.path("img", paste0("DCC_MIDAS_", asset1, "_", asset2, ".png")))

# -- 9.2. Chuẩn bị dữ liệu cho DCC-MIDAS-X -------------------------------------- #
    dccx_data <- prepare_dcc_x_data(
      estimation_result = res,
      macro_data = shocks_data,
      shock_name = shock_name,
      K_x = K_x
    )

# -- 9.3. Ước lượng mô hình DCC-MIDAS-X ----------------------------------------- #
  # Thiết lập giá trị khởi tạo cho các tham số
  par_init <- c(a = 0.03, b = 0.9, m = 0, theta = -0.05, w = 1.1)
  
  # Khai báo hàm điều kiện ràng buộc
  ineq_fun <- function(params, res_data, X_data) {
    return(params[1] + params[2]) 
  }
  
  # Tối ưu hóa bằng solnp
  opt_results <- solnp(
    pars = par_init,                        # Giá trị khởi tạo
    fun = dcc_midas_x_loglik,               # Hàm log-likelihood
    ineqfun = ineq_fun, 
    ineqLB = c(0),                          # a + b > 0
    ineqUB = c(0.999),                      # a + b < 1
    LB = c(1e-4, 1e-4, -10, -10, 1.001),    # Ràng buộc dưới (a, b, m, theta, w)
    UB = c(0.2, 0.999, 10, 10, 50),         # Ràng buộc trên (a, b, m, theta, w)
    res_data = dccx_data$res_matrix,        # Ma trận phần dư
    X_data = dccx_data$X_mat                # Ma trận biến vĩ mô (X)
  )
  
  # Báo cáo kết quả
  report_dcc_midas_results(
    opt_res = opt_results, 
    res_matrix = dccx_data$res_matrix, 
    x_matrix = dccx_data$X_mat,
    model_params = list(
      asset1 = asset1,
      asset2 = asset2,
      shock_name = shock_name,
      K = K,
      N_c = N_c,
      K_c = K_c,
      K_x = K_x,
      univ_model = univ_model,
      distribution = distribution,
      lag_fun = lag_fun
    )
  )