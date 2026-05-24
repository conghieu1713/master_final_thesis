# ============================================================================== #
# ------------                                                       ----------- #
# ------------              ƯỚC LƯỢNG MÔ HÌNH DCC-MIDAS              ----------- #
# ------------                 (PHIÊN BẢN TÙY CHỈNH)                 ----------- #
# ------------                                                       ----------- #
# ============================================================================== #

# Tác giả: Công Hiếu, NGUYỄN (524102110660)
# Học viên Cao học, Đại học Kinh tế TP. Hồ Chí Minh


# ------------------------------------------------------------------------------ #
# -- 1. LÀM SẠCH MÔI TRƯỜNG VÀ KHAI BÁO THƯ VIỆN ------------------------------- #
# ------------------------------------------------------------------------------ #
# Làm sạch môi trường làm việc
  rm(list = ls(all.names = TRUE))
  graphics.off()
  library(xts)     # Khai báo một thư viện để lệnh gỡ bỏ thư viện không gặp lỗi
  invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), 
                 detach, character.only = TRUE, unload = TRUE, force = TRUE))

# Khai báo và tải các thư viện cần thiết
  library(Rsolnp)
  library(xts)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  

# ------------------------------------------------------------------------------ #
# -- 2. THIẾT LẬP CÁC THAM SỐ CHO MÔ HÌNH -------------------------------------- #
# ------------------------------------------------------------------------------ #
# Lựa chọn kịch bản làm mềm return trong GARCH-MIDAS
  do_arma_approach <- 2 # 1. ARMA(1,1), 2. auto.arima

# Cấu hình các cặp tài sản và tham số MIDAS tương ứng
# Bạn có thể tùy chỉnh N_c và K_c cho từng cặp tại đây
asset_pairs_config <- list(
  list(pair = c("VNIndex", "XAUUSD"), N_c = 21, K_c = 756),
  list(pair = c("VNIndex", "GC_F"),   N_c = 21, K_c = 756),
  list(pair = c("VNIndex", "SJC"),    N_c = 21, K_c = 504)
)

# Thiết lập giá trị khởi tạo cho thuật toán tối ưu hóa
# params = c(a, b, w2)
  init_params <- c(a = 0.03, b = 0.7, w2 = 2.0)

# Thư mục làm việc
  dir <- "D:/Post-graduate programme/Final Thesis/Model"
  setwd(dir)

# Tệp phần dư chuẩn hóa
  if (do_arma_approach == 1) {
     file_path <- "arma11/standardized_residuals.csv"
     con_vol_file_path <- "arma11/in_sample_volatility.csv"
  } else if (do_arma_approach == 2) {
     file_path <- "autoarima/standardized_residuals.csv"
     con_vol_file_path <- "autoarima/in_sample_volatility.csv"
  }

# Thiết lập ràng buộc (Bounds) cho các tham số
# Ràng buộc dưới: a > 0, b > 0, w2 > 1
  lower_bounds <- c(0.001, 0.001, 1.001)
# Ràng buộc trên: a < 1, b < 1, w2 < 50
  upper_bounds <- c(0.999, 0.999, 50.0)


# ------------------------------------------------------------------------------ #
# -- 3. ĐỊNH NGHĨA CÁC HÀM CHỨC NĂNG ------------------------------------------- #
# ------------------------------------------------------------------------------ #

# -- 3.1. HÀM NHẬP VÀ TIỀN XỬ LÝ DỮ LIỆU --------------------------------------- #
load_and_prepare_data <- function(file_path, asset1, asset2) {
  cat(">> Đang đọc dữ liệu từ:", file_path, "\n")
  data <- read.csv(file_path, header = TRUE)

  # Kiểm tra xem các cột tài sản có tồn tại không
  if (!asset1 %in% colnames(data) || !asset2 %in% colnames(data)) {
    stop(sprintf("Lỗi: Tên tài sản '%s' hoặc '%s' không tồn tại trong tệp dữ liệu.", asset1, asset2))
  }

  # Chuyển đổi dữ liệu thành đối tượng xts và tạo ma trận phần dư Z_t
  z_asset1 <- xts(data[[asset1]], order.by = as.Date(data$Date))
  colnames(z_asset1) <- asset1
  z_asset2 <- xts(data[[asset2]], order.by = as.Date(data$Date))
  colnames(z_asset2) <- asset2
  z_matrix <- na.omit(as.matrix(cbind(z_asset1, z_asset2)))
  
  cat(sprintf(">> Đã tạo ma trận phần dư Z_t cho %s và %s với %d quan sát.\n", asset1, asset2, nrow(z_matrix)))
  return(z_matrix)
}

# -- 3.2. HÀM TÍNH TRỌNG SỐ BETA MIDAS ----------------------------------------- #
beta_weight <- function(K, w2) {
  k <- 1:K
  w <- (1 - k/K)^(w2 - 1)
  return(w / sum(w))
}

# -- 3.3. HÀM ĐỐI LOG-LIKELIHOOD CHO DCC-MIDAS ----------------------------------- #
# Mô tả: Tính toán giá trị âm của hàm log-likelihood cho mô hình DCC-MIDAS
#        Hàm này được sử dụng làm hàm mục tiêu cho thuật toán tối ưu hóa
#
# Tham số:
#   theta_params: Vector các tham số cần ước lượng (a, b, w2)
#   Z: Ma trận phần dư chuẩn hóa (T x 2)
#   N_c: Cửa sổ tương quan cục bộ (ngày)
#   K_c: Số độ trễ cho thành phần dài hạn (ngày)
#   Trả về: Giá trị đối của log-likelihood
dcc_midas_likelihood <- function(theta_params, Z, N_c, K_c) {
  
  a  <- theta_params[1]
  b  <- theta_params[2]
  w2 <- theta_params[3]
  
  T_obs <- nrow(Z)
  LL <- 0
  
  # Khởi tạo ma trận Q bằng ma trận tương quan mẫu
  Q <- cor(Z) 
  
  # 1. Tính chuỗi tương quan cục bộ (Local Correlation - c_t)
  c_t <- numeric(T_obs)
  for (t in (N_c + 1):T_obs) {
    z1 <- Z[(t - N_c):(t - 1), 1]
    z2 <- Z[(t - N_c):(t - 1), 2]
    
    # Tính tương quan, thay thế bằng 0 nếu có lỗi (ví dụ: phương sai bằng 0)
    c_local <- suppressWarnings(cor(z1, z2))
    if(is.na(c_local)) c_local <- 0 
    c_t[t] <- c_local
  }
  
  # 2. Tính tương quan dài hạn (rho_long) qua hàm MIDAS
  weights <- beta_weight(K_c, w2)
  rho_long <- numeric(T_obs)
  
  # Điểm bắt đầu để có đủ dữ liệu trễ cho cả c_t và rho_long
  start_t <- N_c + K_c + 1
  
  for (t in start_t:T_obs) {
    # Tích chập tương quan cục bộ với trọng số Beta
    rho_long[t] <- sum(weights * c_t[(t-1):(t-K_c)])
    
    # Giới hạn giá trị của rho_long để đảm bảo ma trận tương quan xác định dương
    rho_long[t] <- max(min(rho_long[t], 0.999), -0.999)
  }
  
  # 3. Vòng lặp chính để cập nhật Q_t và tính log-likelihood
  for (t in start_t:T_obs) {
    
    # Ma trận tương quan dài hạn tại thời điểm t
    R_bar <- matrix(c(1, rho_long[t], rho_long[t], 1), nrow = 2, ncol = 2)
    
    # Cập nhật hiệp phương sai giả Q_t theo phương trình DCC
    z_lag <- matrix(Z[t-1, ], ncol = 1)
    Q <- (1 - a - b) * R_bar + a * (z_lag %*% t(z_lag)) + b * Q
    
    # Chuẩn hóa Q_t thành ma trận tương quan R_t
    Q_diag_inv <- diag(1 / sqrt(diag(Q))) 
    R_t <- Q_diag_inv %*% Q %*% Q_diag_inv
    
    # Lấy phần dư tại thời điểm t
    z_t <- matrix(Z[t, ], ncol = 1)
    
    # Tính định thức của R_t, phạt hàm mục tiêu nếu ma trận lỗi
    det_R <- det(R_t)
    if (det_R <= 1e-10 || is.na(det_R)) return(1e10) # Trả về giá trị phạt lớn
    
    # Cộng dồn vào tổng Log-Likelihood
    term1 <- log(det_R)
    term2 <- t(z_t) %*% solve(R_t) %*% z_t
    LL <- LL - 0.5 * (term1 + term2)
  }
  
  # Trả về giá trị âm của Log-Likelihood (vì solnp tìm cực tiểu)
  return(as.numeric(-LL))
}

# -- 3.4. HÀM RÀNG BUỘC BẤT PHƯƠNG TRÌNH --------------------------------------- #
# Mô tả: Định nghĩa hàm ràng buộc bất phương trình cho solnp (a + b < 1).
ineq_fun <- function(theta_params, Z, N_c, K_c) {
  return(theta_params[1] + theta_params[2])
}

# -- 3.5. HÀM ƯỚC LƯỢNG MÔ HÌNH ------------------------------------------------ #
estimate_dcc_midas <- function(z_matrix, init_params, lower_bounds, upper_bounds, N_c, K_c) {
  cat("\n>> Bắt đầu quá trình tối ưu hóa DCC-MIDAS. Vui lòng đợi...\n")
  
  optim_results <- solnp(
    pars = init_params, 
    fun = dcc_midas_likelihood, 
    ineqfun = ineq_fun, 
    ineqLB = 0.001,
    ineqUB = 0.999,
    LB = lower_bounds, 
    UB = upper_bounds, 
    Z = z_matrix, 
    N_c = N_c, 
    K_c = K_c
  )
  
  cat(">> Quá trình tối ưu hóa hoàn tất!\n")
  return(optim_results)
}

# -- 3.6. HÀM BÁO CÁO KẾT QUẢ -------------------------------------------------- #
report_dcc_results <- function(optim_results, z_matrix, N_c, K_c, asset1, asset2) {
  
  # Hàm phụ trợ tạo sao ý nghĩa thống kê
  get_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return(" ")
  }
  
  # 1. Trích xuất tham số và tính toán thống kê
  estimates <- optim_results$pars
  n_pars <- length(estimates)
  n_obs <- nrow(z_matrix)

  cat(">> Đang tính toán ma trận Hessian (numDeriv) để tìm sai số chuẩn...\n")
  hessian_matrix <- numDeriv::hessian(
    func = dcc_midas_likelihood,
    x = estimates,
    Z = z_matrix, 
    N_c = N_c, 
    K_c = K_c
  )
  
  cov_matrix <- tryCatch(
    solve(hessian_matrix),
    error = function(e) {
      warning("Ma trận Hessian không thể nghịch đảo. Sai số chuẩn có thể không chính xác. Lỗi: ", e$message)
      return(matrix(NA, nrow = n_pars, ncol = n_pars))
    }
  )
  
  std_errors <- suppressWarnings(sqrt(diag(cov_matrix)))
  t_values <- estimates / std_errors
  p_values <- 2 * (1 - pnorm(abs(t_values)))
  sig_stars <- sapply(p_values, get_stars)
  
  # 2. Tạo bảng kết quả
  dcc_results_table <- data.frame(
    "Estimate"   = round(estimates, 6),
    "Std. Error" = round(std_errors, 6),
    "t value"    = round(t_values, 4),
    "Pr(>|t|)"   = signif(p_values, 4),
    "Sig."       = sig_stars,
    check.names  = FALSE
  )
  rownames(dcc_results_table) <- c("alpha (a)", "beta (b)", "omega (w2)")
  
  # 3. In bảng kết quả ra Console
  cat("\n========================================================================\n")
  cat("       KẾT QUẢ ƯỚC LƯỢNG: MÔ HÌNH DCC-MIDAS (CUSTOM OPTIMIZATION)      \n")
  cat("========================================================================\n")
  cat(sprintf(" Cặp tài sản: %s - %s\n", asset1, asset2))
  cat("------------------------------------------------------------------------\n")
  print(dcc_results_table, row.names = TRUE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  # Tính toán và in các chỉ số khác
  loglik <- -tail(optim_results$values, 1)
  aic <- 2 * n_pars - 2 * loglik
  bic <- n_pars * log(n_obs) - 2 * loglik
  
  cat(sprintf("Log-Likelihood: %0.4f\n", loglik))
  cat(sprintf("AIC           : %0.4f\n", aic))
  cat(sprintf("BIC           : %0.4f\n", bic))
  cat(sprintf("Số quan sát (T): %d\n", n_obs))
  cat(sprintf("Tham số mô hình: N_c = %d, K_c = %d\n", N_c, K_c))
  cat("========================================================================\n")
  
  return(invisible(dcc_results_table))
}

# -- 3.7. HÀM VẼ BIỂU ĐỒ KẾT QUẢ TÙY CHỈNH ------------------------------------- #
plot_dcc_custom_results <- function(optim_results, z_matrix, con_vol_file_path, asset1, asset2, N_c, K_c, filename = NULL) {
  
  # ====================================================
  # 1. TÍNH TOÁN TƯƠNG QUAN ĐỘNG (THEO test_con_corr.r)
  # ====================================================
  cat("\n>> Đang tính toán tương quan động ngắn hạn và dài hạn...\n")
  a_opt  <- optim_results$pars[1]
  b_opt  <- optim_results$pars[2]
  w2_opt <- optim_results$pars[3]
  
  T_obs <- nrow(z_matrix)
  
  rho_dynamic  <- rep(NA, T_obs) # Tương quan động ngắn hạn
  rho_long_run <- rep(NA, T_obs) # Tương quan dài hạn (MIDAS)
  
  # Tính chuỗi c_t (Tương quan cục bộ)
  c_t <- numeric(T_obs)
  for (t in (N_c + 1):T_obs) {
    z1 <- z_matrix[(t - N_c):(t - 1), 1]
    z2 <- z_matrix[(t - N_c):(t - 1), 2]
    c_local <- suppressWarnings(cor(z1, z2))
    if(is.na(c_local)) c_local <- 0 
    c_t[t] <- c_local
  }
  
  # Tính Tương quan dài hạn (rho_long_run)
  weights <- beta_weight(K_c, w2_opt)
  start_t <- N_c + K_c + 1
  
  for (t in start_t:T_obs) {
    rho_long_run[t] <- sum(weights * c_t[(t-1):(t-K_c)])
    rho_long_run[t] <- max(min(rho_long_run[t], 0.999), -0.999)
  }
  
  # Tính Tương quan động ngắn hạn (rho_dynamic)
  Q <- cor(z_matrix) # Khởi tạo Q
  
  for (t in start_t:T_obs) {
    R_bar <- matrix(c(1, rho_long_run[t], rho_long_run[t], 1), nrow = 2, ncol = 2)
    z_lag <- matrix(z_matrix[t-1, ], ncol = 1)
    Q <- (1 - a_opt - b_opt) * R_bar + a_opt * (z_lag %*% t(z_lag)) + b_opt * Q
    Q_diag_inv <- diag(1 / sqrt(diag(Q))) 
    R_t <- Q_diag_inv %*% Q %*% Q_diag_inv
    rho_dynamic[t] <- R_t[1, 2]
  }
  
  # ====================================================
  # 2. TẢI BIẾN ĐỘNG CÓ ĐIỀU KIỆN TỪ TỆP CSV
  # ====================================================
  cat(">> Đang tải dữ liệu biến động có điều kiện từ:", con_vol_file_path, "\n")
  vol_data <- read.csv(con_vol_file_path) %>%
    mutate(Date = as.Date(Date)) %>%
    select(Date, Vol_Asset1 = all_of(asset1), Vol_Asset2 = all_of(asset2))
  
  # ====================================================
  # 3. TẠO DATA FRAME ĐỂ VẼ BIỂU ĐỒ
  # ====================================================
  plot_df <- data.frame(
    Date = as.Date(rownames(z_matrix)),
    SR_Corr = as.numeric(rho_dynamic),
    LR_Corr = as.numeric(rho_long_run)
  ) %>%
    left_join(vol_data, by = "Date")
  
  # ====================================================
  # 4. LÀM SẠCH DỮ LIỆU - LOẠI BỎ GIAI ĐOẠN "BURN-IN"
  # ====================================================
  burn_in_index <- which(plot_df$SR_Corr != 0 & !is.na(plot_df$SR_Corr))[1]
  if(is.na(burn_in_index)) burn_in_index <- 1
  plot_df_clean <- plot_df %>% slice(burn_in_index:nrow(plot_df))
  
  # ====================================================
  # 5. THIẾT LẬP THEME VÀ VẼ BIỂU ĐỒ
  # ====================================================
  theme1 <- theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.major = element_line(linetype = "dotted", color = "darkgray"),
      panel.grid.minor.x = element_line(linetype = "dotted", color = "darkgray"),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(r = 5)),
      axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1, margin = margin(t = 5)),
      axis.ticks.length = unit(-0.15, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  # Biểu đồ A: Biến động tài sản 1
  p1 <- ggplot(plot_df_clean, aes(x = Date, y = Vol_Asset1)) +
    geom_line(color = "black", linewidth = 0.75) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "2 years", date_minor_breaks = "1 year") +
    labs(title = asset1) + theme1 + ylim(0,8)
  
  # Panel B: Tương quan động
  p2 <- ggplot(plot_df_clean, aes(x = Date)) +
    geom_line(aes(y = SR_Corr), color = "black", alpha = 0.75, linewidth = 0.6) +
    geom_line(aes(y = LR_Corr), color = "red", linewidth = 1) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "2 years", date_minor_breaks = "1 year") +
    labs(title = paste(asset1, "-", asset2)) + theme1 + ylim(-0.3, 0.3)
  
  # Panel C: Biến động tài sản 2
  p3 <- ggplot(plot_df_clean, aes(x = Date, y = Vol_Asset2)) +
    geom_line(color = "black", linewidth = 0.75) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "2 years", date_minor_breaks = "1 year") +
    labs(title = asset2) + theme1 + ylim(0,8)
  
  # ====================================================
  # 6. KẾT HỢP VÀ XUẤT BIỂU ĐỒ
  # ====================================================
  des <- c(
    area(1,1),
    area(1,2),
    area(2,2)
  )
  final_plot <- p1 + p2 + p3 + plot_layout(design = des)
  
  print(final_plot)
  
  if (!is.null(filename)) {
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename, plot = final_plot, width = 12, height = 8, dpi = 300)
    cat(sprintf(">> Đã lưu biểu đồ tại: %s\n", filename))
  }
  
  return(invisible(final_plot))
}


# ------------------------------------------------------------------------------ #
# -- 4. THỰC THI MÔ HÌNH CHO NHIỀU CẶP TÀI SẢN --------------------------------- #
# ------------------------------------------------------------------------------ #

# Vòng lặp qua từng cấu hình cặp tài sản
for (config in asset_pairs_config) {
  asset1 <- config$pair[1]
  asset2 <- config$pair[2]
  N_c    <- config$N_c
  K_c    <- config$K_c
  
  cat(paste("\n\n========================================================================"))
  cat(paste("\n BẮT ĐẦU XỬ LÝ CẶP TÀI SẢN:", asset1, "-", asset2, "\n"))
  cat(sprintf(">> Tham số được sử dụng: N_c = %d, K_c = %d\n", N_c, K_c))
  cat(paste("========================================================================\n"))

  # 1. Tải và chuẩn bị dữ liệu
  z_matrix <- load_and_prepare_data(file_path = file_path, asset1 = asset1, asset2 = asset2)

  # 2. Ước lượng mô hình
  optim_results <- estimate_dcc_midas(
    z_matrix = z_matrix,
    init_params = init_params,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    N_c = N_c,
    K_c = K_c
  )

  # 3. Báo cáo kết quả
  report_dcc_results(
    optim_results = optim_results,
    z_matrix = z_matrix,
    N_c = N_c,
    K_c = K_c,
    asset1 = asset1,
    asset2 = asset2
  )

  # 4. Vẽ biểu đồ kết quả
  plot_dcc_custom_results(
    optim_results = optim_results,
    z_matrix = z_matrix,
    con_vol_file_path = con_vol_file_path,
    asset1 = asset1,
    asset2 = asset2,
    N_c = N_c,
    K_c = K_c,
    filename = file.path("img", paste0("Custom_DCC_MIDAS_", asset1, "_", asset2, "_", do_arma_approach, ".png"))
  )
  
  cat(paste("\n HOÀN TẤT XỬ LÝ CẶP TÀI SẢN:", asset1, "-", asset2, "\n"))
}