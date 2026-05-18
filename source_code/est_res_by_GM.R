# ============================================================================== #
# ------------                                                       ----------- #
# ------------    ƯỚC LƯỢNG PHẦN DƯ & BIẾN ĐỘNG TỪ GARCH-MIDAS       ----------- #
# ------------                                                       ----------- #
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
  library(xts)
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

  # ============================================================================ #
  # -- CẤU HÌNH MÔ HÌNH ------------------------------------------------------- #
  # ============================================================================ #
    # CHỈ ĐỊNH CẤU HÌNH MÔ HÌNH
    model_spec = "GM" # Ví dụ: "GM", "GM2M", "DAGM", "DAGMM2M", "GMX", "DAGMX"
    skew_spec = "NO" # "NO" (đối xứng), "YES" (GJR)
    lag_fun_spec = "Beta" # "Beta" hoặc "Almon"
    
    # LỌC TRƯỚC LỢI SUẤT (PRE-FILTERING)
    do_arma_spec = TRUE # Đặt là TRUE nếu bạn muốn loại bỏ tự tương quan trong chuỗi lợi suất bằng ARMA

    # LỰA CHỌN PHƯƠNG PHÁP LỌC TRƯỚC LỢI SUẤT
    do_arma_approach = 2 # 1. ARMA(1,1); 2. auto.arima

    # CHỈ ĐỊNH BIẾN GARCH-X (HỒI QUY TẦN SỐ CAO)
    # Đặt là NULL nếu mô hình không có thành phần X (ví dụ: GM, DAGM)
    garchx_sheet_spec = NULL # Ví dụ: "rv". Đặt là NULL nếu không dùng GARCH-X.
    garchx_col_spec = NULL   # Ví dụ: "VNIndex". Chỉ định tên cột của biến X.
    
    # CHỈ ĐỊNH CÁC BIẾN MIDAS
    # Biến MIDAS thứ nhất (bắt buộc)
    midas1_sheet_spec = "rv" # Sheet chứa các biến MIDAS
    midas1_freq_spec = "monthly" # Tần suất của biến MIDAS 1
    
    # Biến MIDAS thứ hai (tùy chọn, cho các mô hình 2 thành phần như GM2M, DAGM2M)
    # Đặt là NULL nếu dùng mô hình 1 thành phần (ví dụ: DAGM, GMX, ...)    
    midas2_sheet_spec = NULL 
    midas2_col_spec = NULL # Chỉ định tên cột của biến MIDAS 2.
    midas2_freq_spec = NULL # Tần suất của biến MIDAS 2

    # CÁC THAM SỐ KHÁC
    K_m_1 <- 36 # Độ trễ MIDAS cho yếu tố thứ nhất
    
    # Độ trễ MIDAS cho yếu tố thứ hai (chỉ dùng cho mô hình 2 thành phần)
    # Đặt là NULL nếu không dùng yếu tố thứ hai
    K_m_2 <- NULL 
    dist <- "norm" # Phân phối (norm, std, ged)

    # ĐỊNH NGHĨA CÁC TÀI SẢN CẦN CHẠY VÀ BIẾN MIDAS TƯƠNG ỨNG
    assets_to_process <- list(
      list(asset = "VNIndex", midas_col = "RV_VNIndex"),
      list(asset = "XAUUSD", midas_col = "RV_XAUUSD"),
      list(asset = "GC_F", midas_col = "RV_GC_F"),
      list(asset = "SJC", midas_col = "RV_SJC")
    )


# ------------------------------------------------------------------------------ #
# -- 3. KHAI BÁO THƯ VIỆN ------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(rumidas)
  library(xts)
  library(FinTS)
  library(forecast)

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 4. HÀM CHUẨN BỊ DỮ LIỆU --------------------------------------------------- #
# ------------------------------------------------------------------------------ #
  prepare_ugm_data <- function(asset, rate_data, K_1, midas1_data, midas1_col, midas1_freq, do_arma = FALSE, do_arma_approach = 1, garchx_data = NULL, garchx_col = NULL, K_2 = NULL, midas2_data = NULL, midas2_col = NULL, midas2_freq = NULL) {
    # --- 1. Tính toán ngày bắt đầu an toàn cho dữ liệu tần suất cao (Dynamic Buffer) ---
    calculate_hf_start_date <- function(midas_date, K_val, freq) {
      buffer_period <- switch(freq,
                             "monthly" = months(K_val + 1),
                             "quarterly" = months(3 * (K_val + 1)),
                             "yearly" = years(K_val + 1),
                             months(K_val + 1)) # Mặc định
      return(midas_date %m+% buffer_period)
    }

    midas1_xts_full <- na.omit(xts(midas1_data[[midas1_col]], order.by = midas1_data$Date))
    hf_start_date <- calculate_hf_start_date(start(midas1_xts_full), K_1, midas1_freq)

    midas2_xts_full <- NULL
    if (!is.null(midas2_data) && !is.null(midas2_col) && !is.null(K_2)) {
      midas2_xts_full <- na.omit(xts(midas2_data[[midas2_col]], order.by = midas2_data$Date))
      hf_start_date2 <- calculate_hf_start_date(start(midas2_xts_full), K_2, midas2_freq)
      hf_start_date <- max(hf_start_date, hf_start_date2) # Lấy ngày muộn hơn để đảm bảo cả 2 đều đủ dữ liệu
    }

    # --- 2. Lọc và chuẩn bị dữ liệu tần suất cao ---
    daily_ret_base <- xts(rate_data[[asset]], order.by = rate_data$Date)
    colnames(daily_ret_base) <- asset
    daily_ret_base <- window(daily_ret_base, start = hf_start_date) # Áp dụng bộ đệm động
    daily_ret_base <- na.omit(daily_ret_base) # Đảm bảo không có NA

    # Loại bỏ trung bình (demean) và tùy chọn lọc ARMA
    daily_ret_base <- daily_ret_base - mean(daily_ret_base, na.rm = TRUE)
    
    if (do_arma) {
      if (do_arma_approach == 1) {
        cat(">> Áp dụng bộ lọc ARMA(1,1) cố định cho chuỗi lợi suất.\n")
        arma_fit <- forecast::Arima(daily_ret_base, order = c(1, 0, 1), include.mean = TRUE)
        cat(">> Đã áp dụng mô hình ARMA(1,1).\n")
      } else if (do_arma_approach == 2) {
        cat(">> Tự động tìm và áp dụng bộ lọc ARMA cho chuỗi lợi suất (stationary=TRUE).\n")
        arma_fit <- forecast::auto.arima(daily_ret_base, stationary = TRUE, trace = FALSE, allowdrift = FALSE)
        selected_order <- arma_fit$arma[c(1, 6, 2)] # p, d, q
        cat(sprintf(">> Đã chọn mô hình ARIMA(%d,%d,%d).\n", selected_order[1], selected_order[2], selected_order[3]))
      } else {
        stop("Lỗi: do_arma_approach phải là 1 hoặc 2.")
      }
      daily_ret_base <- xts(as.numeric(residuals(arma_fit)), order.by = index(daily_ret_base))
      colnames(daily_ret_base) <- asset
      cat(">> Lọc ARMA hoàn tất.\n")
    }

    # Chuẩn bị danh sách kết quả trả về
    result_list <- list()

    # Xử lý biến GARCH-X (ind_ret) nếu có
    if (!is.null(garchx_data) && !is.null(garchx_col)) {
      garchx_xts <- xts(garchx_data[[garchx_col]], order.by = garchx_data$Date)
      daily_merged <- merge.xts(daily_ret_base, garchx_xts, join = "inner")
      result_list$daily_ret <- na.omit(daily_merged[,1])
      result_list$ind_ret <- na.omit(daily_merged[,2]) # Đây là đối tượng xts
      colnames(result_list$ind_ret) <- garchx_col
      cat(sprintf(">> Đã thêm biến GARCH-X '%s'.\n", garchx_col))
    } else {
      result_list$daily_ret <- na.omit(daily_ret_base)
    }

    # --- 3. Tạo ma trận MIDAS với dữ liệu tần suất cao đã được lọc ---
    result_list$mv_m <- mv_into_mat(result_list$daily_ret, midas1_xts_full, K_1, midas1_freq)
    
    if (!is.null(midas2_xts_full) && !is.null(K_2)) {
      result_list$mv_m_2 <- mv_into_mat(result_list$daily_ret, midas2_xts_full, K_2, midas2_freq)
      cat(sprintf(">> Đã thêm biến MIDAS thứ hai '%s'.\n", midas2_col))
    }
    
    cat(sprintf("\n>> Đã chuẩn bị dữ liệu cho tài sản '%s' với %d quan sát hàng ngày.\n", asset, nrow(result_list$daily_ret)))
    
    return(result_list)
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 5. HÀM ƯỚC LƯỢNG MÔ HÌNH univ_GM ------------------------------------------ #
# ------------------------------------------------------------------------------ #
  run_ugm_estimation <- function(prepared_data, asset, model_type, skew_type, K_1, distribution, lag_function, K_2 = NULL) {
    
    cat(sprintf("\n[ƯỚC LƯỢNG] Mô hình %s cho %s...\n", model_type, asset))
    cat(sprintf("Cấu hình: Skew=%s, K_1=%d, Distribution=%s, LagFun=%s\n", skew_type, K_1, distribution, lag_function))
    
    # Xây dựng danh sách các tham số cho hàm ugmfit
    args <- list(
      model = model_type,
      skew = skew_type,
      lag_fun = lag_function,
      daily_ret = prepared_data$daily_ret,
      mv_m = prepared_data$mv_m,
      K = K_1,
      distribution = distribution
    )
    
    # Tự động thêm thành phần GARCH-X nếu tồn tại
    if ("ind_ret" %in% names(prepared_data)) {
      args$ind_ret <- prepared_data$ind_ret
    }

    # Tự động thêm thành phần MIDAS thứ hai nếu tồn tại trong dữ liệu đã chuẩn bị
    if ("mv_m_2" %in% names(prepared_data) && !is.null(K_2)) {
      args$mv_m_2 <- prepared_data$mv_m_2
      args$K_2 <- K_2
    }
    
    # Ước lượng mô hình bằng cách gọi hàm ugmfit với các tham số đã được xây dựng
    fit <- tryCatch({
      do.call(ugmfit, args)
    }, error = function(e) {
      cat(">> [LỖI ƯỚC LƯỢNG]:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(fit)) {
      cat(">> [ƯỚC LƯỢNG] Thành công!\n")
    }
    
    return(fit)
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 6. HÀM TỔNG HỢP VÀ BÁO CÁO KẾT QUẢ ---------------------------------------- #
# ------------------------------------------------------------------------------ #
  report_ugm_results <- function(fit, model_params) {
    if (is.null(fit)) {
      cat("\nKhông có kết quả để báo cáo do lỗi ước lượng.\n")
      return(invisible(NULL))
    }
    
    # Trích xuất thông tin
    coefs <- fit$rob_coef_mat
    ll <- fit$loglik
    n_obs <- fit$obs
    n_param <- length(fit$est_pars)
    aic <- fit$inf_criteria[[1]]
    bic <- fit$inf_criteria[[2]]
    mse <- fit$loss_in_s[[1]]
    qlike <- fit$loss_in_s[[2]]
    
    # Chuyển ma trận hệ số thành data frame để thêm cột ý nghĩa
    coefs_df <- as.data.frame(coefs)
    
    # Thêm cột ý nghĩa thống kê ("Sig.")
    p_values <- coefs_df[, 4] # Cột thứ 4 là Pr(>|t|)
    coefs_df$`Sig.` <- ifelse(p_values < 0.01, "***",
                              ifelse(p_values < 0.05, "**",
                                     ifelse(p_values < 0.1, "*", "")))
    
    # In bảng kết quả
    cat("\n======================================================================\n")
    cat(sprintf(" KẾT QUẢ ƯỚC LƯỢNG: MÔ HÌNH %s\n", model_params$model))
    cat("----------------------------------------------------------------------\n")
    cat(sprintf(" Tài sản: %s\n", model_params$asset))
    cat("----------------------------------------------------------------------\n")
    
    config_string <- sprintf(" Cấu hình: Model=%s, Skew=%s, K_1=%d", model_params$model, model_params$skew, model_params$K_1)
    if (!is.null(model_params$K_2)) {
      config_string <- paste0(config_string, sprintf(", K_2=%d", model_params$K_2))
    }
    config_string <- paste0(config_string, sprintf(", Distribution=%s, LagFun=%s\n", model_params$distribution, model_params$lag_fun))
    cat(config_string)
    
    if (!is.null(model_params$garchx_sheet)) {
        cat(sprintf(" Biến GARCH-X: '%s' (từ sheet '%s')\n", model_params$garchx_col, model_params$garchx_sheet))
    }
    cat(sprintf(" Biến MIDAS 1: '%s' (từ sheet '%s', tần suất: %s)\n", 
                model_params$midas1_col, model_params$midas1_sheet, model_params$midas1_freq))
    
    if (!is.null(model_params$midas2_sheet)) {
      cat(sprintf(" Biến MIDAS 2: '%s' (từ sheet '%s', tần suất: %s)\n", 
                  model_params$midas2_col, model_params$midas2_sheet, model_params$midas2_freq))
    }
    cat("======================================================================\n")
    print(coefs_df)
    cat("----------------------------------------------------------------------\n")
    cat(sprintf(" Log-Likelihood: %0.4f\n", ll))
    cat(sprintf(" AIC           : %0.4f\n", aic))
    cat(sprintf(" BIC           : %0.4f\n", bic))
    cat(sprintf(" MSE (%%)       : %0.4f\n", mse))
    cat(sprintf(" Q-Likelihood  : %0.4f\n", qlike))
    cat(sprintf(" Số quan sát   : %d\n", n_obs))
    cat("======================================================================\n")
    cat(" Ghi chú: *** p<0.01, ** p<0.05, * p<0.1. Bảng trên hiển thị sai số chuẩn đã được điều chỉnh (robust standard errors).\n")
    
    return(invisible(fit))
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 7. NHẬP DỮ LIỆU ----------------------------------------------------------- #
# ------------------------------------------------------------------------------ #
    rate_data <- read_excel(file_path, sheet = "rate", guess_max = 100000) %>% mutate(Date = as.Date(Date))
    
    garchx_data <- NULL
    if (!is.null(garchx_sheet_spec)) {
        garchx_data <- read_excel(file_path, sheet = garchx_sheet_spec, guess_max = 100000) %>% mutate(Date = as.Date(Date))
    }
    midas1_data <- read_excel(file_path, sheet = midas1_sheet_spec, guess_max = 100000) %>% mutate(Date = as.Date(Date))
    
    midas2_data <- NULL
    if (!is.null(midas2_sheet_spec)) {
      midas2_data <- read_excel(file_path, sheet = midas2_sheet_spec, guess_max = 100000) %>% mutate(Date = as.Date(Date))
    }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 8. THỰC THI HÀNG LOẠT VÀ LƯU KẾT QUẢ -------------------------------------- #
# ------------------------------------------------------------------------------ #
    # Khởi tạo danh sách để lưu trữ phần dư và biến động
    all_residuals <- list()
    all_volatilities <- list()
    all_long_run_volatilities <- list()

    # Bắt đầu vòng lặp qua từng tài sản đã định nghĩa trong `assets_to_process`
    for (item in assets_to_process) {
      asset_name <- item$asset
      midas1_col_spec <- item$midas_col
      
      cat("\n\n======================================================================\n")
      cat(sprintf("   BẮT ĐẦU XỬ LÝ TÀI SẢN: %s (MIDAS: %s)\n", asset_name, midas1_col_spec))
      cat("======================================================================\n")

      # 1. Chuẩn bị dữ liệu
      prepared_model_data <- prepare_ugm_data(
        asset = asset_name,
        rate_data = rate_data,
        K_1 = K_m_1,
        midas1_data = midas1_data,
        midas1_col = midas1_col_spec,
        midas1_freq = midas1_freq_spec,
        do_arma = do_arma_spec,
        do_arma_approach = do_arma_approach,
        garchx_data = garchx_data,
        garchx_col = garchx_col_spec,
        K_2 = K_m_2,
        midas2_data = midas2_data,
        midas2_col = midas2_col_spec,
        midas2_freq = midas2_freq_spec
      )
      
      # 2. Ước lượng mô hình
      ugm_fit <- run_ugm_estimation(
        prepared_data = prepared_model_data,
        asset = asset_name,
        model_type = model_spec,
        skew_type = skew_spec,
        K_1 = K_m_1,
        distribution = dist,
        lag_function = lag_fun_spec,
        K_2 = K_m_2
      )
      
      # 3. Báo cáo kết quả
      report_ugm_results(
        fit = ugm_fit,
        model_params = list(asset = asset_name, model = model_spec, skew = skew_spec, lag_fun = lag_fun_spec,
                            garchx_sheet = garchx_sheet_spec, garchx_col = garchx_col_spec,
                            midas1_sheet = midas1_sheet_spec, midas1_col = midas1_col_spec, midas1_freq = midas1_freq_spec,
                            midas2_sheet = midas2_sheet_spec, midas2_col = midas2_col_spec, midas2_freq = midas2_freq_spec,
                            K_1 = K_m_1, K_2 = K_m_2, distribution = dist)
      )

      # 4. Trích xuất, kiểm định và lưu kết quả
      if (!is.null(ugm_fit)) {
        cat("\n--- KIỂM ĐỊNH CHẨN ĐOÁN TRÊN PHẦN DƯ CHUẨN HÓA ---\n")
        std_residuals <- prepared_model_data$daily_ret / ugm_fit$est_vol_in_s
        std_residuals <- std_residuals[is.finite(std_residuals)]
        colnames(std_residuals) <- asset_name
        
        lb_test <- Box.test(std_residuals, lag = 12, type = "Ljung-Box")
        cat("\n>>> Kiểm định Ljung-Box (tự tương quan) cho phần dư (lag=12):\n")
        print(lb_test)
        
        arch_test <- ArchTest(as.numeric(std_residuals), lags = 12)
        cat("\n>>> Kiểm định ARCH-LM cho phần dư (lag=12):\n")
        print(arch_test)
        
        # Thêm phần dư vào danh sách tổng
        all_residuals[[asset_name]] <- std_residuals
        cat(sprintf("\n>> Đã trích xuất và lưu trữ phần dư chuẩn hóa cho '%s'.\n", asset_name))
        
        # Trích xuất và lưu trữ biến động trong mẫu (in-sample volatility)
        in_sample_vol <- ugm_fit$est_vol_in_s
        colnames(in_sample_vol) <- asset_name
        all_volatilities[[asset_name]] <- in_sample_vol
        cat(sprintf(">> Đã trích xuất và lưu trữ biến động trong mẫu cho '%s'.\n", asset_name))
        
        # Trích xuất và lưu trữ biến động dài hạn trong mẫu (in-sample long-run volatility)
        long_run_vol <- ugm_fit$est_lr_in_s
        colnames(long_run_vol) <- asset_name
        all_long_run_volatilities[[asset_name]] <- long_run_vol
        cat(sprintf(">> Đã trích xuất và lưu trữ biến động dài hạn trong mẫu cho '%s'.\n", asset_name))
        
      } else {
        cat(sprintf("\n>> Bỏ qua lưu kết quả cho '%s' do lỗi ước lượng.\n", asset_name))
      }
    }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 9. TỔNG HỢP VÀ XUẤT KẾT QUẢ RA TỆP CSV ------------------------------------ #
# ------------------------------------------------------------------------------ #
    if (length(all_residuals) > 0) {
      # Xác định thư mục con dựa trên phương pháp ARMA
      if (do_arma_spec) {
        sub_dir <- ifelse(do_arma_approach == 1, "arma11", "autoarima")
      } else {
        sub_dir <- "no_arma" # Hoặc một tên mặc định khác nếu không lọc ARMA
      }
      
      # Tạo thư mục nếu nó chưa tồn tại
      dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
      
      cat("\n\n======================================================================\n")
      cat("   TỔNG HỢP VÀ LƯU CÁC PHẦN DƯ CHUẨN HÓA\n")
      cat("======================================================================\n")
      
      # Hợp nhất tất cả các chuỗi thời gian phần dư vào một đối tượng xts duy nhất
      merged_residuals_xts <- do.call(merge.xts, all_residuals)
      # Chuyển đổi đối tượng xts thành data.frame để lưu
      residuals_df <- data.frame(Date = index(merged_residuals_xts), coredata(merged_residuals_xts))
      # Xác định tên tệp và lưu
      output_filename <- file.path(sub_dir, "standardized_residuals.csv")
      write.csv(residuals_df, output_filename, row.names = FALSE, na = "")
      
      cat(sprintf(">> Đã lưu thành công %d quan sát của %d tài sản vào tệp '%s'.\n", 
                  nrow(residuals_df), ncol(residuals_df) - 1, output_filename))
    } else {
      cat("\nKhông có phần dư nào được tạo để lưu vào tệp CSV.\n")
    }

    if (length(all_volatilities) > 0) {
      # Xác định thư mục con dựa trên phương pháp ARMA
      if (do_arma_spec) {
        sub_dir <- ifelse(do_arma_approach == 1, "arma11", "autoarima")
      } else {
        sub_dir <- "no_arma" # Hoặc một tên mặc định khác nếu không lọc ARMA
      }
      
      # Tạo thư mục nếu nó chưa tồn tại
      dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
      
      cat("\n\n======================================================================\n")
      cat("   TỔNG HỢP VÀ LƯU BIẾN ĐỘNG TRONG MẪU (IN-SAMPLE VOLATILITY)\n")
      cat("======================================================================\n")
      
      # Hợp nhất tất cả các chuỗi thời gian biến động vào một đối tượng xts duy nhất
      merged_volatilities_xts <- do.call(merge.xts, all_volatilities)
      # Chuyển đổi đối tượng xts thành data.frame để lưu
      volatilities_df <- data.frame(Date = index(merged_volatilities_xts), coredata(merged_volatilities_xts))
      # Xác định tên tệp và lưu
      output_filename_vol <- file.path(sub_dir, "in_sample_volatility.csv")
      write.csv(volatilities_df, output_filename_vol, row.names = FALSE, na = "")
      
      cat(sprintf(">> Đã lưu thành công %d quan sát biến động của %d tài sản vào tệp '%s'.\n", 
                  nrow(volatilities_df), ncol(volatilities_df) - 1, output_filename_vol))
    } else {
      cat("\nKhông có dữ liệu biến động nào được tạo để lưu vào tệp CSV.\n")
    }

    if (length(all_long_run_volatilities) > 0) {
      # Xác định thư mục con dựa trên phương pháp ARMA
      if (do_arma_spec) {
        sub_dir <- ifelse(do_arma_approach == 1, "arma11", "autoarima")
      } else {
        sub_dir <- "no_arma" # Hoặc một tên mặc định khác nếu không lọc ARMA
      }
      
      # Tạo thư mục nếu nó chưa tồn tại
      dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
      
      cat("\n\n======================================================================\n")
      cat("   TỔNG HỢP VÀ LƯU BIẾN ĐỘNG DÀI HẠN TRONG MẪU (IN-SAMPLE LONG-RUN VOLATILITY)\n")
      cat("======================================================================\n")
      
      # Hợp nhất tất cả các chuỗi thời gian biến động vào một đối tượng xts duy nhất
      merged_lr_volatilities_xts <- do.call(merge.xts, all_long_run_volatilities)
      # Chuyển đổi đối tượng xts thành data.frame để lưu
      lr_volatilities_df <- data.frame(Date = index(merged_lr_volatilities_xts), coredata(merged_lr_volatilities_xts))
      # Xác định tên tệp và lưu
      output_filename_lr_vol <- file.path(sub_dir, "in_sample_long_run_volatility.csv")
      write.csv(lr_volatilities_df, output_filename_lr_vol, row.names = FALSE, na = "")
      
      cat(sprintf(">> Đã lưu thành công %d quan sát biến động dài hạn của %d tài sản vào tệp '%s'.\n", 
                  nrow(lr_volatilities_df), ncol(lr_volatilities_df) - 1, output_filename_lr_vol))
    } else {
      cat("\nKhông có dữ liệu biến động dài hạn nào được tạo để lưu vào tệp CSV.\n")
    }
