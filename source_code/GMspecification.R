# ============================================================================== #
# ------------                                                       ----------- #
# ------------       LỰA CHỌN ĐỘ TRỄ CHO univariate GARCH-MIDAS      ----------- #
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
  # Khai báo một thư viện để lệnh gỡ bỏ thư viện không gặp lỗi
  library(xts)
  # Gỡ bỏ các thư viện đang hoạt động để tránh xung đột
  invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), 
                   detach, character.only = TRUE, unload = TRUE, force = TRUE))

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 2. KHAI BÁO THƯ VIỆN ------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(rumidas)
  library(tidyr)
  library(forecast)
  library(xts)

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 3. HÀM CHUẨN BỊ DỮ LIỆU CHO MỖI TÀI SẢN ----------------------------------- #
# ------------------------------------------------------------------------------ #
  prepare_asset_data <- function(asset_name, rv_name, rate_sheet, rv_sheet, K_max, midas_freq, do_arma = FALSE, do_arma_approach = 1) {
    # --- 1. Tính toán ngày bắt đầu an toàn cho dữ liệu tần suất cao (Dynamic Buffer) ---
    calculate_hf_start_date <- function(midas_date, K_val, freq) {
      buffer_period <- switch(freq,
                             "monthly" = months(K_val + 1),
                             "quarterly" = months(3 * (K_val + 1)),
                             "yearly" = years(K_val + 1),
                             months(K_val + 1)) # Mặc định là tháng
      return(midas_date %m+% buffer_period)
    }

    rv_xts_full <- na.omit(xts(rv_sheet[[rv_name]], order.by = as.Date(rv_sheet$Date)))
    colnames(rv_xts_full) <- rv_name
    hf_start_date <- calculate_hf_start_date(start(rv_xts_full), K_max, midas_freq)

    # --- 2. Lọc và chuẩn bị dữ liệu tần suất cao ---
    daily_ret_base <- xts(rate_sheet[[asset_name]], order.by = as.Date(rate_sheet$Date))
    colnames(daily_ret_base) <- asset_name
    daily_ret_base <- window(daily_ret_base, start = hf_start_date) # Áp dụng bộ đệm
    daily_ret_base <- na.omit(daily_ret_base) # Đảm bảo không có NA

    # Trừ đi trung bình của chuỗi lợi suất (bước tiền xử lý phổ biến)
    daily_ret_base <- daily_ret_base - mean(daily_ret_base, na.rm = TRUE)
    
    # Tùy chọn lọc ARMA
    if (do_arma) {
      cat(sprintf("     >> Áp dụng bộ lọc ARMA cho '%s'...\n", asset_name))
      if (do_arma_approach == 1) {
        arma_fit <- forecast::Arima(daily_ret_base, order = c(1, 0, 1), include.mean = TRUE)
        cat("     >> Đã áp dụng mô hình ARMA(1,1).\n")
      } else if (do_arma_approach == 2) {
        arma_fit <- forecast::auto.arima(daily_ret_base, stationary = TRUE, trace = FALSE, allowdrift = FALSE)
        selected_order <- arma_fit$arma[c(1, 6, 2)] # p, d, q
        cat(sprintf("     >> Đã chọn mô hình ARIMA(%d,%d,%d).\n", selected_order[1], selected_order[2], selected_order[3]))
      } else {
        stop("Lỗi: do_arma_approach phải là 1 hoặc 2.")
      }
      daily_ret_base <- xts(as.numeric(residuals(arma_fit)), order.by = index(daily_ret_base))
      colnames(daily_ret_base) <- asset_name
      cat("     >> Lọc ARMA hoàn tất.\n")
    }

    # --- 3. Đồng bộ hóa dữ liệu tần suất cao và thấp ---
    merged_data <- na.omit(merge.xts(daily_ret_base, rv_xts_full))

    return(list(
      daily_ret = merged_data[, 1],
      midas_var = merged_data[, 2],
      obs_count = nrow(merged_data)
    ))
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 4. HÀM ƯỚC LƯỢNG MÔ HÌNH VÀ TRÍCH XUẤT KẾT QUẢ ---------------------------- #
# ------------------------------------------------------------------------------ #
  fit_and_extract_results <- function(daily_ret, mv_m, asset, asym_label, K_val, lag_name) {
    skew_spec <- ifelse(asym_label == "Có yếu tố bất đối xứng", "YES", "NO")

    # Ước lượng mô hình bằng ugmfit từ gói rumidas
    fit <- tryCatch({
      ugmfit(model = "GM", skew = skew_spec, lag_fun = "Beta",
             daily_ret = daily_ret, mv_m = mv_m, K = K_val, distribution = "norm")
    }, error = function(e) {
      cat("     Lỗi khi xử lý", asset, "với K =", K_val, "và Skew =", skew_spec, ":", e$message, "\n")
      return(NULL)
    })

    if (is.null(fit)) {
      return(NULL)
    }

    # Trích xuất LogLikelihood, AIC và BIC từ đối tượng fit của ugmfit
    ll <- fit$loglik
    aic <- fit$inf_criteria[[1]]
    bic <- fit$inf_criteria[[2]]

    # Trả về kết quả dạng data frame
    return(data.frame(
      Asset = asset,
      Specification = asym_label,
      Lag_Length = lag_name,
      LL = round(ll, 2),
      AIC = round(aic, 2),
      BIC = round(bic, 2)
    ))
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 5. HÀM TỔNG HỢP, XỬ LÝ VÀ XUẤT KẾT QUẢ ------------------------------------ #
# ------------------------------------------------------------------------------ #
  process_and_export_specifications <- function(asset_list, lags_K, lag_names, rate_sheet, rv_sheet, out_filename, do_arma, do_arma_approach) {
    all_results <- list()
    
    # Xác định độ trễ lớn nhất để tính bộ đệm
    K_max <- max(lags_K)
    midas_freq <- "monthly" # Tần suất của biến MIDAS (RV)
    
    # Bắt đầu xử lý...
    for (asset in asset_list) {
      cat("\nĐang xử lý tài sản:", asset, "...\n")
      
      # Tạo tên biến realized volatility (rv) tương ứng một cách tự động
      rv_asset <- paste0("RV_", asset)
      # Chuẩn bị dữ liệu với bộ đệm động và tùy chọn lọc ARMA
      prepared_data <- prepare_asset_data(asset, rv_asset, rate_sheet, rv_sheet, K_max, midas_freq, do_arma, do_arma_approach)
      
      if (nrow(prepared_data$daily_ret) == 0) {
        cat("     Không có dữ liệu cho tài sản:", asset, ". Bỏ qua.\n")
        next
      }
      
      # Vòng lặp qua các đặc tả mô hình
      for (asym_label in c("Không có yếu tố bất đối xứng", "Có yếu tố bất đối xứng")) {
        for (i in seq_along(lags_K)) {
          K_val <- lags_K[i]
          lag_name <- lag_names[i]
          
          # Tạo ma trận MIDAS cho K hiện tại.
          # Dựa trên các tệp khác, tần suất của biến MIDAS (RV) là hàng tháng ("monthly").
          mv_m_matrix <- mv_into_mat(prepared_data$daily_ret, prepared_data$midas_var, K_val, "monthly")

          # Ước lượng và trích xuất kết quả
          result <- fit_and_extract_results(prepared_data$daily_ret, mv_m_matrix, asset, asym_label, K_val, lag_name)
          if (!is.null(result)) {
            all_results[[length(all_results) + 1]] <- result
          }
        }
      }
    }
    
    # Gộp thành bảng kết quả cuối cùng
    final_table <- bind_rows(all_results)
    
    # In kết quả ra màn hình console
    cat("\n# ------------------------------------------------------------------------------ #\n")
    cat("# -- KẾT QUẢ TỪ FILE:", toupper(basename(out_filename)), "\n")
    cat("# ------------------------------------------------------------------------------ #\n")
    print(final_table)
    cat("\n")
    
    # Đảm bảo thư mục lưu trữ tồn tại và xuất file CSV
    dir.create(dirname(out_filename), showWarnings = FALSE, recursive = TRUE)
    write.csv(final_table, out_filename, row.names = FALSE)
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 6. THIẾT LẬP VÀ NHẬP DỮ LIỆU ---------------------------------------------- #
# ------------------------------------------------------------------------------ #
  # Đường dẫn thư mục làm việc
  dir <- "D:/Post-graduate programme/Final Thesis/Model"   
  setwd(dir)
  
  # Đọc dữ liệu (chỉ đọc một lần để tối ưu)
  rate_data <- read_excel("official.xlsx", sheet = "rate") #%>% mutate(Date = as.Date(Date)) %>% filter(!is.na(Date))
  rv_data <- read_excel("official.xlsx", sheet = "rv") #%>% mutate(Date = as.Date(Date)) %>% filter(!is.na(Date))

  # LỌC TRƯỚC LỢI SUẤT (PRE-FILTERING) - Tương tự est_res_by_GM.R
  # Đặt là FALSE để giữ nguyên kịch bản gốc (không lọc)
  do_arma_spec = FALSE # Đặt là TRUE nếu bạn muốn loại bỏ tự tương quan trong chuỗi lợi suất bằng ARMA
  do_arma_approach = 2 # 1. ARMA(1,1); 2. auto.arima

  # Danh sách tài sản cần xử lý
  # asset_list <- c("VNIndex", "XAUUSD", "GC_F", "SJC", "BTC", "ETH", "BNB")
  asset_list <- c("VNIndex", "XAUUSD", "GC_F", "SJC")

  # Thiết lập các độ trễ cho mô hình MIDAS
  lags_K <- c(12, 24, 36)
  lag_names <- c("1 năm", "2 năm", "3 năm")

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 7. THỰC THI VÀ LƯU KẾT QUẢ ------------------------------------------------ #
# ------------------------------------------------------------------------------ #
  process_and_export_specifications(
    asset_list = asset_list,
    lags_K = lags_K,
    lag_names = lag_names,
    rate_sheet = rate_data,
    rv_sheet = rv_data,
    out_filename = file.path("DS Results", "GMspecification.csv"),
    do_arma = do_arma_spec,
    do_arma_approach = do_arma_approach
  )

  