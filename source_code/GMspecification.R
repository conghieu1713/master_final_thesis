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
  library(mfGARCH)
  library(tidyr)

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 3. HÀM CHUẨN BỊ DỮ LIỆU CHO MỖI TÀI SẢN ----------------------------------- #
# ------------------------------------------------------------------------------ #
  prepare_asset_data <- function(asset_name, rate_sheet, rv_sheet) {
    df_asset <- rate_sheet %>%
      select(Date, y = all_of(asset_name)) %>%
      inner_join(rv_sheet %>% select(Date, x = all_of(asset_name)), by = "Date") %>%
      rename(date = Date) %>%
      mutate(date = as.Date(date), year_month = floor_date(date, "month")) %>%
      drop_na() %>%
      arrange(date)
    return(df_asset)
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 4. HÀM ƯỚC LƯỢNG MÔ HÌNH VÀ TRÍCH XUẤT KẾT QUẢ ---------------------------- #
# ------------------------------------------------------------------------------ #
  fit_and_extract_results <- function(df_asset, asset, asym_label, K_val, lag_name) {
    is_asym <- ifelse(asym_label == "Có yếu tố bất đối xứng", TRUE, FALSE)
    
    # Ước lượng mô hình...
    fit <- tryCatch({
      fit_mfgarch(data = df_asset, y = "y", x = "x",
                  low.freq = "year_month", var.ratio.freq = "year_month",
                  K = K_val, gamma = is_asym)
    }, error = function(e) {
      cat("     Lỗi khi xử lý", asset, "với K =", K_val, "và Asym =", is_asym, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(fit)) {
      return(NULL)
    }
    
    # Tính toán LogLikelihood, AIC và BIC
    ll <- fit$llh
    k <- length(fit$par)
    aic <- 2 * k - 2 * ll
    bic <- fit$bic
    
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
  process_and_export_specifications <- function(asset_list, lags_K, lag_names, rate_sheet, rv_sheet, out_filename) {
    all_results <- list()
    
    # Bắt đầu xử lý...
    for (asset in asset_list) {
      cat("\nĐang xử lý tài sản:", asset, "...\n")
      
      # Chuẩn bị dữ liệu
      df_asset <- prepare_asset_data(asset, rate_sheet, rv_sheet)
      
      if (nrow(df_asset) == 0) {
        cat("     Không có dữ liệu cho tài sản:", asset, ". Bỏ qua.\n")
        next
      }
      
      # Vòng lặp qua các đặc tả mô hình
      for (asym_label in c("Không có yếu tố bất đối xứng", "Có yếu tố bất đối xứng")) {
        for (i in seq_along(lags_K)) {
          K_val <- lags_K[i]
          lag_name <- lag_names[i]
          
          # Ước lượng và trích xuất kết quả
          result <- fit_and_extract_results(df_asset, asset, asym_label, K_val, lag_name)
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
  rate_data <- read_excel("official.xlsx", sheet = "rate")
  rv_data <- read_excel("official.xlsx", sheet = "rv")

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
    out_filename = file.path("DS Results", "GMspecification.csv")
  )

  