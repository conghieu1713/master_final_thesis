# ============================================================================== #
# ------------                                                       ----------- #
# ------------                     THỐNG KÊ MÔ TẢ                    ----------- #
# ------------                                                       ----------- #  
# ============================================================================== #
                                                                                

# Tác giả: Công Hiếu, NGUYỄN (524102110660)
# Học viên Cao học, Trường Đại học Kinh tế TP.HCM


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
  library(moments)
  library(tseries)
  library(FinTS)
  library(dplyr)
  library(xts)

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 3. HÀM GÁN DẤU SAO Ý NGHĨA THỐNG KÊ --------------------------------------- #
# ------------------------------------------------------------------------------ #
  add_stars <- function(stat, p_val) {
    stars <- ifelse(p_val <= 0.01, "***",
                    ifelse(p_val <= 0.05, "**",
                           ifelse(p_val <= 0.1, "*", "")))
    return(paste0(sprintf("%.3f", stat), stars))
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 4. HÀM TÍNH TOÁN THỐNG KÊ MÔ TẢ ------------------------------------------- #
# ------------------------------------------------------------------------------ #
  calc_stats <- function(x) {
    # Loại bỏ giá trị khuyết (NA)
      x <- na.omit(x) 
    # Bỏ qua nếu dữ liệu có ít hơn 13 quan sát
      if(length(x) < 13) return(NULL)
      
    # Thống kê mô tả cơ bản
      mean_val <- mean(x)
      sd_val   <- sd(x)
      min_val  <- min(x)
      max_val  <- max(x)
      
    # Bỏ qua các kiểm định nếu phương sai bằng 0 (ngăn lỗi chia cho 0)
      if(sd_val == 0) return(NULL)
      
      skew_val <- skewness(x)
      kurt_val <- kurtosis(x)
      
    # Kiểm định JB (Phân phối chuẩn)
      jb <- jarque.bera.test(x)
      jb_str <- add_stars(jb$statistic, jb$p.value)
      
    # Các kiểm định nghiệm đơn vị và tính dừng (Ẩn cảnh báo p-value)
      suppressWarnings({
        # Kiểm định ADF
        adf <- adf.test(x)
        adf_str <- add_stars(adf$statistic, adf$p.value)
        
        # Kiểm định PP
        pp <- pp.test(x)
        pp_str <- add_stars(pp$statistic, pp$p.value)
        
        # Kiểm định KPSS
        kpss <- kpss.test(x)
        kpss_str <- add_stars(kpss$statistic, kpss$p.value)
      })
      
    # Kiểm định ARCH(12) (Phương sai sai số thay đổi)
      arch <- ArchTest(x, lags = 12)
      arch_str <- add_stars(arch$statistic, arch$p.value)
      
    # Kiểm định Ljung-Box LB(12) (Tự tương quan)
      lb <- Box.test(x, lag = 12, type = "Ljung-Box")
      lb_str <- add_stars(lb$statistic, lb$p.value)

    # Kiểm định Ljung-Box LB(12) trên chuỗi thời gian bình phương
      lb2 <- Box.test(x^2, lag = 12, type = "Ljung-Box")
      lb2_str <- add_stars(lb2$statistic, lb2$p.value)



    # Trả về kết quả dạng Data Frame
    return(data.frame(
      Mean = sprintf("%.3f", mean_val),
      `Std. Dev` = sprintf("%.3f", sd_val),
      Min = sprintf("%.3f", min_val),
      Max = sprintf("%.3f", max_val),
      Skew = sprintf("%.3f", skew_val),
      Kurtosis = sprintf("%.3f", kurt_val),
      `JB test` = jb_str,
      `ADF test` = adf_str,
      `PP test` = pp_str,
      `KPSS test` = kpss_str,
      `ARCH (12)` = arch_str,
      `LB (12)` = lb_str,
      `LB2 (12)` = lb2_str,
      check.names = FALSE
    ))
  }
 
# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 5. HÀM TỔNG HỢP VÀ XUẤT DỮ LIỆU (TỐI ƯU HÓA QUY TRÌNH) -------------------- #
# ------------------------------------------------------------------------------ #
  process_and_export <- function(data, cols, out_filename) {
    # Đang xử lý tính toán...
    results_list <- lapply(data[cols], function(col) calc_stats(as.numeric(col)))
    
    # Loại bỏ các cột NULL (nếu có lỗi phương sai = 0 hoặc thiếu dữ liệu)
    results_list <- Filter(Negate(is.null), results_list)
    
    # Gộp thành bảng
    final_table <- do.call(rbind, results_list)
    
    # Thêm cột Tên chuỗi (Return Series)
    final_table <- cbind(`Return Series` = rownames(final_table), final_table)
    rownames(final_table) <- NULL
    
    # In kết quả ra màn hình console với tiêu đề rõ ràng
    cat("\n# ------------------------------------------------------------------------------ #\n")
    cat("# -- KẾT QUẢ TỪ FILE:", toupper(basename(out_filename)), "\n")
    cat("# ------------------------------------------------------------------------------ #\n")
    print(final_table)
    cat("\n")
    
    # Đảm bảo thư mục lưu trữ tồn tại trước khi xuất file
    dir.create(dirname(out_filename), showWarnings = FALSE, recursive = TRUE)
    write.csv(final_table, out_filename, row.names = FALSE)
  }
  
# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 6. NHẬP DỮ LIỆU ----------------------------------------------------------- #
# ------------------------------------------------------------------------------ #
  # Đường dẫn thư mục làm việc (Sử dụng / để tránh lỗi đường dẫn)
    dir <- "D:/Post-graduate programme/Final Thesis/Model"   
    setwd(dir)
    
  # Đọc dữ liệu tần suất cao
    rate <- read_excel("official.xlsx", sheet = "rate", guess_max = 100000)
    rate$Date <- as.Date(rate$Date)
    
  # Đọc dữ liệu tần suất thấp
    shocks <- read_excel("official.xlsx", sheet = "shocks_lg", guess_max = 100000)
    shocks$Date <- as.Date(shocks$Date)

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 7. KẾT QUẢ THỐNG KÊ MÔ TẢ CHO CHUỖI TẦN SUẤT CAO -------------------------- #
# ------------------------------------------------------------------------------ #
    # Thiết lập các biến cần kiểm định
      #cols_to_test_rate <- c("VNIndex", "XAUUSD", "GC_F", "SJC", "BTC", "ETH", "BNB", "BRENT", "WTI")
      cols_to_test_rate <- c("VNIndex", "XAUUSD", "GC_F", "SJC")
      
    # Gọi hàm xử lý và xuất file
      process_and_export(rate, cols_to_test_rate, file.path("DS Results", "Rates DS.csv"))

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 8. KẾT QUẢ THỐNG KÊ MÔ TẢ CHO CHUỖI TẦN SUẤT THẤP ------------------------- #
# ------------------------------------------------------------------------------ #
    # Thiết lập các biến cần kiểm định
      #cols_to_test_shocks <- c("GPR", "GPRT", "GPRA", "GPR_VIE", "GEPU_current", "GEPU_ppp", "WUI", "WPUI", "WSI", "WTUI")
      cols_to_test_shocks <- c("GPR", "GPRT", "GPRA")
      
    # Gọi hàm xử lý và xuất file
      process_and_export(shocks, cols_to_test_shocks, file.path("DS Results", "Shocks DS.csv"))



      

      