# ============================================================================== #
# ------------                                                       ----------- #
# ------------                  MA TRẬN TƯƠNG QUAN                   ----------- #
# ------------                                                       ----------- #
# ============================================================================== #

# Tác giả: Công Hiếu, NGUYỄN (524102110660)
# Học viên Cao học, Đại học Kinh tế TP.Hồ Chí Minh


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

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 3. HÀM TÍNH TOÁN MA TRẬN TƯƠNG QUAN --------------------------------------- #
# ------------------------------------------------------------------------------ #
  cor_with_stars <- function(data) {
    cols <- colnames(data)
    n <- length(cols)
    
    # Tạo một ma trận rỗng để lưu kết quả
    res_matrix <- matrix("", n, n)
    rownames(res_matrix) <- cols
    colnames(res_matrix) <- cols
    
    for (i in 1:n) {
      for (j in 1:i) { # Chỉ lặp qua tam giác dưới và đường chéo
        if (i == j) {
          res_matrix[i, j] <- "1.000"
        } else {
          # Tạo một data frame tạm thời cho cặp cột i, j và loại bỏ giá trị NA
          pair_data <- na.omit(data.frame(x = data[[i]], y = data[[j]]))
          
          if (nrow(pair_data) > 2) {
            # Tính toán hệ số tương quan Pearson và p-value
            test <- cor.test(pair_data$x, pair_data$y, method = "pearson")
            r <- test$estimate
            p <- test$p.value
            
            # Thêm dấu sao ý nghĩa thống kê
            stars <- ifelse(p < 0.01, "***", 
                            ifelse(p < 0.05, "**", 
                                   ifelse(p < 0.1, "*", "")))
            
            # Định dạng kết quả
            res_matrix[i, j] <- paste0(sprintf("%.3f", r), stars)
          } else {
            # Trả về NA nếu dữ liệu không đủ
            res_matrix[i, j] <- "NA"
          }
        }
      }
    }
    return(as.data.frame(res_matrix))
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 4. HÀM TỔNG HỢP VÀ XUẤT DỮ LIỆU ------------------------------------------- #
# ------------------------------------------------------------------------------ #
  process_and_export_correlation <- function(data, cols, out_filename) {
    # Chọn các cột cần thiết và chuyển đổi sang dạng số
    data_subset <- data[, cols]
    data_subset <- as.data.frame(sapply(data_subset, as.numeric))
    
    # Tính toán ma trận tương quan
    correlation_matrix <- cor_with_stars(data_subset)
    
    # In kết quả ra màn hình console với tiêu đề rõ ràng
    cat("\n# ------------------------------------------------------------------------------ #\n")
    cat("# -- KẾT QUẢ TỪ FILE:", toupper(basename(out_filename)), "\n")
    cat("# ------------------------------------------------------------------------------ #\n")
    print(correlation_matrix)
    cat("\n")
    
    # Đảm bảo thư mục lưu trữ tồn tại trước khi xuất file
    dir.create(dirname(out_filename), showWarnings = FALSE, recursive = TRUE)
    write.csv(correlation_matrix, out_filename)
  }

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 5. NHẬP DỮ LIỆU ----------------------------------------------------------- #
# ------------------------------------------------------------------------------ #
  # Đường dẫn thư mục làm việc (Sử dụng / để tránh lỗi đường dẫn)
    dir <- "D:/Post-graduate programme/Final Thesis/Model"   
    setwd(dir)
    
  # Đọc dữ liệu tần suất cao
    rate <- read_excel("official.xlsx", sheet = "rate", guess_max = 100000)
    
  # Đọc dữ liệu tần suất thấp
    shocks <- read_excel("official.xlsx", sheet = "shocks", guess_max = 100000)

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 6. KẾT QUẢ MA TRẬN TƯƠNG QUAN CHO CHUỖI TẦN SUẤT CAO ---------------------- #
# ------------------------------------------------------------------------------ #
    # Thiết lập các biến cần tính toán
      cols_to_test_rate <- c("VNIndex", "XAUUSD", "GC_F", "SJC")
      
    # Gọi hàm xử lý và xuất file
      process_and_export_correlation(rate, cols_to_test_rate, file.path("DS Results", "CorMat HF.csv"))

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# -- 7. KẾT QUẢ MA TRẬN TƯƠNG QUAN CHO CHUỖI TẦN SUẤT THẤP --------------------- #
# ------------------------------------------------------------------------------ #
    # Thiết lập các biến cần tính toán
      cols_to_test_shocks <- c("GPR", "GPRT", "GPRA")
      
    # Gọi hàm xử lý và xuất file
      process_and_export_correlation(shocks, cols_to_test_shocks, file.path("DS Results", "CorMat LF.csv"))
