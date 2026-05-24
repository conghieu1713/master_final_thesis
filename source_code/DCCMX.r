# ============================================================================== #
# ------------                                                       ----------- #
# ------------              ƯỚC LƯỢNG MÔ HÌNH DCC-MIDAS-X            ----------- #
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
library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------ #
# -- 2. THIẾT LẬP CÁC THAM SỐ CHO MÔ HÌNH -------------------------------------- #
# ------------------------------------------------------------------------------ #
# Lựa chọn kịch bản làm mềm return trong GARCH-MIDAS
do_arma_approach <- 1 # 1. ARMA(1,1), 2. auto.arima

# Cấu hình các cặp tài sản và tham số K_lags tương ứng
asset_pairs_config <- list(
  list(pair = c("VNIndex", "XAUUSD"), K_lags = 36),
  list(pair = c("VNIndex", "GC_F"),   K_lags = 36),
  list(pair = c("VNIndex", "SJC"),    K_lags = 24)
)

# Cấu hình các biến ngoại sinh X (không bao gồm K_lags)
exogenous_vars_config <- list(
  list(var = "GPR"),
  list(var = "GPRT"),
  list(var = "GPRA")
)

# Giá trị khởi tạo cho thuật toán tối ưu hóa: a, b, m, theta, w2
init_params <- c(0.05, 0.80, 0.1, -0.05, 2.0)

# Ràng buộc dưới (Lower bounds): a > 0, b > 0, m, theta (tự do), w2 > 1
lower_bounds <- c(0.001, 0.001, -10, -10, 1.001)

# Ràng buộc trên (Upper bounds): a < 1, b < 1, m, theta, w2 < 50
upper_bounds <- c(0.999, 0.999,  10,  10, 50.0)

# Thư mục làm việc và đường dẫn tệp
dir <- "D:/Post-graduate programme/Final Thesis/Model"
setwd(dir)

exogenous_file_path <- "official.xlsx"
exogenous_sheet_name <- "shocks_lg"

if (do_arma_approach == 1) {
  residuals_file_path <- "arma11/standardized_residuals.csv"
} else if (do_arma_approach == 2) {
  residuals_file_path <- "autoarima/standardized_residuals.csv"
}

# ------------------------------------------------------------------------------ #
# -- 3. ĐỊNH NGHĨA CÁC HÀM CHỨC NĂNG ------------------------------------------- #
# ------------------------------------------------------------------------------ #

# -- 3.1. HÀM TẢI VÀ KHỚP DỮ LIỆU (PHIÊN BẢN ĐA TẦN SUẤT) ---------------- #
load_and_align_data <- function(z_file_path, asset1, asset2, x_file_path, x_sheet, x_var) {
  cat(">> Đang đọc phần dư chuẩn hóa từ:", z_file_path, "\n")
  z_data <- read.csv(z_file_path, header = TRUE)
  z_xts <- xts(z_data[, c(asset1, asset2)], order.by = as.Date(z_data$Date))
  z_xts <- na.omit(z_xts)
  colnames(z_xts) <- c(asset1, asset2)
  
  cat(sprintf(">> Đang đọc biến ngoại sinh '%s' từ tệp '%s', sheet '%s'\n", x_var, basename(x_file_path), x_sheet))
  x_data_monthly <- read_excel(x_file_path, sheet = x_sheet)
  x_xts_monthly <- xts(x_data_monthly[[x_var]], order.by = as.Date(x_data_monthly$Date))
  colnames(x_xts_monthly) <- x_var
  
  # Trích xuất định dạng "Năm-Tháng" để tìm điểm giao thoa
  z_ym <- format(index(z_xts), "%Y-%m")
  x_ym <- format(index(x_xts_monthly), "%Y-%m")
  
  # Tìm các tháng chung để đồng bộ
  common_ym <- intersect(unique(z_ym), unique(x_ym))
  
  # Lọc dữ liệu
  z_xts_aligned <- z_xts[z_ym %in% common_ym]
  x_xts_aligned <- x_xts_monthly[x_ym %in% common_ym]
  
  z_matrix <- as.matrix(z_xts_aligned)
  x_monthly <- as.numeric(x_xts_aligned)
  
  # Tạo vector ánh xạ: Ngày thứ t thuộc tháng thứ mấy trong mảng x_monthly
  z_ym_aligned <- format(index(z_xts_aligned), "%Y-%m")
  x_ym_aligned <- format(index(x_xts_aligned), "%Y-%m")
  month_map <- match(z_ym_aligned, x_ym_aligned)
  
  cat(sprintf(">> Đã đồng bộ: %d ngày giao dịch, tương ứng với %d tháng.\n", nrow(z_matrix), length(x_monthly)))
  return(list(z_matrix = z_matrix, x_monthly = x_monthly, month_map = month_map))
}

# -- 3.2. HÀM TÍNH TRỌNG SỐ BETA MIDAS ----------------------------------------- #
beta_weight <- function(K, w2) {
  k <- 1:K
  w <- (1 - k/K)^(w2 - 1)
  return(w / sum(w))
}

# -- 3.3. HÀM ĐỐI LOG-LIKELIHOOD CHO DCC-MIDAS-X -------------------------------- #
dcc_midas_x_likelihood <- function(theta_params, Z, X_monthly, month_map, K) {
  a     <- theta_params[1]
  b     <- theta_params[2]
  m     <- theta_params[3]
  theta <- theta_params[4]
  w2    <- theta_params[5]
  
  T_obs <- nrow(Z)
  N_months <- length(X_monthly)
  LL <- 0
  
  # 1. Tìm ngày t đầu tiên có đủ dữ liệu (thuộc tháng thứ K+1)
  start_month <- K + 1 # Bắt đầu từ tháng K+1 để có đủ K tháng quá khứ
  start_t <- which(month_map >= start_month)[1]
  if (is.na(start_t)) return(1e10) # Trả về phạt nếu không đủ dữ liệu

  # 2. Khởi tạo ma trận Q bằng ma trận tương quan mẫu
  Q <- cor(Z)
  # Q <- cor(Z[1:(start_t - 1), ])
  weights <- beta_weight(K, w2)
  
  # 3. Tính toán tương quan dài hạn ở tần suất THÁNG (tau)
  m_tau <- rep(0, N_months)
  
  for (tau in start_month:N_months) {
    X_lags <- X_monthly[seq(tau - 1, tau - K, by = -1)]
    X_midas <- sum(weights * X_lags)
    m_tau[tau] <- m + theta * X_midas
  }
  
  # 4. Vòng lặp chính ở tần suất NGÀY (t)
  for (t in start_t:T_obs) {
    tau_t <- month_map[t] # Lấy index tháng hiện tại của ngày t
    rho_long_t <- tanh(m_tau[tau_t]) # Rút giá trị tương quan dài hạn tương ứng
    
    R_bar <- matrix(c(1, rho_long_t, rho_long_t, 1), nrow = 2, ncol = 2)
    z_lag <- matrix(Z[t-1, ], ncol = 1)
    Q <- (1 - a - b) * R_bar + a * (z_lag %*% t(z_lag)) + b * Q
    
    Q_diag_inv <- diag(1 / sqrt(diag(Q))) 
    R_t <- Q_diag_inv %*% Q %*% Q_diag_inv
    
    z_t <- matrix(Z[t, ], ncol = 1)
    
    det_R <- det(R_t)
    if (is.na(det_R) || det_R <= 1e-10) return(1e10)
    
    term1 <- log(det_R)
    term2 <- as.numeric(t(z_t) %*% solve(R_t) %*% z_t)
    term3 <- sum(z_t^2) # Tương đương as.numeric(t(z_t) %*% z_t)
    
    LL <- LL - 0.5 * (term1 + term2 - term3)
  }
  
  return(as.numeric(-LL))
}

# -- 3.4. HÀM RÀNG BUỘC BẤT PHƯƠNG TRÌNH --------------------------------------- #
ineq_fun <- function(theta_params, Z, X_monthly, month_map, K) {
  return(theta_params[1] + theta_params[2])
}

# -- 3.5. HÀM ƯỚC LƯỢNG MÔ HÌNH ------------------------------------------------ #
estimate_dcc_midas_x <- function(z_matrix, x_monthly, month_map, init_params, lower_bounds, upper_bounds, K_lags) {
  cat("\n>> Bắt đầu quá trình tối ưu hóa DCC-MIDAS-X...\n")
  optim_results <- solnp(
    pars = init_params, fun = dcc_midas_x_likelihood, ineqfun = ineq_fun, 
    ineqLB = 0.001, ineqUB = 0.999, LB = lower_bounds, UB = upper_bounds, 
    Z = z_matrix, X_monthly = x_monthly, month_map = month_map, K = K_lags
  )
  return(optim_results)
}

# -- 3.6. HÀM BÁO CÁO KẾT QUẢ -------------------------------------------------- #
report_dcc_x_results <- function(optim_results, z_matrix, x_monthly, month_map, K_lags, asset1, asset2, x_var) {
  
  get_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return(" ")
  }
  
  estimates <- optim_results$pars
  n_pars <- length(estimates)
  n_obs <- nrow(z_matrix)

  cat(">> Đang tính toán ma trận Hessian (numDeriv) để tìm sai số chuẩn...\n")
  # CHÚ Ý: Truyền đúng biến x_monthly và month_map
  hessian_matrix <- numDeriv::hessian(
    func = dcc_midas_x_likelihood, x = estimates,
    Z = z_matrix, X_monthly = x_monthly, month_map = month_map, K = K_lags
  )
  
  cov_matrix <- tryCatch(solve(hessian_matrix), error = function(e) {
    warning("Ma trận Hessian không thể nghịch đảo. Sai số chuẩn có thể không chính xác. Lỗi: ", e$message)
    return(matrix(NA, nrow = n_pars, ncol = n_pars))
  })
  
  std_errors <- suppressWarnings(sqrt(diag(cov_matrix)))
  t_values <- estimates / std_errors
  p_values <- 2 * (1 - pnorm(abs(t_values)))
  sig_stars <- sapply(p_values, get_stars)
  
  results_table <- data.frame(
    "Estimate"   = round(estimates, 6),
    "Std. Error" = round(std_errors, 6),
    "t value"    = round(t_values, 4),
    "Pr(>|t|)"   = signif(p_values, 4),
    "Sig."       = sig_stars,
    check.names  = FALSE
  )
  rownames(results_table) <- c("alpha (a)", "beta (b)", "m", "theta", "omega (w2)")
  
  cat("\n========================================================================\n")
  cat("      KẾT QUẢ ƯỚC LƯỢNG: MÔ HÌNH DCC-MIDAS-X (CUSTOM OPTIMIZATION)     \n")
  cat("========================================================================\n")
  cat(sprintf(" Cặp tài sản: %s - %s\n", asset1, asset2))
  cat(sprintf(" Biến ngoại sinh (X): %s\n", x_var))
  cat("------------------------------------------------------------------------\n")
  print(results_table, row.names = TRUE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  loglik <- -tail(optim_results$values, 1)
  aic <- 2 * n_pars - 2 * loglik
  bic <- n_pars * log(n_obs) - 2 * loglik
  
  cat(sprintf("Log-Likelihood: %0.4f\n", loglik))
  cat(sprintf("AIC           : %0.4f\n", aic))
  cat(sprintf("BIC           : %0.4f\n", bic))
  cat(sprintf("Số quan sát (T): %d\n", n_obs))
  cat(sprintf("Tham số MIDAS: K_lags = %d tháng\n", K_lags))
  cat("========================================================================\n")
  
  return(invisible(results_table))
}

# ------------------------------------------------------------------------------ #
# -- 4. THỰC THI MÔ HÌNH CHO NHIỀU CẶP TÀI SẢN VÀ BIẾN NGOẠI SINH -------------- #
# ------------------------------------------------------------------------------ #

# Vòng lặp qua từng biến ngoại sinh
for (x_config in exogenous_vars_config) {
  x_var  <- x_config$var

  # Vòng lặp qua từng cặp tài sản
  for (asset_config in asset_pairs_config) {
    asset1 <- asset_config$pair[1]
    asset2 <- asset_config$pair[2]
    K_lags <- asset_config$K_lags # Lấy K_lags từ cấu hình của cặp tài sản
    
    cat("\n\n========================================================================")
    cat(sprintf("\n BẮT ĐẦU XỬ LÝ: Cặp [%s-%s] với biến X [%s]\n", asset1, asset2, x_var))
    cat(sprintf(">> Tham số MIDAS được sử dụng: K_lags = %d\n", K_lags))
    cat("========================================================================\n")

    # 1. Tải và chuẩn bị dữ liệu
    # Sử dụng try-catch để bỏ qua nếu có lỗi (ví dụ: tệp không tồn tại)
    data_list <- tryCatch({
      load_and_align_data(
        z_file_path = residuals_file_path,
        asset1 = asset1,
        asset2 = asset2,
        x_file_path = exogenous_file_path,
        x_sheet = exogenous_sheet_name,
        x_var = x_var
      )
    }, error = function(e) {
      cat("LỖI khi tải dữ liệu:", e$message, "\n")
      return(NULL)
    })

    if (is.null(data_list)) {
      cat(">> Bỏ qua cặp này do lỗi tải dữ liệu.\n")
      next # Chuyển sang vòng lặp tiếp theo
    }

    # 2. Ước lượng mô hình
    optim_results <- estimate_dcc_midas_x(
      z_matrix = data_list$z_matrix,
      x_monthly = data_list$x_monthly,
      month_map = data_list$month_map,
      init_params = init_params, lower_bounds = lower_bounds,
      upper_bounds = upper_bounds, K_lags = K_lags
    )

    # 3. Báo cáo kết quả
    report_dcc_x_results(
      optim_results = optim_results,
      z_matrix = data_list$z_matrix,
      x_monthly = data_list$x_monthly,     # truyền biến mới
      month_map = data_list$month_map,     # truyền biến mới
      K_lags = K_lags,
      asset1 = asset1,
      asset2 = asset2,
      x_var = x_var
    )
    
    # (Tùy chọn) Bạn có thể thêm hàm vẽ biểu đồ ở đây nếu muốn
    # plot_dcc_x_results(...)
    
    cat(sprintf("\n HOÀN TẤT XỬ LÝ: Cặp [%s-%s] với biến X [%s]\n", asset1, asset2, x_var))
  }
}
