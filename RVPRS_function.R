# Load packages

if (!require(seqminer)) install.packages(seqminer)
if (!require(data.table)) install.packages(data.table)
if (!require(tibble)) install.packages(tibble)
if (!require(dplyr)) install.packages(dplyr)
if (!require(Matrix)) install.packages(Matrix)
if (!require(sparsesvd)) install.packages(sparsesvd)
if (!require(stringr)) install.packages(stringr)

library(seqminer, quietly = TRUE)
library(data.table, quietly = TRUE)
library(tibble, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(sparsesvd, quietly = TRUE)
library(stringr, quietly = TRUE)

# Load PLINK file
load_plink <- function(plink_prefix) {
    plink_obj <- seqminer::openPlink(plink_prefix)
    marker_index <- seq(nrow(plink_obj$bim))
    sample_index <- seq(nrow(plink_obj$fam))
    plink_matrix <- seqminer::readPlinkToMatrixByIndex(plink_prefix, sample_index, marker_index)
    colnames(plink_matrix) <- plink_obj$bim$V2

    return(plink_matrix)
}

# Read Group file and split variants by functional annotations
read_groupfile <- function(groupfile_name, gene_name) {
    groupfile <- file(groupfile_name, "r")
    line <- 0
    var <- NULL
    anno <- NULL

    while (TRUE) {
        line <- line + 1
        marker_group_line <- readLines(groupfile, n = 1)

        if (length(marker_group_line) == 0) {
            break
        }

        marker_group_line_list <- strsplit(marker_group_line, split = c(" +", "\t"))[[1]]

        if (marker_group_line_list[1] == gene_name) {
            if (marker_group_line_list[2] == "var") {
                var <- marker_group_line_list
            } else {
                anno <- marker_group_line_list
            }
        }
    }

    lof_idx <- which(anno == "lof")
    mis_idx <- which(anno == "missense")
    syn_idx <- which(anno == "synonymous")

    lof_var <- var[lof_idx]
    mis_var <- var[mis_idx]
    syn_var <- var[syn_idx]

    out <- list(lof_var, mis_var, syn_var)
    close(groupfile)
    return(out)
}

get_range <- function(v) {
    # input are vectors of variants by functional annotation
    start_pos <- Inf
    end_pos <- 0
    for (i in 1:3) {
        if (length(v[[i]]) > 0) {
            pos <- as.numeric(str_split_fixed(v[[i]], ":", 4)[,2]) # Modified in 0.3
            sub_start_pos <- min(pos)
            sub_end_pos <- max(pos)
            print(sub_start_pos)
            print(sub_end_pos)
            if (start_pos > sub_start_pos) {
                start_pos <- sub_start_pos
            }
            if (end_pos < sub_end_pos) {
                end_pos <- sub_end_pos
            }
        }
    }
    return (list(start_pos, end_pos))
}

# Convert genotypes for missing / alt-first PLINK files
convert_missing <- function(plink_matrix, ref_first=FALSE) {
    if (ref_first == FALSE) {   # alt-first
        # Convert missing genotypes to reference-homozygotes
        plink_matrix[which(plink_matrix == -9)] <- 2
        # Flip genotypes
        plink_matrix <- 2 - plink_matrix
    } else {
        plink_matrix[which(plink_matrix == -9)] <- 0
    }

    return(plink_matrix)
}


# function for counting MAF of input genotype matrix (used in split_plink_matrix)
count_maf <- function(geno_mat) {
    mac_vec <- rep(0, ncol(geno_mat))
    for (i in seq_len(ncol(geno_mat))) {
        # n0 <- length(which(geno_mat[, i]) == 0)
        n1 <- length(which(geno_mat[, i] == 1))
        n2 <- length(which(geno_mat[, i] == 2))

        mac <- n1 + 2 * n2
        mac_vec[i] <- mac
    }
    maf_vec <- mac_vec / (2 * nrow(geno_mat))
    common_list <- which(maf_vec > 0.01)

    out <- list(mac_vec, common_list)
    return(out)
}


# Split genotype matrix by function annotation (lof / mis / syn)
split_plink_matrix <- function(plink_matrix, lof_var, mis_var, syn_var, mac_threshold = 10) {
    lof_mat <- plink_matrix[, lof_var, drop = F]
    mis_mat <- plink_matrix[, mis_var, drop = F]
    syn_mat <- plink_matrix[, syn_var, drop = F]

    # count MAF for each functional annotation
    mac_lof <- count_maf(lof_mat)[[1]]
    common_list_lof <- count_maf(lof_mat)[[2]]

    mac_mis <- count_maf(mis_mat)[[1]]
    common_list_mis <- count_maf(mis_mat)[[2]]

    mac_syn <- count_maf(syn_mat)[[1]]
    common_list_syn <- count_maf(syn_mat)[[2]]

    # collapsing columns with MAC < threshold
    collapse_lof <- which(mac_lof < mac_threshold)
    collapse_mis <- which(mac_mis < mac_threshold)
    collapse_syn <- which(mac_syn < mac_threshold)

    if(length(collapse_lof) > 0) {
        rowsum_lof <- rowSums(lof_mat[, collapse_lof, drop = F])
        rowsum_lof[which(rowsum_lof > 1)] <- 1
        lof_mat_collapsed <- cbind(lof_mat[, -c(collapse_lof, common_list_lof), drop = F], rowsum_lof)
    } else if ((length(common_list_lof) > 0)) {
        lof_mat_collapsed <- lof_mat[, -c(common_list_lof), drop = F]
    } else {
        lof_mat_collapsed <- lof_mat
    }

    if(length(collapse_mis) > 0) {
        rowsum_mis <- rowSums(mis_mat[, collapse_mis, drop = F])
        rowsum_mis[which(rowsum_mis > 1)] <- 1
        mis_mat_collapsed <- cbind(mis_mat[, -c(collapse_mis, common_list_mis), drop = F], rowsum_mis)
    } else if ((length(common_list_mis) > 0)) {
        mis_mat_collapsed <- mis_mat[, -c(common_list_mis), drop = F]
    } else {
        mis_mat_collapsed <- mis_mat
    }

    if(length(collapse_syn) > 0) {
        rowsum_syn <- rowSums(syn_mat[, collapse_syn, drop = F])
        rowsum_syn[which(rowsum_syn > 1)] <- 1
        syn_mat_collapsed <- cbind(syn_mat[, -c(collapse_syn, common_list_syn), drop = F], rowsum_syn)
    } else if ((length(common_list_syn) > 0)) {
        syn_mat_collapsed <- syn_mat[, -c(common_list_syn), drop = F]
    } else {
        syn_mat_collapsed <- syn_mat
    }

    out <- list(lof_mat_collapsed, mis_mat_collapsed, syn_mat_collapsed)
    return(out)
}

# Read phenotype file
read_pheno <- function(pheno_file, pheno_code, iid_col = "f.eid") {
    pheno <- data.table::fread(pheno_file, quote = "")
    out <- subset(pheno, select = c(iid_col, pheno_code))
    return(out)
}


calc_log_lik <- function(delta, S, UtY, Y_UUtY) {
    k <- length(S)
    n <- nrow(Y_UUtY)

    log_lik1 <- 0
    for (i in 1:k) {
        log_lik1 <- log_lik1 + (UtY[i, ])^2 / (S[i] + delta)
    }

    log_lik2 <- 1 / delta * sum((Y_UUtY)^2)

    out <- -0.5 * (n * log(2 * pi) + sum(log(S + delta)) + (n - k) * log(delta)
                   + n + n * log(1 / n * (log_lik1 + log_lik2)))

    return(as.numeric(out))
}

calc_post_beta <- function(K, G, delta, S, UtY, U) {
    K_sparse <- as(K, "dgCMatrix")
    if (length(S) == 1) {
        S <- as.matrix(S)
    }

    out <- K_sparse %*% t(G) %*% U %*% diag(1 / (S + delta)) %*% (UtY)
    return(out)
}

# Run FaST-LMM to obtain posterior beta
fast_lmm <- function(G, Y) {
    # if sum(G) == 0, let effect size = 0
    print("Estimating beta using FaST-LMM")
    print(dim(G))
    G[is.na(G)] <- 0
    print(sum(G))
    if (sum(G, na.rm = T) == 0) {
        return (list(as.matrix(0), 0, 1e6))
    }
    Y <- as.matrix(Y)
    K <- diag(1, nrow = ncol(G))
    L <- chol(K)
    W <- G %*% L
    W_sparse <- as(W, "dgCMatrix")
    svd_mat <- sparsesvd::sparsesvd(W_sparse)

    U <- svd_mat$u
    S <- (svd_mat$d)^2

    UtY <- t(U) %*% Y

    # Y_UUtY = Y - UUtY
    Y_UUtY <- Y - U %*% UtY

    opt <- optim(par = 1, fn = calc_log_lik, S = S, UtY = UtY, Y_UUtY = Y_UUtY,
                 method = c("Brent"), lower = 0, upper = 1e6, control = list(fnscale = -1))
    opt_delta <- opt$par

    tr_GtG <- sum(diag(t(G) %*% G %*% K))

    post_beta <- calc_post_beta(K, G, opt_delta, S, UtY, U)

    return (list(post_beta, tr_GtG, opt_delta))
}


calc_gene_effect_size <- function(G, lof_ncol, post_beta, Y) {
    beta_lof <- post_beta[1:lof_ncol]
    beta_lof <- abs(beta_lof)
    lof_prs <- G[, 1:lof_ncol, drop = F] %*% beta_lof
    lof_prs_norm <- (lof_prs - mean(lof_prs)) / sd(lof_prs)
    m1 <- lm(Y ~ lof_prs_norm[, 1])

    return(m1$coefficients[2])
}

mom_estimator_marginal <- function(G, y) {
    n <- length(y)
    G <- as(G, "dgCMatrix")
    G[is.na(G)] <- 0
    # tr(G Sigma G^T G Sigma G^T) = sum((G Sigma G^T )^2) = sum(G^T G Sigma)^2
    Sigma <- diag(1, nrow = ncol(G))

    system.time({
        t1 <- sum((crossprod(G) %*% Sigma)^2)
    })
    t2 <- sum(diag(crossprod(G) %*% Sigma))
    A <- matrix(c(t1, t2, t2, n), ncol = 2)
    c1 <- as.numeric(t(y) %*% G %*% Sigma %*% t(G) %*% y)
    c2 <- sum(y^2)
    b <- matrix(c(c1, c2), ncol = 1)
    var_comp <- solve(A) %*% b
    print(var_comp)
    # h2_mom_marginal <- var_comp[1, 1] / sum(var_comp)
    return (var_comp)
}

mom_estimator_joint <- function(G1, G2, G3, y) {
    n <- length(y)

    G1 <- as(G1, "dgCMatrix")
    G1[is.na(G1)] <- 0
    G2 <- as(G2, "dgCMatrix")
    G2[is.na(G2)] <- 0
    G3 <- as(G3, "dgCMatrix")
    G3[is.na(G3)] <- 0
    Sigma1 <- diag(1, nrow = ncol(G1)) # L1 %*% t(L1)
    Sigma2 <- diag(1, nrow = ncol(G2)) # L2 %*% t(L2)
    Sigma3 <- diag(1, nrow = ncol(G3)) # L3 %*% t(L3)

    L1 <- chol(Sigma1)
    L2 <- chol(Sigma2)
    L3 <- chol(Sigma3)

    t11 <- sum((crossprod(G1) %*% Sigma1)^2)
    t22 <- sum((crossprod(G2) %*% Sigma2)^2)
    t33 <- sum((crossprod(G3) %*% Sigma3)^2)

    t12 <- sum((t(G1 %*% L1) %*% (G2 %*% L2))^2)
    t13 <- sum((t(G1 %*% L1) %*% (G3 %*% L3))^2)
    t23 <- sum((t(G2 %*% L2) %*% (G3 %*% L3))^2)

    t14 <- sum(diag(t(G1) %*% G1 %*% Sigma1))
    t24 <- sum(diag(t(G2) %*% G2 %*% Sigma2))
    t34 <- sum(diag(t(G3) %*% G3 %*% Sigma3))

    A <- matrix(c(t11, t12, t13, t14, t12, t22, t23, t24, t13, t23, t33, t34, t14, t24, t34, n), ncol = 4)

    c1 <- as.numeric(t(y) %*% G1 %*% Sigma1 %*% t(G1) %*% y)
    c2 <- as.numeric(t(y) %*% G2 %*% Sigma2 %*% t(G2) %*% y)
    c3 <- as.numeric(t(y) %*% G3 %*% Sigma3 %*% t(G3) %*% y)
    c4 <- sum(y^2)
    b <- matrix(c(c1, c2, c3, c4), ncol = 1)

    var_comp <- solve(A) %*% b
    print(var_comp)
    #  h2_mom_joint <- c(var_comp[1, 1], var_comp[2, 1], var_comp[3, 1]) / sum(var_comp)
    return (var_comp)
}

# calculate_single_blup <- function(G, delta, Sigma, y) {
#     # Henderson Mixed Model Equation: beta = (G^T G + delta * Sigma^-1)^-1 * G^T y
#     G <- as(G, "dgCMatrix")
#     beta <- solve(t(G) %*% G + delta * solve(Sigma)) %*% t(G) %*% y
#     return (beta)
# }

calculate_joint_blup <- function(G1, G2, G3, tau1, tau2, tau3, Sigma1, Sigma2, Sigma3, psi, y) {
    G1 <- as(G1, "dgCMatrix")
    G2 <- as(G2, "dgCMatrix")
    G3 <- as(G3, "dgCMatrix")
    G <- cbind(G1, G2, G3)
    G[is.na(G)] <- 0

    Sigma <- bdiag(tau1 * Sigma1, tau2 * Sigma2, tau3 * Sigma3)
    beta <- solve(t(G) %*% G / psi + solve(Sigma)) %*% t(G) %*% y / psi
    return (beta)
}
