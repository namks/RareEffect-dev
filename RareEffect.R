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

read_matrix_by_one_marker <- function(objGeno, var_list, sampleID) {
    n_samples <- length(sampleID)
    mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
    if (length(var_list) == 0) {
        print("No variants in the list")
        return(mat)
    }
    for (i in 1:length(var_list)) {
        t_GVec <- rep(0, n_samples)
        idx <- which(var_list[i] == objGeno$markerInfo$ID)
        SAIGE::Unified_getOneMarker(t_genoType = objGeno$genoType,
            t_gIndex_prev = objGeno$markerInfo$genoIndex_prev[idx],
            t_gIndex = objGeno$markerInfo$genoIndex[idx],
            t_ref = "2",
            t_alt = "1",
            t_marker = objGeno$markerInfo$ID[idx],
            t_pd = objGeno$markerInfo$POS[idx],
            t_chr = toString(objGeno$markerInfo$CHROM[idx]),
            t_altFreq = 0,
            t_altCounts = 0,
            t_missingRate = 0,
            t_imputeInfo = 0,
            t_isOutputIndexForMissing = TRUE,
            t_indexForMissing = 0,
            t_isOnlyOutputNonZero = FALSE,
            t_indexForNonZero = 0,
            t_GVec = t_GVec,
            t_isImputation = FALSE
        )

	t_GVec[t_GVec < 0] <- 0     # Convert missing to zero
        t_GVec_sp <- as(t_GVec, "sparseVector")
        t_GVec_sp_mat <- as(t_GVec_sp, "Matrix")

        mat <- cbind(mat, t_GVec_sp_mat)
    }

    colnames(mat) <- var_list
    rownames(mat) <- sampleID

    return(mat)
}

collapse_matrix <- function(objGeno, var_list, sampleID, modglmm, macThreshold = 10) {
    n_samples <- length(sampleID)
    mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
    if (length(var_list) == 0) {
        print("No variants in the list")
        return(mat)
    }
    mat <- read_matrix_by_one_marker(objGeno, var_list, sampleID)
    mat <- mat[which(rownames(mat) %in% modglmm$sampleID), ]
    MAF <- colSums(mat) / (2 * nrow(mat))
    idx <- which(((MAF < 0.01) | (MAF > 0.99)) & ((MAF < 1) & (MAF > 0)))
    mat <- mat[, idx]
    mat_rare <- mat[, which(colSums(mat) >= macThreshold)]
    mat_UR <- mat[, which(colSums(mat) < macThreshold)]
    UR_rowsum <- rowSums(mat_UR)
    UR_rowsum[which(UR_rowsum > 1)] <- 1
    mat_UR_collapsed <- cbind(mat_rare, UR_rowsum)
    return(mat_UR_collapsed)
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

    GtG <- t(G) %*% G
    tr_GtG <- sum(diag(GtG %*% K))

    post_beta <- calc_post_beta(K, G, opt_delta, S, UtY, U)

    return (list(post_beta, tr_GtG, opt_delta, GtG))
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

run_RareEffect <- function(rdaFile, chrom, geneName, groupFile, traitType, bedFile, bimFile, famFile, macThreshold, collapseLoF, collapsemis, collapsesyn, outputPrefix) {
    # Load SAIGE step 1 results
    load(rdaFile)

    if (traitType == "binary") {
        # modglmm$residuals <- modglmm$Y - modglmm$linear.predictors
        v <- modglmm$fitted.values * (1 - modglmm$fitted.values)    # v_i = mu_i * (1 - mu_i)
        modglmm$residuals <- sqrt(v) * (1 / (modglmm$fitted.values * (1 - modglmm$fitted.values)) * (modglmm$y - modglmm$fitted.values))
    }
    sigma_sq <- var(modglmm$residuals)
    # n_samples <- length(modglmm$residuals)

    # Set PLINK object
    bim <- fread(bimFile)
    fam <- fread(famFile)
    sampleID <- as.character(fam$V2)
    n_samples <- length(sampleID)

    objGeno <- SAIGE::setGenoInput(bgenFile = "",
            bgenFileIndex = "",
            vcfFile = "",
            vcfFileIndex = "",
            vcfField = "",
            savFile = "",
            savFileIndex = "",
            sampleFile = "",
            bedFile=bedFile,
            bimFile=bimFile,
            famFile=famFile,
            idstoIncludeFile = "",
            rangestoIncludeFile = "",
            chrom = "",
            AlleleOrder = "alt-first",
            sampleInModel = sampleID
        )

    print("Analysis started")

    var_by_func_anno <- read_groupfile(groupFile, geneName)

    # LoF
    if (length(var_by_func_anno[[1]]) == 0) {
        print("LoF variant does not exist.")
    }
    # missense
    if (length(var_by_func_anno[[2]]) == 0) {
        print("missense variant does not exist.")
    }
    # synonymous
    if (length(var_by_func_anno[[3]]) == 0) {
        print("synonymous variant does not exist.")
    }

    # Remove variants not in plink file
    for (i in 1:3) {
        var_by_func_anno[[i]] <- var_by_func_anno[[i]][which(var_by_func_anno[[i]] %in% objGeno$markerInfo$ID)]
    }
    print(str(var_by_func_anno))

    # Read genotype matrix
    if (collapseLoF) {
        lof_mat_collapsed <- collapse_matrix(objGeno, var_by_func_anno[[1]], sampleID, modglmm, macThreshold = 0)
        lof_mat_collapsed <- Matrix::Matrix(rowSums(lof_mat_collapsed), ncol = 1, sparse = TRUE, dimnames = list(rownames(lof_mat_collapsed), NULL))
    } else {
        lof_mat_collapsed <- collapse_matrix(objGeno, var_by_func_anno[[1]], sampleID, modglmm, macThreshold)
    }
    mis_mat_collapsed <- collapse_matrix(objGeno, var_by_func_anno[[2]], sampleID, modglmm, macThreshold)
    syn_mat_collapsed <- collapse_matrix(objGeno, var_by_func_anno[[3]], sampleID, modglmm, macThreshold)

    # Set column name
    if (ncol(lof_mat_collapsed) > 0) {
        colnames(lof_mat_collapsed)[ncol(lof_mat_collapsed)] <- "lof_UR"
    }
    if (ncol(mis_mat_collapsed) > 0) {
        colnames(mis_mat_collapsed)[ncol(mis_mat_collapsed)] <- "mis_UR"
    }
    if (ncol(syn_mat_collapsed) > 0) {
        colnames(syn_mat_collapsed)[ncol(syn_mat_collapsed)] <- "syn_UR"
    }

    # Make genotype matrix
    nonempty_mat <- list()
    if (ncol(lof_mat_collapsed) > 0) nonempty_mat <- c(nonempty_mat, list(lof_mat_collapsed))
    if (ncol(mis_mat_collapsed) > 0) nonempty_mat <- c(nonempty_mat, list(mis_mat_collapsed))
    if (ncol(syn_mat_collapsed) > 0) nonempty_mat <- c(nonempty_mat, list(syn_mat_collapsed))

    G <- do.call(cbind, nonempty_mat)
    lof_ncol <- ncol(lof_mat_collapsed)
    mis_ncol <- ncol(mis_mat_collapsed)
    syn_ncol <- ncol(syn_mat_collapsed)

    print("Dimension of G")
    print(dim(G))

    # Obtain residual vector and genotype matrix with the same order
    y_tilde <- cbind(modglmm$sampleID, modglmm$residuals)
    y_tilde <- y_tilde[which(y_tilde[,1] %in% sampleID),]
    # G <- as.matrix(G)
    G_reordered <- G[match(y_tilde[,1], rownames(G)),]
    if (traitType == "binary") {
        vG_reordered <- as.vector(sqrt(v)) * G_reordered    # Sigma_e^(-1/2) G
    }
    n_samples <- nrow(G_reordered)
    print("Dimension of G_reordered")
    print(dim(G_reordered))

    post_beta_lof <- NULL
    post_beta_mis <- NULL
    post_beta_syn <- NULL

    if (lof_ncol > 0) {
        if (traitType == "binary") {
            fast_lmm_lof <- fast_lmm(G = vG_reordered[,c(1:lof_ncol), drop = F], Y = as.numeric(y_tilde[,2]))
        } else {
            fast_lmm_lof <- fast_lmm(G = G_reordered[,c(1:lof_ncol), drop = F], Y = as.numeric(y_tilde[,2]))
        }
        post_beta_lof <- fast_lmm_lof[[1]]
        tr_GtG_lof <- fast_lmm_lof[[2]]
        delta_lof <- fast_lmm_lof[[3]]
        GtG_lof <- fast_lmm_lof[[4]]
        # h2_lof <- (sigma_sq / delta_lof) * tr_GtG_lof / ((sigma_sq / delta_lof) * tr_GtG_lof + sigma_sq * n_samples)
        tau_lof <- as.numeric(sigma_sq / delta_lof)

        # Adjust variance component tau
        if (traitType == "binary") {
            tau_lof_mom_marginal <- mom_estimator_marginal(G = vG_reordered[,c(1:lof_ncol), drop = F], y = as.numeric(y_tilde[,2]))
            tau_mom_joint <- mom_estimator_joint(
                G1 = vG_reordered[,c(1:lof_ncol), drop = F],
                G2 = vG_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F],
                G3 = vG_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F],
                y = as.numeric(y_tilde[,2])
            )
        } else {
            tau_lof_mom_marginal <- mom_estimator_marginal(G = G_reordered[,c(1:lof_ncol), drop = F], y = as.numeric(y_tilde[,2]))
            tau_mom_joint <- mom_estimator_joint(
                G1 = G_reordered[,c(1:lof_ncol), drop = F],
                G2 = G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F],
                G3 = G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F],
                y = as.numeric(y_tilde[,2])
            )
        }
        tau_lof_adj <- tau_lof * tau_mom_joint[1] / tau_lof_mom_marginal[1]
        if (traitType == "binary") {
            h2_lof_adj <- max(tau_lof_adj * tr_GtG_lof / (tau_lof_adj * tr_GtG_lof + sigma_sq * sum(1/v)), 0)
        } else {
            h2_lof_adj <- max(tau_lof_adj * tr_GtG_lof / (tau_lof_adj * tr_GtG_lof + sigma_sq * n_samples), 0)
        }
    }

    if (mis_ncol > 0) {
        if (traitType == "binary") {
            fast_lmm_mis <- fast_lmm(G = vG_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F], Y = as.numeric(y_tilde[,2]))
        } else {
            fast_lmm_mis <- fast_lmm(G = G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F], Y = as.numeric(y_tilde[,2]))
        }
        post_beta_mis <- fast_lmm_mis[[1]]
        tr_GtG_mis <- fast_lmm_mis[[2]]
        delta_mis <- fast_lmm_mis[[3]]
        GtG_mis <- fast_lmm_mis[[4]]
        # h2_mis <- (sigma_sq / delta_mis) * tr_GtG_mis / ((sigma_sq / delta_mis) * tr_GtG_mis + sigma_sq * n_samples)
        tau_mis <- as.numeric(sigma_sq / delta_mis)

        # Adjust variance component tau
        if (traitType == "binary") {
            tau_mis_mom_marginal <- mom_estimator_marginal(G = vG_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F], y = as.numeric(y_tilde[,2]))
        } else {
            tau_mis_mom_marginal <- mom_estimator_marginal(G = G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F], y = as.numeric(y_tilde[,2]))
        }
        tau_mis_adj <- tau_mis * tau_mom_joint[2] / tau_mis_mom_marginal[1]
        if (traitType == "binary") {
            h2_mis_adj <- max(tau_mis_adj * tr_GtG_mis / (tau_mis_adj * tr_GtG_mis + sigma_sq * sum(1/v)), 0)
        } else {
            h2_mis_adj <- max(tau_mis_adj * tr_GtG_mis / (tau_mis_adj * tr_GtG_mis + sigma_sq * n_samples), 0)
        }
    }

    if (syn_ncol > 0) {
        if (traitType == "binary") {
            fast_lmm_syn <- fast_lmm(G = vG_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F], Y = as.numeric(y_tilde[,2]))
        } else {
            fast_lmm_syn <- fast_lmm(G = G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F], Y = as.numeric(y_tilde[,2]))
        }
        post_beta_syn <- fast_lmm_syn[[1]]
        tr_GtG_syn <- fast_lmm_syn[[2]]
        delta_syn <- fast_lmm_syn[[3]]
        GtG_syn <- fast_lmm_syn[[4]]
        # h2_syn <- (sigma_sq / delta_syn) * tr_GtG_syn / ((sigma_sq / delta_syn) * tr_GtG_syn + sigma_sq * n_samples)
        tau_syn <- as.numeric(sigma_sq / delta_syn)

        # Adjust variance component tau
        if (traitType == "binary") {
            tau_syn_mom_marginal <- mom_estimator_marginal(G = vG_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F], y = as.numeric(y_tilde[,2]))
        } else {
            tau_syn_mom_marginal <- mom_estimator_marginal(G = G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F], y = as.numeric(y_tilde[,2]))
        }
        tau_syn_adj <- tau_syn * tau_mom_joint[3] / tau_syn_mom_marginal[1]
        if (traitType == "binary") {
            h2_syn_adj <- max(tau_syn_adj * tr_GtG_syn / (tau_syn_adj * tr_GtG_syn + sigma_sq * sum(1/v)), 0)
        } else {
            h2_syn_adj <- max(tau_syn_adj * tr_GtG_syn / (tau_syn_adj * tr_GtG_syn + sigma_sq * n_samples), 0)
        }
    }

    # Obtain effect size jointly
    if ((tau_lof_adj > 0) & (tau_mis_adj > 0) & (tau_syn_adj > 0)) {
        if (traitType == "binary") {
            post_beta <- calculate_joint_blup(
                G1 = vG_reordered[,c(1:lof_ncol), drop = F],
                G2 = vG_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F],
                G3 = vG_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F],
                tau1 = tau_lof_adj,
                tau2 = tau_mis_adj,
                tau3 = tau_syn_adj,
                Sigma1 = diag(1, ncol(vG_reordered[,c(1:lof_ncol), drop = F])),
                Sigma2 = diag(1, ncol(vG_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F])),
                Sigma3 = diag(1, ncol(vG_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F])),
                psi = as.numeric(sigma_sq),
                y = as.numeric(y_tilde[,2])
            )
        } else {
            post_beta <- calculate_joint_blup(
                G1 = G_reordered[,c(1:lof_ncol), drop = F],
                G2 = G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F],
                G3 = G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F],
                tau1 = tau_lof_adj,
                tau2 = tau_mis_adj,
                tau3 = tau_syn_adj,
                Sigma1 = diag(1, ncol(G_reordered[,c(1:lof_ncol), drop = F])),
                Sigma2 = diag(1, ncol(G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F])),
                Sigma3 = diag(1, ncol(G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F])),
                psi = as.numeric(sigma_sq),
                y = as.numeric(y_tilde[,2])
            )
        }
    } else {
        # If not, calculate beta marginally
        post_beta <- as.vector(rbind(post_beta_lof, post_beta_mis, post_beta_syn))
    }

    post_beta <- as.vector(post_beta)

    # Obtain prediction error variance (PEV)
    diag(GtG_lof) <- diag(GtG_lof) + as.numeric(sigma_sq) / tau_lof_adj
    diag(GtG_mis) <- diag(GtG_mis) + as.numeric(sigma_sq) / tau_mis_adj
    diag(GtG_syn) <- diag(GtG_syn) + as.numeric(sigma_sq) / tau_syn_adj

    PEV_lof <- diag(solve(GtG_lof))
    PEV_mis <- diag(solve(GtG_mis))
    PEV_syn <- diag(solve(GtG_syn))
    PEV_all <- c(PEV_lof, PEV_mis, PEV_syn)

    # Apply Firth bias correction for binary phenotype
    if (traitType == "binary") {
        y_binary <- cbind(modglmm$sampleID, modglmm$y)
        offset1 <- modglmm$linear.predictors - modglmm$coefficients[1]
        l2.var = 1
        maxit = 50

        # LoF
        G.lof.sp <- as(G_reordered[,c(1:lof_ncol), drop = F], "sparseMatrix")
        nMarker.lof <- ncol(G.lof.sp)
        out_single_wL2_lof_sparse <- Run_Firth_MultiVar_Single(G.lof.sp, modglmm$obj.noK, as.numeric(y_binary[,2]), offset1, nMarker.lof, l2.var=1/(2*tau_lof_adj), Is.Fast=FALSE, Is.Sparse=TRUE)[,2]
        print(out_single_wL2_lof_sparse)

        # mis
        G.mis.sp <- as(G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F], "sparseMatrix")
        nMarker.mis <- ncol(G.mis.sp)
        out_single_wL2_mis_sparse <- Run_Firth_MultiVar_Single(G.mis.sp, modglmm$obj.noK, as.numeric(y_binary[,2]), offset1, nMarker.mis, l2.var=1/(2*tau_mis_adj), Is.Fast=FALSE, Is.Sparse=TRUE)[,2]
        print(out_single_wL2_mis_sparse)

        # syn
        G.syn.sp <- as(G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F], "sparseMatrix")
        nMarker.syn <- ncol(G.syn.sp)
        out_single_wL2_syn_sparse <- Run_Firth_MultiVar_Single(G.syn.sp, modglmm$obj.noK, as.numeric(y_binary[,2]), offset1, nMarker.syn, l2.var=1/(2*tau_syn_adj), Is.Fast=FALSE, Is.Sparse=TRUE)[,2]
        print(out_single_wL2_syn_sparse)

        beta_firth <- c(out_single_wL2_lof_sparse, out_single_wL2_mis_sparse, out_single_wL2_syn_sparse)
        effect <- ifelse(abs(post_beta) < log(2), post_beta, beta_firth)
    } else {
        effect <- post_beta
    }

    effect <- as.vector(effect)
    # output related to single-variant effect size
    variant <- colnames(G)

    if (lof_ncol == 1) {
        sgn <- sign(post_beta[1])
    } else {
        MAC <- colSums(G_reordered[,c(1:(lof_ncol - 1)), drop = F])
        sgn <- sign(sum(MAC * post_beta[c(1:(lof_ncol - 1))]))
    }

    h2_all <- sum(h2_lof_adj, h2_mis_adj, h2_syn_adj) * sgn
    h2 <- c(h2_lof_adj, h2_mis_adj, h2_syn_adj, h2_all)
    group <- c("LoF", "mis", "syn", "all")

    tau_lof_out <- c(tau_lof, as.numeric(tau_lof_mom_marginal[1, 1]), as.numeric(tau_mom_joint[1, 1]))
    tau_mis_out <- c(tau_mis, as.numeric(tau_mis_mom_marginal[1, 1]), as.numeric(tau_mom_joint[2, 1]))
    tau_syn_out <- c(tau_syn, as.numeric(tau_syn_mom_marginal[1, 1]), as.numeric(tau_mom_joint[3, 1]))
    tau_out <- rbind(tau_lof_out, tau_mis_out, tau_syn_out)

    effect_out <- as.data.frame(cbind(variant, effect, PEV_all))
    h2_out <- rbind(group, h2)

    effect_outname <- paste0(outputPrefix, "_effect.txt")
    h2_outname <- paste0(outputPrefix, "_h2.txt")
    tau_outname <- paste0(outputPrefix, "_tau.txt")

    print(effect_out)
    print(h2_out)
    write.table(effect_out, effect_outname, row.names=F, quote=F)
    write.table(h2_out, h2_outname, row.names=F, col.names=F, quote=F)
    write.table(tau_out, tau_outname, row.names=F, col.names=F, quote=F)

    print("Analysis completed.")

}
