# For test
source("~/code/RareEffect-dev/RareEffect.R")
source("~/code/RareEffect-dev/Firth.R")
library(SAIGE, lib.loc = "~/utils/SAIGE")

# Example input
rdaFile <- "~/data/RareEffect/output/step1_train/step1_train_f.30760.0.0.rda"   # HDL step 1 result
chrom <- 11     # APOC3 gene (chr11)
bedFile <- "/media/leelabsg-storage1/kisung/WES_470k/merged_by_chr11.norm.bed"
bimFile <- "/media/leelabsg-storage1/kisung/WES_470k/merged_by_chr11.norm.bim"
famFile <- "/media/leelabsg-storage1/kisung/WES_470k/merged_by_chr11.norm.fam"
groupFile <- "/media/leelabsg-storage1/DATA/UKBB/WES/groupfile/LOFTEE/loftee_groupfile_chr11.withQC.txt"
traitType <- "quantitative"
geneName <- "APOC3"     # APOC3 gene
macThreshold <- 10
collapseLoF <- TRUE    # Collapse LoF variants into one super-variant (default: FALSE)
outputPrefix <- "test"


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
system.time({
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
    # write.table(effect_out, effect_outname, row.names=F, quote=F)
    # write.table(h2_out, h2_outname, row.names=F, col.names=F, quote=F)
    # write.table(tau_out, tau_outname, row.names=F, col.names=F, quote=F)
})
print("Analysis completed.")
