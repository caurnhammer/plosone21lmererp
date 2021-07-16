##
# Christoph Aurnhamer <aurnhammer@coli.uni-saarland.de>
# 
# Plot ASBCD lmerERP data
##

require(data.table)
require(ggplot2)
require(gridExtra)
require(grid)

### SESSION PREPARATION
# load plotting functions into workspace
source("plot_lmerERP.r")
source("benjamini-hochberg.r")

## Function definitions
# produce single midline plot. 
lmerERPplot <- function(
    data,
    var= NULL,
    yunit="",
    title="",
    ylims=NULL,
    subject_avg = TRUE,
    ci = TRUE,
    mode = "",
    nointercept = TRUE
) {
    if (mode == "coef") {
        coefs_cols <- colnames(data)[grep("coef", colnames(data))]
        se_cols <- colnames(data)[grep("SE_", colnames(data))]
        cols <- c("Condition", "Item", "Subject", "Timestamp", "TrialNum", "Electrode", coefs_cols, se_cols)
        data <- data[,..cols]
        data1 <- melt(data, id.vars=c("Condition", "Item", "Subject", "Timestamp", "TrialNum", "Electrode"), measure.vars=coefs_cols, variable.name="Coefficient", value.name="coefval")
        data1$seval <- melt(data, id.vars=c("Condition", "Item", "Subject", "Timestamp", "TrialNum", "Electrode"), measure.vars=se_cols, variable.name="SE", value.name="seval")$seval
        data1 <- dcast(data1, formula = Condition + Item + Subject + Timestamp + TrialNum + Coefficient ~ Electrode, value.var=c("coefval", "seval"))
        plot_midline(data=data1, file="plots/coefficients.pdf", title="lmerERP Coefficients", yunit="Intercept + coefficient", subject_avg=subject_avg, ci=ci, ylims=ylims, grouping="Coefficient")
    } else if (mode == "zval") {
        zval_cols <- colnames(data)[grep("zval", colnames(data))]
        sig_cols <- colnames(data)[grep("sig", colnames(data))]
        if (nointercept == TRUE){
            zval_cols <- setdiff(zval_cols, "zval_1")
            sig_cols <- setdiff(sig_cols, "sig_pval_1")
        }
        cols <- c("Condition", "Item", "Subject", "Timestamp", "TrialNum", "Electrode", zval_cols, sig_cols)
        data <- data[,..cols]
        data1 <- melt(data, id.vars=c("Condition", "Item", "Subject", "Timestamp", "TrialNum", "Electrode"), measure.vars=zval_cols, variable.name="zvalue", value.name="zval")
        data1$sig <- melt(data, id.vars=c("Condition", "Item", "Subject", "Timestamp", "TrialNum", "Electrode"), measure.vars=sig_cols, variable.name="signif", value.name="sig")$sig
        data1 <- dcast(data1, formula = Condition + Item + Subject + Timestamp + TrialNum + zvalue ~ Electrode, value.var=c("zval", "sig"))
        plot_midline(data=data1, file="plots/zvalues.pdf", title="lmerERP Effect Sizes", yunit="Z-value", subject_avg=FALSE, ci=FALSE, ylims=ylims, grouping="zvalue")
    } else if (mode == "A_estimates") {
        srp_range <- range(data$z_srp)
        data <- data[, lapply(.SD, mean), by=list(Timestamp, Electrode, Subject), .SDcols=c("1", "srp", "s_1", "s_srp", "i_1", "i_srp")]

        data$est_max <- (data$"1" + data$i_1 + data$s_1) + (data$srp + data$s_srp + data$i_srp) * srp_range[[1]]
        data$est_avg <- (data$"1" + data$i_1 + data$s_1) + (data$srp + data$s_srp + data$i_srp) * 0
        data$est_1SD <- (data$"1" + data$i_1 + data$s_1) + (data$srp + data$s_srp + data$i_srp) * 1
        data$est_min <- (data$"1" + data$i_1 + data$s_1) + (data$srp + data$s_srp + data$i_srp) * srp_range[[2]]

        est_cols <- colnames(data)[grep("est_", colnames(data))]
        data_m <- melt(data, id.vars=c("Timestamp", "Electrode", "Subject"), measure.vars=est_cols, variable.name="logCloze", value.name="estimates")
        data_c <- dcast(data_m, formula = Timestamp + logCloze + Subject ~ Electrode, value.var="estimates")
        plot_midline(data=data_c, file="plots/A.pdf", title=title, yunit=yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping="logCloze")
    } else {
        data <- data[,c("Condition", "Item", "Subject", "Timestamp", "TrialNum", "Electrode", ..var)]
        data <- dcast(data, formula = Condition + Item + Subject + Timestamp + TrialNum ~ Electrode, value.var=var)
        plot_midline(data=data, file=paste0("plots/", var, ".pdf"), title=title, yunit=yunit, subject_avg=subject_avg, ci=ci, ylims=ylims)
    }
}

# produce midline plot for observed, lmerERP, residual, coefficients, z-values
produce_lmer_plots <- function(path){
    # LMER ERP DATA
    lmerERP <- fread(path)

    # Observed data
    lmerERPplot(lmerERP, "EEG", yunit=paste0("Amplitude (", "\u03BC", "Volt)"), title="Observed ERP", ylims=c(8.5, -7), ci=TRUE)

    # Estimated Waveforms
    estimates <- colnames(lmerERP)[which(grepl("est_", colnames(lmerERP)))]
    for (x in estimates) {
        lmerERPplot(lmerERP, x, yunit=paste0("Amplitude (", "\u03BC", "Volt)"), title="Estimated ERPs", ylims=c(8.5, -7), ci=TRUE)    
    }
    
    # Residuals
    residuals <- colnames(lmerERP)[which(grepl("res_", colnames(lmerERP)))]
    for (x in residuals) {
        lmerERPplot(lmerERP, x, yunit=paste0("Amplitude (", "\u03BC", "Volt)"), title="Residuals: Noun Association", ci=TRUE, ylims=c(4, -4))    
    }
    
    # Coefficients
    #lmerERPplot(lmerERP, mode="coef", ci=FALSE, subject_avg=FALSE, ylims=c(15.5, -6))
    lmerERPplot(lmerERP, mode="coef", ci=FALSE, subject_avg=FALSE, ylims=c(7.5, -6))

    # Apply multiple comparisons correction
    lmerERP <- bh_apply(lmerERP, alpha=0.05, time_windows=list(c(300, 500), c(600, 1000)))
    
    lmerERP <- fread(path)
    source("benjamini-hochberg.r")
    lmerERP <- bh_apply(lmerERP, alpha=0.05, time_windows=list(c(350, 450), c(600, 800)))
    # Z-values
    lmerERPplot(lmerERP, mode="zval", ylims=c(9.5, -9.5))

    # condition A estimates
    condA <- fread("data/lmerERP_Asrp.csv")
    lmerERPplot(condA, mode="A_estimates", ylims=c(10, -10), ci=TRUE, subject_avg=TRUE, title="Estimated ERPs: Log(Cloze) Levels")

    condA <- bh_apply(condA, alpha=0.05, time_windows=list(c(350, 450), c(600, 800)))
    lmerERPplot(condA, mode="zval", ylims=c(9.5,-9.5))
}

