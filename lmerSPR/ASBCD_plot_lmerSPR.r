source("/Users/chr/Prog/R/lmerERP/lmerSPR/plot_lmerSPR.R", encoding = "UTF-8")

produce_spr_plots <- function(
    path
) {
    lmerSPR <- fread(path)

    # Observed Data
    plot_lmerSPR(lmerSPR, "logRT", yunit="logRT", title="Observed RTs") #, ylims=c(5.7, 5.9)

    # Estimated RTs
    estimates <- colnames(lmerSPR)[which(grepl("est_", colnames(lmerSPR)))]
    for (x in estimates) {
         plot_lmerSPR(lmerSPR, x, yunit="logRT", title="Estimated RTs")  #, ylims=c(5.7, 5.9)  
    }

    # Residuals
    residuals <- colnames(lmerSPR)[which(grepl("res_", colnames(lmerSPR)))]
    for (x in residuals) {
        plot_lmerSPR(lmerSPR, x, yunit="logRT", title="Residuals: Noun Association + log(Cloze)")  # , ylims=c(0.06, -0.06)
    }

    # Coefficients
    coefs_cols <- colnames(lmerSPR)[grep("coef", colnames(lmerSPR))]
    SE_cols <- colnames(lmerSPR)[grep("SE_", colnames(lmerSPR))]
    cols <- c("Condition", "Item", "Subject", "Region", coefs_cols, SE_cols)
    data <- lmerSPR[,..cols]
    data1 <- melt(data, id.vars=c("Condition", "Region"), measure.vars=coefs_cols, variable.name="Coefficient", value.name="coefficients")
    data1$SE <- melt(data, id.vars=c("Condition", "Region"), measure.vars=SE_cols, variable.name="Coefficient", value.name="SE")$SE
    data1 <- data1[, lapply(.SD, mean), by=list(Region, Coefficient), .SDcols=c("coefficients", "SE")]
    plot_lmerSPR(data1, "coefficients", yunit="SPR Coefficients", title="Coefficients", grouping="Coefficient", ylims=c(5.7,5.85)) 

    # Z-values
    zval_cols <- setdiff(colnames(lmerSPR)[grep("zval_", colnames(lmerSPR))], "zval_1")
    pval_cols <- setdiff(colnames(lmerSPR)[grep("pval_", colnames(lmerSPR))], "pval_1")
    cols <- c("Condition", "Item", "Subject", "Region", zval_cols, pval_cols)
    data <- lmerSPR[,..cols]
    data1 <- melt(data, id.vars=c("Condition", "Region"), measure.vars=zval_cols, variable.name="Zvalue", value.name="zvalue")
    data1$pvalue <- melt(data, id.vars=c("Condition", "Region"), measure.vars=pval_cols, variable.name="Pvalue", value.name="pvalue")$pvalue
    data1 <- data1[, lapply(.SD, mean), by=list(Region, Zvalue), .SDcols=c("zvalue", "pvalue")]
    plot_lmerSPR(data1, "zvalue", yunit="Z-values", title="SPR Effect Size", ylims=c(-0.5, 4), grouping="Zvalue") 

    # Condition A predictions
    dt <- fread("data/lmerSPR_Asrp.csv")
    srp_range <- range(dt$z_srp)
    dt <- dt[, lapply(.SD, mean), by=list(Region, Subject), .SDcols=c("1", "srp", "i_1", "i_srp", "s_1", "s_srp")]
    
    # To Do : column mult to get subject variation in there (random effects)
    dt$est_max <- (dt$"1" + dt$i_1 + dt$s_1) + (dt$srp + dt$s_srp + dt$i_srp) * srp_range[[1]]
    dt$est_avg <- (dt$"1" + dt$i_1 + dt$s_1) + (dt$srp + dt$s_srp + dt$i_srp) * 0
    dt$est_1SD <- (dt$"1" + dt$i_1 + dt$s_1) + (dt$srp + dt$s_srp + dt$i_srp)  * 1
    dt$est_min <- (dt$"1" + dt$i_1 + dt$s_1) + (dt$srp + dt$s_srp + dt$i_srp)  * srp_range[[2]]

    est_cols <- colnames(dt)[grep("est_", colnames(dt))]
    dt_m <- melt(dt, id.vars=c("Region", "Subject"), measure.vars=est_cols, variable.name="estimate", value.name="logRT")
    #dt_c <- dcast(dt_m, formula = Region + logCloze ~ Electrode, value.var="estimates")
    plot_lmerSPR(dt_m, "logRT", yunit="logRT", title="Estimated RTs", grouping="estimate", ylims=c(5.7,5.85)) #, ylims=c(-0.5, 4)) 
}
