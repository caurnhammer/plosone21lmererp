library(ggplot2)

se <- function(
    x,
    na.rm = FALSE
){
    if (na.rm == TRUE) {
        sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
    } else {
        sd(x)/sqrt(length(x))
    }
}

plot_lmerSPR <- function(
    data,
    DV,
    yunit,
    title, 
    ylims = NULL,
    grouping = "Condition"
) {
    if (!(DV %in% c("coefficients", "zvalue"))){
        pm <- aggregate(data[[DV]] ~ Region + data[[grouping]] + Subject, data, FUN=mean)
        colnames(pm)[c(2,4)] <- c("group", "DV")
        
        plusminus <- aggregate(DV ~ Region + group, pm, FUN=mean)
        plusminus$SE <- aggregate(DV ~ Region + group, pm, FUN=se)[,3]
        plusminus <- data.table(plusminus)
    } else if (DV == "zvalue") {
        plusminus <- data
        colnames(plusminus)[c(2,3)] <- c("group", "DV")
        
        plusminus$sig <- plusminus$pvalue < 0.05
        df <- plusminus[,c("Region", "group", "pvalue", "sig")]
        df$posit <- rep(seq(ylims[1]-0.9, ylims[1]-0.5, length=length(unique(plusminus$group))), each=length(unique(plusminus$Region)))
        df$sig <- factor(df$sig, levels=c("TRUE", "FALSE"), labels=c("sign", "insign"))
        df$Region <- plusminus$Region
    } else {
        plusminus <- data
        colnames(plusminus)[c(2,3)] <- c("group", "DV")
    }
    plusminus$Region <- plyr::mapvalues(plusminus$Region, c("critical -1", "critical", "postcritical", "postpostcritical"), c("Pre-Critical", "Critical", "Spillover", "Post-Spillover"))
    plusminus$Region <- factor(plusminus$Region, levels=c("Pre-Critical", "Critical", "Spillover", "Post-Spillover"))
    if (DV == "zvalue"){df$Region <- plusminus$Region}
 
    p <- ggplot(plusminus, aes(x=Region, y=DV, color=group, group=group)) + geom_point(size=2.5, shape="cross") + geom_line(size=0.5)
    p <- p + theme_minimal() 
    if (!(DV %in% c("zvalue", "coefficients"))){
        p <- p + geom_errorbar(aes(ymin=DV-SE, ymax=DV+SE), width=.1, size=0.3)
    }
    if (DV == "coefficients") { 
        p <- p + geom_errorbar(aes(ymin=DV-SE, ymax=DV+SE), width=.1, size=0.3)
        #p <- p + scale_color_manual(name="Coefficients", values=c("#000000", "#00FFFF", "#E349F6", "#FFA500"), labels=c("Intercept", "Noun Association", "log(Cloze)"))
        p <- p + scale_color_manual(name="Coefficients", values=c("#000000", "#E349F6"), labels=c("Intercept", "log(Cloze)"))
    } else if (DV == "zvalue") {
        p <- p + geom_hline(yintercept=0, linetype="dashed")
        p <- p + scale_color_manual(name="Z-value", values=c("#00FFFF", "#E349F6", "#FFA500"), labels=c("Noun Association", "log(Cloze)"))
        p <- p + geom_point(data=df, aes(x=Region, y=posit, shape=sig), size=2.5) 
        p <- p + scale_shape_manual(values=c(20, 32), name="P-values", labels=c("Significant", "Nonsignificant"))
    } else if (grouping == "estimate"){
        p <- p + scale_color_manual(name="log(Cloze)", labels=c("Maximum", "Average", "1 SD", "Minimum"), values=c("#ff0000", "#000000", "#E349F6", "#495cf6"))
    }
    else { # RTs, Residuals
        p <- p + scale_color_manual(name="Condition", labels=c("A: A+E+", "B: A-E+", "C: A+E-", "D: A-E-"), values=c("#000000", "#BB5566", "#004488", "#DDAA33"))
        #p <- p + scale_color_manual(name="Condition", labels=c("A: A+E+", "C: A+E-"), values=c("#000000", "#004488"))
        #p <- p + scale_color_manual(name="Condition", labels=c("C: A+E-", "D: A-E-"), values=c("#004488", "#DDAA33"))
    }

    #if (is.vector(ylims) == TRUE) { p <- p + ylim(ylims[1]-1, ylims[2]) }
    if (is.vector(ylims) == TRUE) { p <- p + ylim(ylims[1], ylims[2]) }
    p <- p + theme(plot.title = element_text(size = 8),
                    axis.text.x= element_text(size=7),
                    legend.position="bottom", 
                    legend.text=element_text(size=5), legend.title=element_text(size=4), 
                    legend.box="vertical", legend.spacing.y=unit(-0.2, "cm"), 
                    legend.margin=margin(0,0,0,0),
                    legend.box.margin=margin(-10,-10,-10,-50))
    p <- p + labs(x="Region", y=yunit, title=title)
    p

    file <- paste0("plots/", DV, ".pdf")
    ggsave(file, p, width=3, height=3)
}
