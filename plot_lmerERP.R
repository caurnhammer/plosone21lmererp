##
# Christoph Aurnhammer <aurnhammer@coli.uni-saarland.de>
#
# EEG plotting options for (lme)rERPs
##

# compute standard error
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

# Return only the legend of an ggplot object (s/o to some person on the internet)
get_legend<-function(
        a.gplot
){
        tmp <- ggplot_gtable(ggplot_build(a.gplot + theme(legend.box="vertical", legend.spacing.y = unit(0.005, 'inch'))))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

# For a single Electrode, plot per-grouping mean (Bootstrapped confidence intervals optional. Slow.)
plot_grandavg_ci <- function(
        dt,                         # subset of data
        ttl,                        # title
        yunit = paste0("Amplitude (", "\u03BC", "Volt\u29"),   # y-axis unit 
        subject_avg = TRUE,         # first average across subjects
        ci = FALSE,                 # add bootrs. conf. intervals
        ylims = NULL,               # optional length=2 vector of ylimits (order descrending)
        grouping = "Condition",     # provide grouping within which data is splitted for plotting
        tws = list(c(350, 450), c(600, 800))
) {     
        if (subject_avg == TRUE) {
            colnames(dt)[ncol(dt)] <- "V2"
            dt <- dt[, lapply(.SD, mean), by=c(grouping, "Subject", "Timestamp"), .SDcols="V2"]
            if ((subject_avg == TRUE & ci == FALSE) == TRUE) {
                dt <- dt[, lapply(.SD, mean), by=c(grouping, "Timestamp"), .SDcols="V2"]
            }
        }

        if (ci == TRUE & grouping != "Coefficient") {
            sedf <- dt[, lapply(.SD, se), by=c(grouping, "Timestamp"), .SDcols="V2"]
            colnames(sedf)[3] <- "SE"
            sedf$SE = sedf$SE * 2
            dt <- merge(dt, sedf, on=c(grouping, "Timestamp"))
            dt <- dt[, lapply(.SD, mean), by=c(grouping, "Timestamp"), .SDcols=c("V2", "SE")]
        }

        if (subject_avg == FALSE & ci == FALSE) {
            if (grouping == "zvalue"){
                colnames(dt)[c(6,7)] <- c("V2", "V3")
                dt <- dt[, lapply(.SD, mean), by=c(grouping, "Timestamp"), .SDcols=c("V2", "V3")]
                colnames(dt)[4] <- "sig"
                df <- dt[,c("zvalue", "Timestamp", "sig")]
                df$posit <- rep(seq(ylims[1]-2, ylims[1], length=length(unique(df$zvalue))), length(unique(df$Timestamp)))   
                df$sig <- factor(df$sig, levels=c("1", "0"), labels=c("sign", "insign"))
            } else if (grouping == "logCloze") {
                colnames(dt)[3] <- "V2"
            } else {
                colnames(dt)[c(6,7)] <- c("V2", "V3")
                dt <- dt[, lapply(.SD, mean), by=c(grouping, "Timestamp"), .SDcols=c("V2", "V3")]            
            }
        }
        
        if (grouping == "Coefficient") {
            plt <- ggplot(dt, aes(x = Timestamp, y=V2, color=eval(parse(text=paste0("dt$", grouping))), fill=eval(parse(text=paste0("dt$", grouping)))))
        } else {
            plt <- ggplot(dt, aes(x = Timestamp, y = V2, color=eval(parse(text=paste0("dt$", grouping))), fill=eval(parse(text=paste0("dt$", grouping)))))
        }
        plt <- plt + scale_y_reverse()
        if (is.vector(ylims) == TRUE) {plt <- plt + ylim(ylims[1], ylims[2])}
        plt <- plt + geom_line()
        plt <- plt + geom_hline(yintercept=0, linetype="dashed")
        plt <- plt + geom_vline(xintercept=0, linetype="dashed")
        if (ci == TRUE) {plt <- plt + geom_ribbon(aes(x=Timestamp, ymax=V2+SE, ymin=V2-SE), alpha=0.20, color=NA)
        } else if (grouping == "Coefficient") {plt <- plt + geom_ribbon(aes(x=Timestamp, ymax=V2+V3, ymin=V2-V3), alpha=0.20, color=NA)}
        plt <- plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = 0.1, linetype = "solid"), 
                           panel.grid.major = element_line(size = 0.3, linetype = "solid", color = "#A9A9A9"), 
                           panel.grid.minor = element_line(size = 0.15, linetype = "solid", color = "#A9A9A9"),
                           legend.position = "top")
        plt <- plt + labs(y=yunit, x= "Time (ms)", title = ttl)
        if (grouping == "Condition") {
            plt <- plt + scale_color_manual(name=grouping, labels=c("A: A+E+", "B: A-E+", "C: A+E-", "D: A-E-"), values=c("#000000", "#BB5566", "#004488", "#DDAA33"))
            plt <- plt + scale_fill_manual(name=grouping, labels=c("A: A+E+", "B: A-E+", "C: A+E-", "D: A-E-"), values=c("#000000", "#BB5566", "#004488", "#DDAA33"))
            #plt <- plt + scale_color_manual(name=grouping, labels=c("C: A+E-", "D: A-E-"), values=c("#004488", "#DDAA33"))
            #plt <- plt + scale_fill_manual(name=grouping, labels=c("C: A+E-", "D: A-E-"), values=c("#004488", "#DDAA33"))
            #plt <- plt + scale_color_manual(name=grouping, labels=c("A: A+E+", "C: A+E-"), values=c("#000000", "#004488"))
            #plt <- plt + scale_fill_manual(name=grouping, labels=c("A: A+E+", "C: A+E-"), values=c("#000000", "#004488"))
        } else if (grouping == "Coefficient") {
            plt <- plt + scale_color_manual(name=grouping, labels=c("Intercept", "Noun Association", "log(Cloze)"), values=c("#000000", "#00FFFF", "#E349F6", "#FFA500"))
            plt <- plt + scale_fill_manual(name=grouping, labels=c("Intercept", "Noun Association", "log(Cloze)"), values=c("#000000", "#00FFFF", "#E349F6", "#FFA500"))
            #plt <- plt + scale_color_manual(name=grouping, labels=c("Intercept", "log(Cloze)"), values=c("#000000", "#E349F6", "#FFA500"))
            #plt <- plt + scale_fill_manual(name=grouping, labels=c("Intercept", "log(Cloze)"), values=c("#000000", "#E349F6", "#FFA500"))
        } else if (grouping == "zvalue") {
            plt <- plt + scale_color_manual(name="Predictor", labels=c("Noun Association", "log(Cloze)"), values=c("#00FFFF", "#E349F6", "#FFA500"))
            plt <- plt + scale_fill_manual(name="Predictor", labels=c("Noun Association", "log(Cloze)"), values=c("#00FFFF", "#E349F6", "#FFA500"))
            #plt <- plt + scale_color_manual(name="Predictor", labels=c("log(Cloze)"), values=c("#E349F6", "#FFA500"))
            #plt <- plt + scale_fill_manual(name="Predictor", labels=c("log(Cloze)"), values=c("#E349F6", "#FFA500"))
            plt <- plt + geom_point(data=df, aes(x=Timestamp, y=posit, shape=sig)) 
            plt <- plt + scale_shape_manual(values=c(20, 32), name="Corrected p-values", labels=c("Significant", "Nonsignificant"))
            plt <- plt + annotate("rect", xmin = tws[1][[1]][1], xmax = tws[1][[1]][2], ymin = ylims[1], ymax = ylims[2], alpha = .15)
            plt <- plt + annotate("rect", xmin = tws[2][[1]][1], xmax = tws[2][[1]][2], ymin = ylims[1], ymax = ylims[2], alpha = .15) 
        } else {
            plt <- plt + scale_color_manual(name="log(Cloze)", labels=c("Maximum", "Average", "1 SD", "Minimum"), values=c("#ff0000", "#000000", "#E349F6", "#495cf6"))
            plt <- plt + scale_fill_manual(name="log(Cloze)", labels=c("Maximum", "Average", "1 SD", "Minimum"), values=c("#ff0000", "#000000", "#E349F6", "#495cf6"))
        }

        plt
}

# Plot midline electrodes.
plot_midline <- function(
        data,                     # input data
        file = FALSE,             # where to store. If FALSE, display.
        title = "Midline ERPs",   # Add title
        yunit = paste0("Amplitude (", "\u03BC", "Volt\u29"), # y-axis label
        subject_avg = TRUE,       # Compute subject-average before plotting (affects width of CIs)
        ci = FALSE,               # include conf interval
        ylims = NULL,             # custom ylims to fix axis scale across electrodes (order descending)
        grouping = "Condition"    # provide grouping within which data is splitted for plotting
) {
        data[,grouping] <- as.factor(eval(parse(text = paste0("data$", grouping))))
        if (grouping == "logCloze") {
          cols <- c(grouping, "Subject", "Timestamp")
        } else {
          cols <- c(grouping, "Item", "Subject", "Timestamp", "TrialNum")
        }
        
        # Make individual plots
        if (grouping == "zvalue") {
            Fzplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, c("zval_Fz", "sig_Fz")]), "Fz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
            Czplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, c("zval_Cz", "sig_Cz")]), "Cz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
            Pzplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, c("zval_Pz", "sig_Pz")]), "Pz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
        } else if (grouping == "Coefficient") {
            Fzplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, c("coefval_Fz", "seval_Fz")]), "Fz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
            Czplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, c("coefval_Cz", "seval_Cz")]), "Cz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
            Pzplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, c("coefval_Pz", "seval_Pz")]), "Pz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
        } else {
            Fzplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, "Fz"]), "Fz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
            Czplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, "Cz"]), "Cz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
            Pzplt <- plot_grandavg_ci(cbind(data[, ..cols], data[, "Pz"]), "Pz", yunit = yunit, subject_avg=subject_avg, ci=ci, ylims=ylims, grouping=grouping)
        }
        
        # Get the legend
        legend <- get_legend(Fzplt)

        # Arrange
        gg <- arrangeGrob(arrangeGrob(Fzplt + theme(legend.position = "none"),
                                 Czplt + theme(legend.position = "none"),
                                 Pzplt + theme(legend.position = "none"), layout_matrix=matrix(seq_len(3 * 1))), legend, heights=c(10, 1.2), top=textGrob(title))
        if (file != FALSE) {
           ggsave(file, gg, device=cairo_pdf, width = 4, height = 7) 
        } else {
           gg
        }
}
