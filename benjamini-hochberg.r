# Benjamini-Hochberg procedure

bh_apply <- function(data, alpha=0.05, time_windows=list(c(300, 500), c(600, 1000))) {
        tws <- time_windows
        dt <- data[, lapply(.SD, mean), by=list(Timestamp, Electrode), .SDcols=colnames(data)[grep("pval", colnames(data))]]
        for (i in colnames(dt)[grep("pval", colnames(dt))]) {
                dt[, paste0("sig_", i)] <- rep(FALSE, nrow(dt))
                for (j in tws) {                        
                        signi <- (p.adjust(dt[Timestamp >= j[1] & Timestamp <= j[2], ..i][[1]], method='fdr') < alpha)
                        indi <- which(dt$Timestamp >= j[1] & dt$Timestamp <= j[2])
                        dt[indi, paste0("sig_", i) := signi]
                }
        }
        
        cols <- c("Timestamp", "Electrode", colnames(dt)[grep("sig", colnames(dt))])
        merge(data, dt[,..cols], by=c("Timestamp", "Electrode"))
}

bh_apply_bonf <- function(data, alpha=0.05, time_windows=list(c(300, 500), c(600, 1000))) {
        tws <- time_windows
        dt <- data[, lapply(.SD, mean), by=list(Timestamp, Electrode), .SDcols=colnames(data)[grep("pval", colnames(data))]]
        for (i in colnames(dt)[grep("pval", colnames(dt))]) {
                dt[, paste0("sig_", i)] <- rep(FALSE, nrow(dt))
                for (j in tws) {                        
                        signi <- (p.adjust(dt[Timestamp >= j[1] & Timestamp <= j[2], ..i][[1]], method='none') < alpha)
                        indi <- which(dt$Timestamp >= j[1] & dt$Timestamp <= j[2])
                        dt[indi, paste0("sig_", i) := signi]
                }
        }
        
        cols <- c("Timestamp", "Electrode", colnames(dt)[grep("sig", colnames(dt))])
        merge(data, dt[,..cols], by=c("Timestamp", "Electrode"))
}

bh_apply_sepelec <- function(data, alpha=0.05, time_windows=list(c(300, 500), c(600, 1000))) {
        tws <- time_windows
        elecs <- c("Fz", "Cz", "Pz")
        dt <- data[, lapply(.SD, mean), by=list(Timestamp, Electrode), .SDcols=colnames(data)[grep("pval", colnames(data))]]
        for (i in colnames(dt)[grep("pval", colnames(dt))]) {
                dt[, paste0("sig_", i)] <- rep(FALSE, nrow(dt))
                for (j in tws) {
                        for (e in elecs) {
                                signi <- (p.adjust(dt[Timestamp >= j[1] & Timestamp <= j[2] & Electrode == e, ..i][[1]], method='fdr') < alpha)
                                indi <- which(dt$Timestamp >= j[1] & dt$Timestamp <= j[2] & dt$Electrode == e)
                                dt[indi, paste0("sig_", i) := signi]
                        }                        
                }
        }
        
        cols <- c("Timestamp", "Electrode", colnames(dt)[grep("sig", colnames(dt))])
        merge(data, dt[,..cols], by=c("Timestamp", "Electrode"))
}


bh_apply_poolall <- function(data, alpha=0.05, time_windows=list(c(300, 500), c(600, 1000))) {
        tws <- time_windows
        elecs <- c("Fz", "Cz", "Pz")
        dt <- data[, lapply(.SD, mean), by=list(Timestamp, Electrode), .SDcols=colnames(data)[grep("pval", colnames(data))]]
        for (i in colnames(dt)[grep("pval", colnames(dt))]) {
                dt[, paste0("sig_", i)] <- rep(FALSE, nrow(dt))
                        signi <- (p.adjust(dt[(Timestamp >= 300 & Timestamp <= 500) | (Timestamp >= 600 & Timestamp <= 1000), ..i][[1]], method='fdr') < alpha)
                        indi <- which((dt$Timestamp >= 300 & dt$Timestamp <= 500) | (dt$Timestamp >= 600 & dt$Timestamp <= 1000))
                        dt[indi, paste0("sig_", i) := signi]
        }
        
        cols <- c("Timestamp", "Electrode", colnames(dt)[grep("sig", colnames(dt))])
        merge(data, dt[,..cols], by=c("Timestamp", "Electrode"))
}

bh_none <- function(data, alpha=0.05, time_windows=list(c(300, 500), c(600, 1000))) {
        tws <- time_windows
        elecs <- c("Fz", "Cz", "Pz")
        dt <- data[, lapply(.SD, mean), by=list(Timestamp, Electrode), .SDcols=colnames(data)[grep("pval", colnames(data))]]
        for (i in colnames(dt)[grep("pval", colnames(dt))]) {
                dt[, paste0("sig_", i)] <- rep(FALSE, nrow(dt))
                signi <- dt[(dt$Timestamp >= 300 & dt$Timestamp <= 500) | (dt$Timestamp >= 600 & dt$Timestamp <= 1000),..i] < alpha
                indi <- which((dt$Timestamp >= 300 & dt$Timestamp <= 500) | (dt$Timestamp >= 600 & dt$Timestamp <= 1000))
                dt[indi, paste0("sig_", i) := signi]
        }

        cols <- c("Timestamp", "Electrode", colnames(dt)[grep("sig", colnames(dt))])
        merge(data, dt[,..cols], by=c("Timestamp", "Electrode"))
}

bh_identify_cutoff <- function(pvals, alpha = 0.05)
{
        pvals <- sort(pvals)     # sort ascending
        m     <- length(pvals)   # number of hypotheses
        k     <- 1               # k iterator
        for (k in seq(1:m)){
                if (pvals[k] >= k/m * alpha) {
                        break()
                }           
        }
        # while (pvals[k] < (k / m) * alpha)
        #         k <- k + 1
        return(pvals[k])
}

bh_plot <- function(pvals, alpha = 0.05)
{
        pvals <- sort(pvals)    # sort ascending
        m     <- length(pvals)  # number of hypothesis 
        ks    <- seq(1, length(pvals), 1)
        plot(pvals ~ ks)
        abline(0, alpha / m)
}
