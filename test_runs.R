library(icesSAG)
source("R/Estimate_Blim.R")

name_data <- "ICESStocks_2021_eqsim_corrected_fishlife.Rdata" # all stocks already there! Could also use SAG directly
#name_data <- unzip(zipfile = "data/ICES_stock_objects_WKREF1_2021.zip", list = TRUE)[,"Name"] # could also use SAG directly


tmpdir <- tempdir()

for (i in 1:length(name_data)){
  
  name_file <- load(unzip("data/ICES_stock_objects_WKREF1_2021.zip", name_data[i], exdir=tmpdir))
  data <- eval(parse(text=name_file))
  
  pdf(paste0("Blim_estim_", name_file, "_2021.pdf"))
  
  for (k in 1:length(data)){
    
    stk <- data[[k]]
    
    if (name(stk)=="cod.27.1-2coastN") name(stk) <- "cod.27.1-2.coastN"
    
    # SAG database
    if (name(stk) %in% c("reb.27.1-2", "reg.27.1-2")) Asstkey <- findAssessmentKey(stock = name(stk), year = 2020) else Asstkey <- findAssessmentKey(stock = name(stk), year = 2021)
    ICES_Blim <- getFishStockReferencePoints(Asstkey)[[1]]["Blim"]
    
    # Estimation of Blim with different methods
    R <- rec(stk)[drop=TRUE]
    S <- ssb(stk)[drop=TRUE]
    years <- as.numeric(names(R))
    Rage <- as.numeric(rownames(rec(stk)))
    par(mfrow=c(2,2), mar=c(2,1,2,1), oma=c(1,1,2,1))
    type1 <- Blim_estim(R,S,Rage,years, is.spasmodic=TRUE, doplot=TRUE)
    if(is.na(type1)) mtext(attr(type1, "code"), side=3, line=0.5, cex=0.6) else mtext(paste0(substr(attr(type1, "code"), 1, 6), ", Blim = ", round(type1)), side=3, line=0.5)
    text(x=par("usr")[2]*0.7, y=par("usr")[4]*0.9, labels = "is.spasmodic=TRUE", cex=0.8)
    type2.6 <- Blim_estim(R,S,Rage,years, doplot=TRUE, log=FALSE)
    text(x=par("usr")[2]*0.9, y=par("usr")[4]*0.9, labels = "types=2:6", cex=0.8)
    if(is.na(type2.6)) mtext(attr(type2.6, "code"), side=3, line=0.5, cex=0.6) else mtext(paste0(substr(attr(type2.6, "code"), 1, 6), ", Blim = ", round(type2.6)), side=3, line=0.5)
    type2 <- Blim_estim(R,S,Rage,years, types=2, doplot=TRUE, log=FALSE)
    text(x=par("usr")[2]*0.9, y=par("usr")[4]*0.9, labels = "types=2", cex=0.8)
    if(is.na(type2)) mtext(attr(type2, "code"), side=3, line=0.5, cex=0.6) else mtext(paste0(substr(attr(type2, "code"), 1, 6), ", Blim = ", round(type2)), side=3, line=0.5)
    type3.6 <- Blim_estim(R,S,Rage,years, types=3:6, doplot=TRUE)
    text(x=par("usr")[2]*0.9, y=par("usr")[4]*0.9, labels = "types=3:6", cex=0.8)
    if(is.na(type3.6)) mtext(attr(type3.6, "code"), side=3, line=0.5, cex=0.6) else mtext(paste0(substr(attr(type3.6, "code"), 1, 6), ", Blim = ", round(type3.6)), side=3, line=0.5)
    mtext(paste0(name(stk), ", ICES 2021 Blim = ", ICES_Blim), side=3, line=0.5, cex=1.2, outer=TRUE)
    
    rm(stk, Asstkey, ICES_Blim); gc()
  }
  
  dev.off()
  
}
