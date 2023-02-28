###########################################
#' Estimation of species richness in a community

ChaoSpecies <- function(data, datatype = c("abundance","abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"), k = 10, conf = 0.95)
{
  if (is.matrix(data) == T || is.data.frame(data) == T){
    #if (ncol(data) != 1 & nrow(data) != 1)
    #stop("Error: The data format is wrong.")
    if(datatype != "incidence_raw"){
      if (ncol(data) == 1){
        data <- data[, 1]
      } else {
        data <- data[1, ]
      }
    } else{
      t <- ncol(data)
      dat <- rowSums(data)
      dat <- as.integer(dat)
      t_infreq <- sum(colSums(data[which(dat<k),])>=1)
      data <- dat
      data <- c(t_infreq, t , data)
    }
    
  }
  
  if(datatype == "abundance_freq_count"){
    data <- as.integer(data)
    length4b <- length(data)
    data <- rep(data[seq(1,length4b,2)],data[seq(2,length4b,2)])
    names(data) <- paste("x",1:length(data),sep="")
    datatype <- "abundance"
  }
  if (datatype == "incidence_freq_count"){
    t <- as.integer(data[1])
    data <- data[-c(1)]
    data <- as.integer(data)
    lengthdat <- length(data)
    data <- rep(data[seq(1,lengthdat,2)],data[seq(2,lengthdat,2)])
    data <- c(t,data)
    names(data) <- c("T", paste("y",1:(length(data)-1),sep=""))
    datatype <- "incidence_freq"
  }
  method <- "all"
  if (k != round(k) || k < 0)
    stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
  if (is.numeric(conf) == FALSE || conf > 1 || conf < 0)
    stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  #  data <- as.numeric(round(data))
  
  if (datatype == "abundance"){
    f <- function(i, data){length(data[which(data == i)])}
    if (f(1, data) == sum(data)){
      stop("Error: The information of data is not enough.")}
    z <- (list(Basic_data_information = basicAbun(data, k)[[1]], Rare_species_group = RareSpeciesGroup(data, k),
               Species_table = round(SpecAbunOut(data, method, k, conf), 3)))
  } else if (datatype == "incidence_raw"){
    dat <- data[-1]; Q <- function(i, data){length(data[which(data == i)])}
    if (Q(1, dat) == sum(dat)){
      stop("Error: The information of data is not enough.")}
    z <- (list(Basic_data_information = basicInci(data[-1], k)[[1]], Infreq_species_group = InfreqSpeciesGroup(data[-1], k),
               Species_table = round(SpecInciOut_raw(data, method, k, conf),3)))
  } else if (datatype == "incidence_freq"){
    dat <- data[-1];
    Q <- function(i, data){length(data[which(data == i)])}
    if (Q(1, dat) == sum(dat)){
      stop("Error: The information of data is not enough.")}
    z <- (list(Basic_data_information = basicInci(data, k)[[1]], Infreq_species_group = InfreqSpeciesGroup(data, k),
               Species_table = round(SpecInciOut(data, method, k, conf),3)))
  } 
  else{
    stop("Error: The data type is wrong.")
  }
  class(z) <- c("ChaoSpecies")
  z
}


###########################################
#' Estimation of the number of shared species between two communities/assemblages

ChaoShared <-
  function(data, datatype = c("abundance", "incidence_freq", "incidence_raw"), units,
           se = TRUE, nboot = 200, conf = 0.95) {
    
    method <- "all"
    if (se == TRUE) {
      if (nboot < 1)
        nboot <- 1
      if (nboot == 1)
        cat("Warning: When \"nboot\" =" ,nboot, ", the bootstrap s.e. and confidence interval can't be calculated.",
            "\n\n")
    }
    
    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
      cat("Warning: \"conf\"(confidence level) must be a numerical value between 0 and 1, e.g. 0.95.",
          "\n")
      cat("          We use \"conf\" = 0.95 to calculate!",
          "\n\n")
      conf <- 0.95
    }
    
    datatype <- match.arg(datatype)
    if (datatype == "abundance") {
      if(class(data)=="list"){data <-cbind(data[[1]],data[[2]]) }
      x1 <- data[, 1]
      x2 <- data[, 2]
      Basic <- BasicFun(x1, x2, nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
      output <- ChaoShared.Ind(x1, x2, method, nboot, conf, se)
      colnames(output) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    }
    if (datatype == "incidence_freq") {
      if(class(data)=="list"){data <-cbind(data[[1]],data[[2]]) }
      y1 <- data[, 1]
      y2 <- data[, 2]
      Basic <- BasicFun(y1, y2, B=nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
      output <- ChaoShared.Sam(y1, y2, method, conf, se)
      colnames(output) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    }
    if (datatype=="incidence_raw"){
      t = units
      if(ncol(data) != sum(t)) stop("Number of columns does not euqal to the sum of key in sampling units")
      dat <- matrix(0, ncol = length(t), nrow = nrow(data))
      n <- 0
      for(i in 1:length(t)){
        dat[, i] <- as.integer(rowSums(data[,(n+1):(n+t[i])] ) )
        n <- n+t[i]
      }
      t <- as.integer(t)
      dat <- apply(dat, MARGIN = 2, as.integer)
      dat <- data.frame(rbind(t, dat),row.names = NULL)
      y1 <- dat[,1]
      y2 <- dat[,2]
      datatype = "incidence_freq"
      Basic <- BasicFun(y1, y2, B=nboot, datatype)
      output <- ChaoShared.Sam(y1, y2, method, conf, se)
      colnames(output) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    }
    out <- list(Basic_data_information=Basic,
                Estimation_results=output)
    class(out) <- c("ChaoShared")
    return(out)
  }

###########################################
#' Estimation of species diversity (Hill numbers)


Diversity=function(data, datatype=c("abundance","abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"), q=NULL)
{
  if (is.matrix(data) == T || is.data.frame(data) == T){
    #if (ncol(data) != 1 & nrow(data) != 1)
    #stop("Error: The data format is wrong.")
    if(datatype != "incidence_raw"){
      if (ncol(data) == 1){
        data <- data[, 1]
      } else {
        data <- as.vector(data[1, ])
      }
    } else{
      t <- ncol(data)
      dat <- rowSums(data)
      dat <- as.integer(dat)
      data <- c(t , dat)
    }
    
  }
  X <- data
  if(datatype == "abundance_freq_count"){
    data <- as.integer(data)
    length4b <- length(data)
    data <- rep(data[seq(1,length4b,2)],data[seq(2,length4b,2)])
    names(data) <- paste("x",1:length(data),sep="")
    datatype <- "abundance"
    X <- data
  }
  if(datatype=="abundance"){
    type="abundance"
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    
    BASIC.DATA <- matrix(round(c(sum(X), sum(X>0), 1-sum(X==1)/sum(X), CV.Ind(X)),3), ncol = 1)
    nickname <- matrix(c("n", "D", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("    Sample size", "    Number of observed species",
                              "    Estimated sample coverage",
                              "    Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    
    table0 <- matrix(0,5,4)
    table0[1,]=c(Chao1(X)[-5])
    table0[2,]=c(Chao1_bc(X))
    table0[3,]=round(SpecAbuniChao1(X, k=10, conf=0.95)[1,],1)
    table0[4,]=round(c(SpecAbunAce(X)),1)
    table0[5,]=round(c(SpecAbunAce1(X)),1)
    colnames(table0) <- c("Estimate", "s.e.", paste(Chao1(X)[5]*100,"%Lower", sep=""), paste(Chao1(X)[5]*100,"%Upper", sep=""))
    rownames(table0) <- c("    Chao1 (Chao, 1984)","    Chao1-bc ", "    iChao1","    ACE (Chao & Lee, 1992)",
                          "    ACE-1 (Chao & Lee, 1992)")
    
    SHANNON=Shannon_index(X)
    table1=round(SHANNON[c(1:5),],3)
    table1=table1[-2,]              ##2016.05.09
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Jackknife",
    #                      " Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("     MLE","     Jackknife",
                          "     Chao & Shen","     Chao et al. (2013)")
    
    table1_exp=round(SHANNON[c(6:10),],3)
    table1_exp=table1_exp[-2,]      ##2016.05.09
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Jackknife",
    #                         " Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("     MLE","     Jackknife",
                              "     Chao & Shen","     Chao et al. (2013)")
    
    table2=round(Simpson_index(X)[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("     MVUE","     MLE")
    
    table2_recip=round(Simpson_index(X)[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("     MVUE","     MLE")
    
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=NULL, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=q, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2
    
    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- nrow(Hill)
    rownames(Hill) <- paste("    ",1:q_hill)
    z <- list("datatype"= type,"Basic_data"=BASIC.DATA,"Species_richness"=table0,
              "Shannon_index"=table1,"Shannon_diversity"=table1_exp,
              "Simpson_index"=table2,"Simpson_diversity"=table2_recip,
              "Hill_numbers"= Hill)
  }
  if(datatype == "incidence_freq_count"){
    t <- as.integer(data[1])
    data <- data[-c(1)]
    data <- as.integer(data)
    lengthdat <- length(data)
    data <- rep(data[seq(1,lengthdat,2)],data[seq(2,lengthdat,2)])
    data <- c(t,data)
    names(data) <- c("T", paste("y",1:(length(data)-1),sep=""))
    datatype <- "incidence_freq"
    X <- data
  }
  if(datatype=="incidence_freq"){
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    type="incidence"
    U<-sum(X[-1])
    D<-sum(X[-1]>0)
    T<-X[1]
    C<-Chat.Sam(X,T)
    CV_squre<-max( D/C*T/(T-1)*sum(X[-1]*(X[-1]-1))/U^2-1, 0)
    CV<-CV_squre^0.5
    BASIC.DATA <- matrix(round(c(D,T,U, C, CV),3), ncol = 1)
    nickname <- matrix(c("D", "T","U", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("    Number of observed species", "    Number of Sampling units","    Total number of incidences",
                              "    Estimated sample coverage",
                              "    Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    #BASIC.DATA <- basicInci(X, k=10)[[1]]
    ############################################################
    table0=SpecInci(X, k=10, conf=0.95)
    rownames(table0) <- c("    Chao2 (Chao, 1987)","    Chao2-bc ", "    iChao2","    ICE (Lee & Chao, 1994)",
                          "    ICE-1 (Lee & Chao, 1994)")
    SHANNON=Shannon_Inci_index(X)
    table1=round(SHANNON[c(1,4),],3)
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("     MLE","     Chao et al. (2013)")
    table1_exp=round(SHANNON[c(5,8),],3)
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("     MLE","     Chao et al. (2013)")
    
    SIMPSON=Simpson_Inci_index(X)
    table2=round(SIMPSON[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("     MVUE","     MLE")
    
    table2_recip=round(SIMPSON[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("     MVUE","     MLE")
    
    
    ############################################################
    #Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", from=0, to=3, interval=0.25, B=50, conf=0.95))
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence_freq", q=NULL, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence_freq", q=q, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2
    
    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- nrow(Hill)
    rownames(Hill) <- paste("    ",1:q_hill)
    #z <- list("BASIC.DATA"=BASIC.DATA,"HILL.NUMBERS"= Hill)
    z <- list("datatype"= type,"Basic_data"=BASIC.DATA,"Species_richness"=table0,
              "Shannon_index"=table1,"Shannon_diversity"=table1_exp,
              "Simpson_index"=table2,"Simpson_diversity"=table2_recip,
              "Hill_numbers"= Hill)
  }
  if(datatype=="incidence_raw"){
    type="incidence"
    datatype = "incidence_freq"
    U<-sum(X[-1])
    D<-sum(X[-1]>0)
    T<-X[1]
    C<-Chat.Sam(X,T)
    CV_squre<-max( D/C*T/(T-1)*sum(X[-1]*(X[-1]-1))/U^2-1, 0)
    CV<-CV_squre^0.5
    BASIC.DATA <- matrix(round(c(D,T,U, C, CV),3), ncol = 1)
    nickname <- matrix(c("D", "T","U", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("    Number of observed species", "    Number of Sampling units","    Total number of incidences",
                              "    Estimated sample coverage",
                              "    Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    ############################################################
    table0=SpecInci(X, k=10, conf=0.95)
    rownames(table0) <- c("    Chao2 (Chao, 1987)","    Chao2-bc ", "    iChao2","    ICE (Lee & Chao, 1994)",
                          "    ICE-1 (Lee & Chao, 1994)")
    SHANNON=Shannon_Inci_index(X)
    table1=round(SHANNON[c(1,4),],3)
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("     MLE","     Chao et al. (2013)")
    table1_exp=round(SHANNON[c(5,8),],3)
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("     MLE","     Chao et al. (2013)")
    
    SIMPSON=Simpson_Inci_index(X)
    table2=round(SIMPSON[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("     MVUE","     MLE")
    
    table2_recip=round(SIMPSON[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("     MVUE","     MLE")
    
    
    ############################################################
    #Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", from=0, to=3, interval=0.25, B=50, conf=0.95))
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=NULL, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=q, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2
    
    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- nrow(Hill)
    rownames(Hill) <- paste("   ",1:q_hill)
    #z <- list("BASIC.DATA"=BASIC.DATA,"HILL.NUMBERS"= Hill)
    
    z <- list("datatype"= type,"Basic_data"=BASIC.DATA,"Species_richness"=table0,
              "Shannon_index"=table1,"Shannon_diversity"=table1_exp,
              "Simpson_index"=table2,"Simpson_diversity"=table2_recip,
              "Hill_numbers"= Hill)
  }
  class(z) <- c("spadeDiv")
  return(z)
}

###########################################
#' Estimation of two-assemblage similarity measures

SimilarityPair=function(X, datatype = c("abundance","incidence_freq", "incidence_raw"), units,nboot=200)
{ 
  
  if(datatype=="abundance")
  {
    if(class(X)=="list"){X <- do.call(cbind,X)}
    type <- "abundance"
    info1 <- c("S.total"=sum(rowSums(X)>0), "n1"=sum(X[,1]), "n2"=sum(X[,2]), 
               "D1"=sum(X[,1]>0), "D2"=sum(X[,2]>0), "D12"=sum(X[,1]>0 & X[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("f[11]"=sum(X[,1]==1 & X[,2]==1), 
               "f[1+]"=sum(X[,1]==1 & X[,2]>0), "f[+1]"=sum(X[,1]>0 & X[,2]==1),
               "f[2+]"=sum(X[,1]==2 & X[,2]>0), "f[+2]"=sum(X[,1]>0 & X[,2]==2),"f[22]"=sum(X[,1]==2 & X[,2]==2))
    info <- c(info1, info2)
    ################################################################2016.07.08-(P.L.Lin)
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    temp <- list()
    weight <- c(sum(X[,1])/(sum(X[,1])+sum(X[,2])), sum(X[,2])/(sum(X[,1])+sum(X[,2])))
    weight <- - sum(weight*log(weight)) / log(2)
    mat <- Jaccard_Sorensen_Abundance_equ(datatype,X[, 1],X[, 2], nboot)[, c(1, 2)]
    mat <- cbind(mat, mat[, 1]-1.96*mat[, 2],  mat[, 1]+1.96*mat[, 2])
    MLE_Jaccard <- mat[1, ]
    Est_Jaccard <- mat[2, ]
    MLE_Sorensen <- mat[3, ]
    Est_Sorensen <- mat[4, ]
    mat2 <-  Two_Horn_equ(X[,1], X[,2], method="all", weight="unequal", nboot = nboot)
    MLE_Ee_Horn <- mat2$mle
    MLE_Ee_Horn <- plus_CI(c(MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_Horn <- mat2$est
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1],Est_Ee_Horn[2]))
    mat3 <- Two_BC_equ(X[, 1],X[, 2], datatype="abundance", nboot)
    MLE_Ee_Braycurtis <- mat3$mle
    Est_Ee_Braycurtis <- mat3$est
    mat4 <- SimilarityTwo(X,2,nboot,method="unequal weight")
    MLE_Ee_C22 <- mat4$CqN[1, ]
    Est_Ee_C22 <- mat4$CqN[2, ]
    MLE_Ee_U22 <- mat4$UqN[1, ]
    Est_Ee_U22 <- mat4$UqN[2, ]
    mat5 <- Two_Horn_equ(X[,1], X[,2], method="all", weight="equal", nboot = nboot)
    MLE_ew_Horn <- mat5$mle
    Est_ew_Horn <- mat5$est
    mat6 <- SimilarityTwo(X,2,nboot,method="equal weight")
    MLE_ew_C22 <- mat6$CqN[1, ]
    Est_ew_C22 <- mat6$CqN[2, ]
    MLE_ew_U22 <- mat6$UqN[1, ]
    Est_ew_U22 <- mat6$UqN[2, ]
    #MLE_ew_Braycurtis <- plus_CI(MLE_Braycurtis_equ(X[,1],X[,2],w1=0.5))
    #Est_ew_Braycurtis <- plus_CI(KH_Braycurtis_equ(X[,1],X[,2],w1=0.5))
    MLE_ew_ChaoSoresen <- mat[11,]
    Est_ew_ChaoSoresen <- mat[12, ]
    MLE_ew_ChaoJaccard <- mat[9, ]
    Est_ew_ChaoJaccard <- mat[10, ]
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22, MLE_ew_ChaoJaccard, MLE_ew_ChaoSoresen)
    rownames(temp[[2]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[3]] <- t(as.matrix(MLE_Ee_Horn))
    rownames(temp[[3]]) <- c("Horn size weighted(q=1)")  
    temp[[4]] <- rbind(MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[4]]) <- c("C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")  
    temp[[5]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[5]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[6]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22, Est_ew_ChaoJaccard, Est_ew_ChaoSoresen)
    rownames(temp[[6]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[7]] <- t(as.matrix(Est_Ee_Horn))
    rownames(temp[[7]]) <- c("Horn size weighted(q=1)")  
    temp[[8]] <- rbind(Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    temp <- lapply(temp, FUN = function(x){
      colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") 
      return(x)
    })
    rownames(temp[[8]]) <- c("C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")  
    z <- list("datatype"=type,"info"=info, "Empirical_richness"=temp[[1]], "Empirical_relative"=temp[[2]], "Empirical_WtRelative"=temp[[3]],
              "Empirical_absolute"=temp[[4]], "estimated_richness"=temp[[5]], "estimated_relative"=temp[[6]], "estimated_WtRelative"=temp[[7]], "estimated_absolute"=temp[[8]]) 
  }      
  ##---------------------------------------------------------------
  if(datatype=="incidence_raw"){
    data <- X
    t = units
    if(ncol(data) != sum(t)) stop("Number of columns does not euqal to the sum of key in sampling units")
    dat <- matrix(0, ncol = length(t), nrow = nrow(data))
    n <- 0
    for(i in 1:length(t)){
      dat[, i] <- as.integer(rowSums(data[,(n+1):(n+t[i])] ) )
      n <- n+t[i]
    }
    t <- as.integer(t)
    dat <- apply(dat, MARGIN = 2, as.integer)
    dat <- data.frame(rbind(t, dat),row.names = NULL)
    y1 <- dat[,1]
    y2 <- dat[,2]
    X <- cbind(y1, y2)
    type <- "incidence_freq"
    X <- as.data.frame(X)
  }
  if(datatype=="incidence_freq") type <- "incidence_freq" 
  if(datatype=="incidence_freq" | type == "incidence_freq")
  { 
    if(class(X)=="list"){X <- do.call(cbind,X)}
    no.assemblage=length(X[1,])
    Y=X[-1,]  
    type="incidence"
    info1 <- c("S.total"=sum(rowSums(Y)>0), "T1"=X[1,1], "T2"=X[1,2], "U1"=sum(Y[,1]), "U2"=sum(Y[,2]), 
               "D1"=sum(Y[,1]>0), "D2"=sum(Y[,2]>0), "D12"=sum(Y[,1]>0 & Y[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("Q[11]"=sum(Y[,1]==1 & Y[,2]==1), 
               "Q[1+]"=sum(Y[,1]==1 & Y[,2]>0), "Q[+1]"=sum(Y[,1]>0 & Y[,2]==1),
               "Q[2+]"=sum(Y[,1]==2 & Y[,2]>0), "Q[+2]"=sum(Y[,1]>0 & Y[,2]==2),  "Q[22]"=sum(Y[,1]==2 & Y[,2]==2))
    info <- c(info1, info2)
    
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    temp <- list()
    weight <- c(sum(Y[,1])/(sum(Y[,1])+sum(Y[,2])), sum(Y[,2])/(sum(Y[,1])+sum(Y[,2])))
    weight <- - sum(weight*log(weight)) / log(2)
    mat <- Jaccard_Sorensen_Abundance_equ(datatype="incidence",X[, 1],X[, 2], nboot)[, c(1, 2)]
    mat <- cbind(mat, mat[, 1]-1.96*mat[, 2],  mat[, 1]+1.96*mat[, 2])
    MLE_Jaccard <- mat[1, ]
    Est_Jaccard <- mat[2, ]
    MLE_Sorensen <- mat[3, ]
    Est_Sorensen <- mat[4, ]
    mat2 <-  Two_Horn_equ(X[,1], X[,2], datatype = "incidence", method="all", weight="unequal", nboot)
    MLE_Ee_Horn <- mat2$mle
    MLE_Ee_Horn <- plus_CI(c(MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_Horn <- mat2$est
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1],Est_Ee_Horn[2]))
    mat3 <- C2N_ee_se_inc(X, nboot)
    MLE_Ee_C22 <- plus_CI(mat3[1,])
    Est_Ee_C22 <- plus_CI(mat3[3,])
    MLE_Ee_U22 <- plus_CI(mat3[2,])
    Est_Ee_U22 <- plus_CI(mat3[4,])
    mat4 <- Two_Horn_equ(X[,1], X[,2], datatype = "incidence", method="all", weight="equal", nboot)
    MLE_ew_Horn <- mat4$mle
    Est_ew_Horn <- mat4$est
    mat5 <- SimilarityTwo(X, 2, nboot, method="equal weight", datatype="incidence")
    MLE_ew_C22 <- mat5$CqN[1, ]
    Est_ew_C22 <- mat5$CqN[2, ]
    MLE_ew_U22 <- mat5$UqN[1, ]
    Est_ew_U22 <- mat5$UqN[2, ]
    MLE_ew_ChaoSoresen <- mat[11,]
    Est_ew_ChaoSoresen <- mat[12, ]
    MLE_ew_ChaoJaccard <- mat[9, ]
    Est_ew_ChaoJaccard <- mat[10, ]
    mat5 <- Two_BC_equ(X[, 1],X[, 2], datatype="incidence", nboot)
    MLE_Ee_Braycurtis <- mat5$mle
    Est_Ee_Braycurtis <- mat5$est
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22, MLE_ew_ChaoJaccard, MLE_ew_ChaoSoresen)
    rownames(temp[[2]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[3]] <- t(as.matrix(MLE_Ee_Horn))
    rownames(temp[[3]]) <- c("Horn size weighted(q=1)")  
    temp[[4]] <- rbind(MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[4]]) <- c("C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")  
    temp[[5]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[5]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[6]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22, Est_ew_ChaoJaccard, Est_ew_ChaoSoresen)
    rownames(temp[[6]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[7]] <- t(as.matrix(Est_Ee_Horn))
    rownames(temp[[7]]) <- c("Horn size weighted(q=1)")  
    temp[[8]] <- rbind(Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    rownames(temp[[8]]) <- c("C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")  
    temp <- lapply(temp, FUN = function(x){
      colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") 
      return(x)
    })
    z <- list("datatype"=type,"info"=info, "Empirical_richness"=temp[[1]], "Empirical_relative"=temp[[2]], "Empirical_WtRelative"=temp[[3]],
              "Empirical_absolute"=temp[[4]], "estimated_richness"=temp[[5]], "estimated_relative"=temp[[6]], "estimated_WtRelative"=temp[[7]], "estimated_absolute"=temp[[8]]) 
    
  }  
  class(z) <- c("spadeTwo")
  return(z)   
}


###########################################
#' Estimation of multiple-community similarity measures



SimilarityMult=function(X,datatype=c("abundance","incidence_freq", "incidence_raw"),units,q=2,nboot=200,goal="relative")
{ 
  method <- goal
  if(datatype=="abundance"){
    if(class(X)=="list"){X <- do.call(cbind,X)}
    type <- "abundance"
    N <- no.community <- ncol(X)
    temp <- c("N"=ncol(X), "S.total"=sum(rowSums(X)>0))
    n <- apply(X,2,sum)
    D <- apply(X,2,function(x)sum(x>0))
    
    if(N > 2){
      temp1 <- temp2 <- rep(0, N*(N-1)/2)
      k <- 1
      for(i in 1:(N-1)){     
        for(j in (i+1):N){
          temp1[k] <- paste('D',i,j,sep="")
          temp2[k] <- sum(X[,i]>0 & X[,j]>0)
          k <- k + 1
        }
      }
    }
    names(temp2) <- temp1
    names(n) <- paste('n',1:N, sep="")
    names(D) <- paste('D',1:N, sep="")
    info <- c(temp, n, D, temp2)
    if(N == 3) info <- c(temp, n, D, temp2, D123=sum(X[,1]>0 & X[,2]>0 & X[,3]>0))
    info <- c(info, nboot=nboot)
    ################################################################2016.07.11-(P.L.Lin)
    temp <- list()
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    n <- apply(X = X, MARGIN = 2, FUN = sum)
    weight <- n/sum(n)
    weight <- - sum(weight*log(weight)) / log(N)
    mat <- SimilarityMul(X, 0, nboot, method ="unequal weight")
    MLE_Jaccard <- mat$UqN[1, ]
    Est_Jaccard <- mat$UqN[2, ]
    MLE_Sorensen <- mat$CqN[1, ]
    Est_Sorensen <- mat$CqN[2, ]
    mat2 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("unequal"))
    MLE_Ee_Horn <- mat2$mle
    Est_Ee_Horn <- mat2$est
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1], Est_Ee_Horn[2]))
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1], MLE_Ee_Horn[2]))
    mat3 <- BC_equ(X, datatype="abundance", nboot)
    MLE_Ee_Braycurtis <- mat3$mle
    Est_Ee_Braycurtis <- mat3$est
    mat4 <- SimilarityMul(X,2,nboot,method="unequal weight")
    MLE_Ee_C22 <- mat4$CqN[1, ]
    Est_Ee_C22 <- mat4$CqN[2, ]
    MLE_Ee_U22 <- mat4$UqN[1, ]
    Est_Ee_U22 <- mat4$UqN[2, ]
    mat5 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("equal"))
    MLE_ew_Horn <- mat5$mle
    Est_ew_Horn <- mat5$est
    mat6 <- SimilarityMul(X,2,nboot,method="equal weight")
    MLE_ew_C22 <- mat6$CqN[1, ]
    Est_ew_C22 <- mat6$CqN[2, ]
    MLE_ew_U22 <- mat6$UqN[1, ]
    Est_ew_U22 <- mat6$UqN[2, ]
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22)
    rownames(temp[[2]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")
    temp[[3]] <- t(as.matrix(MLE_Ee_Horn))
    rownames(temp[[3]]) <- c("Horn size weighted(q=1)")  
    temp[[4]] <- rbind(MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[4]]) <- c("C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")  
    temp[[5]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[5]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[6]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22)
    rownames(temp[[6]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")
    temp[[7]] <- t(as.matrix(Est_Ee_Horn))
    rownames(temp[[7]]) <- c("Horn size weighted(q=1)")  
    temp[[8]] <- rbind(Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    rownames(temp[[8]]) <- c("C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")    
    temp <- lapply(temp, FUN = function(x){
      colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") 
      return(x)
    })
    if(q == 0){
      temp_PC <- rep(0, N*(N-1)/2)
      C02=matrix(0,choose(no.community,2),4)
      U02=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(1,N,N)
      C_SM_2=matrix(1,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          mat <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal effort')
          C02[k,] <- mat[1, ]
          U02[k,] <- mat[2, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C02[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- U02[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C02"=C02, "U02"=U02)
      C_SM <- list("C02"=C_SM_1, "U02"=C_SM_2)
    }
    if(q == 1 & method=="relative"){
      temp_PC <- rep(0, N*(N-1)/2)
      C12=matrix(0,choose(no.community,2),4)
      Horn=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(1,N,N)
      C_SM_2=matrix(1,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          C12[k,] <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal weight')[1, ]
          Horn[k,] <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal effort')[2, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C12[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- Horn[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C12"=C12, "Horn"=Horn)
      C_SM <- list("C12"=C_SM_1, "Horn"=C_SM_2)
    }
    if(q == 1 & method=="absolute"){
      temp_PC <- rep(0, N*(N-1)/2)
      k=1
      C_SM_1=matrix(1,N,N)
      C12=matrix(0,choose(no.community,2),4)
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          C12[k,] <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal effort')[1, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C12[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C12"=C12)
      C_SM <- list("C12"=C_SM_1)
    }
    if(q == 2){
      temp_PC <- rep(0, N*(N-1)/2)
      if(method=="absolute") method2 <- 'equal effort'
      if(method=="relative") method2 <- 'equal weight'
      C22=matrix(0,choose(no.community,2),4)
      U22=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(1,N,N)
      C_SM_2=matrix(1,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          mat <- Cq2_est_equ(X[,c(i,j)], q, nboot, method=method2)
          C22[k,] <- mat[1, ]
          U22[k,] <- mat[2, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C22[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- U22[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C22"=C22, "U22"=U22)
      C_SM <- list("C22"=C_SM_1,"U22"=C_SM_2)
    }
    Cqn_PC <- lapply(Cqn_PC, function(x){
      colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") ; rownames(x) <- temp_PC
      return(x)
    })
    z <- list("datatype"=datatype,"info"=info, "Empirical_richness"=temp[[1]], "Empirical_relative"=temp[[2]], "Empirical_WtRelative"=temp[[3]],
              "Empirical_absolute"=temp[[4]], "estimated_richness"=temp[[5]], "estimated_relative"=temp[[6]], "estimated_WtRelative"=temp[[7]], "estimated_absolute"=temp[[8]], "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "goal"=method, "q"=q)
  }
  if(datatype == "incidence_raw"){
    data <- X
    t = units
    if(ncol(data) != sum(t)) stop("Number of columns does not euqal to the sum of key in sampling units")
    dat <- matrix(0, ncol = length(t), nrow = nrow(data))
    n <- 0
    for(i in 1:length(t)){
      dat[, i] <- as.integer(rowSums(data[,(n+1):(n+t[i])] ) )
      n <- n+t[i]
    }
    t <- as.integer(t)
    dat <- apply(dat, MARGIN = 2, as.integer)
    X <- data.frame(rbind(t, dat),row.names = NULL)
    if(ncol(X) <= 2) stop("Multiple Commumity measures is only for the data which has three community or more")
    type = "incidence_freq"
  }
  if(datatype=="incidence_freq") type <- "incidence_freq" 
  if(datatype=="incidence_freq" | type == "incidence_freq"){
    if(class(X)=="list"){X <- do.call(cbind,X)}
    type <- "incidence"
    Y <- X
    X <- X[-1,]
    t <- as.vector(Y[1,])
    N <- no.community <- ncol(X)
    temp <- c("N"=ncol(X), "S.total"=sum(rowSums(X)>0))
    n <- apply(X,2,sum)
    D <- apply(X,2,function(x)sum(x>0))
    if(N > 2){
      temp1 <- temp2 <- rep(0, N*(N-1)/2)
      k <- 1
      for(i in 1:(N-1)){     
        for(j in (i+1):N){
          temp1[k] <- paste('D',i,j,sep="")
          temp2[k] <- sum(X[,i]>0 & X[,j]>0)
          k <- k + 1
        }
      }
    }
    names(temp2) <- temp1
    names(t) <- paste('T',1:N, sep="")
    names(n) <- paste('u',1:N, sep="")
    names(D) <- paste('D',1:N, sep="")
    info <- c(temp, t, n, D, temp2)
    if(N == 3) info <- c(temp, t, n, D, temp2, D123=sum(X[,1]>0 & X[,2]>0 & X[,3]>0))
    info <- unlist(c(info, nboot=nboot))
    ################################################################2016.07.20-(P.L.Lin)
    temp <- list()
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    n <- apply(X = X, MARGIN = 2, FUN = sum)
    weight <- n/sum(n)
    weight <- - sum(weight*log(weight)) / log(N)
    mat <- SimilarityMul(Y, 0, nboot, method ="unequal weight", datatype="incidence")
    MLE_Jaccard <- mat$UqN[1, ]
    Est_Jaccard <- mat$UqN[2, ]
    MLE_Sorensen <- mat$CqN[1, ]
    Est_Sorensen <- mat$CqN[2, ]
    mat2 <- Horn_Multi_equ(Y, datatype="incidence", nboot, method=c("unequal"))
    MLE_Ee_Horn <- mat2$mle
    Est_Ee_Horn <- mat2$est
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1], Est_Ee_Horn[2]))
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1], MLE_Ee_Horn[2]))
    mat3 <- BC_equ(Y, datatype="incidence", nboot)
    MLE_Ee_Braycurtis <- mat3$mle
    Est_Ee_Braycurtis <- mat3$est
    mat4 <- C2N_ee_se_inc(Y, nboot)
    MLE_Ee_C22 <- plus_CI(mat4[1, ])
    Est_Ee_C22 <- plus_CI(mat4[3, ])
    MLE_Ee_U22 <- plus_CI(mat4[2, ])
    Est_Ee_U22 <- plus_CI(mat4[4, ])
    mat5 <- Horn_Multi_equ(Y, datatype="incidence", nboot, method=c("equal"))
    MLE_ew_Horn <- mat5$mle
    Est_ew_Horn <- mat5$est
    mat6 <- SimilarityMul(Y, 2, nboot, datatype = "incidence", method="equal weight")
    MLE_ew_C22 <- mat6$CqN[1, ]
    Est_ew_C22 <- mat6$CqN[2, ]
    MLE_ew_U22 <- mat6$UqN[1, ]
    Est_ew_U22 <- mat6$UqN[2, ]
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22)
    rownames(temp[[2]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")  
    temp[[3]] <- t(as.matrix(MLE_Ee_Horn))
    rownames(temp[[3]]) <- c("Horn size weighted(q=1)")  
    temp[[4]] <- rbind(MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[4]]) <- c("C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")  
    temp[[5]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[5]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[6]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22)
    rownames(temp[[6]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")   
    temp[[7]] <- t(as.matrix(Est_Ee_Horn))
    rownames(temp[[7]]) <- c("Horn size weighted(q=1)")  
    temp[[8]] <- rbind(Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    rownames(temp[[8]]) <- c("C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")    
    temp <- lapply(temp, FUN = function(x){
      colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") 
      return(x)
    })
    if(q == 0){
      temp_PC <- rep(0, N*(N-1)/2)
      C02=matrix(0,choose(no.community,2),4)
      U02=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(1,N,N)
      C_SM_2=matrix(1,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          mat <- Cq2_est_equ(Y[,c(i,j)], q, nboot,datatype="incidence", method='equal effort')
          C02[k,] <- mat[1, ]
          U02[k,] <- mat[2, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C02[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- U02[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C02"=C02, "U02"=U02)
      C_SM <- list("C02"=C_SM_1, "U02"=C_SM_2)
    }
    if(q == 1 & method=="relative"){
      temp_PC <- rep(0, N*(N-1)/2)
      C12=matrix(0,choose(no.community,2),4)
      Horn=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(1,N,N)
      C_SM_2=matrix(1,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          C12[k,] <- Cq2_est_equ(Y[,c(i,j)], q, nboot,datatype="incidence", method='equal weight')[1, ]
          Horn[k,] <- Cq2_est_equ(Y[,c(i,j)], q, nboot,datatype="incidence", method='equal effort')[2, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C12[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- Horn[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C12"=C12, "Horn"=Horn)
      C_SM <- list("C12"=C_SM_1, "Horn"=C_SM_2)
    }
    if(q == 1 & method=="absolute"){
      temp_PC <- rep(0, N*(N-1)/2)
      k=1
      C_SM_1=matrix(1,N,N)
      C12=matrix(0,choose(no.community,2),4)
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          C12[k,] <- Cq2_est_equ(Y[,c(i,j)], q, nboot,datatype="incidence", method='equal effort')[1, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C12[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C12"=C12)
      C_SM <- list("C12"=C_SM_1)
    }
    if(q == 2){
      temp_PC <- rep(0, N*(N-1)/2)
      if(method=="absolute") method2 <- 'equal effort'
      if(method=="relative") method2 <- 'equal weight'
      C22=matrix(0,choose(no.community,2),4)
      U22=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(1,N,N)
      C_SM_2=matrix(1,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          mat <- Cq2_est_equ(Y[,c(i,j)], q, nboot,datatype="incidence", method=method2)
          C22[k,] <- mat[1, ]
          U22[k,] <- mat[2, ]
          temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C22[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- U22[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C22"=C22, "U22"=U22)
      C_SM <- list("C22"=C_SM_1,"U22"=C_SM_2)
    }
    Cqn_PC <- lapply(Cqn_PC, function(x){
      colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") ; rownames(x) <- temp_PC
      return(x)
    })
    z <- list("datatype"=datatype,"info"=info, "Empirical_richness"=temp[[1]], "Empirical_relative"=temp[[2]], "Empirical_WtRelative"=temp[[3]],
              "Empirical_absolute"=temp[[4]], "estimated_richness"=temp[[5]], "estimated_relative"=temp[[6]], "estimated_WtRelative"=temp[[7]], "estimated_absolute"=temp[[8]], "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "goal"=method, "q"=q)
  }
  class(z) <- c("spadeMult")
  z
}


