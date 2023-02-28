###TEST SCRIPT

Sim <- SimilarityPair(recoverytsf_list[[1]][[7]],"incidence_freq", nboot = 200)

X <- recoverytsf_list[[1]][[2]]

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
##################

#mat <- Jaccard_Sorensen_Abundance_equ(datatype="incidence",X[, 1],X[, 2], nboot)[, c(1, 2)]
# Jaccard_Sorensen_Abundance_equ=function(datatype = c("abundance", "incidence"),X1,X2,boot)

X1 = X[, 1]
X2 = X[, 2]
boot = nboot

shared.species.hat=PanEstFun.Sam(X1,X2)
alpha.species.hat=SpecInciChao2(X1, k=10, conf=0.95)[1,1]+SpecInciChao2(X2, k=10, conf=0.95)[1,1]
Esti.Jaccard=shared.species.hat/(alpha.species.hat-shared.species.hat)
Esti.Sorensen=2*shared.species.hat/alpha.species.hat
w=X1[1];z=X2[1]
X1=X1[-1];X2=X2[-1]

n1=sum(X1);n2=sum(X2)
if(datatype == "abundance"){w=n1;z=n2}
I=which(X1>0 & X2>0)

MLE.Jaccard=length(I)/sum(X1+X2>0)
MLE.Sorensen=2*length(I)/(sum(X1>0)+sum(X2>0))

MLE.Lennon=length(I)/(min(sum(X1>0),sum(X2>0)))
MLE.Bray_Curtis=sum(sapply(I,function(I) 2*min(X1[I],X2[I])))/(n1+n2) 
Morisita_Horn=2*sum(X1[I]/n1*X2[I]/n2)/sum((X1/n1)^2+(X2/n2)^2)
Morisita_Original=2*sum(X1[I]/n1*X2[I]/n2)/sum( X1*(X1-1)/n1/(n1-1) +X2*(X2-1)/n2/(n2-1))  

U_tilde=sum(X1[I])/n1;V_tilde=sum(X2[I])/n2
fplus1=sum(X1>0 & X2==1);fplus2=sum(X1>0 & X2==2)
fplus2=ifelse(fplus2>0,fplus2,1)
f1plus=sum(X1==1 & X2>0);f2plus=sum(X1==2 & X2>0)
f2plus=ifelse(f2plus>0,f2plus,1)
U_hat=U_tilde+(z-1)/z*fplus1/2/fplus2*sum(X1[I][X2[I]==1])/n1
U_hat=ifelse(U_hat<=1,U_hat,1)
V_hat=V_tilde+(w-1)/w*f1plus/2/f2plus*sum(X2[I][X1[I]==1])/n2
V_hat=ifelse(V_hat<=1,V_hat,1)
JAu=U_tilde*V_tilde/(U_tilde+V_tilde-U_tilde*V_tilde)
JAa=U_hat*V_hat/(U_hat+V_hat-U_hat*V_hat)
SAu=2*U_tilde*V_tilde/(U_tilde+V_tilde)
SAa=2*U_hat*V_hat/(U_hat+V_hat)

p <- Two_com_correct_obspi_Inc(X1=c(w,X1),X2=c(z,X2))
p1 <- p[,1];p2 <- p[,2]


boot.Jaccard=rep(0,boot)
boot.Esti.Jaccard=rep(0,boot)
boot.Sorensen=rep(0,boot)
boot.Esti.Sorensen=rep(0,boot)
boot.Lennon=rep(0,boot)
boot.Bray_Curtis=rep(0,boot)
boot.Morisita_Horn=rep(0,boot)
boot.Morisita_Original=rep(0,boot)


boot.U_hat=rep(0,boot)
boot.V_hat=rep(0,boot)
boot.JAu=rep(0,boot)
boot.JAa=rep(0,boot)
boot.SAu=rep(0,boot)
boot.SAa=rep(0,boot)
for(h in 1:boot)
{
  if(datatype == "abundance")
  {
    boot.X1=rmultinom(1,w,p1)
    boot.X2=rmultinom(1,z,p2)
    boot.n1=sum( boot.X1);boot.n2=sum(boot.X2)
    boot.shared.species.hat=PanEstFun(boot.X1,boot.X2)
    boot.alpha.species.hat=SpecAbunChao1(boot.X1, k=10, conf=0.95)[1,1]+SpecAbunChao1(boot.X2, k=10, conf=0.95)[1,1]
    boot.Esti.Jaccard[h]=boot.shared.species.hat/(boot.alpha.species.hat-boot.shared.species.hat)
    boot.Esti.Sorensen[h]=2*boot.shared.species.hat/boot.alpha.species.hat
  }
  if(datatype == "incidence")
  {
    boot.X1=sapply(1:length(p1), function(i)  rbinom(1, w, p1[i]))
    boot.X2=sapply(1:length(p2), function(i)  rbinom(1, z, p2[i]))
    boot.shared.species.hat=PanEstFun.Sam(c(w,boot.X1), c(z,boot.X2) )
    boot.alpha.species.hat=SpecInciChao2(c(w,boot.X1), k=10, conf=0.95)[1,1]+SpecInciChao2(c(z,boot.X2), k=10, conf=0.95)[1,1]
    boot.Esti.Jaccard[h]=boot.shared.species.hat/(boot.alpha.species.hat-boot.shared.species.hat)
    boot.Esti.Sorensen[h]=2*boot.shared.species.hat/boot.alpha.species.hat
    n1=sum(boot.X1)
    n2=sum(boot.X2)
  }
  I=which( boot.X1>0 &  boot.X2>0)
  
  boot.Jaccard[h]= length(I)/sum(boot.X1+boot.X2>0)
  boot.Sorensen[h]=2*length(I)/(sum(boot.X1>0)+sum(boot.X2>0))
  boot.Lennon[h]=length(I)/(min(sum(boot.X1>0),sum(boot.X2>0)))
  boot.Bray_Curtis[h]=sum(sapply(I,function(I) 2*min(boot.X1[I],boot.X2[I])))/sum(boot.X1+boot.X2)
  boot.Morisita_Horn[h]=2*sum(boot.X1[I]/n1*boot.X2[I]/n2)/sum((boot.X1/n1)^2+(boot.X2/n2)^2)
  boot.Morisita_Original[h]=2*sum(boot.X1[I]/n1*boot.X2[I]/n2)/sum( boot.X1*(boot.X1-1)/n1/(n1-1) +boot.X2*(boot.X2-1)/n2/(n2-1))
  
  boot.U_tilde=sum( boot.X1[I])/n1
  boot.V_tilde=sum( boot.X2[I])/n2
  fplus1=sum(boot.X1>0 & boot.X2==1);fplus2=sum(boot.X1>0 & boot.X2==2)
  fplus2=ifelse(fplus2>0,fplus2,1)
  f1plus=sum(boot.X1==1 & boot.X2>0);f2plus=sum(boot.X1==2 & boot.X2>0)
  f2plus=ifelse(f2plus>0,f2plus,1)
  boot.U_hat[h]=boot.U_tilde+(z-1)/z*fplus1/2/fplus2*sum(boot.X1[I][boot.X2[I]==1])/n1
  boot.U_hat[h]=ifelse(boot.U_hat[h]<=1,boot.U_hat[h],1)
  boot.V_hat[h]=boot.V_tilde+(w-1)/w*f1plus/2/f2plus*sum(boot.X2[I][boot.X1[I]==1])/n2
  boot.V_hat[h]=ifelse(boot.V_hat[h]<=1,boot.V_hat[h],1)
  boot.JAu[h]=boot.U_tilde*boot.V_tilde/(boot.U_tilde+boot.V_tilde-boot.U_tilde*boot.V_tilde)
  boot.JAa[h]=boot.U_hat[h]*boot.V_hat[h]/(boot.U_hat[h]+boot.V_hat[h]-boot.U_hat[h]*boot.V_hat[h])
  boot.SAu[h]=2*boot.U_tilde*boot.V_tilde/(boot.U_tilde+boot.V_tilde)
  boot.SAa[h]=2*boot.U_hat[h]*boot.V_hat[h]/(boot.U_hat[h]+boot.V_hat[h])
}
a=matrix(0,12,6)
a[1,]=c(min(MLE.Jaccard,1),sd(boot.Jaccard),rep(0,4))
a[2,]=c(Esti.Jaccard,sd(boot.Esti.Jaccard),rep(0,4))
a[3,]=c(min(MLE.Sorensen,1),sd(boot.Sorensen),rep(0,4))
a[4,]=c(Esti.Sorensen,sd(boot.Esti.Sorensen),rep(0,4))
a[5,]=c(MLE.Lennon,sd(boot.Lennon),rep(0,4))
a[6,]=c(min(MLE.Bray_Curtis,1),sd(boot.Bray_Curtis),rep(0,4))
a[7,]=c(min(Morisita_Horn,1),sd(boot.Morisita_Horn),rep(0,4))
a[8,]=c(Morisita_Original,sd(boot.Morisita_Original),rep(0,4))
a[9,]=c(min(JAu,1),sd(boot.JAu),U_tilde,V_tilde,rep(0,2))
a[10,]=c(min(JAa,1),sd(boot.JAa),U_hat,sd(boot.U_hat),V_hat,sd(boot.V_hat))
a[11,]=c(min(SAu,1),sd(boot.SAu),U_tilde,V_tilde,rep(0,2))
a[12,]=c(min(SAa,1),sd(boot.SAa),U_hat,sd(boot.U_hat),V_hat,sd(boot.V_hat))
round(a,4)
}


###########

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



SimilarityPair=function(X, datatype = c("abundance","incidence_freq", "incidence_raw"), units,nboot=200)
{ 
  
  
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
