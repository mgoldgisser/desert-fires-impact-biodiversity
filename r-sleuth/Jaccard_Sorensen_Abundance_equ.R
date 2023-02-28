Two_com_correct_obspi_Inc=function(X1,X2)
{
  X1 = X[,1]
  X2 = X[,2]
  X1=c(w,X1)
  X2=c(z,X2)
  
  x1=X1; x2=X2
  t1=X1[1]; X1=X1[-1]
  t2=X2[1]; X2=X2[-1]
  Q11=sum(X1==1)
  Q12=sum(X1==2)
  Q21=sum(X2==1)
  Q22=sum(X2==2)
  C1=Chat.Ind_Inc(x1,t1)
  C2=Chat.Ind_Inc(x2,t2)
  
  PP1=correct_obspi_Inc(x1)
  PP2=correct_obspi_Inc(x2)
  D12=which(X1>0 & X2>0)
  
  Q0hat_1=ceiling( (t1-1)/t1* ifelse(Q12 == 0,  Q11 * (Q11 - 1) / 2,  Q11 ^ 2/ 2 / Q12)   )   
  Q0hat_2=ceiling( (t2-1)/t2* ifelse(Q22 == 0,  Q21 * (Q21 - 1) / 2,  Q21 ^ 2/ 2 / Q22)   )
  #-----------------------------------------------------------------------------
  
  r1=which(X1>0 & X2==0)
  Q.1=length(which(X1>0 & X2==1))
  Q.2=length(which(X1>0 & X2==2))
  Q.0=ceiling((t2-1)/t2* ifelse(Q.2>0,Q.1^2/2/Q.2,Q.1*(Q.1-1)/2))
  #------------------------------------------------------------------------------
  r2=which(X1==0 & X2>0)
  Q1.=length(which(X1==1 & X2>0))
  Q2.=length(which(X1==2 & X2>0))
  Q0.=ceiling((t1-1)/t1*ifelse(Q2.>0,Q1.^2/2/Q2.,Q1.*(Q1.-1)/2))
  #------------------------------------------------------------------------------
  t11=length(which(X1==1 & X2==1))
  t22=length(which(X1==2 & X2==2))
  Q00hat=ceiling((t1-1)/t1*(t2-1)/t2* ifelse(t22 == 0,  t11 * (t11 - 1) / 4,  t11 ^ 2/ 4 / t22)   )
  #------------------------------------------------------------------------------
  temp1=max(length(r1),Q.0)-length(r1)+Q0.+Q00hat
  temp2=max(length(r2),Q0.)-length(r2)+Q.0+Q00hat
  p1_us_sum=min(U2_equ_Inc(x1,x2),1-C1)
  p1_us=p1_us_sum/temp1
  p2_us_sum=min(U2_equ_Inc(x2,x1),1-C2)
  p2_us=p2_us_sum/temp2
  if(Q0hat_1-temp1>0){p0_1= (1-C1-p1_us_sum)/(Q0hat_1-temp1) }
  if(Q0hat_1-temp1<=0){p0_1=0} 
  if(Q0hat_2-temp2>0){p0_2= (1-C2-p2_us_sum)/(Q0hat_2-temp2) }
  if(Q0hat_2-temp2<=0){p0_2=0}
  #------------------------------------------------------------------------------
  P1=PP1[D12]
  P2=PP2[D12]
  if(length(r1)> Q.0)
  {
    P1=c(P1,PP1[r1])
    Y=c(rep(p2_us, Q.0), rep(0,length(r1)-Q.0))
    P2=c(P2,sample(Y,length(Y)) )
  }
  if(length(r1)< Q.0)
  {
    P1=c(P1,PP1[r1],rep( p1_us,Q.0-length(r1)))
    P2=c(P2,rep(p2_us, Q.0) )
  }
  #----------------------------------------------------------------------------   
  if(length(r2)> Q0.)
  {
    Y=c(rep(p1_us,Q0.),rep(0,length(r2)- Q0.))
    P1=c(P1,sample(Y,length(Y)))
    P2=c(P2,PP2[r2] )
  }
  if(length(r2)< Q0.)
  {
    P1=c(P1,rep(p1_us,Q0.))
    P2=c(P2,PP2[r2],rep( p2_us,Q0.-length(r2)) )
  }
  P1=c(P1,rep( p1_us,Q00hat))
  P2=c(P2,rep( p2_us,Q00hat))
  P1=c(P1, rep(   p0_1,max( Q0hat_1-temp1,0)) , rep(    0  ,max( Q0hat_2-temp2,0))         )
  P2=c(P2, rep(     0 ,max( Q0hat_1-temp1,0)) , rep(   p0_2,max( Q0hat_2-temp2,0))         ) 
  #------------------------------------------------------------------------------------
  a=cbind(P1,P2)
  return(a)
}

Jaccard_Sorensen_Abundance_equ=function(datatype = c("abundance", "incidence"),X1,X2,boot)
{
  if(datatype == "incidence")
  {
    shared.species.hat=PanEstFun.Sam(X1,X2)
    alpha.species.hat=SpecInciChao2(X1, k=10, conf=0.95)[1,1]+SpecInciChao2(X2, k=10, conf=0.95)[1,1]
    Esti.Jaccard=shared.species.hat/(alpha.species.hat-shared.species.hat)
    Esti.Sorensen=2*shared.species.hat/alpha.species.hat
    w=X1[1];z=X2[1]
    X1=X1[-1];X2=X2[-1]
  }
  n1=sum(X1);n2=sum(X2)
  if(datatype == "abundance"){w=n1;z=n2}
  I=which(X1>0 & X2>0)
  
  MLE.Jaccard=length(I)/sum(X1+X2>0)
  MLE.Sorensen=2*length(I)/(sum(X1>0)+sum(X2>0))
  #############################################
  if(datatype == "abundance")
  {
    shared.species.hat=PanEstFun(X1,X2)
    alpha.species.hat=SpecAbunChao1(X1, k=10, conf=0.95)[1,1]+SpecAbunChao1(X2, k=10, conf=0.95)[1,1]
    Esti.Jaccard=shared.species.hat/(alpha.species.hat-shared.species.hat)
    Esti.Sorensen=2*shared.species.hat/alpha.species.hat
  }
  #############################################
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
  
  if(datatype=="incidence"){
    p <- Two_com_correct_obspi_Inc(X1=c(w,X1),X2=c(z,X2))
    p1 <- p[,1];p2 <- p[,2]
  }else{
    p <- Two_com_correct_obspi(X1,X2)
    p1 <- p[,1] ; p2 <- p[,2]
  }
  
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
