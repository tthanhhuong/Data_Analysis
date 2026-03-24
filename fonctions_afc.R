Codisjc<-function(X)
{
#       codage disjonctif complet de la matrice des variables qualitatives X
#
#          Entree
#
# X matrice des variables qualitatives
#
#          Sorties
#
# nbmodal vecteur des nombres de modalites de chaque variable
# nbmodaltot nombre total de modalites
# U matrice du codage
# B matrice de Burt associee a U
#
#          Fonctions exterieures appelees
#    Aucune
#
#                    J.F. Durand et R. Sabatier
#                         MAJ 03/04/97
#
   X<-as.matrix(X)
   n<-nrow(X)
   p<-ncol(X)
   U<-NULL
   nbmodal<-NULL
   nommod<-NULL
   ylab<-NULL
   nbmodaltot<-as.numeric(0)
   for(j in 1:p)
                {
                 nbmodal[j]<-max(X[,j])
                 nommod<-paste(dimnames(X)[[2]][j],(1:nbmodal[j]),sep="")
                 nbmodaltot<-nbmodaltot+nbmodal[j]
                 V<-matrix(0,n,nbmodal[j])
                 for(i in 1:n) V[i,X[i,j]]<-1
                 U<-cbind(U,V)
                 ylab<-c(ylab,nommod)
                }
   B<-t(U)%*%U
   dimnames(U)<-list(dimnames(X)[[1]],ylab)
   dimnames(B)<-list(ylab,ylab)
   return(list(nbmodal=nbmodal,nbmodaltot=nbmodaltot,U=U,B=B))
}

acpxqd<-function(X,Q=1,D=1,centrer=T,cor=T,k=3,impres=T,graph=T,Xl=1,Xc=1,aideus=1,aideva=1)
{
#
#                DSD d'un triplet (X,Q,D)
#
#    Entrees
#
# X matrice des variables
# Q metrique des u.s. : si Q=1 alors Q=Idn, si Q est un vecteur alors Q =diag(vecteur), sinon Q=matrice
# D metrique des variables, si D=1 alors 1/nIdn, si D est un vecteur alors D=diag(vecteur), sinon D=matrice 
# centrer=F on ne centre pas X, sinon centrer=T
# cor=F on ne reduit pas, centrage et reduction si cor=T
# k nombre de composantes retenues, par defaut 3
# impres impression des resultats si T, si F rien
# graph trace des graphiques si T, si F pas de trace
# Xl matrice eventuelle des u.s. supplememtaires
# Xc matrice eventuelle des variables supplememtaires
# aideus=1 aides a l'interpretation pour les u.s., CTR et COS, pour les k premiers axes, sinon aideus=0
# aideva=1 aides a l'interpretation pour les variables, CTR et COS, pour les k premiers axes, sinon aideva=0
#
#    Sorties
#
# nomfichX nom du fichier X
# Xini tableau initial
# X tableau apres eventuel centrage et/ou reduction
# Q metrique utilisee pour les us
# D metrique utilisee pour les variables
# inertiaX inertie du triplet
# k rang de la DSD retenue
# valp vecteur de toutes les valeurs propres
# Fa matrice des k premiers facteurs principaux de la DSD
# A matrice des k premiers axes principaux de la DSD
# C matrice des k premieres composantes principales de la DSD
# CRTus et COSus contributions absolues et relatives des u.s. pour les k premieres composantes
# CTRva et COSva contributions absolues et relatives des variables pour les k premiers axes
# As matrice des k premiers axes pour les eventuelles variables supplementaires
# Cs matrice des k premieres composantes pour les eventuelles u.s. supplementaires
#
#    Fonctions exterieures appelees
#
#  Dcentred
#
#                               R.Sabatier
#                              MAJ 17/09/2000
# lecture des donnees
   nomfichX<-deparse(substitute(X))
   if(!is.matrix(X))stop("X n'est pas une matrice")
   X<-as.matrix(X)
   n<-nrow(X)
   p<-ncol(X)
  
if(is.null(dimnames(X)))dimnames(X)<-list(paste("i",1:n,sep=""),paste("v",1:p,sep=""))
if(length(dimnames(X)[[1]])==0)dimnames(X)[[1]]<-paste("i",1:n,sep="")
if(length(dimnames(X)[[2]])==0)dimnames(X)[[2]]<-paste("v",1:p,sep="")
   Xini<-as.matrix(X)
   if(p<2)return()
   if(k>p)return()
   if(impres==F)graph<-F
   if(centrer==F)cor<-F
   if(cor==T)centrer<-T
    p2<-p*p
    n2<-n*n
    As<-NULL
    Cs<-NULL
    CTRus<-NULL
    COSus<-NULL
    CTRva<-NULL
    COSva<-NULL
# calcul de la metrique Q
   if(length(Q)==1) Q<-diag(p)
   if(length(Q)==p) Q<-diag(Q,nrow=p)
   if(length(Q)==p2) Q<-as.matrix(Q)
# calcul de la metrique D
    if(length(D)==1) D<-diag(rep(1/n,n),nrow=n)
    if(length(D)==n) D<-diag(D,nrow=n)
    if(length(D)==n2) D<-as.matrix(D)
# verifications u.s. et variables supplementaires
     if(length(Xl)!=1) 
                      {
                       isup<-1
                       Xl<-as.matrix(Xl)
                       nns<-nrow(Xl)
                       if(ncol(Xl)!=ncol(X))return()
if(is.null(dimnames(Xl)))dimnames(Xl)<-list(paste("is",1:nns,sep=""),paste("v",1:p,sep=""))
if(length(dimnames(Xl)[[1]])==0)dimnames(Xl)[[1]]<-paste("is",1:nns,sep="")
if(length(dimnames(Xl)[[2]])==0)dimnames(Xl)[[2]]<-paste("v",1:p,sep="")
                       }
                       else isup<-0
     if(length(Xc)!=1) 
                        {
                        vsup<-1
                        Xc<-as.matrix(Xc)
                        pps<-ncol(Xc)
                        if(nrow(Xc)!=nrow(X))return()
if(is.null(dimnames(Xc)))dimnames(Xc)<-list(paste("i",1:n,sep=""),paste("vs",1:pps,sep=""))
if(length(dimnames(Xc)[[1]])==0)dimnames(Xc)[[1]]<-paste("i",1:n,sep="")
if(length(dimnames(Xc)[[2]])==0)dimnames(Xc)[[2]]<-paste("vs",1:pps,sep="")
                        }
                        else vsup<-0
# centrage et/ou reduction
   if(centrer==T)
                {
                 centrage<-Dcentred(X,D=D)
                 X<-centrage$Xc
                 meanX<-matrix(centrage$moy)
                 varX<-matrix(centrage$var)
                 dimnames(meanX)<-list(dimnames(X)[[2]],"moy")
                 dimnames(varX)<-list(dimnames(X)[[2]],"var")
                 if(isup==1)Xl<-sweep(Xl,2,meanX)
                 if(vsup==1)
                            {
                             centragevs<-Dcentred(Xc,D=D)
                             Xc<-centragevs$Xc
                            }           
                 }
   if(cor==T)
             {
              X<-centrage$Xcr
              if(isup==1)Xl<-sweep(Xl,2,sqrt(varX),FUN="/")
              if(vsup==1)Xc<-centragevs$Xcr
             }
   if(impres)
             {
              cat("            - D.S.D. d'un triplet (X,Q,D) -\n")
              cat("            -------------------------------\n")
              if(isup==1)
              cat("               avec u.s. supplementaires\n")
              if(vsup==1)
              cat("               avec variables supplementaires\n")
              cat("__________________________________________________________________\n")
              if(cor==T)
                        {
                         cat("                donnees centrees et reduites\n")
                             cat("                moyenne des variables de X\n")
                             print(t(meanX))
                             cat("                variance des variables de X\n")
                             print(t(varX))
                            }
              if(centrer==T && cor==F)
                            {
                             cat("                donnees centrees\n")
                             cat("                moyenne des variables de X\n")
                             print(t(meanX))
                             cat("                variance des variables de X\n")
                             print(t(varX))
                            }
             cat("__________________________________________________________________\n")
              }
# initialisation
   valp<-NULL                      # vecteur des valeurs propres
   A<-NULL                         # matrice des facteurs principaux 
   C<-NULL                         # matrice des composantes principales
# decomposition de Choleski de Q
   Qhalf<-chol(round(Q,6))
#   Qhalf<-Qhalf[,order(attr(Qhalf,"pivot"))]
   Y<-X%*%t(Qhalf)
# diagonalisation
   diago<-t(Y)%*%D%*%Y
   eg<-eigen(diago,symmetric=T)
   valp<-eg$values
   for(i in 1:ncol(diago))if(valp[i]<0)valp[i]<-0
   inertiaX<-sum(valp)
   valp2<-sum(valp^2)
# impression de l'histogramme des valeurs propres
   if(impres)
             {
              cat(paste("Inertie totale de (X,Q,D) =",format(inertiaX),"\n"))   
              cat("histogramme des valeurs propres (o ou n) ?\n")
              plth<-scan("",character(),1)
              if(length(plth)==0)plth<-"n"
              if(plth=="o" || plth=="O")
                                        {
                                         par(mfrow=c(1,1))
                                          barplot(valp,space=2,names=paste("v. p.",1:p))
                                         title("valeurs propres")
                                         }
#  choix du nombre d'axes si non donne en argument              
              cat("__________________________________________________________________\n")
              valtab<-matrix(0,p,4)
              dimnames(valtab)<-list(format(1:p),c("val.pro.","% inert.","% cumul."," RV"))
              valtab[,1]<-round(valp,digits=4)
              valtab[,2]<-round(valp/inertiaX*100,digits=2)
              for(i in 1:p)
                           {
                            valtab[i,3]<-sum(valtab[1:i,2])
                            valtab[i,4]<-sqrt(sum(valp[1:i]^2)/valp2)
                            }
              print(valtab)
              if(k==0)
                      {
                       repeat
                             {
                              cat("Combien d'axes voulez-vous ? (<=",p,")\n")
                              k<-scan("",n=1)
                              if(k!=0)break
                             }
                      }
              } 
# calcul des composantes principales
   Cent<-Y%*%eg$vectors
   C<-Cent[,1:k]
   if(k==1)C<-as.matrix(C)
   dimnames(C)<-list(dimnames(X)[[1]],paste("c",format(1:k),sep=""))
# calcul des axes principaux
   Aent<-matrix(NA,nrow=ncol(X),ncol=p)
   if(length(Q)!=1)fac<-solve(Qhalf)%*%eg$vectors
                   else fac<-eg$vectors
   for(i in 1:p) Aent[,i]<-(fac[,i])*sqrt(valp[i])
   A<-Aent[,1:k]
   dimnames(A)<-list(dimnames(X)[[2]],paste("a",format(1:k),sep=""))
# calcul des facteurs principaux
   Fa<-matrix(NA,nrow=ncol(X),ncol=k)
   Fa<-Q%*%A
   dimnames(Fa)<-list(dimnames(X)[[2]],paste("f",format(1:k),sep=""))
# calcul des coordonnees u.s. et variables supplementaires
   if(isup==1)
              {
               Cs<-Xl%*%Q%*%A
               for(i in 1:k) Cs[,i]<-(Cs[,i])/sqrt(valp[i])
               dimnames(Cs)<-list(dimnames(Xl)[[1]],paste("c",format(1:k),sep=""))
              }
   if(vsup==1)
              {
               As<-t(Xc)%*%D%*%C
               for(i in 1:k) As[,i]<-(As[,i])/sqrt(valp[i])
               dimnames(As)<-list(dimnames(Xc)[[2]],paste("a",format(1:k),sep=""))
              }
# calcul et impression des aides a l'interpretation
    if(aideus==1)
                   {
                        CTRus<-matrix(NA,nrow=n,ncol=k)
                        COSus<-CTRus
                        for(i in 1:n)for(j in 1:k)CTRus[i,j]<-round(D[i,i]*C[i,j]^2/valp[j]*10000,digits=0)
                        dimnames(CTRus)<-list(dimnames(X)[[1]],paste("CTR",format(1:k),sep=""))
                        for(i in 1:n)
                                 {
                                  ss<-0
                                  for(j in 1:p)ss<-ss+Cent[i,j]^2
                                  for(j in 1:k)COSus[i,j]<-round(Cent[i,j]^2/ss*10000,digits=0)
                                  }
                         dimnames(COSus)<-list(dimnames(X)[[1]],paste("COS",format(1:k),sep=""))
                    }
    if(aideva==1)
                    {  
                        CTRva<-matrix(NA,nrow=p,ncol=k)
                        COSva<-CTRva
                        for(i in 1:p)for(j in 1:k)CTRva[i,j]<-round(Q[i,i]*A[i,j]^2/valp[j]*10000,digits=0)
                        dimnames(CTRva)<-list(dimnames(X)[[2]],paste("CTR",format(1:k),sep=""))
                        for(i in 1:p)
                                 {
                                  ss<-0
                                  for(j in 1:p)ss<-ss+Aent[i,j]^2
                                  for(j in 1:k)COSva[i,j]<-round(Aent[i,j]^2/ss*10000,digits=0)
                                 }
                        dimnames(COSva)<-list(dimnames(X)[[2]],paste("COS",format(1:k),sep=""))
                     }
   if(impres)
           {
             if(aideus==1)
                {
                 cat("_______________________________________________________________\n")
                 cat("aides a l'interpretation pour les u.s. (o/n) ?\n")
                 plta<-scan("",character(),1)
                 if(length(plta)==0)plta<-"n"
                 if(plta=="o" || plth=="O")
                               {
                                cat("Contributions absolues des",n,"u.s. pour les",k,"premieres composantes\n")
                                print(CTRus)
                                cat("Contributions relative des",n,"u.s. pour les",k,"premieres composantes\n")
                                print(COSus)
                                }
                }
              if(aideva==1)
                {
                 cat("_______________________________________________________________\n")
                 cat("aides a l'interpretation pour les variables (o/n) ?\n")
                 plta<-scan("",character(),1)
                 if(length(plta)==0)plta<-"n"
                 if(plta=="o" || plth=="O")
                               {
                                cat("Contributions absolues des",p,"variables pour les",k,"premiers axes\n")
                                print(CTRva)
                                cat("Contributions relative des",p,"variables pour les",k,"premiers axes\n")
                                print(COSva)
                                }
                }
            }
# trace des eventuels graphiques
   if(graph==T)
            {
             cat("_______________________________________________________________\n")
             repeat
                   {
                   cat("graphique pour les u.s. (o/n) ?\n")
                   pltc<-scan("",character(),1)
                   if((length(pltc)==0)|(pltc=="n"))break
                                       else
                                           {
                                            cat("axe horizontal (<=",k,") ?\n")
                                            pltch<-scan("",numeric(),1)
                                            cat("axe vertical (<=",k,") ?\n")
                                            pltcv<-scan("",numeric(),1)
                                            par(mfrow=c(1,1),pty="s")
                                            axespar<-c(pltch,pltcv)
                                            if(isup==1)Ctot<-rbind(C,Cs) else Ctot<-C
                                            pltx<-Ctot[,pltch]
                                            plty<-Ctot[,pltcv]
                                            plot(x=pltx,y=plty,xlab=paste("c",axespar[1]," ",round(valp[axespar[1]],
                                                  digits=4),"(",round(valp[axespar[1]]/inertiaX*100,digits=2),"%)"),
                                                  ylab=paste("c",axespar[2]," ",round(valp[axespar[2]],
                                                  digits=4),"(",round(valp[axespar[2]]/inertiaX*100,digits=2),"%)"),
                                                  type="n")
                                            abline(h=0)
                                            abline(v=0)
                                            text(x=pltx,y=plty,dimnames(Ctot)[[1]])
                                           }
                   }
            cat("_______________________________________________________________\n")
             repeat
                   {
                   if(cor==F)
                            {
                   cat("graphique pour les variables (o/n) ?\n")
                   plta<-scan("",character(),1)
                   if((length(plta)==0)|(plta=="n"))break
                                       else
                                           {
                                            cat("axe horizontal (<=",k,") ?\n")
                                            pltah<-scan("",numeric(),1)
                                            cat("axe vertical (<=",k,") ?\n")
                                            pltav<-scan("",numeric(),1)
                                            par(mfrow=c(1,1),pty="s")
                                            axespar<-c(pltah,pltav)
                                            if(vsup==1)Atot<-rbind(A,As) else Atot<-A
                                            pltx<-Atot[,pltah]
                                            plty<-Atot[,pltav]
                                            plot(x=pltx,y=plty,xlab=paste("a",axespar[1]," ",round(valp[axespar[1]],
                                                  digits=4),"(",round(valp[axespar[1]]/inertiaX*100,digits=2),"%)"),
                                                  ylab=paste("a",axespar[2]," ",round(valp[axespar[2]],
                                                  digits=4),"(",round(valp[axespar[2]]/inertiaX*100,digits=2),"%)"),
                                                  type="n")
                                            text(x=pltx,y=plty,dimnames(Atot)[[1]])
                                          }
                              }
                        else
                              {
                   cat("cercle des correlations (o/n) ?\n")
                   plta<-scan("",character(),1)
                   if((length(plta)==0)|(plta=="n"))break
                                       else
                                           {
                                            cat("axe horizontal (<=",k,") ?\n")
                                            pltah<-scan("",numeric(),1)
                                            cat("axe vertical (<=",k,") ?\n")
                                            pltav<-scan("",numeric(),1)
                                            par(mfrow=c(1,1),pty="s")
                                            axespar<-c(pltah,pltav)
                                            theta<-seq(0,20,.05)
                                            x<-cos(theta)
                                            y<-sin(theta)
                                            if(vsup==1)Atot<-rbind(A,As) else Atot<-A
                                            pltx<-Atot[,pltah]
                                            plty<-Atot[,pltav]
                                            plot(x,y,type="l",xlab=paste("a",axespar[1]," ",round(valp[axespar[1]],
                                                  digits=4),"(",round(valp[axespar[1]]/inertiaX*100,digits=2),"%)"),
                                                  ylab=paste("a",axespar[2]," ",round(valp[axespar[2]],
                                                  digits=4),"(",round(valp[axespar[2]]/inertiaX*100,digits=2),"%)")
)
                                            abline(h=0)
                                            abline(v=0)
                                            text(x=pltx,y=plty,dimnames(Atot)[[1]])
                                          }  
                              }                
                   }
# fin des graphiques
            }
   return(list(nomfichX=nomfichX,Xini=Xini,X=X,Q=Q,D=D,inertiaX=inertiaX,k=k,valp=valp,Fa=Fa,A=A,C=C,CTRus=CTRus,COSus=COSus,CTRva=CTRva,COSva=COSva,As=As,Cs=Cs))
}

afc<-function(N,DI=1,DJ=1,M=1,k=3,ns=F,impres=T,graph=T)
{
#
#       AFC et AFC generalisee par rapport a un modele avec metriques lignes 
#                       et colonnes quelconques
#
#    Entrees
#
# N tableau de contingence
# DI metrique des lignes, si absent c'est la marge
# DJ metrique des colonnes, si absent c'est la marge
# M modele (matrice meme dim que N) si c'est une AFC par rapport a un modele, 
#                                   sinon M=1
# k nombre de composantes retenues, par defaut 3
# ns=T si AFC non symetrique, si ns=F AFC normale
# impres impression des resultats si T, si F rien
# graph trace des graphiques si T, si F pas de trace
#
#    Sorties
#
# nomfichX nom du fichier X
# Nini tableau initial
# N tableau apres eventuelle transformation
# DI metrique ligne
# DJ metrique colonne
# inertia inertie du triplet
# k rang de la DSD retenue
# valp vecteur de toutes les valeurs propres
# A matrice des k premiers axes de la DSD
# C matrice des k premieres composantes principales de la DSD
# CRTus et COSus contributions absolues et relatives des u.s. pour les k premieres composantes
# CTRva et COSva contributions absolues et relatives des variables pour les k premiers axes
#
#    Fonctions exterieures appelees
#
#       acpxqd, Dcentred
#
#                               R.Sabatier
#                              MAJ 03/04/97
#
# lecture des donnees
   nomfichN<-deparse(substitute(N))
   N<-as.matrix(N)
   M<-as.matrix(M)
   Nini<-N
   I<-nrow(N)
   J<-ncol(N)
   IJ<-I*J
   if(I<=2)return()
   if(J<=2)return()
   if(length(M)!=1 && ns==T)stop("conflits de parametres entre M et ns")
   if(impres==F)graph<-F
   n<-sum(N)
   P<-N/n
   I1<-matrix(1,nrow=I,ncol=1)
   J1<-matrix(1,nrow=J,ncol=1)
# metriques si AFC non symetrique
   if(ns==T)
            {
             DJ<-diag(J)
             DI<-diag(as.vector(P%*%J1),nrow=I)
            }
    else
            {
#          calcul de la metrique DI
            if(length(DI)==1) DI<-diag(as.vector(P%*%J1),nrow=I)
            if(length(DI)==I) DI<-diag(DI,nrow=I)
#          calcul de la metrique DJ
            if(length(DJ)==1) DJ<-diag(as.vector(t(I1)%*%P),nrow=J)
            if(length(DJ)==J) DJ<-diag(DJ,nrow=J)
            }
   if(impres)
             { 
              cat("_______________________________________________________________\n")
              if(ns==F)cat("            - AFC de ",nomfichN," -\n")
              if(ns==T)cat("            - AFC non symetrique de ",nomfichN," -\n")
              if(length(M)==IJ)cat("              par rapport a un modele\n")
             }
# calcul de X, mis dans N
   if(length(M)==1)
                   {
                    if(ns==F) for(i in 1:I)for(j in 1:J)
                              N[i,j]<-P[i,j]/(DJ[j,j]*DI[i,i])-1
                    if(ns==T) for(i in 1:I)for(j in 1:J)
                      N[i,j]<-(P[i,j]/DI[i,i])-diag(t(I1)%*%P,nrow=J)[j,j]
                   }
   if(length(M)==IJ)for(i in 1:I)for(j in 1:J)
            N[i,j]<-(P[i,j]-M[i,j]/n)/(DJ[j,j]*DI[i,i])
   res<-acpxqd(N,Q=DJ,D=DI,centrer=F,k=k,impres=F,graph=F)
   inertia<-res$inertiaX
   chi2<-n*inertia
   ddl<-(I-1)*(J-1)
   valp<-res$valp
   C<-res$C
   k<-res$k
   A<-res$A
   CTRus<-res$CTRus
   COSus<-res$COSus
   CTRva<-res$CTRva
   COSva<-res$COSva
   if(impres==T)
             {
              cat("_______________________________________________________________\n")
              cat(paste("Inertie totale du triplet    =",format(inertia),"\n"))
              cat(paste("nombre de modalites lignes   =",I,"\n"))
              cat(paste("nombre de modalites colonnes =",J,"\n"))
              cat(paste("effectif du tableau          =",n,"\n"))
              if(ns==F)cat(paste("Chi-deux d'independance      =",format(chi2)," avec ",ddl,
                            " d.d.l. , p =",format(pchisq(chi2,ddl))),"\n")
              cat("_______________________________________________________________\n")
              cat("histogramme des valeurs propres (o ou n) ?\n")
              plth<-scan("",character(),1)
              if(length(plth)==0)plth<-"n"
              if(plth=="o" || plth=="O")
                                        {              
                                         pp<-length(valp)
                                         valtab<-matrix(0,pp,3)
                                         dimnames(valtab)<-list(format(1:pp),c("val.pro.","% inert.","% cumul."))
                                         valtab[,1]<-round(valp,digits=5)
                                         valtab[,2]<-round(valp/inertia*100,digits=2)
                                         for(i in 1:pp)valtab[i,3]<-sum(valtab[1:i,2])
                                         print(valtab)
                                         par(mfrow=c(1,1))
                                         barplot(valp,space=2,names=paste("v. p.",1:J))
                                         title("valeurs propres")
                                         }
 
              }
# calcul et impression des aides a l'interpretation
   if((impres==T)|length(D)==I)
                {
                 cat("_______________________________________________________________\n")
                 cat("aides a l'interpretation pour les modalites lignes (o/n) ?\n")
                 plta<-scan("",character(),1)
                 if(length(plta)==0)plta<-"n"
                 if(plta=="o" || plth=="O")
                               {
                                cat("Contributions absolues des",I,"modalites pour les",k,"premieres composantes\n")
                                print(CTRus)
                                cat("\nContributions relative des",I,"modalites pour les",k,"premieres composantes\n")
                                print(COSus)
                                }
                }
   if((impres==T)|length(D)==J)
                {
                 cat("_______________________________________________________________\n")
                 cat("aides a l'interpretation pour les modalites colonnes (o/n) ?\n")
                 plta<-scan("",character(),1)
                 if(length(plta)==0)plta<-"n"
                 if(plta=="o" || plth=="O")
                               {
                                cat("Contributions absolues des",J,"modalites pour les",k,"premiers axes\n")
                                print(CTRva)
                                cat("\nContributions relative des",J,"modalites pour les",k,"premiers axes\n")
                                print(COSva)
                                }
                }
# trace des eventuels graphiques
   if(graph==T)
            {
             cat("_______________________________________________________________\n")
             repeat
                   {
                   cat("graphique pour les composantes (lignes) (o/n) ?\n")
                   pltc<-scan("",character(),1)
                   if((length(pltc)==0)|(pltc=="n"))break
                                       else
                                           {
                                            cat("axe horizontal (<=",k,") ?\n")
                                            pltch<-scan("",numeric(),1)
                                            cat("axe vertical (<=",k,") ?\n")
                                            pltcv<-scan("",numeric(),1)
                                            par(mfrow=c(1,1),pty="s")
                                            axespar<-c(pltch,pltcv)
                                            pltx<-C[,pltch]
                                            plty<-C[,pltcv]
                                            plot(x=pltx,y=plty,xlab=paste("c",axespar[1]," ",round(valp[axespar[1]],
                                                  digits=4),"(",round(valp[axespar[1]]/inertia*100,digits=2),"%)"),
                                                  ylab=paste("c",axespar[2]," ",round(valp[axespar[2]],
                                                  digits=4),"(",round(valp[axespar[2]]/inertia*100,digits=2),"%)"),
                                                  type="n")
                                            abline(h=0)
                                            abline(v=0)
                                            text(x=pltx,y=plty,dimnames(C)[[1]])
                                           }
                   }
             cat("_______________________________________________________________\n")
             repeat
                   {
                   cat("graphique pour les axes (colonnes) (o/n) ?\n")
                   plta<-scan("",character(),1)
                   if((length(plta)==0)|(plta=="n"))break
                                       else
                                           {
                                            cat("axe horizontal (<=",k,") ?\n")
                                            pltah<-scan("",numeric(),1)
                                            cat("axe vertical (<=",k,") ?\n")
                                            pltav<-scan("",numeric(),1)
                                            par(mfrow=c(1,1),pty="s")
                                            axespar<-c(pltah,pltav)
                                            pltx<-A[,pltah]
                                            plty<-A[,pltav]
                                            plot(x=pltx,y=plty,xlab=paste("a",axespar[1]," ",round(valp[axespar[1]],
                                                  digits=4),"(",round(valp[axespar[1]]/inertia*100,digits=2),"%)"),
                                                  ylab=paste("a",axespar[2]," ",round(valp[axespar[2]],
                                                  digits=4),"(",round(valp[axespar[2]]/inertia*100,digits=2),"%)"),
                                                  type="n")
                                            abline(h=0)
                                            abline(v=0)
                                            text(x=pltx,y=plty,dimnames(A)[[1]])
                                          }
                    }
             cat("_______________________________________________________________\n")
             repeat
                   {
                   cat("representations simultanees (o/n) ?\n")
                   plta<-scan("",character(),1)
                   if((length(plta)==0)|(plta=="n"))break
                                       else
                                           {
                                            cat("axe horizontal (<=",k,") ?\n")
                                            pltah<-scan("",numeric(),1)
                                            cat("axe vertical (<=",k,") ?\n")
                                            pltav<-scan("",numeric(),1)
                                            par(mfrow=c(1,1),pty="s")
                                            axespar<-c(pltah,pltav)
                                            Z<-matrix(0,nrow=I+J,ncol=2)
                                            Z[,1]<-c(C[,axespar[1]],A[,axespar[1]])
                                            Z[,2]<-c(C[,axespar[2]],A[,axespar[2]])
                                            dimnames(Z)<-list(c(dimnames(C)[[1]],dimnames(A)[[1]]),c("1","2"))
                                            plot(Z,xlab=paste(axespar[1]," ",round(valp[axespar[1]],
                                                  digits=4),"(",round(valp[axespar[1]]/inertia*100,digits=2),"%)"),
                                                  ylab=paste(axespar[2]," ",round(valp[axespar[2]],
                                                  digits=4),"(",round(valp[axespar[2]]/inertia*100,digits=2),"%)"),
                                                  type="n")
                                            abline(h=0)
                                            abline(v=0)
                                            text(Z,dimnames(Z)[[1]])
                                          }
                    }
       if(ns==F)
            {
            cat("_______________________________________________________________\n")
             repeat
                   {
                   cat("representations barycentriques (o/n) ?\n")
                   plta<-scan("",character(),1)
                   if((length(plta)==0)|(plta=="n"))break
                                       else
                                           {
                                            cat("lignes au barycentre des colonnes (1) ou l'inverse (2) ?\n")
                                            bar<-scan("",numeric(),1)
                                            cat("axe horizontal (<=",k,") ?\n")
                                            pltah<-scan("",numeric(),1)
                                            cat("axe vertical (<=",k,") ?\n")
                                            pltav<-scan("",numeric(),1)
                                            par(mfrow=c(1,1),pty="s")
                                            axespar<-c(pltah,pltav)
                                            if(bar==1)
                                                      {
                                                       Z<-matrix(0,nrow=I+J,ncol=2)
                                                       Z[,1]<-c(C[,axespar[1]],A[,axespar[1]]/sqrt(valp[axespar[1]]))
                                                       Z[,2]<-c(C[,axespar[2]],A[,axespar[2]]/sqrt(valp[axespar[2]]))
                                                       dimnames(Z)<-list(c(dimnames(C)[[1]],dimnames(A)[[1]]),c("1","2"))
                                                       }
                                            if(bar==2)
                                                      {
                                                       Z<-matrix(0,nrow=I+J,ncol=2)
                                                       Z[,1]<-c(C[,axespar[1]]/sqrt(valp[axespar[1]]),A[,axespar[1]])
                                                       Z[,2]<-c(C[,axespar[2]]/sqrt(valp[axespar[2]]),A[,axespar[2]])
                                                       dimnames(Z)<-list(c(dimnames(C)[[1]],dimnames(A)[[1]]),c("1","2"))
                                                       }
                                            plot(Z,xlab=paste(axespar[1]," ",round(valp[axespar[1]],
                                                  digits=4),"(",round(valp[axespar[1]]/inertia*100,digits=2),"%)"),
                                                  ylab=paste(axespar[2]," ",round(valp[axespar[2]],
                                                  digits=4),"(",round(valp[axespar[2]]/inertia*100,digits=2),"%)"),
                                                  type="n")
                                            abline(h=0)
                                            abline(v=0)
                                            text(Z,dimnames(Z)[[1]])
                                          }
                    }
              }
# fin des graphiques
            }
    return(list(nomfichN=nomfichN,Nini=Nini,N=N,DI=DI,DJ=DJ,inertia=inertia,k=k,valp=valp,A=A,C=C,CTRus=CTRus,COSus=COSus,CTRva=CTRva,COSva=COSva))           
}

Dcentred<-function(X,D=1)
{
#       D-centrage de la matrice X
#
#          Entree
#
# X matrice des variables
# D metrique des poids
#
#          Sorties
#
# Xc matrice deduite de X, D-centree
# Xcr matrice deduite de X, D-centree et reduite
# moy vecteur des D-moyennes
# var vecteur des D-variances
#
#          Fonctions exterieures appelees
#       Aucune.
#
#                       R. Sabatier
#                      MAJ 03/04/97

   X<-as.matrix(X)
   n<-nrow(X)
   p<-ncol(X)
   n2<-n*n
   n1<-matrix(1,nrow=n,ncol=1)
   if(length(D)==1)D<-diag(rep(1/n,n),nrow=n)
   if(length(D)==n)D<-diag(D,nrow=n)
   if(length(D)==n2)D<-as.matrix(D)
   moy<-t(n1)%*%D%*%X
   Xc<-sweep(X,2,moy)
   var<-diag(t(Xc)%*%D%*%Xc)
   ect<-sqrt(var)
   Xcr<-sweep(Xc,2,ect,FUN="/")
  return(list(Xc=Xc,Xcr=Xcr,moy=moy,var=var))
}

star.graph<-function(XY,Gr,D=1,taille=1)
{
#
#          Graphique en etoile
#
#    Entrees
#
# XY matrice des coordonnees des points
# Gr matrice des modalites d'appartenance des points
# D metrique des variables, si D=1 alors 1/nIdn, si D est un vecteur alors D=diag(vecteur), sinon D=matrice 
# taille dimension en 'inches' des noms des classes sur les graphiques
#
#    Sorties
#       Aucune
#
#    Fonctions exterieures appelees
#       Aucune.
#
#                               R.Sabatier
#                               MAJ 03/04/97
#
   nomfichXY<-deparse(substitute(XY))
   XY<-as.matrix(XY)
   n<-nrow(XY)
   if(length(dimnames(XY)[[1]])==0)dimnames(XY)[[1]]<-paste("i",1:n,sep="")
   Gr<-as.matrix(Gr)
   q<-ncol(Gr)
   if(length(dimnames(Gr)[[2]])==0)dimnames(Gr)[[2]]<-paste("Var",1:q,sep="")
   n2<-n*n
   k<-ncol(XY)
   if(k<2)return()
   if(n!=nrow(Gr)) return()
# calcul de la metrique D
   if(length(D)==1) D<-diag(rep(1/n,n),nrow=n)
   if(length(D)==n) D<-diag(D,nrow=n)
   if(length(D)==n2) D<-as.matrix(D)
# verification du graphisme
   if(!exists(".Device",frame=0))
            {
              cat("Initialisez le graphique !!!\n")
              return()
            }
# trace des graphiques
            {
             cat("_______________________________________________________________\n")
             repeat
                   {
                   cat("graphique pour les u.s. (o/n) ?\n")
                   pltc<-scan("",character(),1)
                   if((length(pltc)==0)|(pltc=="n"))break
                   else
            {
# choix des numeros axes et de la variable qualitative                                            
                     cat("axe horizontal (<=",k,") ?\n")
                     pltch<-scan("",numeric(),1)
                     cat("axe vertical (<=",k,") ?\n")
                     pltcv<-scan("",numeric(),1)
                     cat("numero de la var. qual. (<=",q,") ?\n")
                     pltvqual<-scan("",numeric(),1)
# calcul des moyennes par classe
                     nbmodal<-max(Gr[,pltvqual])
                     Cs<-matrix(nrow=nbmodal,ncol=k,0)
                     if(nbmodal<=1) return()
                     poicla<-matrix(ncol=1,nrow=nbmodal,0)
                     for(i in 1:n)poicla[Gr[i,pltvqual]]<-poicla[Gr[i,pltvqual]]+D[i,i]
                     for(i in 1:n)
                          {
                           cl<-Gr[i,pltvqual]
                           poi<-D[i,i]
                           x<-XY[i,pltch]
                           y<-XY[i,pltcv]
                           Cs[cl,pltch]<-Cs[cl,pltch]+x*poi
                           Cs[cl,pltcv]<-Cs[cl,pltcv]+y*poi
                          }
                     for(i in 1:nbmodal)
                          {
                           no<-poicla[i]
                           Cs[i,pltch]<-Cs[i,pltch]/no
                           Cs[i,pltcv]<-Cs[i,pltcv]/no
                          }
                     dimnames(Cs)<-list(c(1:nbmodal),dimnames(XY)[[2]])
                     par(mfrow=c(1,1),pty="s")
                     axespar<-c(pltch,pltcv)
                     Ctot<-rbind(XY,Cs)
                     plot(Ctot[,axespar],xlab=paste("c",axespar[1]),
                     ylab=paste("c",axespar[2]),main=dimnames(Gr)[[2]][pltvqual],
                     type="p",pch=18)
                     for(i in 1:n)
                          {
                           cl<-Gr[i,pltvqual]
                           segments(Cs[cl,pltch],Cs[cl,pltcv],XY[i,pltch],XY[i,pltcv],col=2)
                          }
                     abline(h=0)
                     abline(v=0)
                     points(Cs[,axespar[1]],Cs[,axespar[2]],pch=16,mkh=(taille+0.09),col=0)
                     points(Cs[,axespar[1]],Cs[,axespar[2]],pch=1,mkh=(taille+0.09),col=1)
                     text(Cs[,axespar[1]],Cs[,axespar[2]],dimnames(Cs)[[1]],cex=taille,col=3)
            }
                   }
            }
}

etaACM<-function(X,sortie)
{
# calcul des rapports de correlation des variables qualitatives actives suite à une ACM pour les k axes retenus
#
#
#    Entrees
#
# X               matrice du codage (non disjonctif complet) dont on a réalisé l'ACM apres utilisation de Codisjc
# sortie nom de l'objet suite à l'ACM de U : sortie de la fonction afc
#
#
#    Sorties
#
# sortierap rapports de correlation de chaque variable qualitative pour les k dimensions retenues dand l'ACM 
#
#
#    Fonctions exterieures appelees
#
#       Codisjc
#
#                               R.Sabatier
#                              MAJ 20/11/2018
#


k=sortie$k
truc=Codisjc(X)
nbmodal=truc$nbmodal
nbvar=dim(X)[[2]]
n=dim(X)[[1]]

sortierap=matrix(0,nrow=nbvar,ncol=k)
colnames(sortierap)=c(paste("valp.",1:k,sep=""))
rownames(sortierap)=colnames(X)

# browser()

for(j in 1:k)
{
	sortierap[,j]=tapply(sortie$CTRva[,j],INDEX=rep(1:nbvar,nbmodal),FUN="sum")*nbvar*sortie$valp[j]/10000
}

valpACM=apply(sortierap,2,mean) # => valp de l'AFCM

# fin et return
return(list(k=k,sortierap=sortierap,valpACM=valpACM))           

}