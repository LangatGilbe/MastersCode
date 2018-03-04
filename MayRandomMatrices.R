library("plotrix")
library("Matrix")
library("magic")

###### Allesina an Tang 2012
### Random matrices construction
### Equivalent to Mays random matrix
May_rnd = function(S,C,sg,m){
  s = S*S; nc = round(C*S*S, digits = 0)
  mm = matrix(1,S,S); diag(mm)=0                                       
  while (s>(nc)){                                           
    a <- sample(S,1)                                        
    b <- sample(S,1)                                        
    if (mm[a,b]>0 ){      
      mm[a,b] <- 0                                                 
      s <- s-1                                          
    }
    
  }
  AB = matrix(rnorm((S*S),0,sg), nr=S)
  M = diag(-m,S,S)
  J = (mm*AB) + M
  return(J)
}


S=1000; C=0.25; sg =0.05; m= 1
J = May_rnd(S,C,sg,m)
E = eigen(J)$values
plot(Re(E), Im(E), pch=20, xlim = c(-5,5), ylim = c(-2.5,2.5))
draw.circle(-m,0,(sg*sqrt(C*S)))
abline(h=0, lty=1)
abline(v=0, lty=1)
abline(v=-m, lty=2)

############################################################################################################
### Allesina criteria competition matrix
Competition = function(S,C,sg,m){
  s = S*S; nc = round(C*S*S, digits = 0)
  mm = matrix(1,S,S); diag(mm)=0                                       
  while (s>nc){                                           
    a <- sample(S,1)                                        
    b <- sample(S,1)                                        
    if (mm[a,b]>0 && mm[b,a]>0 ){      
      mm[a,b] <- 0
      mm[b,a] <- 0
    }
    s <- s-1
  }

  AB = - abs(matrix(rnorm((S*S),0,sg), nr=S))
  M = diag(-m,S,S)
  J = (mm*AB) + M
  return(J)
}


S=100; C=0.25; sg =0.01; m= 1;pi=3.142
J = Competition(S,C,sg,m)
E = eigen(J)$values
E1= ((pi-2)/pi); E2= ((pi+2)/pi)
plot(Re(E), Im(E), pch=20, xlim = c(-1.5,0), ylim = c(-0.05,0.05))
#draw.circle(-m,0,(sg*sqrt(c*S))
draw.ellipse(-1,0,b=(sg*sqrt(S*C)*E1),a=(sg*sqrt(S*C)*E2))
abline(h=0, lty=1)
abline(v=0, lty=1)
abline(v=-m, lty=2)




###Random matrices plus dispersal
##Based on Gravel et.al paper 2016

GRAVEL_rnd = function(S,C,sg,m,sites,d){
  
  s = S*S*sites; nc = round(C*S*S*sites, digits = 0)
  #mmm = matrix(1,S,S); diag(mmm)<-0  
  mm=array(rbinom(S*S*sites,1,1),c(S,S,sites))
  
  for (i in 1:S){ for (j in 1: sites){ mm[i,i,j]=0  }}
  
  
  while (s>(nc)){
    for (site in 1:sites){
      
      a <- sample(S,1)                                        
      b <- sample(S,1)
      
      if (mm[a,b,site] >0 ){      
        mm[a,b,site] <- 0
        s <- s-1
      }
      
    }
    
  }
  AB = array(rnorm((S*S*sites),0,sg), c(S,S,sites))
  M = array(0, c(S,S,sites)) ; for (i in 1:S){ for (j in 1: sites){ M[i,i,j]=-1 }}
  Community_matrix = (mm*AB)+ M
  
  Intra_specific_Matrix =diag(nrow=(S*sites))
  
  Interpatch_dispersal = `diag<-`(matrix(0,S,S),d/(sites-1))
  
  A=matrix(rep(t(Interpatch_dispersal),sites),nrow = S*sites,byrow=T)
  Dispersal_matrix = matrix(A,nrow = S*sites,ncol = S*sites)
  
  diag(Dispersal_matrix)<--d
  
  site=1
  Patch=list()
  while(site< sites){
    Patch[site]=list(Community_matrix[,,site],Community_matrix[,,(site+1)])
    site=site+1
    Patch[[sites]]<-Community_matrix[,,sites]
  }
  
  A_diag=do.call(adiag,Patch)
  Global_Community=-Intra_specific_Matrix+Dispersal_matrix +A_diag
  
  
  return(Global_Community)
}


S=100;sites=10;d=8;sg=1;C=0.3;m=2
Global_Community = GRAVEL_rnd(S,C,sg,m,sites,d)
Global_Eigenval= eigen(Global_Community)$values
plot(Re(Global_Eigenval),Im(Global_Eigenval),pch=20,xlim=c(-30,10),ylim = c(-10,10))
abline(h=0,lty=1)
abline(v=0,lty=1)
abline(v=-1,lty=2)
draw.circle(-(m+(sites*d)/(sites-1)),0,sg*sqrt(C*(S-1)*(sites-1)/sites))
draw.circle(-m,0,sg*sqrt(C*(S-1)/sites))


### Competition matrix plus dispersal

Dispersal_compet<-function(S,sites,d,sigma,C,m){
  s = S*S*sites; nc = round(C*S*S*sites, digits = 0)
  Int_mat=array(rbinom(S*S*sites,1,1),c(S,S,sites))
  
  for (i in 1:S){ for (j in 1: sites){ Int_mat[i,i,j]=0  }}
  
  while (s>(nc)){
    for (site in 1:sites){
      
      a <- sample(S,1)                                        
      b <- sample(S,1)
      
      if (Int_mat[a,b,site] >0 ){      
        Int_mat[a,b,site] <- 0
        Int_mat[b,a,site]<-0
      }
      s <- s-1
    }
    
  }
  
  Competition_Strength = -abs(array(rnorm((S*S*sites),0,sigma), c(S,S,sites)))
  M = array(0, c(S,S,sites)) ; for (i in 1:S){ for (j in 1: sites){ M[i,i,j]=-1 }}
  Community_matrix = (Int_mat*Competition_Strength)+ M
  
  Intra_specific_Matrix =diag(nrow=(S*sites))
  
  Interpatch_dispersal = `diag<-`(matrix(0,S,S),d/(sites-1))
  
  A=matrix(rep(t(Interpatch_dispersal),sites),nrow = S*sites,byrow=T)
  Dispersal_matrix = matrix(A,nrow = S*sites,ncol = S*sites)
  
  diag(Dispersal_matrix)<--d
  
  site=1
  Patch=list()
  while(site< sites){
    Patch[site]=list(Community_matrix[,,site],Community_matrix[,,(site+1)])
    site=site+1
    Patch[[sites]]<-Community_matrix[,,sites]
  }
  
  A_diag=do.call(adiag,Patch)
  Competitive_dispersal_Community=-Intra_specific_Matrix+Dispersal_matrix +A_diag
  
  
  return(Competitive_dispersal_Community)
  
}

S=100; C=0.3; sigma =1; m= 1;pi=3.142;d=1;sites=10
Competitive_dispersal_Community = Dispersal_compet(S,sites,d,sigma,C,m)
CompetdispersalEigenvalues = eigen(Competitive_dispersal_Community)$values
E1= ((pi-2)/pi); E2= ((pi+2)/pi)
plot(Re(CompetdispersalEigenvalues), Im(CompetdispersalEigenvalues), pch=20, xlim = c(-30,30), ylim = c(-2.5,2.5))
#draw.circle(-m,0,(sg*sqrt(c*S))
draw.ellipse(-(-m+1+d),0,b=(sigma*sqrt(S*C)*E1),a=(sigma*sqrt(S*C)*E2))
abline(h=0, lty=1)
abline(v=0, lty=1)
abline(v=-m, lty=2)


