#########################################################################################
## Implementattion of Modified Lotka-Volterra model on Euler method
#########################################################################################


##############################################
#Clearing memory and Saving code to folder
##############################################
 
rm(list=ls(all=TRUE)) 
setwd("E:/MSC/Code/Part 2 Code/Trial")
#set.seed(123)
##########################
### LOAD LIBRARIES
##########################
library(tidyverse)
library(BiodiversityR)
library(pracma)
library(vegan)
library(bipartite)
library(rnetcarto)
library(ggthemes)
library(reshape2)
library(magic)



#############################
## FUNCTION DEFINITIONS
#############################

###
# lvm(t,pop,parms)
# Use:	Function to calculate derivative of multispecies Lotka-Volterra equations
# Input:
# 	t: time (not used here, because there is no explicit time dependence)
# 	pop: vector containing the current abundance of all species
# 	parms: 	dummy variable,(used to pass on parameter values), Here, int_mat,str_mat,carry_cap,int_growth are parameter values
#           with: int_mat     - interaction matrix,
#                 str_mat     - competition strength matrix
#                 carry_cap   -  carrying capacity
#                 int_growth  - intrinsic growth rate
# Output:
#	dN: derivative of the Modified Lotka-Volterra equations
NonSpatial_switchingdatalist=c()
Non_switchingdatalist=list()
Elimination_switchingdatalist=list()
Optimization_switchingdatalist=list()
Abundsite_datalist= list()

combinedabundancelist = list()


  
  ###############################################################
  NoSpatial_LVM <- function(t,pop,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,site){ 
    
    dN_LVM=int_growth*pop*((carry_cap-(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])%*%pop)/carry_cap)
    return(dN_LVM)
  }
  
  
  ########################################################################################
  #LVM funtion implementation  of NOn switching on EULER
  ######################
  NonSpatial_euler_meth = function( t_int, y_int, stepsize, t_end,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,site){      
    m = length(y_int)+1
    ## Number of steps and creation of output matrix
    nsteps = ceiling((t_end-t_int)/stepsize)
    Y_out = array(dim=c(nsteps+1,m,sites))
   
    ## loop for implementing the function over nsteps
    for (i in 1:nsteps) {
      
      for(site in 1:sites){
        Y_out[1,,site] = c(t_int, y_int)
      
      Y_out[i+1,,site]= Y_out[i,,site]+ stepsize*c(1, NoSpatial_LVM(Y_out[i,1,site],Y_out[i,2:m,site],Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,site))
      }
    }
    
    return(Y_out)
  }
  
  
  Spatial_LVM <- function(t,pop,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site){ 
    
    dN=int_growth*pop*((carry_cap-(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])%*%pop)/carry_cap) -
        d*Original_pop[,site] + t(Mortality*(Emigration_Prop[site,]%*%t(Original_pop)))
    return(dN)
  }
  
  #####
  #Sampling fuction 
  #sample_g(x)
  
  sample_g <- function(x) {
    if (length(x) <= 1) {
      return(x)
    } else {
      return(sample(x,1))
    }
  }
  
  #########################
  #LVM funtion implementation  of NOn switching on EULER
  ##########################
  Spatial_euler_meth = function( t_int, y_int, stepsize, t_end,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site){      
    
    m = length(y_int)+1
    nsteps = ceiling((t_end-t_int)/stepsize)
    Y_out = array(dim=c(nsteps+1,m,sites))
    
    
    ## Number of steps and creation of output matrix
    # Y_out[1,,1] = c(t_int, y_int)
    # Y_out[1,,2] = c(t_int, y_int)
    # Y_out[1,,3] = c(t_int, y_int)
    
    for (i in 1:nsteps) {
       
      updated_pop= y_int
      
       for(site in 1:sites){
        
        Y_out[1,,site] =c(t_int,y_int)
        ## loop for implementing the function over nsteps
        Y_out[i+1,,site]= Y_out[i,,site]+ stepsize*c(1, Spatial_LVM(Y_out[i,1,site],Y_out[i,2:m,site],Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site))
      
        updated_pop= Y_out[i+1,2:m,site]
        }
    }
    return(Y_out)
  }
  
  #########################################################
  ##LVM funtion implementation  of Elimination switching algorithm on EULER
  # #########################################################
  # 
  Meta_euler_meth_elimination_switch = function(t_int, y_int, stepsize, t_end,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site){
    m = length(y_int)+1
    ## Number of steps and creation of output matrix
    nsteps = ceiling((t_end-t_int)/stepsize)
    Y_out = array(dim=c(nsteps+1,m,sites))
    
    ## loop for implementing the function over nsteps
    
    for (i in 1:nsteps) {
      current_pop=y_int
      
      for(site in 1:sites){
        
        Y_out[1,,site] = c(t_int, y_int)
        
      ## implementing the switching rule
      Meta_Elim_Intstr=(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])* current_pop
      
      #switching_elim=1
      #while(switching_elim<=5){
      repeat{
        IntRowsums=rowSums(Meta_Interaction_Matrix[,,site])
        introwsums_greater1= which(IntRowsums>1,arr.ind = T)
        j_elim=sample_g(introwsums_greater1)
        if (sum(Meta_Interaction_Matrix[,,site][j_elim,])< S){
          break
        }
      }
      
      Vec_elim=Meta_Elim_Intstr[j_elim,]
      # New_Vec_elim = Vec_elim
      # ## Ensuring the maximum value picked is not on the main diagonal
      # k=sort(New_Vec_elim)[(S-1)]
      # Next_max=which(Vec_elim==j)
      # j=sample_g(Next_max)
      MAX=0
      j=1
      for (k in 1:length(Vec_elim)){
        if(k!=j_elim){
          if (MAX < Vec_elim[k]){
            MAX=Vec_elim[k]
            j=k
          }
        }
      }
      
      NI=which(Meta_Interaction_Matrix[,,site][j_elim,]!=1,arr.ind = T)     ## from the sampled row in str_mat check the position of all zeros
      f=sample_g(NI)
      Meta_Interaction_Matrix[,,site][j_elim,c(j,f)]= Meta_Interaction_Matrix[,,site][j_elim,c(f,j)]  ## swapping the value
      
      ##############################################################################
       #switching_elim=switching_elim+1
      #}
      
      ## Calcualting the nestedness of matrix int_mat and appending to NestA at each time step
      # Nesta = nested(Meta_Interaction_Matrix[,,site],method="NODF",rescale = FALSE, normalised = TRUE)
      # NestA <- c(NestA,Nesta)  ## Appending nestedness values at each time step to empty vector
      
      # ### connectance
      # connectance_a=networklevel(Meta_Interaction_Matrix[,,site],index = "connectance")
      # connectance_A = c(connectance_A,connectance_a)
      ##############################################################################
      
      
      Y_out[i+1,,site]= Y_out[i,,site]+ stepsize*c(1, Spatial_LVM(Y_out[i,1,site],Y_out[i,2:m,site],Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site))
      current_pop= Y_out[i+1,2:m,site]
      
      }
    }
    return(Y_out)
    
  }
  
  ## LVM funtion implementation of optimization switching algorithm on EULER
  #################################################################################################
  Meta_euler_meth_optimization_switch = function(t_int, y_int, stepsize, t_end,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site){
    m = length(y_int)+1
    ## Number of steps and creation of output matrix
    nsteps = ceiling((t_end-t_int)/stepsize)
    Y_out = array(dim=c(nsteps+1,m,sites))
    
    ## loop for implementing the function over nsteps
    
    for (i in 1:nsteps) {
      new_pop=y_int
      
      for (site in 1:sites){
        
        Y_out[1,,site] = c(t_int, y_int)
      ## implementing the switching rule
      Opt_Intstr=(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])*new_pop
      
      #switching_opt=1
      #while(switching_opt<=5){
      repeat{
        IntRowsums=rowSums(Meta_Interaction_Matrix[,,site])
        introwsums_greater=which(IntRowsums>2,arr.ind = T)
        if (length(introwsums_greater)>1){
          j_opt= sample_g(introwsums_greater)
        }else{
          j_opt=introwsums_greater
        }
        if (sum(Meta_Interaction_Matrix[,,site][j_opt,])< S){
          break
        }}
    
      Vec_opt= Opt_Intstr[j_opt,]
      j_k=which(Vec_opt!=Vec_opt[j_opt] & Vec_opt!=0,arr.ind = T)
      k_opt=sample_g(j_k)

    
      
      j_m=which(Vec_opt!=Vec_opt[j_opt]& Vec_opt!=Vec_opt[k_opt],arr.ind = T)
      new_removed_ind=c()
      removed_ind=c(j_opt,k_opt)
      enter=TRUE
      while(enter==TRUE){
        m_opt=sample_g(j_m)
        switch_Meta_Interaction_Matrix=Meta_Interaction_Matrix[,,site]
        switch_Meta_Strength_Matrix=Meta_Strength_Matrix[,,site]
        
        Meta_Interaction_Matrix[,,site][j_opt,c(k_opt,m_opt)]=Meta_Interaction_Matrix[,,site][j_opt,c(m_opt,k_opt)]
        Meta_Strength_Matrix[,,site][j_opt,c(k_opt,m_opt)]=Meta_Strength_Matrix[,,site][j_opt,c(m_opt,k_opt)]
        
        dN_switch=int_growth*new_pop*((carry_cap-(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])%*%new_pop)/carry_cap)-
          d*Original_pop[,site] + t(Mortality*(Emigration_Prop[site,]%*%t(Original_pop)))
        
        if ( dN_switch[j_opt] > 0){
             enter=FALSE
          }
        else{
          
          Meta_Interaction_Matrix[,,site]=switch_Meta_Interaction_Matrix
          Meta_Strength_Matrix[,,site]=switch_Meta_Strength_Matrix
          
          new_removed_ind=c(new_removed_ind,m_opt)
  
          j_m=j_m[!(j_m %in% new_removed_ind)]
          
          if (length(j_m)==0){
            enter=FALSE
          }
        }
      }
      ######################################################
      #switching_opt=switching_opt+1
      
      #}
      
      # Calcualting the nestedness of matrix int_mat and appending to NestA at each time step
      # Nestopt = nested(Meta_Interaction_Matrix[,,site],method="NODF",rescale = FALSE, normalised = TRUE)
      # Nest_opt <- c(Nest_opt,Nestopt)  ## Appending nestedness values at each time step to empty vector
      
      
      
      # ### connectance
      # connectance_Opt=networklevel(Meta_Interaction_Matrix[,,site],index = "connectance")
      # connectance_opt = c(connectance_opt,connectance_Opt)
      
      #######################################################################################
      
      Y_out[i+1,,site]= Y_out[i,,site]+ stepsize*c(1, Spatial_LVM(Y_out[i,1,site],Y_out[i,2:m,site],Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site))
      new_pop= Y_out[i+1,2:m,site]
      
      }
    }
    return(Y_out)
    
  }
  # 
  # ####################################################################################################
  
  
  
  S=20
  p=S+1
  sites=3
  #site=1

  int_growth=runif(S,min =0.5,max =0.5)
  carry_cap= runif(S, min=1, max=1)
  pop=runif(S,min =0,max =1)
  
  #Original_pop=matrix(runif(S*sites,min =0,max =1), nrow=S, ncol=sites)
  #for(site in 1:sites){
  #pop=Original_pop[,site]
  #}
  #pop
  
  Original_pop=replicate(3,pop)
  
  d=0.2
  #Emigration_Prop[,-site]*Original_pop[,-site]
  Emigration_Prop=matrix(runif(sites*sites,min =d,max =d), nrow = sites, ncol=sites)
  for (emmirate in 1:sites) {
    Emigration_Prop[emmirate,emmirate]=0
  }
  #Immigration_Prop=matrix(runif(S*sites,min =0.9,max =0.9),nrow = S, ncol = sites)
  #Immigration_Prop=t(Emigration_Prop)
  #Immigration_Prop=abs(rnorm(S,0,0.2))
  Mortality= 0.5 
  # Mortality= matrix(runif(sites*sites,min =0.75,max =0.75), nrow = sites, ncol=sites)
  # for (mortal in 1:sites) {
  #   Mortality[mortal,mortal]=0
  # }
  
  
  steps=10001
  #steps=nsteps+1
  ww=replicate(101,0)
  cc=c(0:100)
  
  ##D for str_mat (competition strength matrix)
  # D = matrix(rnorm(S*S,0,0.2),nrow=S)
  # str_mat=abs(D)
  # for (l in 1:S) {
  #   str_mat[l,l]=1
  # }
  # 
  # Meta_Strength_Matrix= array(replicate(sites,str_mat), dim=c(S,S,sites))
  Meta_Strength_Matrix= array(rep(abs(rnorm(50,0,0.25)),S*S*sites),dim = c(S,S,sites))
  for (j in 1:S){
    for(patch1 in 1:sites){
      Meta_Strength_Matrix[j,j,patch1]=1
    }
  }
  
  len =length(seq(0.1,0.9,0.05))
  con_abund_NonSw = c()
  con_abund_NonSpatialSw=c()
  con_abund_elimination = c()
  con_abund_optimization = c()
  
  #####################
  Meta_max_eigen_NonSw=c();Global_maxeigen_Nonswitch=c();Lead_Eigen_Nonspatial=c()
  Meta_max_eigen_NonSpatialSw=c();Global_maxeigen_Nonspatial=c();Lead_Eigen_Nonswitch=c()
  Meta_max_eigen_elim=c();Global_maxeigen_Elim=c();Lead_Eigen_Elim=c()
  Meta_max_eigen_opt=c();Global_maxeigen_Opt=c() ;Lead_Eigen_Opt=c()
  ###################
  
  ####
  #Appending containers
  ##
  connectance_values=c() ; connect_abundance=c(); connect_NonSpatialabundance=c()
  connect_abund_elim_switch=c();connect_abund_opt_switch=c()
  connect_nest_elim_switch= c();connect_nest_opt_switch= c()
  connect_stab_elim_switch= c();connect_stab_opt_switch= c()
  # connect_mod_elim_switch=c();connect_mod_opt_switch=c()
  connect_nested= c() ;connect_stability= c() #connect_modular = c()
  
  connectA=c()
  
  #### Loop for geneting species dynamics over time for each connectance value
  # z - connectance
  #con_seq - connectance sequence, form o.1 to 0.9 with step of 0.05
  Total_abundance=list()
  con_seq=seq(0.1,0.9,0.05)
  for (z in 1:length(con_seq)){
    # 
    # int_mat=matrix(rbinom(S*S,1,con_seq[z]),S,S)
    # for (j in 1:S) {
    #   int_mat[j,j]=1
    # }
    # 
    # Meta_Interaction_Matrix= array(replicate(sites,int_mat), dim=c(S,S,sites))

    Meta_Interaction_Matrix=array(rep(rbinom(50,1,con_seq[z]),S*S*sites), dim=c(S,S,sites))
    for (meta in 1:S){
      for(patch_int in 1:sites){
        Meta_Interaction_Matrix[meta,meta,patch_int]=1
      }
    }
    
    ##################################################################
    # connect_sim= networklevel(Meta_Interaction_Matrix[,,site],index="connectance")
    # connectA = c(connectA,connect_sim) ### appending the connectance values to connectA
    # connectsim=data.frame(connectA)
    
    
  
    
    
    # #initial_nest=c()
    # initial_nest = nested(Meta_Interaction_Matrix[,,site],method='NODF',rescale=FALSE,normalised = TRUE)
    # connect_nested = append(connect_nested,initial_nest)
    # initial_nest=t(replicate(steps,initial_nest))
    
    ## Figures saved in a folder.
    ## Euler method implemented over time 0-100 with 0.01 time step.
    ## plotting the species dynamics over the time interval.
    ## Non-swithing population dynamics plots for all connectance values.
 #############################################################################################   
    Meta_NonSpatial_Switching=NonSpatial_euler_meth(0,pop, 0.01, 100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,site)
    ################################################################################################################
    
    Meta_NonSpatial_TAbundnance = c()
    
    
    for(patch in 1: sites){
      Total_AbundNonspatial = c()
      for(i in 1:length(Meta_NonSpatial_Switching[,1,1])){
        
        meta_abund_nonspatial=sum(Meta_NonSpatial_Switching[i,2:(S+1), patch])
        Total_AbundNonspatial=c(Total_AbundNonspatial,meta_abund_nonspatial)
      }
      
      Meta_NonSpatial_TAbundnance=data.frame(cbind(Meta_NonSpatial_TAbundnance,Total_AbundNonspatial))
      
    }
    connect_NonSpatialabundance=data.frame(rbind(connect_NonSpatialabundance,Meta_NonSpatial_TAbundnance[steps,]))
    
    
    Meta_Non_Switching=Spatial_euler_meth(0,pop, 0.01, 100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site)
    ################
   #  #######################################################################################################
   
    Meta_Nonswitch_TAbundnance = c()
    
    
    for(patch in 1: sites){
      Total_AbundNonswitch = c()
      for(i in 1:length(Meta_Non_Switching[,1,1])){
        
        meta_abund_nonswitch=sum(Meta_Non_Switching[i,2:(S+1), patch])
        Total_AbundNonswitch=c(Total_AbundNonswitch,meta_abund_nonswitch)
      }
      
      Meta_Nonswitch_TAbundnance=data.frame(cbind(Meta_Nonswitch_TAbundnance,Total_AbundNonswitch))
      
    }
    
    connect_abundance=data.frame(rbind(connect_abundance,Meta_Nonswitch_TAbundnance[steps,]))
    ###################################################################################################
   #  Meta_NonSpatialswitching =transform(Meta_NonSpatial_Switching, T_abundanceNonSpatialsw=rowSums(Meta_NonSpatial_Switching[,2:(S+1),]))
    # Meta_Spatialswitching =transform(Meta_Non_Switching, T_abundanceNonsw=rowSums(Meta_Non_Switching[2:(S+1)]))
    # 
    # 
    # #V2=cbind(Non_Switching,T_abund1)
    # con_abund_NonSpatialSw1 = c()
    # con_abund_NonSpatialSw1=unlist(Meta_NonSpatialswitching[steps,2:(S+1)])
    # con_abund_NonSpatialSw =rbind(con_abund_NonSpatialSw,con_abund_NonSpatialSw1)
    # connect_NonSpatialabundance=append(connect_NonSpatialabundance,Meta_NonSpatialswitching[steps,(S+2)])
    # 
    # 
    # con_abund_NonSw1 = c()
    # con_abund_NonSw1=unlist(Meta_Spatialswitching[steps,2:(S+1)])
    # con_abund_NonSw =rbind(con_abund_NonSw,con_abund_NonSw1)
    # connect_abundance=append(connect_abundance,Meta_Spatialswitching[steps,(S+2)])
    
    
    
    
    # ## Computing Nestedness, Stability, Connectance and Modularity
    # 
    # Y=t(replicate(S,pop))  
    # NestA = c()
    # Nest_int= nested(meta_int_mat[,,site],method = "NODF",rescale = FALSE,normalised = TRUE)
    # NestA= c(NestA,Nest_int)
    # 
    # 
    # connectance_A=c()
    # Init_connectance=networklevel(meta_int_mat[,,site],index = "connectance")
    # connectance_A = c(connectance_A,Init_connectance)
    
    ## Elimination switching population dynamics plots over the time step for every connectance value.
    ## Plots saved in a folder.
    
    ###############################################################################################################
    # ########################################################################################################
    Meta_elimination_switch=Meta_euler_meth_elimination_switch(0, pop, 0.01,100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site)

     ################ 
    # con_abund_elim = c()
    # con_abund_elim=unlist(Meta_elimination_switch[steps,2:(S+1)])
    # con_abund_elimination =rbind(con_abund_elimination,con_abund_elim)
    # #rowSums(con_abund_elimination)
    # 
    # ## computing species total abundance and appending on the population dynamics data frame.
    #################
    Meta_Elimination_TAbundnance = c()
     
     
       for(patch in 1: sites){
         Total_AbundElim = c()
         for(i in 1:length(Meta_elimination_switch[,1,1])){
           
           meta_abund_elimation=sum(Meta_elimination_switch[i,2:(S+1), patch])
           Total_AbundElim=c(Total_AbundElim,meta_abund_elimation)
           }
         
         Meta_Elimination_TAbundnance=data.frame(cbind(Meta_Elimination_TAbundnance,Total_AbundElim))
    
          }
    
    connect_abund_elim_switch=data.frame(rbind(connect_abund_elim_switch, Meta_Elimination_TAbundnance[steps,]))
    ########################################################################
    #Meta_SpatialElimination=transform(Meta_elimination_switch, T_abundanceElim=rowSums(Meta_elimination_switch[2:(S+1),,]))
    # 
    # connect_abund_elim_switch=append(connect_abund_elim_switch, Meta_SpatialElimination[steps,(S+3)])
    # connect_nest_elim_switch= append(connect_nest_elim_switch,mean( Meta_SpatialElimination[,(S+2)]))
    # #connect_stab_elim_switch= append(connect_stab_elim_switch,mean(W_elimination[,(S+3)]))
    # #connect_mod_elim_switch = append(connect_mod_elim_switch,mean(W_elimination[,S+5]))
    # 
    
    # connectance_opt=c()
    # Init_opt_connectance=networklevel(meta_int_mat[,,site],index = "connectance")
    # connectance_opt = c(connectance_opt,Init_opt_connectance)
    # 
    # ######################################################################################################
    # ## Optimization switching population dynamics plots over the time step for every connectance value.
    # ## Plots saved in a folder.
    Meta_optimization_switch=Meta_euler_meth_optimization_switch(0, pop, 0.01,100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site)
    ########################################
    
    Meta_Optimization_TAbundnance = c()
    
    
    for(patch in 1: sites){
      Total_AbundOpt = c()
      for(i in 1:length(Meta_optimization_switch[,1,1])){
        
        meta_abund_Opt=sum(Meta_optimization_switch[i,2:(S+1), patch])
        Total_AbundOpt=c(Total_AbundOpt,meta_abund_Opt)
      }
      
      Meta_Optimization_TAbundnance=data.frame(cbind(Meta_Optimization_TAbundnance,Total_AbundOpt))
      
    }
    
    connect_abund_opt_switch=data.frame(rbind(connect_abund_opt_switch,Meta_Optimization_TAbundnance[steps,]))
    
    #Total abunance over time step
    
    Total_abund= data.frame(cbind(Meta_Non_Switching[,1,1],Meta_NonSpatial_TAbundnance, Meta_Nonswitch_TAbundnance,Meta_Elimination_TAbundnance, Meta_Optimization_TAbundnance))
    Total_Abund=plyr::rename(Total_abund,c(Meta_Non_Switching...1..1. ="time",Total_AbundNonspatial ="TotalAbundance_NonSpatialS1",Total_AbundNonspatial.1 ="TotalAbundance_NonSpatialS2",Total_AbundNonspatial.2 ="TotalAbundance_NonSpatialS3",
                                           Total_AbundNonswitch="TotalAbundance_NonSwitchS1",Total_AbundNonswitch.1="TotalAbundance_NonSwitchS2",Total_AbundNonswitch.2="TotalAbundance_NonSwitchS3",
                                           Total_AbundElim ="Total_AbundanceElim_S1", Total_AbundElim.1 ="Total_AbundanceElim_S2",Total_AbundElim.2 ="Total_AbundanceElim_S3",
                                           Total_AbundOpt = "Total_AbundanceOpt_S1",Total_AbundOpt.1 = "Total_AbundanceOpt_S2",Total_AbundOpt.2 = "Total_AbundanceOpt_S3"))
    
    for (site in 1:sites){
      Total_Abund_site = cbind(Total_Abund[1],Total_Abund[(site+1)],Total_Abund[(site+4)],Total_Abund[(site+7)],Total_Abund[(site+10)])
      # 
      Total_Abund_site_melted=melt(Total_Abund_site,id.vars = "time")
      ggplot(data=Total_Abund_site_melted,aes(x=time,y=value, group=variable))+
        #geom_hline(aes(yintercept=0),size=1,colour='red',linetype="dotted")+
        geom_line(aes(color=variable),size=1)+
        scale_color_discrete(breaks=c("TotalAbundance_NonSpatialS1","TotalAbundance_NonSwitchS1","Total_AbundanceElim_S1","Total_AbundanceOpt_S1"),
                             labels=c("Non_Spatial","Non-switch","Elimination", "Optimization"))+
        theme_few()+

        labs(y=expression("Re"~(lambda)),x="Time",aspect.ratio = 2)+
        theme(legend.position = c(0.85,0.15),legend.title = element_blank(), aspect.ratio = 1)+
        theme(axis.title=element_text(size=12,face="bold"),
              axis.text =element_text(face = "bold",size=12))
      ggsave(paste0("TotalAbundance_",toString(c(site,z)),".pdf"),width = 5,height = 5)
      
    }
    
    
    ##Stability computatation
    
    Modified_Spatial_LVM_stab <- function(pop){ 
      
      dN_Spatial_stab=int_growth*pop*((carry_cap-(IM*SM)%*%pop)/carry_cap) -
        d*Original_pop[,site] + t(Mortality*(Emigration_Prop[site,]%*%t(Original_pop)))
      
      return(dN_Spatial_stab)
    }
    
    Modified_NoSpatial_LVM <- function(pop){ 
      
      dN_NonSpatial=int_growth*pop*((carry_cap-(IM*SM)%*%pop)/carry_cap)
      return(dN_NonSpatial)
    }
    
    max_eigen_NonSw=c();max_eigen_elim=c();max_eigen_NonSpatialSw=c()
    max_eigen_opt=c();Max_eigen_NonSW=c();Max_eigen_NonSpatialSw=c()
    Max_eigen_elim=c();Max_eigen_opt=c()
    Jacobian_Nonswitch = list();Jacobian_Nonspatial=list()
    Jacobian_Elim = list();Jacobian_Opt=list()
    
    for (site in 1:sites){
    
    IM=Meta_Interaction_Matrix[,,site]
    SM =Meta_Strength_Matrix[,,site]
    ######################################################################################
    VSpatialNonSW=Meta_Non_Switching[steps,2:(S+1),site]
    VNonSpatialSW=Meta_NonSpatial_Switching[steps,2:(S+1),site]
    
    V_elim_sw=Meta_elimination_switch[steps,2:(S+1),site]
    V_opt_sw=Meta_optimization_switch[steps,2:(S+1),site]
    
    Jac_NonSW=jacobian( Modified_Spatial_LVM_stab,VSpatialNonSW)
    Jacobian_Nonswitch[[site]] = jacobian( Modified_Spatial_LVM_stab,VSpatialNonSW)
    
    Jac_NonSpatialSW=jacobian(Modified_NoSpatial_LVM,VNonSpatialSW)
    Jacobian_Nonspatial[[site]]= jacobian(Modified_NoSpatial_LVM,VNonSpatialSW)
    
    Jac_Meta_elimination_switch=jacobian(Modified_Spatial_LVM_stab,V_elim_sw)
    Jacobian_Elim[[site]]=jacobian(Modified_Spatial_LVM_stab,V_elim_sw)
    
    Jac_optimization_Switch=jacobian(Modified_Spatial_LVM_stab,V_opt_sw)
    Jacobian_Opt[[site]]=jacobian(Modified_Spatial_LVM_stab,V_opt_sw)
    
    eigen_NonSW=eigen(Jac_NonSW)$values
    eigen_NonSpatialSW=eigen(Jac_NonSpatialSW)$values
    eigen_elim=eigen(Jac_Meta_elimination_switch)$values
    eigen_opt=eigen(Jac_optimization_Switch)$values
    
    Max_eigen_NonSW=data.frame(max(Re(eigen_NonSW)))
    Max_eigen_NonSpatialSw=data.frame(max(Re(eigen_NonSpatialSW)))
    Max_eigen_elim =data.frame(max(Re(eigen_elim)))
    Max_eigen_opt =data.frame(max(Re(eigen_opt)))
    
    max_eigen_NonSw=append(max_eigen_NonSw,Max_eigen_NonSW)
    max_eigen_NonSpatialSw=append(max_eigen_NonSpatialSw,Max_eigen_NonSpatialSw)
    max_eigen_elim=append(max_eigen_elim,Max_eigen_elim)
    max_eigen_opt=append(max_eigen_opt,Max_eigen_opt)
    
    
    }
    
    Meta_max_eigen_NonSpatialSw = data.frame(cbind(Meta_max_eigen_NonSpatialSw,unlist(max_eigen_NonSpatialSw)))
    Meta_max_eigen_NonSw = data.frame(cbind(Meta_max_eigen_NonSw,unlist(max_eigen_NonSw)))
    Meta_max_eigen_elim = data.frame(cbind(Meta_max_eigen_elim,unlist(max_eigen_elim)))
    Meta_max_eigen_opt = data.frame(cbind(Meta_max_eigen_opt,unlist(max_eigen_opt)))
    
    
    
    ###########################################################################
    ### Commputing stability over time for different values of connectance
    Meta_elim_max_eigen=c()
    #elim_Max_eigen=c()
    Meta_opt_max_eigen=c()
    #opt_Max_eigen=c()
    #Meta_Max_eigen_NonSpatialSWT=c()
    Meta_max_eigen_NonSpatialSwT= c()
    #Meta_Max_eigen_NonSWT=c()
    Meta_max_eigen_NonSWT=c()
    
    for (site in 1:sites){
      
      Max_eigen_NonSWT=c()
      max_eigen_NonSwT=c()
      
      for(nsw in 1:length(Meta_Non_Switching[,1,1])){
      
        VSpatialNonSWT=Meta_Non_Switching[nsw,2:(S+1),site]
        
        Jac_NonSWT=jacobian(Modified_Spatial_LVM_stab,VSpatialNonSWT)
        eigen_NonSWT=eigen(Jac_NonSWT)$values
        Max_eigen_NonSWT=data.frame(max(Re(eigen_NonSWT)))
        max_eigen_NonSwT=append(max_eigen_NonSwT,Max_eigen_NonSWT)
        }
      Meta_max_eigen_NonSWT = data.frame(cbind(Meta_max_eigen_NonSWT,unlist(max_eigen_NonSwT)))
    }
    
    Meta_max_eigen_NonSWT_Stab=plyr::rename(Meta_max_eigen_NonSWT,c(cbind.Meta_max_eigen_NonSWT..unlist.max_eigen_NonSwT..="StabNonSwitchsite1",
                                                                    unlist.max_eigen_NonSwT.="StabNonSwitchsite2",unlist.max_eigen_NonSwT..1="StabNonSwitchsite3"))

  
    for (site in 1:sites){
      
      max_eigen_NonSpatialSwT =c()
      Max_eigen_NonSpatialSWT =c()
      
      for(nonspatialsw in 1:length(Meta_NonSpatial_Switching[,1,site])){
     
         VNonSpatialSWT=Meta_NonSpatial_Switching[nonspatialsw,2:(S+1),site]
         ##########################
         Jac_NonSpatialSWT=jacobian(Modified_NoSpatial_LVM,VNonSpatialSWT)
         eigen_NonSpatialSWT=eigen(Jac_NonSpatialSWT)$values
         Max_eigen_NonSpatialSWT=data.frame(max(Re(eigen_NonSpatialSWT)))
         max_eigen_NonSpatialSwT=append(max_eigen_NonSpatialSwT,Max_eigen_NonSpatialSWT)
         }
      
      Meta_max_eigen_NonSpatialSwT = data.frame(cbind(Meta_max_eigen_NonSpatialSwT,unlist(max_eigen_NonSpatialSwT)))
    }
    
    Meta_max_eigen_NonSpatialSwT_Stab=plyr::rename(Meta_max_eigen_NonSpatialSwT,c(cbind.Meta_max_eigen_NonSpatialSwT..unlist.max_eigen_NonSpatialSwT..="StabNonSpatialsite1",
                                                                                  unlist.max_eigen_NonSpatialSwT.="StabNonSpatialsite2",unlist.max_eigen_NonSpatialSwT..1="StabNonSpatialsite3"))
      
    for (site in 1:sites){
      
      elim_max_eigen=c()
      elim_Max_eigen=c()
      
      
    for(elim_sw in 1:length(Meta_elimination_switch[,1,site])){
      
      V_elim=Meta_elimination_switch[elim_sw,2:(S+1),site]
      ##########################
      Jac_elim=jacobian(Modified_Spatial_LVM_stab,V_elim)
      elim_eigen=eigen(Jac_elim)$values
      elim_Max_eigen=data.frame(max(Re(elim_eigen)))
      elim_max_eigen=append(elim_max_eigen,elim_Max_eigen)
      
    }
      Meta_elim_max_eigen = data.frame(cbind(Meta_elim_max_eigen,unlist(elim_max_eigen)))
    }

    Meta_elim_max_eigen_Stab=plyr::rename(Meta_elim_max_eigen,c(cbind.Meta_elim_max_eigen..unlist.elim_max_eigen..="StabElimsite1",
                                                           unlist.elim_max_eigen.="StabElimsite2",unlist.elim_max_eigen..1="StabElimsite3"))
  
    for (site in 1:sites){
    
      opt_max_eigen =c()
      opt_Max_eigen = c()
   
       for(opt_sw in 1:length(Meta_optimization_switch[,1,site])){
      
         V_opt=Meta_optimization_switch[opt_sw,2:(S+1),site]
         ##########################
         Jac_opt=jacobian(Modified_Spatial_LVM_stab,V_opt)
         opt_eigen=eigen(Jac_opt)$values
         opt_Max_eigen=data.frame(max(Re(opt_eigen)))
         opt_max_eigen=append(opt_max_eigen,opt_Max_eigen)
         }
      
      Meta_opt_max_eigen = data.frame(cbind(Meta_opt_max_eigen, unlist(opt_max_eigen)))
    
    }
    
    
    Meta_opt_max_eigen_Stab=plyr::rename(Meta_opt_max_eigen,c(cbind.Meta_opt_max_eigen..unlist.opt_max_eigen..="StabOptsite1",
                                                                unlist.opt_max_eigen.="StabOptsite2",unlist.opt_max_eigen..1="StabOptsite3"))
    
    Stabb_time= data.frame(cbind(Meta_Non_Switching[,1,1],Meta_max_eigen_NonSpatialSwT_Stab,Meta_max_eigen_NonSWT_Stab,Meta_elim_max_eigen_Stab,Meta_opt_max_eigen_Stab))
    Stab_time=plyr::rename(Stabb_time,c(Meta_Non_Switching...1..1.="time"))
    # Stab_time_site1 = cbind(Stab_time[1],Stab_time[2],Stab_time[5],Stab_time[8],Stab_time[11])
    # Stab_time_site2 = cbind.data.frame(Stab_time$time,Stab_time$StabNonSpatialsite2,Stab_time$StabNonSwitchsite2,Stab_time$StabElimsite2,Stab_time$StabOptsite2)
    # Stab_time_site3 = cbind.data.frame(Stab_time$time,Stab_time$StabNonSpatialsite3,Stab_time$StabNonSwitchsite3,Stab_time$StabElimsite3,Stab_time$StabOptsite3)
    
    for (site in 1:sites){
      Stab_time_site = cbind(Stab_time[1],Stab_time[(site+1)],Stab_time[(site+4)],Stab_time[(site+7)],Stab_time[(site+10)])
    
    # Stab_time_site1_melted=melt(Stab_time_site,id.vars = "time")
    #  ggplot(data=Stab_time_site1_melted,aes(x=time,y=value, group=variable))+
    #    geom_hline(aes(yintercept=0),size=1,colour='red',linetype="dotted")+
    #    geom_line(aes(color=variable),size=1)+
    #    scale_color_discrete(breaks=c("StabNonSpatialsite1","StabNonSwitchsite1","StabElimsite1","StabOptsite1"),
    #                         labels=c("Non_Spatial","Non-switch","Elimination", "Optimization"))+
    #    theme_few()+
    #    
    #    labs(y=expression("Re"~(lambda)),x="Time",aspect.ratio = 2)+
    #    theme(legend.position = c(0.85,0.85),legend.title = element_blank(), aspect.ratio = 1)+
    #    theme(axis.title=element_text(size=12,face="bold"),
    #         axis.text =element_text(face = "bold",size=12))
    # ggsave(paste0("Stabtime_",toString(c(site,z)),".pdf"),width = 5,height = 5)
    # 
    }
    
    ######################################
    #Global stability computation
    Jacobians_Nonspatial = adiag(Jacobian_Nonspatial[[1]],Jacobian_Nonspatial[[2]],Jacobian_Nonspatial[[3]])
    Jacobians_Nonswitch = adiag(Jacobian_Nonswitch[[1]],Jacobian_Nonswitch[[2]],Jacobian_Nonswitch[[3]])
    Jacobians_Elim = adiag(Jacobian_Elim[[1]],Jacobian_Elim[[2]],Jacobian_Elim[[3]])
    Jacobians_Opt = adiag(Jacobian_Opt[[1]],Jacobian_Opt[[2]],Jacobian_Opt[[3]])
    
    Intra_specific_Matrix =diag(nrow=(S*sites))
    Dispersal_matrix=`diag<-`(matrix(0,(S*sites),(S*sites)),-d)
    Interpatch_dispersal = `diag<-`(matrix(0,S,S),d*0.5)
    Dispersal_matrix[(S+1):(S*2),1:(S*2)]=Interpatch_dispersal;Dispersal_matrix[1:S,(S*2+1):(S*3)]=Interpatch_dispersal
    Dispersal_matrix[1:S,(S+1):(S*2)]=Interpatch_dispersal;Dispersal_matrix[(S*2+1):(S*3),1:(S*2)]=Interpatch_dispersal
    Dispersal_matrix[(S+1):(S*2),(S*2+1):(S*3)]=Interpatch_dispersal;Dispersal_matrix[(S*2+1):(S*3),(S+1):(S*2)]=Interpatch_dispersal
    
    Global_Jacobian_Nonspatial=-Intra_specific_Matrix+Dispersal_matrix +Jacobians_Nonspatial
    Global_eigen_Nonspatial = eigen(Global_Jacobian_Nonspatial)$values
    Lead_Eigen_Nonspatial = data.frame(max(Re(Global_eigen_Nonspatial)))
    Global_maxeigen_Nonspatial = append(Global_maxeigen_Nonspatial,Lead_Eigen_Nonspatial)
    
    Global_Jacobian_Nonswitch=-Intra_specific_Matrix+Dispersal_matrix +Jacobians_Nonswitch
    Global_eigen_Nonswitch = eigen(Global_Jacobian_Nonswitch)$values
    Lead_Eigen_Nonswitch = data.frame(max(Re(Global_eigen_Nonswitch)))
    Global_maxeigen_Nonswitch = append(Global_maxeigen_Nonswitch,Lead_Eigen_Nonswitch)
    
    Global_Jacobian_Elim=-Intra_specific_Matrix+Dispersal_matrix +Jacobians_Elim
    Global_eigen_Elim = eigen(Global_Jacobian_Elim)$values
    Lead_Eigen_Elim = data.frame(max(Re(Global_eigen_Elim)))
    Global_maxeigen_Elim = append(Global_maxeigen_Elim,Lead_Eigen_Elim)
    
    Global_Jacobian_Opt=-Intra_specific_Matrix+Dispersal_matrix +Jacobians_Opt
    Global_eigen_Opt = eigen(Global_Jacobian_Opt)$values
    Lead_Eigen_Opt = data.frame(max(Re(Global_eigen_Opt)))
    Global_maxeigen_Opt = append(Global_maxeigen_Opt,Lead_Eigen_Opt)

  }
  
  #Global Productivity For each connectance value
  Global_TotalAbundance_Nonspatial = transform(connect_NonSpatialabundance, GlobalAbundance_Nonspatial=rowSums(connect_NonSpatialabundance[1:sites]))
  Global_TotalAbundance_Nonswitch = transform(connect_abundance,GlobalAbundance_Noswitch=rowSums(connect_abundance[1:sites]))
  Global_TotalAbundance_Elimination =transform(connect_abund_elim_switch,GlobalAbundance_Elim= rowSums(connect_abund_elim_switch[1:sites]))
  Global_TotalAbundance_Optimization = transform(connect_abund_opt_switch,GlobalAbundance_Opt=rowSums(connect_abund_opt_switch[1:sites]))
  
  Global_Productivity = data.frame(cbind(Global_TotalAbundance_Nonspatial,Global_TotalAbundance_Nonswitch,
                                         Global_TotalAbundance_Elimination,Global_TotalAbundance_Optimization))
  
  Local_Stability = data.frame(cbind(stability_Nonspatialsite=t(Meta_max_eigen_NonSpatialSw),stability_Nonswitchsite=t(Meta_max_eigen_NonSw),
                                     stability_Elimsite=t(Meta_max_eigen_elim),stability_Optsite=t(Meta_max_eigen_opt)))
  
  Global_Stability=data.frame(rbind(unlist(Global_maxeigen_Nonspatial),unlist(Global_maxeigen_Nonswitch),unlist(Global_maxeigen_Elim),unlist(Global_maxeigen_Opt)))
  
  
  Abundance_Connectance = data.frame(cbind(Connectance=con_seq ,Global_Productivity,Local_Stability, t(Global_Stability)))
                                           #global_StabNospatial=t(unlist(Global_maxeigen_Nonspatial)),global_StabNonswitch=t(unlist(Global_maxeigen_Nonswitch)),
                                           #global_StabEli=t(unlist(Global_maxeigen_Elim)),Global_StabOpt=t(unlist(Global_maxeigen_Opt))))
  
row.names(Abundance_Connectance)<- NULL



##Productivity per site
for (Patch in 1:sites){

  Productivity_Site = cbind(Abundance_Connectance[1], Abundance_Connectance[c((Patch+1),(Patch+5),(Patch+9),(Patch+13))])
  ProductivitySite_melted=melt(Productivity_Site,id.vars = "Connectance")
  ggplot(data=ProductivitySite_melted,aes(x=Connectance,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2)+
  geom_smooth(method = "loess",se=F,aes(color=variable))+
  theme_few()+
  labs(y="Productivity",x="Connectance")+
  theme(legend.position = "None",aspect.ratio = 1)+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text =element_text(face = "bold",size=12))

  ggsave(paste0("Productivitysite_",toString(Patch),".pdf"),width = 5,height = 5)
  ggsave(paste0("Productivitysite_",toString(Patch),".tiff"),width = 5,height = 5)

}

for (Patches in 1:sites){
  
  LocalStability_Site = cbind(Abundance_Connectance[1], Abundance_Connectance[c((Patches+17),(Patches+20),(Patches+23),(Patches+26))])
  LocalStabilitySite_melted=melt(LocalStability_Site,id.vars = "Connectance")
  ggplot(data=LocalStabilitySite_melted,aes(x=Connectance,y=value, group=variable))+
    geom_hline(aes(yintercept=0),size=1,colour='red',linetype='dotted')+
    geom_point(aes(shape=variable,color=variable),size=2)+
    geom_smooth(method = "loess",se=F,aes(color=variable))+
    theme_few()+
    labs(y=expression("Re"~(lambda)),x="Connectance")+
    theme(legend.position = "None",aspect.ratio = 1)+
    theme(axis.title=element_text(size=12,face="bold"),
          axis.text =element_text(face = "bold",size=12))
  
  ggsave(paste0("localstabilitysite_",toString(Patches),".pdf"),width = 5,height = 5)
  ggsave(paste0("localstabilitysite_",toString(Patches),".tiff"),width = 5,height = 5)
  
}

##Global Productivity 

Productivity= cbind(Abundance_Connectance[1],Abundance_Connectance[c((sites+2),(sites+6),(sites+10),(sites+14))])
Productivity_melted=melt(Productivity,id.vars = "Connectance")
ggplot(data=Productivity_melted,aes(x=Connectance,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2)+
  geom_smooth(method = "loess",se=F,aes(color=variable))+
  theme_few()+
  labs(y="Productivity",x="Connectance")+
  theme(legend.position = "None",aspect.ratio = 1)+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text =element_text(face = "bold",size=12))
ggsave("Productivity.pdf",width = 5,height = 5)
ggsave("Productivity.tiff",width = 5,height = 5)

##Global Stability
GlobalStability= cbind(Abundance_Connectance[1],Abundance_Connectance[30:33])
GlobalStability_melted=melt(GlobalStability,id.vars = "Connectance")
ggplot(data=GlobalStability_melted,aes(x=Connectance,y=value,group=variable))+
  geom_hline(aes(yintercept=0),size=1,colour='red',linetype='dotted')+
  geom_point(aes(shape=variable,color=variable),size=2)+
  geom_smooth(method = "loess",se=F,aes(color=variable))+
  theme_few()+
  labs(y=expression("Re"~(lambda)),x="Connectance")+
  theme(legend.position = "None",aspect.ratio = 1)+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text =element_text(face = "bold",size=12))
ggsave("Global_Stab.pdf",width = 5,height = 5)
ggsave("GLobal_Stab.tiff",width = 5,height = 5)


  save(Meta_NonSpatial_Switching,Meta_Non_Switching,Meta_elimination_switch,Meta_optimization_switch, file="Sites.RData")
  #save(Modified_NoSpatial_LVM,Modified_Spatial_LVM_stab,NoSpatial_LVM, Non_euler_meth, Spatial_LVM,Spatial_euler_meth,Meta_euler_meth_elimination_switch,Meta_euler_meth_optimization_switch,file="Metafunctions.RData")
  