############################################################
#Project:       Master Research project
#Description :  MUltispecies Spatial Competition Lotka-Voterra model Euler simulation functions

# By:     Gilbert K Langat

setwd("E:/MSC/Code Repo/MastersCode")
############################################################
## FUNCTION DEFINITIONS

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
#############################################################

NoSpatial_LVM <- function(t,pop,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,site){ 
  
  dN_LVM=int_growth*pop*((carry_cap-(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])%*%pop)/carry_cap)
  return(dN_LVM)
  }


#############################################################
#LVM funtion implementation  of NOn switching on EULER
#############################################################

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

###########################################################################
# Spatial Lotka-Volterra competition model

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

##############################################################
#LVM funtion implementation  of NOn switching on EULER
#############################################################
Spatial_euler_meth = function( t_int, y_int, stepsize, t_end,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site){      
  
  m = length(y_int)+1
  
  # Number of time steps
  nsteps = ceiling((t_end-t_int)/stepsize)
  
  # dimensions of output matrix
  Y_out = array(dim=c(nsteps+1,m,sites))
  
  # Y_out[1,,1] = c(t_int, y_int)
  # Y_out[1,,2] = c(t_int, y_int)
  # Y_out[1,,3] = c(t_int, y_int)
  
  # loop for implementing the Euler function over nsteps(time steps)
  for (i in 1:nsteps) {
    
    #Updating the population densities at the start of each step
    updated_pop= y_int
    
    #for throught the number of sites(patches)
    
    for(site in 1:sites){
      
      Y_out[1,,site] =c(t_int,y_int)
      
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
    
    #Updating the population after each timestep
    current_pop=y_int
    
    # loop over the number of patches
    for(site in 1:sites){
      
      Y_out[1,,site] = c(t_int, y_int)
      
      ###############################################
      # Elimination switching rule implemention
      #############################################
      
      #Community matrix
      Meta_Elim_Intstr=(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])* current_pop
      
      #switching_elim=1
      #while(switching_elim<=5){
      
      # Ensuring the selected species interacts with more than one species
      
      repeat{
        
        IntRowsums=rowSums(Meta_Interaction_Matrix[,,site])
        introwsums_greater1= which(IntRowsums>1,arr.ind = T)
        j_elim=sample_g(introwsums_greater1)
        if (sum(Meta_Interaction_Matrix[,,site][j_elim,])< S){
          break
        }
        }
      
      Vec_elim=Meta_Elim_Intstr[j_elim,]
      ##########################
      # New_Vec_elim = Vec_elim
      # ## Ensuring the maximum value picked is not on the main diagonal
      # k=sort(New_Vec_elim)[(S-1)]
      # Next_max=which(Vec_elim==j)
      # j=sample_g(Next_max)
      #########################
      
      # Ensuring the maximum value picked is not on the main diagonal
      
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
      
      # Check the position of all zeros and sample 1 element
      
      NI=which(Meta_Interaction_Matrix[,,site][j_elim,]!=1,arr.ind = T)     
      f=sample_g(NI)
      
      ## swapping between the zero selected the original value
      
      Meta_Interaction_Matrix[,,site][j_elim,c(j,f)]= Meta_Interaction_Matrix[,,site][j_elim,c(f,j)]
      
      #switching_elim=switching_elim+1
      
      #} # End of repeated loop of switches
      
      # Euler iteration
      
      Y_out[i+1,,site]= Y_out[i,,site]+ stepsize*c(1, Spatial_LVM(Y_out[i,1,site],Y_out[i,2:m,site],Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site))
      current_pop= Y_out[i+1,2:m,site]
      
    
      }
  
    }
  
  return(Y_out)
  

  }# End of function

###############################################################################
# Optimization switching euler simulation

##############################################################################

Meta_euler_meth_optimization_switch = function(t_int, y_int, stepsize, t_end,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site){
  
  m = length(y_int)+1
  
  # Number of time steps
  
  nsteps = ceiling((t_end-t_int)/stepsize)
  
  # Output matrix dimensions
  Y_out = array(dim=c(nsteps+1,m,sites))
  
  #Updating species population at end of each time step
  new_pop=y_int
  
  # loop for implementing the function over nsteps(time steps)
  
  for (i in 1:nsteps) {
  
    for (site in 1:sites){
      
      Y_out[1,,site] = c(t_int, y_int)
      
      ############################################
      # Optimization switching rule implementation at each time step
      ############################################
      
      #community matrix
      
      Opt_Intstr=(Meta_Interaction_Matrix[,,site]*Meta_Strength_Matrix[,,site])*new_pop
      
      # switch repeatedly a number of times at each time step
      
      #switching_opt=1
      #while(switching_opt<=5){
      
      # Select species interaction with more than one partner
      
      repeat{
        IntRowsums=rowSums(Meta_Interaction_Matrix[,,site])
        introwsums_greater=which(IntRowsums>2,arr.ind = T)
        if (length(introwsums_greater)>1){
          j_opt= sample_g(introwsums_greater)
          }
        else{
          j_opt=introwsums_greater
          
          }
        
        if (sum(Meta_Interaction_Matrix[,,site][j_opt,])< S){
          break
        
        }
        }
      
      #the selected vector
      Vec_opt= Opt_Intstr[j_opt,]
      
      # Choose the  2 non-interacting partners (zeros in int_mat) & sample_g one
      
      j_k=which(Vec_opt!=Vec_opt[j_opt] & Vec_opt!=0,arr.ind = T)
      k_opt=sample_g(j_k)
      j_m=which(Vec_opt!=Vec_opt[j_opt]& Vec_opt!=Vec_opt[k_opt],arr.ind = T)
      
      
      new_removed_ind=c()
      removed_ind=c(j_opt,k_opt)
      
      # Checking if switching increase the species growth
      
      enter=TRUE
      while(enter==TRUE){
        m_opt=sample_g(j_m)
        
        #Updating the interaction and strength matrices
        
        switch_Meta_Interaction_Matrix=Meta_Interaction_Matrix[,,site]
        switch_Meta_Strength_Matrix=Meta_Strength_Matrix[,,site]
        
        #Swapping the selected partners 
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
        
          }#End of else loop
      
        } #End of while loop
      
      
      #switching_opt=switching_opt+1
      
      #} #End of repeated switching

      # Euler iteration over nsteps 
      
      Y_out[i+1,,site]= Y_out[i,,site]+ stepsize*c(1, Spatial_LVM(Y_out[i,1,site],Y_out[i,2:m,site],Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site))
      new_pop= Y_out[i+1,2:m,site]
      
    
      }
  
    } # End of site loop
  
  return(Y_out)
  

  }# End of time iteration loop

#################################################
#Stability computation function 
#################################################

#Spatial case

Modified_Spatial_LVM_stab <- function(pop){ 
  
  dN_Spatial_stab=int_growth*pop*((carry_cap-(IM*SM)%*%pop)/carry_cap) -
    d*Original_pop[,site] + t(Mortality*(Emigration_Prop[site,]%*%t(Original_pop)))
  
  return(dN_Spatial_stab)

  }

#Non-spatial case
Modified_NoSpatial_LVM <- function(pop){ 
  
  dN_NonSpatial=int_growth*pop*((carry_cap-(IM*SM)%*%pop)/carry_cap)
  return(dN_NonSpatial)

  }

######################################################################################################
# SAVING THE FUNCTIONS 
save(NoSpatial_LVM,NonSpatial_euler_meth,Spatial_LVM,Spatial_euler_meth,Meta_euler_meth_elimination_switch,Meta_euler_meth_optimization_switch,
     Modified_NoSpatial_LVM,Modified_Spatial_LVM_stab,sample_g,file="Spatialmodelfuctions.RData")

