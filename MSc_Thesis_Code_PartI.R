#########################################################################################
#CODE FOR RESEARCH MASTERS THESIS
#Project: Implementattion of Modified Lotka-Volterra model on Euler method
#By: Gilbert K Langat



##############################################
#Clearing memory and Saving code to folder
##############################################

rm(list=ls(all=TRUE)) 


################################################################################
## IMPORTING FUNCTIONs 
################################################################################

source("E:/MSC/Code Repo/MastersCode/Part1_Functions.R")

setwd("E:/MSC/Code Repo/MastersCode")

#set.seed(123)
##########################
# LOAD LIBRARIES

library(BiodiversityR) # for nestedness
library(pracma) #Numerical computations 
library(vegan)  # nestedness
library(bipartite)  # Community structuring 
library(tidyverse)
library(ggthemes)  # Data visualization 
library(reshape2)  # for organizing data



################################################################################
#Appending containers
#

CombinedAbund_datalist=list()   ;CombinedCompet_datalist1=list()
CombinedCompet_datalist2=list() ;CompetCombineddatalist=list()
Overall_Speciesranking=list()   ;con_abund_NonSw = c()
con_abund_elimination = c()     ;con_abund_optimization = c()

len =length(seq(0.1,0.9,0.05))

################################################################################

#Loop over different values of S(number of species)
#S number of species

##############################################
species=c(10,20,50)

for (S in species){
      
      p=S+1
      
      int_growth=runif(S,min =0.5,max =0.5)
      carry_cap= runif(S, min=1, max=1)
      pop=runif(S,min =0,max =1)
      
      steps=10001
      ww=replicate(101,0)
      cc=c(0:100)
    
      Abundsigma_datalist=list()    ;Compet_sigma_datalist1=list()
      Compet_sigma_datalist2=list() ;Species_rankdatalist=list()

 
      con_abund_elimination=c();con_abund_optimization=c();con_abund_NonSw=c()
      connectance_values=c() ; connect_abundance=c()
      connect_abund_elim_switch=c();connect_abund_opt_switch=c()
      connect_nest_elim_switch= c();connect_nest_opt_switch= c()
      connect_stab_elim_switch= c();connect_stab_opt_switch= c()
    
      connect_nested= c() ;connect_stability= c() 
      max_eigen_NonSw=c();max_eigen_elim=c()
      max_eigen_opt=c();Max_eigen_NonSW=c()
      Max_eigen_elim=c();Max_eigen_opt=c()
      connectA=c()         
      connect_sim=c()
      
      ##################################################
      # Loop over different values of sigma(Standard deviation)
      # SD - standard deviation
      #*std_dvtn - SD sequence, from 0.1 to 0.5 with step of 0.2
      ##################################################
      
      std_dvtn=seq(0.1,0.5,0.2)
      for (SD in 1:length(std_dvtn)){
        
            D = matrix(rnorm(S*S,0,std_dvtn[SD]),nrow=S)
            
            str_mat=abs(D)
            for (l in 1:S) {
              str_mat[l,l]=1
            }
            con_abund_elimination=c();con_abund_optimization=c();con_abund_NonSw=c()
            connectance_values=c() ; connect_abundance=c()
            connect_abund_elim_switch=c();connect_abund_opt_switch=c()
            connect_nest_elim_switch= c();connect_nest_opt_switch= c()
            connect_stab_elim_switch= c();connect_stab_opt_switch= c()
            
            connect_nested= c() ;connect_stability= c() 
            max_eigen_NonSw=c();max_eigen_elim=c()
            max_eigen_opt=c();Max_eigen_NonSW=c()
            Max_eigen_elim=c();Max_eigen_opt=c()
            connectA=c()         
            connect_sim=c()
            
            #############################################################################
            # Loop for geneting species dynamics over time for each connectance value
            # z - connectance
            #con_seq - connectance sequence, from 0.1 to 0.9 with step of 0.05
            ############################################################################
            
            con_seq=seq(0.1,0.9,0.1)
          
            for (z in con_seq){
              int_mat=matrix(rbinom(S*S,1,z),S,S)
              ind<-lower.tri(int_mat)
              int_mat[ind]<-t(int_mat)[ind]
              #mat[lower.tri(mat,diag = F)]<-rbinom((((S*(S+1))/2)-S),1,z)
              #upper_mat<-t(mat)
              #int_mat <-upper_mat+mat
              #int_mat=matrix(rbinom(S*S,1,z),S,S)
              for (j in 1:S) {
                int_mat[j,j]=1
                }
              connect_sim= networklevel(int_mat,index="connectance")
              connectA = append(connectA,connect_sim) ### appending the connectance values to connectA
              connectsim=data.frame(connectA)
                      
                        
                        initial_nest=c()
                        initial_nest = nested(int_mat,method='NODF',rescale=FALSE,normalised = TRUE)
                        connect_nested = append(connect_nested,initial_nest)
                        initial_nest=t(replicate(steps,initial_nest))
                        
                        # Figures saved in a folder.
                        # Euler method implemented over time 0-100 with 0.01 time step.
                        # plotting the species dynamics over the time interval.
                        # Non-swithing population dynamics plots for all connectance values.
                        
                        Non_Switching=NonSwitch_Euler(0,pop, 0.01, 100,int_mat,str_mat,carry_cap,int_growth)
                    
                        #########################################################################################    
                        # None_Switch= cbind(Non_Switching[1],Non_Switching[2:p])
                        # None_Switch_melted=melt(None_Switch,id.vars = "X1")
                        # ggplot(data=None_Switch_melted,aes(x=X1,y=value, group=variable))+
                        #   #geom_point(aes(shape=variable),size=1)+
                        #   #scale_shape_manual(values = c(1,16,18))+
                        #   geom_line(aes(linetype=variable),size=1)+
                        #   #scale_linetype_manual(values = c("solid","longdash","dotted"))+
                        #   theme_few()+
                        #   labs(y="Abundance",x="Time",aspect.ratio = 2)+
                        #   theme(legend.position = "None",aspect.ratio = 1)+
                        #   theme(axis.title=element_text(size=12,face="bold"),
                        #         axis.text =element_text(face = "bold",size=12))
                        # ggsave(paste0("Noneswitchdynamics_",toString(z),".pdf"),width = 5,height = 5)
                        
                        # T_abund1=c()
                        # for(i in 1:length(Non_Switching[,1])){
                        #   abund=sum(Non_Switching[i,2:p])
                        #   T_abund1=c(T_abund1,abund)
                        # }
                    
                        #########################################################################################    
                        V2 =transform(Non_Switching, T_abundanceNonsw=rowSums(Non_Switching[2:(S+1)]))
                        #V2=cbind(Non_Switching,T_abund1)
                        
                        con_abund_NonSw1 = c()
                        con_abund_NonSw1=unlist(Non_Switching[steps,2:(S+1)])
                        con_abund_NonSw =rbind(con_abund_NonSw,con_abund_NonSw1)
                        connect_abundance=append(connect_abundance,V2[steps,(S+2)])
                        
                        ## Computing Nestedness, Stability, Connectance and Modularity
                        
                        Y=t(replicate(S,pop))  
                        NestA = c()
                        Nest_int= nested(int_mat,method = "NODF",rescale = FALSE,normalised = TRUE)
                        NestA= c(NestA,Nest_int)


                        
                        connectance_A=c()
                        Init_connectance=networklevel(int_mat,index = "connectance")
                        connectance_A = c(connectance_A,Init_connectance)
                        
                        ##################################################################
                        # Elimination switching population dynamics plots over the time step for every connectance value.
                        # Plots saved in a folder.
                        ##################################################################
                        
                        elimination_switch=Elimination_switch_Euler(0, pop, 0.01,100,int_mat, str_mat,NestA,connectance_A,carry_cap,int_growth)
                       ########################################################################################
                         # Eliminate_Switch= cbind(elimination_switch[1],elimination_switch[2:p])
                        # Eliminate_Switch_melted=melt(Eliminate_Switch,id.vars = "X1")
                        # ggplot(data=Eliminate_Switch_melted,aes(x=X1,y=value, group=variable))+
                        #   #geom_point(aes(shape=variable),size=1)+
                        #   #scale_shape_manual(values = c(1,16,18))+
                        #   geom_line(aes(linetype=variable),size=1)+
                        #   #scale_linetype_manual(values = c("solid","longdash","dotted"))+
                        #   theme_few()+
                        #   labs(y="Productivity",x="Time",aspect.ratio = 2)+
                        #   theme(legend.position = "None",aspect.ratio = 1)+
                        #   theme(axis.title=element_text(size=12,face="bold"),
                        #         axis.text =element_text(face = "bold",size=12))
                        # ggsave(paste0("Elimdynamics_",toString(z),".pdf"),width = 5,height = 5)
                        # 
                        
                       ######################################################################################## 
                        
                        
                        con_abund_elim = c()
                        con_abund_elim=unlist(elimination_switch[steps,2:(S+1)])
                        con_abund_elimination =rbind(con_abund_elimination,con_abund_elim)
                        
                        
                        # computing species total abundance and appending on the population dynamics data frame.
                       
                        W_elimination=transform(elimination_switch, T_abundanceElim=rowSums(elimination_switch[2:(S+1)]))
                        
                        connect_abund_elim_switch=append(connect_abund_elim_switch,W_elimination[steps,(S+4)])
                        connect_nest_elim_switch= append(connect_nest_elim_switch,mean(W_elimination[,(S+2)]))
                        
                        # Computing Nestedness, Stability, Connectance and Modularity
                        
                        Y=t(replicate(S,pop))
                        Nest_opt = c()
                        Nestopt= nested(int_mat,method = "NODF",rescale = FALSE,normalised = TRUE)
                        Nest_opt= c(Nest_opt,Nestopt)


                        connectance_opt=c()
                        Init_opt_connectance=networklevel(int_mat,index = "connectance")
                        connectance_opt = c(connectance_opt,Init_opt_connectance)
                
                        
                        # Optimization switching population dynamics plots over the time step for every connectance value.
                        # Plots saved in a folder. 
                        
                        optimization_switch=Optimization_switch_Euler(0, pop, 0.01,100,int_mat, str_mat,NestA,connectance_A,carry_cap,int_growth) 
                        ###########################################################################################
                         # Opt_Switch= cbind(optimization_switch[1],optimization_switch[2:p])
                        # Opt_Switch_melted=melt(Opt_Switch,id.vars = "X1")
                        # ggplot(data=Opt_Switch_melted,aes(x=X1,y=value, group=variable))+
                        #   #geom_point(aes(shape=variable),size=1)+
                        #   #scale_shape_manual(values = c(1,16,18))+
                        #   geom_line(aes(linetype=variable),size=1)+
                        #   #scale_linetype_manual(values = c("solid","longdash","dotted"))+
                        #   theme_few()+
                        #   labs(y="Productivity",x="Time",aspect.ratio = 2)+
                        #   theme(legend.position = "None",aspect.ratio = 1)+
                        #   theme(axis.title=element_text(size=12,face="bold"),
                        #         axis.text =element_text(face = "bold",size=12))
                        # ggsave(paste0("Optdynamics_",toString(z),".pdf"),width = 5,height = 5)
                      
                        ###########################################################################################  
                        
                        con_abund_opt = c()
                        con_abund_opt=unlist(optimization_switch[steps,2:(S+1)])
                        con_abund_optimization =rbind(con_abund_optimization,con_abund_opt)
                        
                        
                        # computing species total abundance and appending on the population dynamics data frame.
                        
                        W_optimization=transform(optimization_switch, T_abundanceOpt=rowSums(optimization_switch[2:(S+1)]))
                        
                        connect_abund_opt_switch=append(connect_abund_opt_switch,W_optimization[steps,(S+4)])
                        connect_nest_opt_switch= append(connect_nest_opt_switch,mean(W_optimization[,(S+2)]))
                        
                       
                        # Computing the maximum real eigenivalue of the community matrix to measure community resilience.
                    
                        IM=int_mat
                        SM = str_mat
                        VNonSW=c(do.call("cbind",Non_Switching[steps,2:(S+1)]))
                        V_elim_sw=c(do.call("cbind",elimination_switch[steps,2:(S+1)]))
                        V_opt_sw=c(do.call("cbind",optimization_switch[steps,2:(S+1)]))
                      
                      
                        Jac_NonSW=jacobian(LVM_stability,VNonSW)
                        Jac_elimination_Switch=jacobian(LVM_stability,V_elim_sw)
                        Jac_optimization_Switch=jacobian(LVM_stability,V_opt_sw)
                        eigen_NonSW=eigen(Jac_NonSW)$values
                        eigen_elim=eigen(Jac_elimination_Switch)$values
                        eigen_opt=eigen(Jac_optimization_Switch)$values
                        Max_eigen_NonSW=data.frame(max(Re(eigen_NonSW)))
                        Max_eigen_elim =data.frame(max(Re(eigen_elim)))
                        Max_eigen_opt =data.frame(max(Re(eigen_opt)))
                        max_eigen_NonSw=append(max_eigen_NonSw,Max_eigen_NonSW)
                        max_eigen_elim=append(max_eigen_elim,Max_eigen_elim)
                        max_eigen_opt=append(max_eigen_opt,Max_eigen_opt)
                       
                        elim_max_eigen=c()
                        elim_Max_eigen=c()
                        opt_max_eigen=c()
                        opt_Max_eigen=c()
                        Max_eigen_NonSWT=c()
                        max_eigen_NonSwT=c()
                        
                        # Commputing stability over time for different values of connectance
                        # NOn- switching stability (Compute the maximum lead eigenvalue)
                        
                        for(nsw in 1:length(Non_Switching[,1])){
                          VNonSWT=c(do.call("cbind",Non_Switching[nsw,2:(S+1)]))
                          
                          Jac_NonSWT=jacobian(LVM_stability,VNonSWT)
                          eigen_NonSWT=eigen(Jac_NonSWT)$values
                          Max_eigen_NonSWT=data.frame(max(Re(eigen_NonSWT)))
                          max_eigen_NonSwT=append(max_eigen_NonSwT,Max_eigen_NonSWT)
                        }
                      
                        # Elimination switching stability
                        
                        for(elim_sw in 1:length(elimination_switch[,1])){
                          V_elim=c(do.call("cbind",elimination_switch[elim_sw,2:(S+1)]))
                          
                          Jac_elim=jacobian(LVM_stability,V_elim)
                          elim_eigen=eigen(Jac_elim)$values
                          elim_Max_eigen=data.frame(max(Re(elim_eigen)))
                          elim_max_eigen=append(elim_max_eigen,elim_Max_eigen)
                        }
                      
                        # Optimization switching stability
                        
                        for(opt_sw in 1:length(optimization_switch[,1])){
                          V_opt=c(do.call("cbind",optimization_switch[opt_sw,2:(S+1)]))
                          
                          Jac_opt=jacobian(LVM_stability,V_opt)
                          opt_eigen=eigen(Jac_opt)$values
                          opt_Max_eigen=data.frame(max(Re(opt_eigen)))
                          opt_max_eigen=append(opt_max_eigen,opt_Max_eigen)
                        }
                      
                      
                        # Stability values for the switching forms stored in a dataframe
                        
                        Stabb_time= data.frame(cbind(elimination_switch[1],unlist(max_eigen_NonSwT),unlist(elim_max_eigen),unlist(opt_max_eigen)))
                        #######################################################################################
                        Stab_time=plyr::rename(Stabb_time,c(X1="time",unlist.max_eigen_NonSwT.="max_eigen_NonSwT",unlist.elim_max_eigen.="elim_max_eigen",unlist.opt_max_eigen.="opt_max_eigen"))
                        Stab_time_melted=melt(Stab_time,id.vars = "time")
                        # ggplot(data=Stab_time_melted,aes(x=time,y=value, group=variable))+
                        #   #geom_point(aes(shape=variable),size=1)+
                        #   #scale_shape_manual(values = c(1,16,18))+
                        #   geom_line(aes(linetype=variable),size=1)+
                        #   scale_linetype_manual(values = c("solid","longdash","dotted"))+
                        #   scale_linetype_manual(values=c("solid","dashed", "dotted"))+
                        #   geom_hline(aes(yintercept=0),size=1,colour='red')+
                        #   theme_few()+
                        #   labs(y=expression("Re"~(lambda)),x="Time",aspect.ratio = 1)+
                        #   theme(legend.position = "None",aspect.ratio = 1)+
                        #   theme(axis.title=element_text(size=10,face="bold"),
                        #         axis.text =element_text(face = "bold",colour = "black",size=10))
                        # ggsave(paste0("Stabtime_",toString(c(S,SD,z)),".pdf"),width = 5,height = 5)
                        #######################################################################################
                      
            
                        }# End of Connectance Loop
            
            
            # Coonectance, total abundance, Stability, Nestedness, and modularity values for each connectance value saved in one dataframe
            
              Abund_connectance1=do.call(rbind, Map(data.frame,Connectance=con_seq, Abundance_Nonswitch=rowSums(con_abund_NonSw),
                                                    Abund_elim_switch=rowSums(con_abund_elimination),Abund_opt_switch=rowSums(con_abund_optimization),
                                                    Nested_Nonswitch=connect_nested,Nest_elim_switch=connect_nest_elim_switch,Nest_opt_switch=connect_nest_opt_switch,
                                                    Stability_NonSW=max_eigen_NonSw,Stability_Elim_sw=max_eigen_elim,Stability_Opt_sw=max_eigen_opt))
              
              Abund_connectance=sort.data.frame(Abund_connectance1)
             
              Abund_connectance$SD = SD
              Abundsigma_datalist[[SD]]=Abund_connectance
              Abundsigma_data=bind_rows(Abundsigma_datalist)
  
              #############################################################################################
            
              # Connectance versus Total Abundance|Nestedness|Stability plots
              #Productivity Plot
              ############################################################################################# 
              
             Productivity= cbind(Abund_connectance[1],Abund_connectance[2:4])
              Productivity_melted=melt(Productivity,id.vars = "Connectance")
              ggplot(data=Productivity_melted,aes(x=Connectance,y=value, group=variable))+
                geom_point(aes(shape=variable,color=variable),size=2.5)+
                #scale_shape_manual(values = c(1,2,5))+
                geom_smooth(method = "loess",span=0.65,se=F,aes(color=variable),size=0.75)+
                scale_color_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                                     labels=c("Non-switch","Elimination", "Optimization"))+
                scale_shape_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                                     labels=c("Non-switch","Elimination", "Optimization"))+
                theme_few()+
                labs(y="Productivity",x="Connectance",aspect.ratio = 1)+
                theme(legend.position = c(0.8,0.85), legend.title=element_blank(),aspect.ratio = 1)+
                theme(axis.title=element_text(face="plain"),
                      axis.text =element_text(face = "plain",colour = "black"))
              ggsave(paste0("Productivity_",toString(c(S,SD)),".pdf"),width = 7.5,height = 7.5)
              ggsave(paste0("Productivity_",toString(c(S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
            
              #Nestedness Plot
              
              Nestedness= cbind(Abund_connectance[1],Abund_connectance[5:7])
              Nestedness_melted=melt(Nestedness,id.vars = "Connectance")
              ggplot(data=Nestedness_melted,aes(x=Connectance,y=value, group=variable))+
                geom_point(aes(shape=variable,color=variable),size=2.5)+
                #scale_shape_manual(values = c(1,2,5))+
                geom_smooth(method = "loess",span=0.65,se=F,aes(color=variable),size=0.75)+
                #geom_smooth(method = "loess",se=F,col='black',aes(linetype=variable),size=0.5)+
                scale_color_discrete(breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                                     labels=c("Non-switch","Elimination", "Optimization"))+
                scale_shape_discrete(breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                                     labels=c("Non-switch","Elimination", "Optimization"))+
                theme_few()+
                labs(y="Nestedness",x="Connectance",aspect.ratio = 1)+
                theme(legend.position = c(0.15,0.85),legend.title = element_blank(),aspect.ratio = 1)+
                theme(axis.title=element_text(face="plain"),
                      axis.text =element_text(face = "plain",colour = "black"))
              
              # Saving the plots to Folder (pdf & Tiff files )
              
              ggsave(paste0("Nestedness_",toString(c(S,SD)),".pdf"),width = 7.5,height = 7.5)
              ggsave(paste0("Nestedness_",toString(c(S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
            
              # Community Resilience .
            
              Local_Stab= cbind(Abund_connectance[1],Abund_connectance[8:10])
              Local_Stab_melted=melt(Local_Stab,id.vars = "Connectance")
              ggplot(data=Local_Stab_melted,aes(x=Connectance,y=value, group=variable))+
                geom_hline(aes(yintercept=0),size=0.5,colour='red',linetype="dotted")+
                geom_point(aes(shape=variable,color=variable),size=2.5)+
                #scale_shape_manual(values = c(1,2,5))+
                geom_smooth(method = "loess",se=F,aes(color=variable),size=0.75)+
                #geom_smooth(method = "loess",se=F,col='black',aes(linetype=variable),size=0.5)+
                scale_color_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                                     labels=c("Non-switch","Elimination", "Optimization"))+
                scale_shape_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                                     labels=c("Non-switch","Elimination", "Optimization"))+
                theme_few()+
                labs(y=expression("Re"~(lambda)),x="Connectance",aspect.ratio = 1)+
                theme(legend.position = c(0.85,0.15),legend.title = element_blank(), aspect.ratio = 1)+
                theme(axis.title=element_text(face="plain"),
                      axis.text =element_text(face = "plain",colour = "black"))
              
              # Saving the plots to Folder (pdf & Tiff files )
              
              ggsave(paste0("LocalStability_",toString(c(S,SD)),".pdf"),width = 7.5,height = 7.5)
              ggsave(paste0("LocalStability_",toString(c(S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
              
              
              ############################################################################################
              
              # Species competitiveness computation. For each connectance values stored in a dataframe.
              #Competitiveness computed by taking the sum of the row values of comunity matrix for species competition from species and
              #   column sum for species competion effect to other species.
              ############################################################################################
              #Combinining Nonswitching, elimination and optimization average population abundances at for each connectance value           
              
              Compet_Abunda= do.call(rbind, Map(data.frame, row_wise=rowSums(str_mat), column_wise=colSums(str_mat),NonSw_abundance=colMeans(con_abund_NonSw),
                                               Nswitch_elim_abundance=colMeans(con_abund_elimination),Nswitch_opt_abundance=colMeans(con_abund_optimization) ))
            
              Sorted_Compet_Abunda_Col=Compet_Abunda[order(Compet_Abunda[,"column_wise"]),]
              Sorted_Compet_Abunda_row=Compet_Abunda[order(Compet_Abunda[,"row_wise"]),]
              
             
              
              #Inputing the row and column competitiveness for different standard deviation values at each
              
              Sorted_Compet_Abunda_Col$SD=SD
              Compet_sigma_datalist1[[SD]]=Sorted_Compet_Abunda_Col
              
              Sorted_Compet_Abunda_row$SD=SD
              Compet_sigma_datalist2[[SD]]=Sorted_Compet_Abunda_row
              
              #########################################################################           
              ## Species competitiveness and Pressure plots
     
              ######################################################################### 
              
              Compe_Row= cbind(Sorted_Compet_Abunda_row[1],Sorted_Compet_Abunda_row[3:5])
              Compe_Row_melted=melt(Compe_Row,id.vars = "row_wise")
              ####################################################
              # ggplot(data=Compe_Row_melted,aes(x=row_wise,y=value, group=variable))+
              #   geom_point(aes(shape=variable), size=2.5)+
              #   #scale_shape_manual(values = c(1,2,5))+
              #   geom_smooth(method = "loess",se=F,col='black',aes(linetype=variable),size=0.5)+
              #   scale_linetype_manual(values=c("solid","dashed", "dotted"))+
              #   theme_few()+
              #   labs(y="Abundance",x="Competitiveness",aspect.ratio = 1)+
              #   theme(legend.position = "None",aspect.ratio = 1)+
              #   theme(axis.title=element_text(face="plain"),
              #         axis.text =element_text(face = "plain",colour = "black"))
              # ### saving the plots uniquely. Each plot for corresponding SD and S
              # 
              # ggsave(paste0("CompeRow_",toString(c(S,SD)),".pdf"),width = 7.5,height = 7.5)
              # ggsave(paste0("CompeRow_",toString(c(S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
              # 
              #####################################################################
            
              Compe_Column= cbind(Sorted_Compet_Abunda_row[2],Sorted_Compet_Abunda_row[3:5])
              Compe_Column_melted=melt(Compe_Column,id.vars = "column_wise")
              ##############################################################################
              # ggplot(data=Compe_Column_melted,aes(x=column_wise,y=value, group=variable))+
              #   geom_point(aes(shape=variable), size =2.5)+
              #   #scale_shape_manual(values = c(1,2,5))+
              #   geom_smooth(method = "loess",se=F,col='black',aes(linetype=variable), size=0.5)+
              #   scale_linetype_manual(values=c("solid","dashed", "dotted"))+
              #   theme_few()+
              #   labs(y="Abundance",x="pressure",aspect.ratio = 1)+
              #   theme(legend.position = "None",aspect.ratio = 1)+
              #   theme(axis.title=element_text(face="plain"),
              #         axis.text =element_text(face = "plain",colour = "black"))
              # ggsave(paste0("CompeColumn_",toString(c(S,SD)),".pdf"),width = 7.5,height = 7.5)
              # ggsave(paste0("CompeColumn_",toString(c(S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
              #ggsave("Compe_Column.pdf",width = 5,height = 5)

              ##############################################################################
                           

              ############################################################################################              
              ### Species diversity measure using rank abundance curve
           
              ############################################################################################   
            
              conn_per_abund_NOnSW=data.frame(t(con_abund_NonSw))
              conn_per_abund_SW=data.frame(t(con_abund_elimination))
              conn_per_abund_Opt=data.frame(t(con_abund_optimization))
              
              Species_rankConnect=list()
              
            # Rank abundance curves. And ranking the switch based on non switching ranking
            
              for (C in 1:length(con_seq)) {
              
                rakingNonswitch=order(-conn_per_abund_NOnSW[,C])
                sortedNonswitch=conn_per_abund_NOnSW[,C][rakingNonswitch]
                rankingSwitch=order(-conn_per_abund_SW[,C])
                ranking_Opt_Switch=order(-conn_per_abund_Opt[,C])
                sortedElimSwitch=conn_per_abund_SW[,C][rankingSwitch]
                sortedOptSwitch=conn_per_abund_Opt[,C][ranking_Opt_Switch]
              
              
            
              Species_rank= data.frame(cbind(speciesrank=rank(-sortedNonswitch),sortedNonswitch,sortedElimSwitch,sortedOptSwitch))
              
              Species_rank$C=C
              Species_rankConnect[[C]]=Species_rank
              Species_rank_melted=melt(Species_rank,id.vars = c("speciesrank","C"))
              
              ###############################################           
              ##Rank abundance plot         
              ###################################################################################################       
             
              # ggplot(data=Species_rank_melted,aes(x=speciesrank,y=value, group=variable))+
             
              #    geom_point(aes(shape=variable),size=2.5)+
              #    #scale_shape_manual(values = c(1,2,5))+
             #    geom_line(aes(linetype=variable),size=0.5)+
             #    scale_linetype_manual(values = c("solid","dashed", "dotted"))+
             #    theme_few()+
             #    labs(y="Abundance",x="Species rank")+
             #    theme(legend.position = "None",aspect.ratio = 1)+
             #    theme(axis.title=element_text(face="plain"),
             #          axis.text =element_text(face = "plain",colour = "black"))
             #  ggsave(paste0("ranking_",toString(c(S,SD,C)),".pdf"),width = 7.5,height = 7.5)
             #  ggsave(paste0("ranking_",toString(c(S,SD,C)),".tiff"),width = 7.5,height = 7.5,units = 'in', dpi=300)
              
              ###################################################################################################            
            
              }
            
              Species_rankndata=bind_rows(Species_rankConnect)
            
              Species_rankndata$SD=SD
            
              Species_rankdatalist[[SD]]=Species_rankndata
      
              } # End of standard deviation(SD) loop
      
      ######################################################################
      # Combining the dataframes produced by each value of S uniquely into one dataframe
      # Abundsigma_data - species abundances  
      # Compet_sigma_data1- competitveness along the column
      # Compet_sigma_data2- competitveness along the row
      ######################################################################
      
      Species_rankdata=bind_rows(Species_rankdatalist)
      Species_rankdata$S=S
      Overall_Speciesranking[[S]]=Species_rankdata
      
      Abundsigma_data=bind_rows(Abundsigma_datalist)
      Abundsigma_data$S=S
      CombinedAbund_datalist[[S]]=Abundsigma_data

      Compet_sigma_data1=bind_rows(Compet_sigma_datalist1)
      Compet_sigma_data1$S=S
      CombinedCompet_datalist1[[S]]=Compet_sigma_data1
      
      Compet_sigma_data2=bind_rows(Compet_sigma_datalist2)
      Compet_sigma_data2$S=S
      CombinedCompet_datalist2[[S]]=Compet_sigma_data2
      
      } # End of species abundance (S) loop.


#########################################################################################
# Combining the dataframes into one at end of whole loop
# SpeciesAbunda_data - species abundances  
# SpeciesCompet_data1- competitveness along the column
# SpeciesCompet_data2- competitveness along the row
#Overall_CombinedSpecies_rank - Species ranking 
#########################################################################################

SpeciesAbunda_data=bind_rows(CombinedAbund_datalist)
SpeciesCompet_data1=bind_rows(CombinedCompet_datalist1)
SpeciesCompet_data2=bind_rows(CombinedCompet_datalist2)
Overall_COmbinedSpecies_rank=bind_rows(Overall_Speciesranking)

########################################################################################


SpeciesCommunity=SpeciesAbunda_data
SpeciesCommunity$SD[SpeciesCommunity$SD==1] =std_dvtn[1] 
SpeciesCommunity$SD[SpeciesCommunity$SD==2] =std_dvtn[2]
SpeciesCommunity$SD[SpeciesCommunity$SD==3] =std_dvtn[3]

# Computing May's Criteria (Sigma* Sqrt(SC))

SpeciesCommunity$MaysRule <- with(SpeciesCommunity, SD* sqrt(S*Connectance))


#Productivity Vs May's Criterion

MayProductivity= cbind(log(SpeciesCommunity[13]),SpeciesCommunity[2:4])
MayProductivity_melted=melt(MayProductivity,id.vars = "MaysRule")
pdtn=ggplot(data=MayProductivity_melted,aes(x=MaysRule,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  #scale_shape_manual(values = c(1,2,5))+
  ylim(3,20) +
  # stat_smooth_func(geom = "text",method = "lm",hjust=-3,parse=T)+
  geom_smooth(method = "loess",se=F,aes(color=variable),size=0.75)+
  #geom_smooth(method = "loess",se=F,col='black',aes(linetype=variable), size= 0.75)+
  scale_color_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Productivity",x=expression("ln"(~sigma*sqrt("SC"))),aspect.ratio = 1)+
  theme(legend.position = c(0.15,0.15),legend.title = element_blank(), aspect.ratio = 1)+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste("May'sRule_vs_Productivity.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("May'sRule_vs_Productivity.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


# Stability vs May's Criterion

MayLocal_Stab= cbind(log(SpeciesCommunity[13]),SpeciesCommunity[8:10])
MayLocal_Stab_melted=melt(MayLocal_Stab,id.vars = "MaysRule")
stable=ggplot(data=MayLocal_Stab_melted,aes(x=MaysRule,y=value, group=variable))+
  geom_hline(aes(yintercept=0),size=0.9,colour='red',linetype="dashed")+
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  ylim(-0.5, 0.1)+
  geom_smooth(method = "loess",span=0.65,se=F,aes(color=variable),size=0.75)+
  scale_color_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y=expression("Re"~(lambda)),x=expression("ln"(~sigma*sqrt("SC"))),aspect.ratio = 1)+
  theme(legend.position = c(0.15,0.9),legend.title = element_blank(),aspect.ratio = 1)+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste("May'sRule_vs_stability.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("May'sRule_vs_stability.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)

# Load Libraries

require(gridExtra)
require(cowplot)

#Combining Productivity and Stability plots

MayComparison<- plot_grid(stable+theme(legend.position = "None"),pdtn+theme(legend.position = "None"),
                          labels = "AUTO", align = 'h', label_size = 12, hjust = -0.5, vjust = -0.5,scale=c(1.,0.96))+
  theme(plot.margin = unit(c(1,-0.2,-4.5,0.5), "cm")) 
legend_A <- get_legend(stable+theme(legend.position = "top"))
plot_grid(MayComparison,legend_A,ncol=1, rel_heights = c(2, 2))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste("maycriteria.pdf"), width = 7.5,height = 7.5)
ggsave(paste0("maycriteria.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


# Nestedness vs May's criterion

MayNestedness= cbind(SpeciesCommunity[13],SpeciesCommunity[5:7])
MayNestedness_melted=melt(MayNestedness,id.vars = "MaysRule")
ggplot(data=MayNestedness_melted,aes(x=MaysRule,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  geom_smooth(method = "loess",span=0.6,se=F,aes(color=variable),size=0.75)+
  scale_color_discrete(breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Nestedness",x=expression(~sigma*sqrt("SC")),aspect.ratio = 1)+
  theme(legend.position = c(0.85,0.15),legend.title = element_blank(), aspect.ratio = 1)+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste0("May'sRule_vs_Nestedness.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("May'sRule_vs_Nestedness.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


#################################################################################################
#Combined Community plots for different values of sigma and S (number of species)
# Community plots comparison under different values of  standard deviation and species population
# S along the rows and SD( for community standard deviation) along the columns on the partitions
################################################################################################

# productivity versus connectance plots

CommunityProductivity= cbind(SpeciesCommunity[1],SpeciesCommunity[2:4],SpeciesCommunity[11:12])
CommunityProductivity_melted=melt(CommunityProductivity,id.vars = c("Connectance","SD","S"))
ggplot(data=CommunityProductivity_melted,aes(x=Connectance,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  facet_grid(S~SD, scales = "free",labeller = "label_both")+  ## partitioning S-rows and SD-columns
  ## line of best fit using LOESS 
  geom_smooth(method = "loess",se=F,aes(color=variable),size=0.75)+
  scale_color_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Productivity",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = "bottom",legend.title = element_blank(),aspect.ratio = 1)+
  theme(axis.title=element_text(face="bold"),axis.text =element_text(face = "plain",colour = "black"))

##Saving plots in a folder (pdf & Tiff)

ggsave(paste("productivity.pdf"), width = 7.5,height = 7.5)
ggsave(paste0("Productivity.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300) 


### Community stability versus connectance

CommunityLocal_Stab= cbind(SpeciesCommunity[1],SpeciesCommunity[8:10],SpeciesCommunity[11:12])
CommunityLocal_Stab_melted=melt(CommunityLocal_Stab,id.vars = c("Connectance","SD","S"))
ggplot(data=CommunityLocal_Stab_melted,aes(x=Connectance,y=value, group=variable))+
  geom_hline(aes(yintercept=0),size=0.5,colour='red',linetype="dashed")+  ## dashed line at y=0. 
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  facet_grid(S~SD, scales = "free",labeller = "label_both")+
  ## line of best fit using LOESS 
  geom_smooth(method = "loess",se=F,aes(color=variable),size=0.75)+
  scale_color_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y=expression("Re"~(lambda)),x="Connectance",aspect.ratio = 1)+
  theme(legend.position = "bottom",legend.title = element_blank(),aspect.ratio = 1)+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste("LocalStability.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("Localstability.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


### Nestedness versus connectance

CommunityNestedness= cbind(SpeciesCommunity[1],SpeciesCommunity[5:7],SpeciesCommunity[11:12])
CommunityNestedness_melted=melt(CommunityNestedness,id.vars = c("Connectance","SD","S"))
ggplot(data=CommunityNestedness_melted,aes(x=Connectance,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  facet_grid(S~SD, scales="free",labeller = "label_both")+
  geom_smooth(method = "loess",se=F,aes(color=variable), size=0.75)+
  scale_color_discrete(breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Nestedness",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = "bottom",legend.title = element_blank(),aspect.ratio = 1)+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste0("Nestedness.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("Nestedness.tiff"), width = 7.5, height = 7.5, units = 'in', dpi = 300)


#Competitiveness

SpeciesCCompetitiveness=SpeciesCompet_data2        ##row_wise
SpeciesCCompetitiveness$SD[SpeciesCCompetitiveness$SD==1] =std_dvtn[1]
SpeciesCCompetitiveness$SD[SpeciesCCompetitiveness$SD==2] =std_dvtn[2]
SpeciesCCompetitiveness$SD[SpeciesCCompetitiveness$SD==3] =std_dvtn[3]

SpeciesCCompetitiveness_Row= cbind(SpeciesCCompetitiveness[1],SpeciesCCompetitiveness[3:7])
SpeciesCCompetitiveness_melted=melt(SpeciesCCompetitiveness_Row,id.vars = c("row_wise","SD","S"))
ggplot(data=SpeciesCCompetitiveness_melted,aes(x=row_wise,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  geom_smooth(method = "loess",se=F,aes(color=variable), size=0.75)+
  facet_wrap(S~SD, scales = 'free', nrow = 3)+
  scale_color_discrete(breaks=c("NonSw_abundance","Nswitch_elim_abundance","Nswitch_opt_abundance"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("NonSw_abundance","Nswitch_elim_abundance","Nswitch_opt_abundance"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Abundance",x="Competitiveness",aspect.ratio = 1)+
  theme(legend.position = "bottom",legend.title = element_blank(),aspect.ratio = 1)+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste0("Competitiveness.tiff"), width = 10, height = 10, units = 'in', dpi = 300)
ggsave("Competitiveness.pdf",width = 10,height = 10)


SpeciesCCompetitiveness_Colmn=SpeciesCompet_data1        ##Column_wise
SpeciesCCompetitiveness_Colmn$SD[SpeciesCCompetitiveness_Colmn$SD==1] =std_dvtn[1]
SpeciesCCompetitiveness_Colmn$SD[SpeciesCCompetitiveness_Colmn$SD==2] =std_dvtn[2]
SpeciesCCompetitiveness_Colmn$SD[SpeciesCCompetitiveness_Colmn$SD==3] =std_dvtn[3]

SpeciesCCompetitiveness_Column= cbind(SpeciesCCompetitiveness_Colmn[2],SpeciesCCompetitiveness_Colmn[3:7])
SpeciesCCompetitiveness_Colmn_melted=melt(SpeciesCCompetitiveness_Column,id.vars = c("column_wise","SD","S"))
ggplot(data=SpeciesCCompetitiveness_Colmn_melted,aes(x=column_wise,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable), size=2.5)+
  facet_wrap(S~SD, scales = 'free', nrow = 3)+
  geom_smooth(method = "loess",se=F,aes(color=variable), size=0.75)+
  scale_color_discrete(breaks=c("NonSw_abundance","Nswitch_elim_abundance","Nswitch_opt_abundance"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  scale_shape_discrete(breaks=c("NonSw_abundance","Nswitch_elim_abundance","Nswitch_opt_abundance"),
                       labels=c("Non-switch","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Abundance",x="Pressure",aspect.ratio = 1)+
  theme(legend.position = "bottom",legend.title = element_blank(),aspect.ratio = 1)+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))

#Saving plots in a folder (pdf & Tiff)

ggsave(paste0("Pressure.tiff"), width = 10, height = 10, units = 'in', dpi = 300)
ggsave("Pressure.pdf",width = 10,height = 10)


##########################################################

Overall_COmbinedSpecies_rank$C[Overall_COmbinedSpecies_rank$C==which(con_seq==min(con_seq),arr.ind = T)]=con_seq[which(con_seq==min(con_seq),arr.ind = T)]
Overall_COmbinedSpecies_rank$C[Overall_COmbinedSpecies_rank$C==which(con_seq==median(con_seq),arr.ind = T)]=con_seq[which(con_seq==median(con_seq),arr.ind = T)]
Overall_COmbinedSpecies_rank$C[Overall_COmbinedSpecies_rank$C==which(con_seq==max(con_seq),arr.ind = T)]=con_seq[which(con_seq==max(con_seq),arr.ind = T)]

Overall_COmbinedSpecies_rank$SD[Overall_COmbinedSpecies_rank$SD==1]=std_dvtn[1]
Overall_COmbinedSpecies_rank$SD[Overall_COmbinedSpecies_rank$SD==2]=std_dvtn[2]
Overall_COmbinedSpecies_rank$SD[Overall_COmbinedSpecies_rank$SD==3]=std_dvtn[3]

Overall_COmbinedSpecies_rankmelt=melt(Overall_COmbinedSpecies_rank,id.vars = c("speciesrank","C","SD","S"))

# Species ranking in 20 species community

ggplot(subset(Overall_COmbinedSpecies_rankmelt,C%in%c(0.1,0.5,0.9) & S%in%20),aes(x=speciesrank,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2.5)+
  facet_grid(SD~S+C, scales = "free", labeller="label_both")+
  geom_line(aes(color=variable),size=0.75)+
  scale_shape_discrete(breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                       labels=c("Non-switch","Elimination","Optimization"))+
  scale_color_discrete(breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                       labels=c("Non-switch","Elimination","Optimization"))+
  theme_few()+
  labs(y="Abundance",x="Species rank")+
  theme(legend.position = "bottom",legend.title=element_blank(),aspect.ratio = 1)+
     theme(axis.title=element_text(face="plain"),
           axis.text =element_text(face = "plain",colour = "black"))
   ggsave("Species ranking20.pdf",width = 10,height = 10)
   ggsave(paste0("Species ranking20.tiff"), width = 8, height = 8, units = 'in', dpi = 300)
   
   
# Species ranking in 50 species community
   
ggplot(subset(Overall_COmbinedSpecies_rankmelt,C%in%c(0.1,0.5,0.9) & S%in%50),aes(x=speciesrank,y=value, group=variable))+
     geom_point(aes(shape=variable,color=variable),size=2.5)+
     #scale_shape_manual(values = c(1,2,5))+
     facet_grid(SD~S+C, scales = "free", labeller="label_both")+
     geom_line(aes(color=variable),size=0.75)+
     scale_shape_discrete(breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                          labels=c("Non-switch","Elimination","Optimization"))+
     scale_color_discrete(breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                          labels=c("Non-switch","Elimination","Optimization"))+
     theme_few()+
     labs(y="Abundance",x="Species rank")+
     theme(legend.position = "bottom",legend.title = element_blank(),aspect.ratio = 1)+
     theme(axis.title=element_text(face="plain"),
           axis.text =element_text(face = "plain",colour = "black"))
   ggsave("Species ranking50.pdf",width = 10,height = 10)
   ggsave(paste0("Species ranking50.tiff"), width = 8, height = 8, units = 'in', dpi = 300)

  
   
   ## Saving the data sets 
####################################
   save(SpeciesCommunity,SpeciesAbunda_data,SpeciesCompet_data1,SpeciesCompet_data2,Overall_COmbinedSpecies_rank,Abund_connectance,
        V2, W_optimization, W_elimination, file="Community.RData")
   
###############################################################
   ##AOB
   
   MayProductivity1= cbind(SpeciesCommunity[13],SpeciesCommunity[2:4])
   MayProductivity_melted1=melt(MayProductivity1,id.vars = "MaysRule")
   pdtn1=ggplot(data=MayProductivity_melted1,aes(x=MaysRule,y=value, group=variable))+
     geom_point(aes(shape=variable,color=variable),size=2.5)+
     #scale_shape_manual(values = c(1,2,5))+
     ylim(3,20) +
     # stat_smooth_func(geom = "text",method = "lm",hjust=-3,parse=T)+
     geom_smooth(method = "loess",se=F,aes(color=variable),size=0.75)+
     #geom_smooth(method = "loess",se=F,col='black',aes(linetype=variable), size= 0.75)+
     scale_color_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                          labels=c("Non-switch","Elimination", "Optimization"))+
     scale_shape_discrete(breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                          labels=c("Non-switch","Elimination", "Optimization"))+
     theme_few()+
     labs(y="Productivity",x=expression(~sigma*sqrt("SC")),aspect.ratio = 1)+
     theme(legend.position = c(0.15,0.15),legend.title = element_blank(), aspect.ratio = 1)+
     theme(axis.title=element_text(face="plain"),
           axis.text =element_text(face = "plain",colour = "black"))
   
   #Saving plots in a folder (pdf & Tiff)
   
   ggsave(paste("Maysrule_vs_Productivity.pdf"),width = 7.5,height = 7.5)
   ggsave(paste0("Maysrule_vs_Productivity.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
   
   
   # Stability vs May's Criterion
   
   MayLocal_Stab1= cbind(SpeciesCommunity[13],SpeciesCommunity[8:10])
   MayLocal_Stab_melted1=melt(MayLocal_Stab1,id.vars = "MaysRule")
   stable1=ggplot(data=MayLocal_Stab_melted1,aes(x=MaysRule,y=value, group=variable))+
     geom_hline(aes(yintercept=0),size=0.9,colour='red',linetype="dashed")+
     geom_point(aes(shape=variable,color=variable),size=2.5)+
     ylim(-0.5, 0.1)+
     geom_smooth(method = "loess",span=0.65,se=F,aes(color=variable),size=0.75)+
     scale_color_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                          labels=c("Non-switch","Elimination", "Optimization"))+
     scale_shape_discrete(breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                          labels=c("Non-switch","Elimination", "Optimization"))+
     theme_few()+
     labs(y=expression("Re"~(lambda)),x=expression(~sigma*sqrt("SC")),aspect.ratio = 1)+
     theme(legend.position = c(0.15,0.9),legend.title = element_blank(),aspect.ratio = 1)+
     theme(axis.title=element_text(face="plain"),
           axis.text =element_text(face = "plain",colour = "black"))
   
   #Saving plots in a folder (pdf & Tiff)
   
   ggsave(paste("Maysrule_vs_stability.pdf"),width = 7.5,height = 7.5)
   ggsave(paste0("Maysrule_vs_stability.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
   
   #Combining Productivity and Stability plots
   
   MayComparison1<- plot_grid(stable1+theme(legend.position = "None"),pdtn1+theme(legend.position = "None"),
                             labels = "AUTO", align = 'h', label_size = 12, hjust = -0.5, vjust = -0.5,scale=c(1.,0.96))+
     theme(plot.margin = unit(c(1,-0.2,-4.5,0.5), "cm")) 
   legend_A <- get_legend(stable1+theme(legend.position = "top"))
   plot_grid(MayComparison1,legend_A,ncol=1, rel_heights = c(2, 2))
   
   #Saving plots in a folder (pdf & Tiff)
   
   ggsave(paste("mayscriteria.pdf"), width = 7.5,height = 7.5)
   ggsave(paste0("mayscriteria.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
   
   ####################
   ##
   ##
   ##/** END OF PART ONE 