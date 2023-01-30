# Vector containing the names of the future columns of the the summary table
names_col=c("arch","mean_off", "s1","s2","s3","s4","e12","e13","e14","e23","e24","e34","r12","r23","r34","h12","h13","h14","h23","h24","h34","N","K","f0000","f0001","f0010","f0011","f0100","f0101","f0110","f0111","f1000","f1001","f1010","f1011","f1100","f1101","f1110","f1111","n_rep","p_anc","p_A1","p_B1","p_A2","p_B2","p_parA","p_parB","p_hybA1","p_hybB1","p_all_inc","tL_A1","tL_B1","tL_A2","tL_B2","t_dmi1","t_dmi2","t_anc","t_A1","t_B1","t_A2","t_B2","t_parA","t_parB","t_hybA1","t_hybB1","t_all_inc","t_res_given_hybA1","t_res_given_hybA2","t_both","p_alive","p_alive_anc","p_alive_A1","p_alive_B1","p_alive_A2","p_alive_B2","p_alive_parA","p_alive_parB","p_alive_hybA1","p_alive_hybB1","p_alive_all_inc","tL_alive_A1","tL_alive_B1","tL_alive_A2","tL_alive_B2","t_alive_dmi1","t_alive_dmi2","t_alive_both","NL_A1","NL_B1","NL_A2","NL_B2","NL_alive_A1","NL_alive_B1","NL_alive_A2","NL_alive_B2","N_dmi1","N_dmi2","N_alive_dmi1","N_alive_dmi2","N_both","N_alive_both","t_diff","t_diff_alive","t_diff_alive_hybA1","t_diff_alive_hybB1","N_min_alive")

# Vectors with the architecture names
list_Arch=c("adj_ABAB","adj_ABBA","cro_AABB","cro_ABBA","nes_ABAB","nes_AABB","sin_DMI","other")
list_Arch_legend=c("Adj. ABAB","Adj. ABBA","Cro. AABB","Cro. ABBA","Nes. AABB","Nes. ABAB","Sin. DMI")

arch_label=c(`adj_ABAB`="Adj. ABAB",`adj_ABBA`="Adj. ABBA",`cro_AABB`="Cro. AABB",`cro_ABBA`="Cro. ABBA",`nes_ABAB`="Nes. ABAB",`nes_AABB`="Nes. AABB")

names_col_detail=c("mean_off", "s1","s2","s3","s4","e12","e13","e14","e23","e24","e34","r12","r23","r34","h12","h13","h14","h23","h24","h34","N","K","f0000","f0001","f0010","f0011","f0100","f0101","f0110","f0111","f1000","f1001","f1010","f1011","f1100","f1101","f1110","f1111","hap","alive","t_dmi1","t_dmi2","N_dmi1","N_dmi2","t_res","N_res")

# function to check that all files were properly written and not interupted in the middle (due to partial file transfers)
check_non_finish=function(root){
  list_files=list.files(root,pattern=root) # list all files in the folder 
  for (k in list_files){ 

    # read the rest of the file. each line si a different iterations and eahc lines contains 8 elements: the allele that fixed (0 or 1) and the time of fixation of the alleles for all four loci.
    data=read.table(paste(c(root,k),sep="",collapse="/"),skip=1,sep='')
    if(min(data[,c(2,4,6,8)])==0){print(c(k,sum(rowSums(data[,c(2,4,6,8)]==0))))
    }
  }
}


# main function that reads the different output files form the C++ program and generate a data frame with all the probabilities and mean time
add_dataframe=function(root,dataframe,truncature=0 ,debug=FALSE){
  list_files=list.files(root,pattern=root)# list all files in the folder 
  if(debug){print(root)}
  for (k in list_files){ 
    # print(k)
    # first read only the first line of the file to extract the parameter values
    header=readLines(paste(c(root,k),sep="",collapse="/"),n = 1) 
    # store them temporary in a vector
    temp=strsplit(unlist(strsplit(header,"=")),",")
    vec=c()
    #print(temp)
    for (i in 3:22){
      vec=c(vec,temp[[i]][1])
    }
    vec=c(vec,strsplit(temp[[24]]," ")[[1]][1])
    vec=c(vec,strsplit(temp[[25]]," ")[[1]][1])
    vec=c(vec,strsplit(temp[[25]]," ")[[1]][4:19])
    vec=as.double(vec)
    
    # read the rest of the file. each line is a different iterations and each lines contains 14 elements. The first 8 elements corresponds to the allele that fixed (0 or 1) and the time of fixation of the alleles for all four loci. The 9th element is the maximum population size, the 10th the smallest population size, and the last 4 to the population size when the polymorphism is lost at the given locus.
    data=read.table(paste(c(root,k),sep="",collapse="/"),skip=1,sep='')
    # check that we only read "truncature" first simulations to have always the same precision
    if (truncature>0){
      if (dim(data)[1]>=truncature){
        data=data[1:truncature,]
      } else {
        #print(k)
        next
      }
    }
    # recreate the fixed haplotype
    hap= 2^3*data[,1]+2^2*data[,3]+2*data[,5]+data[,7]+1
    hapcount=hist(hap,plot = FALSE,breaks = seq(0,16))$counts
    hapcount_alive=hist(hap[data$V9>0],plot = FALSE,breaks = seq(0,16))$counts
    # add counts of fixation of 16 possibles haplotypes
    fixtime=colMeans(data)[c(2,4,6,8)]
    fixtime_alive=colMeans(data[data$V9>0,])[c(2,4,6,8)]
    
    # add average time of fixation of each haplotype as well as the mean square time to fixation (simplest way to have easily access to variance)
    t_out=rep(NA,16)
    t_out_alive=rep(NA,16)
    for (i in (1:16)){
      t_out[i]=sum(apply(data[hap==i,1:8],1,max))
      t_out_alive[i]=sum(apply(data[hap==i&data$V9>0,1:8],1,max))
    }
    
    # add time of resolution of the first and second DMIs; we sued the value of the epistatis is parameters to find the matching pairs of incompatible alleles.
    t_dmi=c()
    t_dmi_alive=c()
    N_dmi=c()
    N_dmi_alive=c()
  
    
    # add resolution of the DMI conditioned on the fixation of a hybrid haplotype. There is possibility depending fitness scheme. 
    # last element save the time of resolution of both DMIs, without any conditioning on the outcome.
    t_dmi_ghs=c(0,0,0,0,0,0)
    t_res_alive=0
    N_res=0
    N_res_alive=0
    
    N_min=mean(data$V10[data$V9>0])
    arch=-1
    # e14 and e23 != 0 (ie nes AABB or nes ABAB)
    if(vec[8]!=0 & vec[9]!=0 ){
      t_dmi=c(t_dmi,mean(apply(data[,c(2,8)],1,min))) 
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(2,8)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(11,14)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(11,14)],1,min)))
      t_dmi=c(t_dmi,mean(apply(data[,c(4,6)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(4,6)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(12,13)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(12,13)],1,min)))
      
      
      t_dmi_ghs[1]=mean(apply(rbind(apply(data[hap==4,c(2,8)],1,min),apply(data[hap==4,c(4,6)],1,min)),2,max))
      t_dmi_ghs[2]=mean(apply(rbind(apply(data[hap==6,c(2,8)],1,min),apply(data[hap==6,c(4,6)],1,min)),2,max))
      t_dmi_ghs[5]=mean(apply(rbind(apply(data[hap==11,c(2,8)],1,min),apply(data[hap==11,c(4,6)],1,min)),2,max))
      t_dmi_ghs[6]=mean(apply(rbind(apply(data[hap==13,c(2,8)],1,min),apply(data[hap==13,c(4,6)],1,min)),2,max))
      t_dmi_ghs[7]=mean(apply(rbind(apply(data[,c(2,8)],1,min),apply(data[,c(4,6)],1,min)),2,max))
      t_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(2,8)],1,min),apply(data[data$V9>0,c(4,6)],1,min)),2,max))
      N_res=mean(apply(rbind(apply(data[,c(11,14)],1,min),apply(data[,c(12,13)],1,min)),2,max))
      N_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(11,14)],1,min),apply(data[data$V9>0,c(12,13)],1,min)),2,max))
      t_diff=c(mean(abs(apply(data[,c(2,8)],1,min)-apply(data[,c(4,6)],1,min))),mean(abs(apply(data[data$V9>0,c(2,8)],1,min)-apply(data[data$V9>0,c(4,6)],1,min))))
      # parent: hap6:0101 and hap11:1010 ie nes ABAB 
      if ( vec[33]>0& vec[28]>0 & vec[33]+vec[28]==vec[21]){
        arch=list_Arch[5]
        to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[2], hapcount[3], hapcount[5], hapcount[11],hapcount[6],hapcount[13],hapcount[4],sum(hapcount[c(7,8,10,12,14,15,16)]))
        to_add=c(to_add,fixtime[c(1,4,3,2)])
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[2], (t_out/hapcount)[3], (t_out/hapcount)[5], (t_out/hapcount)[11],(t_out/hapcount)[6],(t_out/hapcount)[13],(t_out/hapcount)[4]),sum(t_out[c(7,8,10,12,14,15,16)])/sum(hapcount[c(7,8,10,12,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(6,1,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[2], hapcount_alive[3], hapcount_alive[5], hapcount_alive[11],hapcount_alive[6],hapcount_alive[13],hapcount_alive[4],sum(hapcount_alive[c(7,8,10,12,14,15,16)])))
        to_add=c(to_add,fixtime_alive[c(1,4,3,2)])
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[c(11,14,13,12)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[c(11,14,13,12)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,mean(abs(apply(data[hap==13&data$V9>0,c(2,8)],1,min)-apply(data[hap==13&data$V9>0,c(4,6)],1,min))))
        to_add=c(to_add,mean(abs(apply(data[hap==4&data$V9>0,c(2,8)],1,min)-apply(data[hap==4&data$V9>0,c(4,6)],1,min))))
        
      } 
      # parent: hap4:0011 and hap13:1100 ie nes AABB
      if (vec[35]>0& vec[26]>0 & vec[35]+vec[26]==vec[21]){
        
        arch=list_Arch[6]
        to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[2], hapcount[5], hapcount[3], hapcount[13],hapcount[4],hapcount[11],hapcount[6],sum(hapcount[c(7,8,10,12,14,15,16)]))
        to_add=c(to_add,fixtime[c(1,4,2,3)])
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[2], (t_out/hapcount)[5], (t_out/hapcount)[3], (t_out/hapcount)[13],(t_out/hapcount)[4],(t_out/hapcount)[11],(t_out/hapcount)[6]),sum(t_out[c(7,8,10,12,14,15,16)])/sum(hapcount[c(7,8,10,12,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(5,2,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[2], hapcount_alive[5], hapcount_alive[3], hapcount_alive[13],hapcount_alive[4],hapcount_alive[11],hapcount_alive[6],sum(hapcount_alive[c(7,8,10,12,14,15,16)])))
        to_add=c(to_add,fixtime_alive[c(1,4,2,3)])
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[c(11,14,12,13)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[c(11,14,12,13)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,mean(abs(apply(data[hap==11&data$V9>0,c(2,8)],1,min)-apply(data[hap==11&data$V9>0,c(4,6)],1,min))))
        to_add=c(to_add,mean(abs(apply(data[hap==6&data$V9>0,c(2,8)],1,min)-apply(data[hap==6&data$V9>0,c(4,6)],1,min))))
        
      }
    }
    # e13 and e24 != 0 (ie croAABB or croABBA)
    if(vec[7]!=0 & vec[10]!=0){
      t_dmi=c(t_dmi,mean(apply(data[,c(2,6)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(2,6)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(11,13)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(11,13)],1,min)))
      t_dmi=c(t_dmi,mean(apply(data[,c(4,8)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(4,8)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(12,14)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(12,14)],1,min)))
      t_dmi_ghs[1]=mean(apply(rbind(apply(data[hap==4,c(2,6)],1,min),apply(data[hap==4,c(4,8)],1,min)),2,max))
      t_dmi_ghs[3]=mean(apply(rbind(apply(data[hap==7,c(2,6)],1,min),apply(data[hap==7,c(4,8)],1,min)),2,max))
      t_dmi_ghs[4]=mean(apply(rbind(apply(data[hap==10,c(2,6)],1,min),apply(data[hap==10,c(4,8)],1,min)),2,max))
      t_dmi_ghs[6]=mean(apply(rbind(apply(data[hap==13,c(2,6)],1,min),apply(data[hap==13,c(4,8)],1,min)),2,max))
      t_dmi_ghs[7]=mean(apply(rbind(apply(data[,c(2,6)],1,min),apply(data[,c(4,8)],1,min)),2,max))
      t_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(2,6)],1,min),apply(data[data$V9>0,c(4,8)],1,min)),2,max))
      N_res=mean(apply(rbind(apply(data[,c(11,13)],1,min),apply(data[,c(12,14)],1,min)),2,max))
      N_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(11,13)],1,min),apply(data[data$V9>0,c(12,14)],1,min)),2,max))
      t_diff=c(mean(abs(apply(data[,c(2,6)],1,min)-apply(data[,c(4,8)],1,min))),mean(abs(apply(data[data$V9>0,c(2,6)],1,min)-apply(data[data$V9>0,c(4,8)],1,min))))
      # parent: hap7:0110 and hap10:1001 ie croABBA
      if ( vec[32]>0& vec[29]>0 & vec[32]+vec[29]==vec[21]){
        arch=list_Arch[4]
        to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[3], hapcount[2], hapcount[5], hapcount[10],hapcount[7],hapcount[13],hapcount[4],sum(hapcount[c(6,8,11,12,14,15,16)]))
        to_add=c(to_add,fixtime[c(1,3,4,2)])
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[3], (t_out/hapcount)[2], (t_out/hapcount)[5], (t_out/hapcount)[10],(t_out/hapcount)[7],(t_out/hapcount)[13],(t_out/hapcount)[4]),sum(t_out[c(6,8,11,12,14,15,16)])/sum(hapcount[c(6,8,11,12,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(6,1,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[3], hapcount_alive[2], hapcount_alive[5], hapcount_alive[10],hapcount_alive[7],hapcount_alive[13],hapcount_alive[4],sum(hapcount_alive[c(6,8,11,12,14,15,16)])))
        to_add=c(to_add,fixtime_alive[c(1,3,4,2)])
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[c(11,13,14,12)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[c(11,13,14,12)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,mean(abs(apply(data[hap==10&data$V9>0,c(2,6)],1,min)-apply(data[hap==10&data$V9>0,c(4,8)],1,min))))
        to_add=c(to_add,mean(abs(apply(data[hap==7&data$V9>0,c(2,6)],1,min)-apply(data[hap==7&data$V9>0,c(4,8)],1,min))))
      } 
      # parent: hap4:0011 and hap13:1100 ie croAABB
      if (vec[35]>0& vec[26]>0 & vec[35]+vec[26]==vec[21]){
        arch=list_Arch[3]
        to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[3], hapcount[5], hapcount[2], hapcount[13],hapcount[4],hapcount[10],hapcount[7],sum(hapcount[c(6,8,11,12,14,15,16)]))
        to_add=c(to_add,fixtime[c(1,3,2,4)])
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[3], (t_out/hapcount)[5], (t_out/hapcount)[2], (t_out/hapcount)[13],(t_out/hapcount)[4],(t_out/hapcount)[10],(t_out/hapcount)[7]),sum(t_out[c(6,8,11,12,14,15,16)])/sum(hapcount[c(6,8,11,12,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(4,3,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[3], hapcount_alive[5], hapcount_alive[2], hapcount_alive[13],hapcount_alive[4],hapcount_alive[10],hapcount_alive[7],sum(hapcount_alive[c(6,8,11,12,14,15,16)])))
        to_add=c(to_add,fixtime_alive[c(1,3,2,4)])
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[c(11,13,12,14)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[c(11,13,12,14)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,mean(abs(apply(data[hap==13&data$V9>0,c(2,6)],1,min)-apply(data[hap==13&data$V9>0,c(4,8)],1,min))))
        to_add=c(to_add,mean(abs(apply(data[hap==4&data$V9>0,c(2,6)],1,min)-apply(data[hap==4&data$V9>0,c(4,8)],1,min))))
      }
    }
    # e12 and e34 != 0 (ie adjABAB or adjABBA)
    if(vec[6]!=0 & vec[11]!=0){
      t_dmi=c(t_dmi,mean(apply(data[,c(2,4)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(2,4)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(11,12)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(11,12)],1,min)))
      t_dmi=c(t_dmi,mean(apply(data[,c(6,8)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(6,8)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(13,14)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(13,14)],1,min)))
      t_dmi_ghs[2]=mean(apply(rbind(apply(data[hap==6,c(2,4)],1,min),apply(data[hap==6,c(6,8)],1,min)),2,max))
      t_dmi_ghs[3]=mean(apply(rbind(apply(data[hap==7,c(2,4)],1,min),apply(data[hap==7,c(6,8)],1,min)),2,max))
      t_dmi_ghs[4]=mean(apply(rbind(apply(data[hap==10,c(2,4)],1,min),apply(data[hap==10,c(6,8)],1,min)),2,max))
      t_dmi_ghs[5]=mean(apply(rbind(apply(data[hap==11,c(2,4)],1,min),apply(data[hap==11,c(6,8)],1,min)),2,max))
      t_dmi_ghs[7]=mean(apply(rbind(apply(data[,c(2,4)],1,min),apply(data[,c(6,8)],1,min)),2,max))
      t_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(2,4)],1,min),apply(data[data$V9>0,c(6,8)],1,min)),2,max))
      N_res=mean(apply(rbind(apply(data[,c(11,12)],1,min),apply(data[,c(13,14)],1,min)),2,max))
      N_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(11,12)],1,min),apply(data[data$V9>0,c(13,14)],1,min)),2,max))
      t_diff=c(mean(abs(apply(data[,c(2,4)],1,min)-apply(data[,c(6,8)],1,min))),mean(abs(apply(data[data$V9>0,c(2,4)],1,min)-apply(data[data$V9>0,c(6,8)],1,min))))
      # parent: hap6:0101 and hap11:1010 ie adjABAB
      if ( vec[33]>0& vec[28]>0 & vec[33]+vec[28]==vec[21]){
        arch=list_Arch[1]
        # merge all the different data in a signle vectore
        to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[5], hapcount[3], hapcount[2], hapcount[11],hapcount[6],hapcount[10],hapcount[7],sum(hapcount[c(4,8,12,13,14,15,16)]))
        to_add=c(to_add,fixtime)
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[5], (t_out/hapcount)[3], (t_out/hapcount)[2], (t_out/hapcount)[11],(t_out/hapcount)[6],(t_out/hapcount)[10],(t_out/hapcount)[7]),sum(t_out[c(4,8,12,13,14,15,16)])/sum(hapcount[c(4,8,12,13,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(4,3,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[5], hapcount_alive[3], hapcount_alive[2], hapcount_alive[11],hapcount_alive[6],hapcount_alive[10],hapcount_alive[7],sum(hapcount_alive[c(4,8,12,13,14,15,16)])))
        to_add=c(to_add,fixtime_alive)
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[seq(11,14)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[seq(11,14)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,mean(abs(apply(data[hap==10&data$V9>0,c(2,4)],1,min)-apply(data[hap==10&data$V9>0,c(6,8)],1,min))))
        to_add=c(to_add,mean(abs(apply(data[hap==7&data$V9>0,c(2,4)],1,min)-apply(data[hap==7&data$V9>0,c(6,8)],1,min))))
      } 
      # parent: hap7:0110 and hap10:1001 ie adjABBA
      if (vec[32]>0& vec[29]>0 & vec[32]+vec[29]==vec[21]){
        arch=list_Arch[2]
        to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[5], hapcount[2], hapcount[3], hapcount[10],hapcount[7],hapcount[11],hapcount[6],sum(hapcount[c(4,8,12,13,14,15,16)]))
        to_add=c(to_add,fixtime[c(1,2,4,3)])
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[5], (t_out/hapcount)[2], (t_out/hapcount)[3], (t_out/hapcount)[10],(t_out/hapcount)[7],(t_out/hapcount)[11],(t_out/hapcount)[6]),sum(t_out[c(4,8,12,13,14,15,16)])/sum(hapcount[c(4,8,12,13,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(5,2,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[5], hapcount_alive[2], hapcount_alive[3], hapcount_alive[10],hapcount_alive[7],hapcount_alive[11],hapcount_alive[6],sum(hapcount_alive[c(4,8,12,13,14,15,16)])))
        to_add=c(to_add,fixtime_alive[c(1,2,4,3)])
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[c(11,12,14,13)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[c(11,12,14,13)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,mean(abs(apply(data[hap==11&data$V9>0,c(2,4)],1,min)-apply(data[hap==11&data$V9>0,c(6,8)],1,min))))
        to_add=c(to_add,mean(abs(apply(data[hap==6&data$V9>0,c(2,4)],1,min)-apply(data[hap==6&data$V9>0,c(6,8)],1,min))))
      }
    }
    # single locus case
    if((vec[6]!=0 & sum(vec[7:11])==0)|(sum(vec[6:11])==0&all(vec[2:3]!=0)&all(vec[4:5]==0))){
      arch=list_Arch[7]
      t_dmi=c(t_dmi,mean(apply(data[,c(2,4)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(2,4)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(11,12)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(11,12)],1,min)))
      t_dmi=c(t_dmi,NA)
      t_dmi_alive=c(t_dmi_alive,NA)
      N_dmi=c(N_dmi,NA)
      N_dmi_alive=c(N_dmi_alive,NA)
      t_dmi_ghs[2]=NA
      t_dmi_ghs[3]=NA
      t_dmi_ghs[4]=NA
      t_dmi_ghs[5]=NA
      t_dmi_ghs[7]=NA
      t_res_alive=NA
      N_res=NA
      N_res_alive=NA
      t_diff=c(NA,NA)
      
#parent 5 and 9

        to_add=c(length(hap),hapcount[1],NA,NA,NA,NA,hapcount[5],hapcount[9],NA,NA,hapcount[13])
        to_add=c(to_add,fixtime)
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[5], (t_out/hapcount)[3], (t_out/hapcount)[2], (t_out/hapcount)[11],(t_out/hapcount)[6],(t_out/hapcount)[10],(t_out/hapcount)[7]),sum(t_out[c(4,8,12,13,14,15,16)])/sum(hapcount[c(4,8,12,13,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(4,3,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],NA,NA,NA,NA, hapcount_alive[5], hapcount_alive[9], NA,NA,hapcount_alive[13]))
        to_add=c(to_add,fixtime_alive)
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[seq(11,14)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[seq(11,14)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,NA)
        to_add=c(to_add,NA)
       
    }
    # no dmis
    if(sum(vec[6:11])==0&all(vec[2:5]!=0)){
      t_dmi=c(t_dmi,mean(apply(data[,c(2,4)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(2,4)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(11,12)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(11,12)],1,min)))
      t_dmi=c(t_dmi,mean(apply(data[,c(6,8)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(6,8)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(13,14)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(13,14)],1,min)))
      t_dmi_ghs[2]=mean(apply(rbind(apply(data[hap==6,c(2,4)],1,min),apply(data[hap==6,c(6,8)],1,min)),2,max))
      t_dmi_ghs[3]=mean(apply(rbind(apply(data[hap==7,c(2,4)],1,min),apply(data[hap==7,c(6,8)],1,min)),2,max))
      t_dmi_ghs[4]=mean(apply(rbind(apply(data[hap==10,c(2,4)],1,min),apply(data[hap==10,c(6,8)],1,min)),2,max))
      t_dmi_ghs[5]=mean(apply(rbind(apply(data[hap==11,c(2,4)],1,min),apply(data[hap==11,c(6,8)],1,min)),2,max))
      t_dmi_ghs[7]=mean(apply(rbind(apply(data[,c(2,4)],1,min),apply(data[,c(6,8)],1,min)),2,max))
      t_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(2,4)],1,min),apply(data[data$V9>0,c(6,8)],1,min)),2,max))
      N_res=mean(apply(rbind(apply(data[,c(11,12)],1,min),apply(data[,c(13,14)],1,min)),2,max))
      N_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(11,12)],1,min),apply(data[data$V9>0,c(13,14)],1,min)),2,max))
      t_diff=c(NA,NA)
      #if ( vec[33]>0& vec[28]>0 & vec[33]+vec[28]==vec[21]){
        arch=list_Arch[1]
        # merge all the different data in a signle vectore
        to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[5], hapcount[3], hapcount[2], hapcount[11],hapcount[6],hapcount[10],hapcount[7],sum(hapcount[c(4,8,12,13,14,15,16)]))
        to_add=c(to_add,fixtime)
        to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[5], (t_out/hapcount)[3], (t_out/hapcount)[2], (t_out/hapcount)[11],(t_out/hapcount)[6],(t_out/hapcount)[10],(t_out/hapcount)[7]),sum(t_out[c(4,8,12,13,14,15,16)])/sum(hapcount[c(4,8,12,13,14,15,16)]))
        to_add=c(to_add,t_dmi)
        to_add=c(to_add,t_dmi_ghs[c(4,3,7)])
        to_add=c(to_add,sum(data$V9>0))
        to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[5], hapcount_alive[3], hapcount_alive[2], hapcount_alive[11],hapcount_alive[6],hapcount_alive[10],hapcount_alive[7],sum(hapcount_alive[c(4,8,12,13,14,15,16)])))
        to_add=c(to_add,fixtime_alive)
        to_add=c(to_add,t_dmi_alive)
        to_add=c(to_add,t_res_alive)
        to_add=c(to_add,colMeans(data)[seq(11,14)])
        to_add=c(to_add,colMeans(data[data$V9>0,])[seq(11,14)])
        to_add=c(to_add,N_dmi)
        to_add=c(to_add,N_dmi_alive)
        to_add=c(to_add,N_res)
        to_add=c(to_add,N_res_alive)
        to_add=c(to_add,t_diff)
        to_add=c(to_add,NA)
        to_add=c(to_add,NA)
      } 
    
    if(sum(vec[6:11])==0&all(vec[3:5]==0)){
      t_dmi=c(t_dmi,mean(apply(data[,c(2,4)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(2,4)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(11,12)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(11,12)],1,min)))
      t_dmi=c(t_dmi,mean(apply(data[,c(6,8)],1,min)))
      t_dmi_alive=c(t_dmi_alive,mean(apply(data[data$V9>0,c(6,8)],1,min)))
      N_dmi=c(N_dmi,mean(apply(data[,c(13,14)],1,min)))
      N_dmi_alive=c(N_dmi_alive,mean(apply(data[data$V9>0,c(13,14)],1,min)))
      t_dmi_ghs[2]=mean(apply(rbind(apply(data[hap==6,c(2,4)],1,min),apply(data[hap==6,c(6,8)],1,min)),2,max))
      t_dmi_ghs[3]=mean(apply(rbind(apply(data[hap==7,c(2,4)],1,min),apply(data[hap==7,c(6,8)],1,min)),2,max))
      t_dmi_ghs[4]=mean(apply(rbind(apply(data[hap==10,c(2,4)],1,min),apply(data[hap==10,c(6,8)],1,min)),2,max))
      t_dmi_ghs[5]=mean(apply(rbind(apply(data[hap==11,c(2,4)],1,min),apply(data[hap==11,c(6,8)],1,min)),2,max))
      t_dmi_ghs[7]=mean(apply(rbind(apply(data[,c(2,4)],1,min),apply(data[,c(6,8)],1,min)),2,max))
      t_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(2,4)],1,min),apply(data[data$V9>0,c(6,8)],1,min)),2,max))
      N_res=mean(apply(rbind(apply(data[,c(11,12)],1,min),apply(data[,c(13,14)],1,min)),2,max))
      N_res_alive=mean(apply(rbind(apply(data[data$V9>0,c(11,12)],1,min),apply(data[data$V9>0,c(13,14)],1,min)),2,max))
      
      #if ( vec[33]>0& vec[28]>0 & vec[33]+vec[28]==vec[21]){
      
      arch=list_Arch[8]
      # merge all the different data in a signle vectore
      to_add=c(length(hap),hapcount[1],hapcount[9], hapcount[5], hapcount[3], hapcount[2], hapcount[11],hapcount[6],hapcount[10],hapcount[7],sum(hapcount[c(4,8,12,13,14,15,16)]))
      to_add=c(to_add,fixtime)
      to_add=c(to_add, c((t_out/hapcount)[1],(t_out/hapcount)[9], (t_out/hapcount)[5], (t_out/hapcount)[3], (t_out/hapcount)[2], (t_out/hapcount)[11],(t_out/hapcount)[6],(t_out/hapcount)[10],(t_out/hapcount)[7]),sum(t_out[c(4,8,12,13,14,15,16)])/sum(hapcount[c(4,8,12,13,14,15,16)]))
      to_add=c(to_add,t_dmi)
      to_add=c(to_add,t_dmi_ghs[c(4,3,7)])
      to_add=c(to_add,sum(data$V9>0))
      to_add=c(to_add,c(hapcount_alive[1],hapcount_alive[9], hapcount_alive[5], hapcount_alive[3], hapcount_alive[2], hapcount_alive[11],hapcount_alive[6],hapcount_alive[10],hapcount_alive[7],sum(hapcount_alive[c(4,8,12,13,14,15,16)])))
      to_add=c(to_add,fixtime_alive)
      to_add=c(to_add,t_dmi_alive)
      to_add=c(to_add,t_res_alive)
      to_add=c(to_add,colMeans(data)[seq(11,14)])
      to_add=c(to_add,colMeans(data[data$V9>0,])[seq(11,14)])
      to_add=c(to_add,N_dmi)
      to_add=c(to_add,N_dmi_alive)
      to_add=c(to_add,N_res)
      to_add=c(to_add,N_res_alive)
      to_add=c(to_add,c(NA,NA,NA,NA))
    } 
    if (length(t_dmi)==1){t_dmi=c(t_dmi,NA)}
    vec=c(arch,vec)
    vec=c(vec,to_add)
    vec=c(vec,N_min)
    #---------------
    #add the vector to the dataframe
    dataframe=rbind(dataframe,vec)
  }
  return(dataframe)
}


# NOT SURE
add_dataframe_detail=function(root,dataframe){
  list_files=list.files(root,pattern=root) # list all files in the folder 
  for (k in list_files){ 
    # first read only the first line of the file to extratx the parameter values
    header=readLines(paste(c(root,k),sep="",collapse="/"),n = 1) 
    # store them temporary in a vector
    temp=strsplit(unlist(strsplit(header,"=")),",")
    vec=c()
    for (i in 3:22){
      vec=c(vec,temp[[i]][1])
    }
    vec=c(vec,strsplit(temp[[24]]," ")[[1]][1])
    vec=c(vec,strsplit(temp[[25]]," ")[[1]][1])
    vec=c(vec,strsplit(temp[[25]]," ")[[1]][4:19])
    vec=as.double(vec)
    
    # read the rest of the file. each line si a different iterations and eahc lines contains 8 elements: the allele that fixed (0 or 1) and the time of fixation of the alleles for all four loci.
    data=read.table(paste(c(root,k),sep="",collapse="/"),skip=1,sep='')
    
    vec2=t(matrix(vec,nrow=length(vec),ncol =1000))
    
    # recreate the fixed haplotype
    hap= 2^3*data[,1]+2^2*data[,3]+2*data[,5]+data[,7]+1
    vec2=cbind(vec2,hap)
    vec2=cbind(vec2,as.integer(data$V9>0))
    
    
    # add time of resolution of the first and second DMIs; we sued the value of the epistatissi parameters to find the matching pairs of incompatible alleles.
    t_dmi=c()
    N_dmi=c()
    if(vec[6]!=0){
      t_dmi=cbind(t_dmi,apply(data[,c(2,4)],1,min))
      N_dmi=cbind(N_dmi,apply(data[,c(11,12)],1,min))
    }
    if(vec[7]!=0){
      t_dmi=cbind(t_dmi,apply(data[,c(2,6)],1,min))
      N_dmi=cbind(N_dmi,apply(data[,c(11,13)],1,min))
    }
    if(vec[8]!=0){
      t_dmi=cbind(t_dmi,apply(data[,c(2,8)],1,min))
      N_dmi=cbind(N_dmi,apply(data[,c(11,14)],1,min))
    }
    if(vec[9]!=0){
      t_dmi=cbind(t_dmi,apply(data[,c(4,6)],1,min))
      N_dmi=cbind(N_dmi,apply(data[,c(12,13)],1,min))
    }
    if(vec[10]!=0){
      t_dmi=cbind(t_dmi,apply(data[,c(4,8)],1,min))
      N_dmi=cbind(N_dmi,apply(data[,c(12,14)],1,min))
    }
    if(vec[11]!=0){
      t_dmi=cbind(t_dmi,apply(data[,c(6,8)],1,min))
      N_dmi=cbind(N_dmi,apply(data[,c(13,14)],1,min))
    }
    if (length(t_dmi)==1){t_dmi=c(t_dmi,NA)}
    
    
    # add resolution of the dmi conditioned on the fixation of a hybrid hyplotype. There is possibility depending fitness scheme. 
    # last element save the tine of resolution of both DMIs, without any conditioning on the outcome.

    t_res=c()
    N_res==c()

    
    if(vec[8]!=0 & vec[9]!=0){
      t_res=apply(rbind(apply(data[,c(2,8)],1,min),apply(data[,c(4,6)],1,min)),2,max)
      N_res=apply(rbind(apply(data[,c(11,14)],1,min),apply(data[,c(12,13)],1,min)),2,max)
    }
    if(vec[7]!=0 & vec[10]!=0){
      t_res=apply(rbind(apply(data[,c(2,6)],1,min),apply(data[,c(4,8)],1,min)),2,max)
      N_res=apply(rbind(apply(data[,c(11,13)],1,min),apply(data[,c(12,14)],1,min)),2,max)
    }
    if(vec[6]!=0 & vec[11]!=0){
      t_res=apply(rbind(apply(data[,c(2,4)],1,min),apply(data[,c(6,8)],1,min)),2,max)
      N_res=apply(rbind(apply(data[,c(11,12)],1,min),apply(data[,c(13,14)],1,min)),2,max)
    }
    
    
    
    # merge all the different data in a signle vectore
    vec2=cbind(vec2,t_dmi)
    vec2=cbind(vec2,N_dmi)
    vec2=cbind(vec2,t_res)
    vec2=cbind(vec2,N_res)
    colnames(vec2)=names_col_detail
    #add the vecotr to the dataframe
    dataframe=rbind(dataframe,vec2)
  }
  return(dataframe)
}

# this function is used to only get the simulations corresponding to the "Adjacent ABAB" architecture
adj_ABAB=function(data=whole_dataset,e12=-0.2,e34=-0.2){return(data$e12==e12&data$e13==0&data$e14==0&data$e23==0&data$e24==0&data$e34==e34&data$f0101>0&data$f1010>0&data$f0101+data$f1010==data$N)}

# this function is used to only get the simulations corresponding to the "Adjacent ABBA" architecture
adj_ABBA=function(data=whole_dataset,e12=-0.2,e34=-0.2){return(data$e12==e12&data$e13==0&data$e14==0&data$e23==0&data$e24==0&data$e34==e34&data$f0110>0&data$f1001>0&data$f0110+data$f1001==data$N)}

# this function is used to only get the simulations corresponding to the "Crossed ABAB" architecture
cro_AABB=function(data=whole_dataset,e13=-0.2,e24=-0.2){return(data$e12==0&data$e13==e13&data$e14==0&data$e23==0&data$e24==e24&data$e34==0&data$f1100>0&data$f0011>0&data$f1100+data$f0011==data$N)}
# this function is used to only get the simulations corresponding to the "Nested AABB" architecture
nest_AABB=function(data=whole_dataset,e14=-0.2,e23=-0.2){return(data$e12==0&data$e13==0&data$e14==e14&data$e23==e23&data$e24==0&data$e34==0&data$f1100>0&data$f0011>0&data$f1100+data$f0011==data$N)}

# this function is used to only get the simulations corresponding to the "Nested ABAB" architecture
nest_ABAB=function(data=whole_dataset,e14=-0.2,e23=-0.2){return(data$e12==0&data$e13==0&data$e14==e14&data$e23==e23&data$e24==0&data$e34==0&data$f1010>0&data$f0101>0&data$f0101+data$f1010==data$N)}

# this function is used to only get the simulations corresponding to the "Crossed ABBA" architecture
cro_ABBA=function(data=whole_dataset,e13=-0.2,e24=-0.2){return(data$e12==0&data$e13==e13&data$e14==0&data$e23==0&data$e24==e24&data$e34==0&data$f0110>0&data$f1001>0&data$f0110+data$f1001==data$N)}

# this function is used to selection specific selection coefficient for all derived alleles.
sel=function(data=whole_dataset,s1=-.001,s2=-.001,s3=-.001,s4=-.001){return(data$s1==s1&data$s2==s2&data$s3==s3&data$s4==s4)}

# this function is used to simualtions where all loci are equidistant.
equidistant=function(data=whole_dataset,r=NA){if(is.na(r)){return(data$r12==data$r23&data$r34==data$r23)}else{return(data$r12==r&data$r23==r&data$r34==r)}} 
#


# function to double check that the architecture have been properly identified
check_arch=function(temp_data){
  if(dim(temp_data)[1]==0){
    return(TRUE)
  }
    if(temp_data$arch[1]=="adj_ABAB"){
      return(all(temp_data$e13==0)&&all(temp_data$e14==0)&&all(temp_data$e23==0)&&all(temp_data$e24==0)&&all(temp_data$f0101>0)&&all(temp_data$f1010>0))
    }
  if(temp_data$arch[1]=="adj_ABBA"){
    return(all(temp_data$e13==0)&&all(temp_data$e14==0)&&all(temp_data$e23==0)&&all(temp_data$e24==0)&&all(temp_data$f0110>0)&&all(temp_data$f1001>0))
  }
  if(temp_data$arch[1]=="cro_AABB"){
    return(all(temp_data$e12==0)&&all(temp_data$e14==0)&&all(temp_data$e23==0)&&all(temp_data$e34==0)&&all(temp_data$f1100>0)&&all(temp_data$f0011>0))
  }
  if(temp_data$arch[1]=="cro_ABBA"){
    return(all(temp_data$e12==0)&&all(temp_data$e14==0)&&all(temp_data$e23==0)&&all(temp_data$e34==0)&&all(temp_data$f1001>0)&&all(temp_data$f0110>0))
  }
  if(temp_data$arch[1]=="nes_ABAB"){
    return(all(temp_data$e12==0)&&all(temp_data$e13==0)&&all(temp_data$e24==0)&&all(temp_data$e34==0)&&all(temp_data$f1010>0)&&all(temp_data$f0101>0))
  }
  if(temp_data$arch[1]=="nes_AABB"){
    return(all(temp_data$e12==0)&&all(temp_data$e13==0)&&all(temp_data$e24==0)&&all(temp_data$e34==0)&&all(temp_data$f1100>0)&&all(temp_data$f0011>0))
  }
  if(temp_data$arch[1]=="sin_DMI"){
    return(all(temp_data$e13==0)&&all(temp_data$e14==0)&&all(temp_data$e23==0)&&all(temp_data$e34==0)&&all(temp_data$s4==0))
  }
  return(FALSE)
}

# function to only keep case with direct selection is equal to -0.001, K=10^6 and mean offspring number to 1.01
default_regime=function(data){return(data$s1==-0.001&data$s2==-0.001&data$s3==-0.001&data$s4==-0.001&data$mean_off==1.01&data$K==1000000)}

# function to make stack plot of all the possible outcome
stack_plot=function(temp_data,ymax=1,title=""){
  plot(-1,xlim=c(20,10000),ylim=c(0,ymax),log="x",ylab="Probability",xlab="Initial pop. size",main=title)
  polygon(c(temp_data$N,rev(temp_data$N)),c(0*temp_data$N,rev(temp_data$p_alive_hybA1/temp_data$p_alive)),col="yellow",border = "NA")
  polygon(c(temp_data$N,rev(temp_data$N)),c(temp_data$p_alive_hybA1/temp_data$p_alive,rev((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1)/temp_data$p_alive)),col="green",border = "NA")
  polygon(c(temp_data$N,rev(temp_data$N)),c((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1)/temp_data$p_alive,rev((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parB)/temp_data$p_alive)),col="blue",border = "NA")
  polygon(c(temp_data$N,rev(temp_data$N)),c((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parB)/temp_data$p_alive,rev((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parA+temp_data$p_alive_parB)/temp_data$p_alive)),col="red",border = "NA")
  polygon(c(temp_data$N,rev(temp_data$N)),c((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parB+temp_data$p_alive_parA)/temp_data$p_alive,rev((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parA+temp_data$p_alive_parB+temp_data$p_alive_B1+temp_data$p_alive_B2)/temp_data$p_alive)),col="cyan",border = "NA")
  polygon(c(temp_data$N,rev(temp_data$N)),c((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parB+temp_data$p_alive_parA+temp_data$p_alive_B1+temp_data$p_alive_B2)/temp_data$p_alive,rev((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parA+temp_data$p_alive_parB+temp_data$p_alive_B1+temp_data$p_alive_B2+temp_data$p_alive_A1+temp_data$p_alive_A2)/temp_data$p_alive)),col="pink",border = "NA")
  polygon(c(temp_data$N,rev(temp_data$N)),c((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parB+temp_data$p_alive_parA+temp_data$p_alive_B1+temp_data$p_alive_B2+temp_data$p_alive_A1+temp_data$p_alive_A2)/temp_data$p_alive,rev((temp_data$p_alive_hybA1+temp_data$p_alive_hybB1+temp_data$p_alive_parA+temp_data$p_alive_parB+temp_data$p_alive_B1+temp_data$p_alive_B2+temp_data$p_alive_A1+temp_data$p_alive_A2+temp_data$p_alive_anc)/temp_data$p_alive)),col="black",border = "NA")
}