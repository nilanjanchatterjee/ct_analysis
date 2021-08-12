
########################################################################################################
# function evaldesign
#
# this function evaluates the performance of a design for the single-season single-species occupancy model
# with constant probabilities of occupancy and detectability. the function generates simulated histories, 
# calculates the corresponding MLEs for psi and p and evaluates estimator performance
#
# e.g. myres<-evaldesign(psi=0.2,p=0.3,s=62,k=4,nits=10000,doprint=1,doplot=1) 
#
# 	function input
#		- psi: assumed value for the probability of occupancy
#		- p: assumed value for the probability of detection
#		- s: number of sites surveyed
#		- k: number of replicated surveys per site
# 		- nits: number of iterations in the simulation
#		- doprint: print results on screen
#		- doplot: plot the distribution of MLEs
# 	function output (myres$)
#		- dist: matrix containing the results of the simulation. Each row contains details for 
#                  each history type summarized by the sufficient statistics (SD,d) 
#	    	  dist[,1] contain SD (number of sites where the species was detected)
#	    	  dist[,2] contains d (total number of detections in the history)
#	    	  dist[,3] contains the probability of obtaining that history type (i.e. pair (SD,d))
#	    	  dist[,4] contains the estimate for the probability of occupancy
#	    	  dist[,5] contains the estimate for the probability of detection
#	      - biaspsi,varpsi,MSEpsi: occupancy estimator bias/variance/MSE (excl empty histories)
#	      - biasp,varp,MSEp: detectability estimator bias/variance/MSE (excl empty histories)
#		- covar: occupancy and detectability estimator covariance (excl empty histories)
#		- critA: sum of the MSEs (excl empty histories)
#		- critB: determinant of the MSE matrix (excl empty histories)
#		- biaspsi_B,varpsi_B,MSEpsi_B: occupancy estimator bias/variance/MSE (excl also boundary estimates) 
#		- ...
#    		- pempty: percentage of empty histories obtained 
#    		- pbound: percentage of histories obtained that produce boundary estimates (i.e. psi=1)
#
########################################################################################################

evaldesign <- function(mat, psi,p,s,k,nits=10000,doprint=TRUE,doplot=TRUE) {

	t1=proc.time()

	mydist<-matrix(NA,0,5)
	jj<-1

	for (ii in 1:nits)   
	{
	# generate a history
	Sp<-rbinom(1,s,psi)
	#hist<-rbind(matrix(rbinom(Sp*k,1,p),Sp,k),matrix(0,s-Sp,k))
  hist <- mat[sample(1:nrow(mat),s, replace = F),sample(ncol(mat), k, replace = F), drop = FALSE]
	# summarize history (sufficient statistics)
	SD<-sum(rowSums(hist)>0)
	d<-sum(hist)

	# count history
	if (jj==1)
	{
		mydist<-rbind(mydist,c(SD,d,1,NA,NA))
	} else 
	{
		temp<-which((mydist[1:nrow(mydist),1]==SD))
		temp2<-which((mydist[temp,2]==d))	
		if (length(temp2)==0)
		{
			mydist<-rbind(mydist,c(SD,d,1,NA,NA))
		} else 
		{ 
			#add one count to this {SD,d} history 
			mydist[temp[temp2],3]<-mydist[temp[temp2],3]+1
		}
	}	
	jj<-jj+1
	}		#end for nits


	# analyze each history type obtained in simulation
	for (ii in 1:nrow(mydist))
	{  
	SD<-mydist[ii,1]
	d<-mydist[ii,2]

	if (SD==0){			  #empty histories
		psihat<-NA
		phat<-NA
	} else 
	{	
            if (((s-SD)/s)<((1-d/(s*k))^k)) {	#boundary estimates
			psihat <- 1
			phat <- d/(s*k)
		} else
		{
			params = c(psi,p)   #init vals for optim
			fitted1 = optim(params,loglikf,SD=SD,d=d,s=s,k=k,lin=0)
			psihat <- 1/(1+exp(-fitted1$par[1]))      
			phat <- 1/(1+exp(-fitted1$par[2])) 
		}
	}

	# store estimates
	mydist[ii,4]<-psihat
	mydist[ii,5]<-phat   
	mydist[ii,3]<-mydist[ii,3]/nits
	}
	
	# MLE properties: distribution removing empty histories
	mydist2<-mydist[mydist[,1]!=0,]	
	mydist2[,3]<-mydist2[,3]/sum(mydist2[,3])
	mymeanpsi<-sum(mydist2[,3]*mydist2[,4],na.rm='true')
	mybiaspsi<-mymeanpsi-psi
	myvarpsi<-sum((mydist2[,4]-mymeanpsi)^2*mydist2[,3])
	myMSEpsi<-myvarpsi+mybiaspsi^2
	mymeanp<-sum(mydist2[,3]*mydist2[,5],na.rm='true')
	mybiasp<-mymeanp-p
	myvarp<-sum((mydist2[,5]-mymeanp)^2*mydist2[,3])
	myMSEp<-myvarp+mybiasp^2
	mycovar<-sum((mydist2[,5]-mymeanp)*(mydist2[,4]-mymeanpsi)*mydist2[,3])
	mycritA<-myMSEpsi+myMSEp
	mycritD<-myMSEpsi*myMSEp-mycovar^2

	# MLE properties: distribution removing also boundary estimates
	mydist3<-mydist2[mydist2[,4]!=1,]	
	mydist3[,3]<-mydist3[,3]/sum(mydist3[,3])
	mymeanpsi_B<-sum(mydist3[,3]*mydist3[,4],na.rm='true')
	mybiaspsi_B<-mymeanpsi_B-psi
	myvarpsi_B<-sum((mydist3[,4]-mymeanpsi_B)^2*mydist3[,3])
	myMSEpsi_B<-myvarpsi_B+mybiaspsi_B^2
	mymeanp_B<-sum(mydist3[,3]*mydist3[,5],na.rm='true')
	mybiasp_B<-mymeanp_B-p
	myvarp_B<-sum((mydist3[,5]-mymeanp_B)^2*mydist3[,3])
	myMSEp_B<-myvarp_B+mybiasp_B^2
	mycovar_B<-sum((mydist3[,5]-mymeanp_B)*(mydist3[,4]-mymeanpsi_B)*mydist3[,3])
	mycritA_B<-myMSEpsi_B+myMSEp_B
	mycritD_B<-myMSEpsi_B*myMSEp_B-mycovar_B^2

	pempty<-100*mydist[mydist[,1]==0,3]
	if (!length(pempty)) pempty=0
	pbound<-100*sum(mydist[mydist[,4]==1,3],na.rm='true')
	if (!length(pbound)) pbound=0

	# compute processing time
	t2=proc.time()

	# print in screen results
	if (doprint) {
	cat("\n--------------------------------------------------------------------------\n",sep = "")
	cat("Evaluation of design K = ",k," S = ",s, " (TS = ",s*k,")\n",sep = "")
	cat("--------------------------------------------------------------------------\n",sep = "")
	cat("estimator performance (excl empty histories)\n",sep = "")
	cat("psi: bias = ",sprintf("%+0.4f",mybiaspsi),"   var = ",sprintf("%+0.4f",myvarpsi),"   MSE = ",sprintf("%+0.4f",myMSEpsi),"\n",sep = "")
	cat("  p: bias = ",sprintf("%+0.4f",mybiasp),"   var = ",sprintf("%+0.4f",myvarp),"   MSE = ",sprintf("%+0.4f",myMSEp),"\n",sep = "")
	cat("    covar = ",sprintf("%+0.4f",mycovar)," critA = ",sprintf("%+0.4f",mycritA)," critD = ",sprintf("%+.3e",mycritD),"\n",sep = "")
	cat("    RMSE_occ = ",sprintf("%+0.4f",sqrt(myMSEpsi))," RMSE_occ_det = ",sprintf("%+0.4f",sqrt(myMSEpsi+myMSEp)) ,"\n",sep = "")
	cat("estimator performance (excl also histories leading to boundary estimates)\n",sep = "")
	cat("psi: bias = ",sprintf("%+0.4f",mybiaspsi_B),"   var = ",sprintf("%+0.4f",myvarpsi_B),"   MSE = ",sprintf("%+0.4f",myMSEpsi_B),"\n",sep = "")
	cat("  p: bias = ",sprintf("%+0.4f",mybiasp_B),"   var = ",sprintf("%+0.4f",myvarp_B),"   MSE = ",sprintf("%+0.4f",myMSEp_B),"\n",sep = "")
	cat("    covar = ",sprintf("%+0.4f",mycovar_B)," critA = ",sprintf("%+0.4f",mycritA_B)," critD = ",sprintf("%+.3e",mycritD_B),"\n",sep = "")
	cat(" empty histories = ",sprintf("%0.1f",pempty),"%\n",sep = "")
	cat(" boundary estimates = ",sprintf("%0.1f",pbound),"%\n",sep = "")
	cat("this took ", (t2-t1)[1],"seconds \n")
	cat("--------------------------------------------------------------------------\n\n",sep = "")
	}

	# write results 
	myres <- list(dist=mydist,biaspsi=mybiaspsi,varpsi=myvarpsi,MSEpsi=myMSEpsi,biasp=mybiasp,varp=myvarp,
			  MSEp=myMSEp,covar=mycovar,critA=mycritA,critD=mycritD,biaspsi_B=mybiaspsi_B,varpsi_B=myvarpsi_B,
			  MSEpsi_B=myMSEpsi_B,biasp_B=mybiasp_B,varp_B=myvarp_B,MSEp_B=myMSEp_B,covar_B=mycovar_B,
			  critA_B=mycritA_B,critD_B=mycritD_B,pempty=pempty,pbound=pbound)

	if (doplot) plotMLEdist(myres$dist,p,psi)

	return(myres)

}
########################################################################################################



########################################################################################################
# function loglikf 
#
# this function computes the likelihood function for the single-season single-species occupancy model
# with constant probabilities of occupancy and detectability, given a history summarized by (SD,d) 
#
# 	function input
#		- params: values of psi and p where the likelihood is evaluated
#		- SD: number of sites where the species was detected
#		- d: total number of detections in the history
#		- s: number of sites surveyed
#		- k: number of replicated surveys per site
# 		- lin: indicates whether the values of psi and p are given in linear of logit domain
# 	function output
#		- loglik: value of the likelihood function
#
########################################################################################################

loglikf <- function(params, SD, d, s, k, lin) {

    if (lin){
        psi         = params[1]   
        p           = params[2]  
    } else {
        psi         = 1/(1+exp(-params[1]))  
        p           = 1/(1+exp(-params[2]))
    }
    
    loglik = -(SD*log(psi)+d*log(p)+(k*SD-d)*log(1-p)+(s-SD)*log((1-psi)+psi*(1-p)^k))
}

########################################################################################################



########################################################################################################
# function plotMLEdist 
#
# this function displays the distribution of the MLEs obtained for the given design and values of 
# occupancy and detectability (single-season single-species occupancy model with constant probabilities)
#
########################################################################################################

plotMLEdist <- function(mydist, p, psi)
{        
	  mcol<-rev(rainbow(200))
	  mcol<-c(mcol[52:200],mcol[1:16])
	  range=cbind(1,length(mcol))	
        x <-mydist[,3]
	  y <-as.integer((range[2]-range[1])*(x-min(x))/(max(x)-min(x)))+1
        plot(mydist[,5],mydist[,4],col=mcol[y],xlim=cbind(0,1),ylim=cbind(0,1),xlab=expression(hat("p")),
		 ylab=expression(hat("psi")),pch=19)
        abline(v=p, col="lightgray")
        abline(h=psi, col="lightgray")
}
########################################################################################################

occmat <-matrix(data = NA, nrow = 7, ncol = 10)
for(i in 7:1)
{
  myres<-evaldesign(psi=0.3383,p=0.0691,s=400,k=5*i+5,nits=1000,doprint=0,doplot=0)
  occmat[i,1]=sqrt(myres$MSEpsi)
  occmat[i,2]= sqrt(myres$MSEpsi+myres$MSEp)
  
  myres1<-evaldesign(psi=0.3383,p=0.0691,s=300,k=5*i+5,nits=1000,doprint=0,doplot=0)
  occmat[i,3]=sqrt(myres1$MSEpsi)
  occmat[i,4]= sqrt(myres1$MSEpsi+myres$MSEp)
  
  myres2<-evaldesign(psi=0.3383,p=0.0691,s=200,k=5*i+5,nits=1000,doprint=0,doplot=0)
  occmat[i,5]=sqrt(myres2$MSEpsi)
  occmat[i,6]= sqrt(myres2$MSEpsi+myres$MSEp)
  
  myres3<-evaldesign(psi=0.3383,p=0.0691,s=100,k=5*i+5,nits=1000,doprint=0,doplot=0)
  occmat[i,7]=sqrt(myres3$MSEpsi)
  occmat[i,8]= sqrt(myres3$MSEpsi+myres$MSEp)
  
  myres4<-evaldesign(psi=0.3383,p=0.0691,s=50,k=5*i+5,nits=1000,doprint=0,doplot=0)
  occmat[i,9]=sqrt(myres4$MSEpsi)
  occmat[i,10]= sqrt(myres4$MSEpsi+myres$MSEp)
  
}
occmat


### Updated code for sub-sampling

library(readxl)
smb <-read_excel("D:/Papers/Occupancy/camera_occu_species.xlsx", sheet = "Sambar")
dim(smb)
mat <-smb[1:400, 4:29]