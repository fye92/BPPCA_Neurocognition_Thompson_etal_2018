#############################################################################################################
#############################################################################################################
## Code for:                                                                                               ##
## Thompson, W. K., Barch, D., Bjork, J., Gonzalez, R., Nagel, B., Nixon, S. J., & Luciana, M. (2018).     ## 
## The Structure of Cognition in 9 and 10 year-old Children and Associations with Problem Behaviors:       ##
## Findings from the ABCD Study’s Baseline Neurocognitive Battery. Developmental Cognitive Neuroscience.   ##
##                                                                                                         ##
## Wes Thompson (wes.stat@gmail.com)                                                                       ##
## Dec 22, 2018                                                                                            ##
#############################################################################################################
#############################################################################################################

####################
####################
## Load libraries ##
####################
####################

library(mvtnorm)
library(tableone)
library(parallel)
library(rstan)
library(loo)
library(gamm4)
library(Hmisc)
library(FactoMineR)
library(nFactors)
library(reshape2)
library(psych)
library(data.table)
library(mice)
library(abind)
library(cvTools)
library(modEvA)

####################################
####################################
## Load and manipulate nda17 data ##
####################################
####################################


## Read Rds file from DEAP (ABCD NDA version 1.1 release)
## scripts for creating the file "nda17_release1.1.Rds" can be found at: https://github.com/ABCD-STUDY/analysis-nda17

setwd('/Users/fengdanye/Documents/Yale Research/ABCD/ABCDStudyNDA_Release_3.0/rds/results')
nda17 = readRDS("nda3.0.Rds")

#########################################
#########################################
## Select & process data for analyses  ##
#########################################
#########################################

## Select demographics
ind_demog = c(which(names(nda17)=="age"),which(names(nda17)=="female"),which(names(nda17)=="race_ethnicity"),
			which(names(nda17)=="high.educ"),which(names(nda17)=="married"),which(names(nda17)=="household.income"))
names(nda17)[ind_demog]
summary(nda17[,ind_demog])

## Select nesting variables
ind_nest = c(which(names(nda17)=="abcd_site"),which(names(nda17)=="rel_family_id"));summary(nda17[,ind_nest])

nda17$abcd_site = as.character(nda17$abcd_site)
nda17$abcd_site[nda17$abcd_site=="site22"] = "site07"
nda17$abcd_site = factor(nda17$abcd_site)

## Select neuropsychological measures
ind_pea_ravlt = c(which(names(nda17)=="pea_ravlt_sd_trial_i_tc"),which(names(nda17)=="pea_ravlt_sd_trial_ii_tc"),
	which(names(nda17)=="pea_ravlt_sd_trial_iii_tc"),which(names(nda17)=="pea_ravlt_sd_trial_iv_tc"),
	which(names(nda17)=="pea_ravlt_sd_trial_v_tc")); names(nda17)[ind_pea_ravlt]; summary(nda17[,ind_pea_ravlt])
nda17$pea_ravlt_ld = apply(nda17[,ind_pea_ravlt],1,sum)

par(mfrow=c(1,2))
hist(nda17$pea_ravlt_ld)
hist(nda17$lmt_scr_perc_correct)

ind_np = c(which(names(nda17)=="nihtbx_picvocab_uncorrected"),which(names(nda17)=="nihtbx_flanker_uncorrected"),
	which(names(nda17)=="nihtbx_list_uncorrected"),which(names(nda17)=="nihtbx_cardsort_uncorrected"),
	which(names(nda17)=="nihtbx_pattern_uncorrected"),which(names(nda17)=="nihtbx_picture_uncorrected"),
	which(names(nda17)=="nihtbx_reading_uncorrected"),which(names(nda17)=="pea_ravlt_ld"),
	which(names(nda17)=="lmt_scr_perc_correct")); names(nda17)[ind_np]

## Select dependent measures 
ind_dv = c(which(names(nda17)=="cbcl_scr_syn_external_r"),
			which(names(nda17)=="cbcl_scr_syn_internal_r"),which(names(nda17)=="cbcl_scr_07_stress_r"))
names(nda17)[ind_dv]

######################
######################
## Rename variables ##
######################
######################

names(nda17)[which(names(nda17)=="age")] = "Age"; nda17$Age = nda17$Age/12
names(nda17)[which(names(nda17)=="female")] = "Female"
names(nda17)[which(names(nda17)=="race_ethnicity")] = "RaceEthnicity"
names(nda17)[which(names(nda17)=="high.educ")] = "HighestParentalEducation"
names(nda17)[which(names(nda17)=="married")] = "HouseholdMaritalStatus"
names(nda17)[which(names(nda17)=="household.income")] = "HouseholdIncome"
names(nda17)[which(names(nda17)=="abcd_site")] = "Site"
names(nda17)[which(names(nda17)=="rel_relationship")] = "Relationship"

names(nda17)[which(names(nda17)=="nihtbx_picvocab_uncorrected")] = "PicVocab"
names(nda17)[which(names(nda17)=="nihtbx_flanker_uncorrected")] = "Flanker"
names(nda17)[which(names(nda17)=="nihtbx_list_uncorrected")] = "List"
names(nda17)[which(names(nda17)=="nihtbx_cardsort_uncorrected")] = "CardSort"
names(nda17)[which(names(nda17)=="nihtbx_pattern_uncorrected")] = "Pattern"
names(nda17)[which(names(nda17)=="nihtbx_picture_uncorrected")] = "Picture"
names(nda17)[which(names(nda17)=="nihtbx_reading_uncorrected")] = "Reading"
names(nda17)[which(names(nda17)=="pea_ravlt_ld")] = "RAVLT"
names(nda17)[which(names(nda17)=="lmt_scr_perc_correct")] = "LMT"
names(nda17)[which(names(nda17)=="pea_wiscv_tss")] = "WISC-V"
names(nda17)[which(names(nda17)=="cbcl_scr_syn_external_r")] = "Externalizing"
names(nda17)[which(names(nda17)=="cbcl_scr_syn_internal_r")] = "Internalizing"
names(nda17)[which(names(nda17)=="cbcl_scr_07_stress_r")] = "Stress"

nda17$compl = "Incomplete"
nda17$compl[complete.cases(nda17[,c(ind_nest,ind_np)])] = "Complete"
table(nda17$compl)

## Create Table 1
vars1 <- c("Age", "Female", "RaceEthnicity", "HighestParentalEducation", "HouseholdMaritalStatus",
            "HouseholdIncome","Relationship","Site")
tab1 <- CreateTableOne(vars = vars1, data = nda17, strata = "compl")
tabAsStringMatrix <- print(tab1, printToggle = FALSE, noSpaces = TRUE)
tab1 = knitr::kable(tabAsStringMatrix)

## Create Table 2
vars2 <- c("PicVocab", "Flanker", "List", "CardSort", "Pattern",
            "Picture","Reading","RAVLT","WISC-V",
            "LMT","Externalizing","Internalizing","Stress")
tab2 <- CreateTableOne(vars = vars2, data = nda17,strata = "compl")
tabAsStringMatrix <- print(tab2, printToggle = FALSE, noSpaces = TRUE)
tab2 = knitr::kable(tabAsStringMatrix)

## Create SM Figures 1 & 2
pdf("figures_and_tables/fig_sm1.pdf")
par(mfrow=c(3,3))
for(p in 1:9){
	hist(scale(nda17[,ind_np[p]]), xlab="", freq=FALSE,main=names(nda17)[ind_np][p])
}
dev.off()

pdf("figures_and_tables/fig_sm2.pdf")
par(mfrow=c(2,2))
for(p in 1:3){
	hist(scale(nda17[,ind_dv[p]]), xlab="", freq=FALSE,main=names(nda17)[ind_dv][p])
}
dev.off()

###################################
###################################
## Subset variables for analyses ##
###################################
###################################

data = nda17[,c(1,which(names(nda17)=="eventname"),ind_nest,ind_demog,which(names(nda17)=="rel_relationship"),which(names(nda17)=="rel_group_id"),ind_np,ind_dv)]
data = data[which(data$eventname=='baseline_year_1_arm_1'),]
names(data);dim(data)
data$src_subject_id = as.character(data$src_subject_id)
names(data)[names(data)=="src_subject_id"] = "pid"
data$site_num = as.numeric(substr(data$Site,5,6))
data$fam_num = 0
ind=0
for(i in sort(unique(data$rel_family_id))){
	ind = ind+1
	data$fam_num[data$rel_family_id==i & !is.na(data$rel_family_id)] = ind
}
data = data[order(data$site_num,data$fam_num,data$rel_group_id),]

####################
####################
## Fit usual PCA  ## 
####################
####################

## PCA on unimputed data 
ind_Y = c(12:20); names(data)[ind_Y] # the nine neurocog measures
Y = as.matrix(scale(data[complete.cases(data[,c(ind_Y)]),ind_Y]))
ev = eigen(cor(Y))
ap = parallel(subject=nrow(Y),var=ncol(Y),rep=100,cent=.05)
nS = nScree(x=ev$values,aparallel=ap$eigen$qevpea)
plotnScree(nS)
ncomp = 3
#y.pca = psych::principal(Y, rotate="promax", nfactors=ncomp, scores=TRUE)
#y.pca$loadings
y.pca = psych::principal(Y, rotate="varimax", nfactors=ncomp, scores=TRUE)
y.pca$loadings

#########################################
#########################################
## Fit 9-variable PCA Model using Stan ##
#########################################
#########################################

###############
## Prep data ##
###############

ind_Y = 12:20; names(data)[ind_Y]
ind_demog = c(5:10); names(data)[ind_demog]
ind_nest = c(3,4); names(data)[ind_nest]

data1 = data[complete.cases(data[,c(ind_Y)]),]; names(data1); dim(data1)
site_num = rep(NA,length(data1$site_num))
for(s in 1:length(unique(data1$site_num))){
	site_num[data1$site_num == unique(data1$site_num)[s]] = s
}
data1$site_num = site_num; rm(site_num)
data1$id_fam = 0
data1$fam_size = 0
ind=0
for(s in 1:length(unique(data1$site_num))){
	data_s = data1[data1$site_num == s, ]
	for(f in 1:length(unique(data_s$rel_family_id))){
		data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = 
			sum(data_s$rel_family_id == unique(data_s$rel_family_id)[f])
		if(sum(data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]])>1){	
			ind=ind+1	
			data_s$id_fam[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = ind
		}	
	}
	data1[data1$site_num == s, ] = data_s
}
data1 = data1[order(data1$site_num,data1$id_fam),]

Site = data1$site_num
Fam = data1$id_fam
ns = length(unique(Site))
ns_s = rep(0,ns)
for(s in 1:ns){
	ns_s[s] = sum(Site==s)
}
nf_s = rep(0,ns)
for(s in 1:ns){
	nf_s[s] = length(unique(Fam[Site==s]))-1
}
nf = sum(nf_s)
Y = (as.matrix(scale(data1[,ind_Y]))); summary(Y)
N = nrow(Y)
P = ncol(Y)

####################
## Run stan model ##
####################

Nsamples = 1000
Nchains = 3

model_file = "bppca.stan"
smod = stan_model(model_file)

D_max = 5
sa.list = list()
log_lik.list = list()
looic.list = list()
for(d in 1:D_max){
	print(d)
	pca_data <- list(Y = Y, N = N, P = P, D = d, Fam =  Fam, Site = Site, ns = ns, nf = nf)
	set.seed(314)
	sa.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")
	#log_lik.list[[d]] <- extract_log_lik(sa.list[[d]])
	log_lik.list[[d]] <- extract(sa.list[[d]],"log_lik_marg")[[1]]
	looic.list[[d]] = loo(log_lik.list[[d]])
	save(sa.list,log_lik.list,looic.list,file="results/bppca_results.RData")
	print("###############################")
}	

#################################
## Model selection using LOOIC ##
#################################

load("results/bppca_results.RData")

looic.obj = compare(looic.list[[1]],looic.list[[2]],looic.list[[3]],looic.list[[4]],looic.list[[5]])
print(looic.obj)

d=3

print(sa.list[[d]], pars=c("Q"), probs=c(.025,.5,.975))
print(sa.list[[d]], pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"), probs=c(.025,.5,.975))

## Create SM Figure 3
jpeg("figures_and_tables/fig_sm3.jpeg")
traceplot(sa.list[[d]], pars=c("sigma_eps","sigma2_a","sigma2_b","sigma_c","sigma_d"))
dev.off()

#######################################
## Extract parameters from the model ##
#######################################

sa = sa.list[[d]]
Q = extract(sa,"Q",permuted=FALSE)
Theta = extract(sa,"Theta",permuted=FALSE)
Lambda = extract(sa,"Lambda",permuted=FALSE)
sigma2_a = extract(sa,"sigma2_a",permuted=FALSE)
sigma2_b = extract(sa,"sigma2_b",permuted=FALSE)
sigma_c = extract(sa,"sigma_c",permuted=FALSE)
sigma_d = extract(sa,"sigma_d",permuted=FALSE)
sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)

par(mfrow=c(2,2))
hist(sigma2_a,main="Var(a)",xlim=c(0,.12),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_a,.5),3)))
hist(sigma2_b,main="Var(b)",xlim=c(.35,.65),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma2_b,.5),3)))
hist(sigma_c^2,main="Var(c)",xlim=c(0,.012),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_c^2,.5),3)))
hist(sigma_d^2,main="Var(d)",xlim=c(.06,.13),freq=FALSE,ylab = "Proportion",xlab=paste("median =",round(quantile(sigma_d^2,.5),3)))

vc_tab = array(0, dim=c(4,3))
vc_tab[1,] = quantile(sigma2_a,c(.025,.5,.975))
vc_tab[2,] = quantile(sigma_c^2,c(.025,.5,.975))
vc_tab[3,] = quantile(sigma2_b,c(.025,.5,.975))
vc_tab[4,] = quantile(sigma_d^2,c(.025,.5,.975))

## Reshape parameters and reorient loadings with PCA rotation
Q_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
R_new = array(0, dim=c(P,P,Nchains*Nsamples/2))
S_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
W_new = array(0, dim=c(d,d,Nchains*Nsamples/2)) 
Theta_old = array(0, dim=c(d,N,Nchains*Nsamples/2))
Theta_new = array(0, dim=c(d,N,Nchains*Nsamples/2))
Lambda_old = array(0, dim=c(P,d,Nchains*Nsamples/2))
Lambda_new = array(0, dim=c(P,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind+ 1
		Theta_old[,,ind] = t(array(Theta[i,j,],dim=c(N,d)))
		Lambda_old[,,ind] = array(Lambda[i,j,],dim=c(P,d))
		Q_new[,,ind] = array(Q[i,j,],dim=c(P,P))
		W_new[,,ind] = diag(eigen(Q_new[,,ind] - diag(rep(sigma_eps[i,j,]^2,P)))$values[1:d]^.5)
		R_new[,,ind] = diag(diag(Q_new[,,ind])^-.5)%*%Q_new[,,ind]%*%diag(diag(Q_new[,,ind])^-.5)
		S_new[,,ind] = diag(eigen(Q_new[,,ind])$values[1:d]^.5)
		Lambda_new[,,ind] = eigen(Q_new[,,ind])$vectors[,1:d]%*%S_new[,,ind]
	}
}

## Fix signs of factors and rotate factor scores
ll = rep(0,d)
ll_mx = rep(0,d)
for(k in 1:d){
	for(j in 1:(Nsamples*Nchains/2)){
		if(max(abs(Lambda_new[,k,j]))>ll_mx[k])
			ll[k] = which(abs(Lambda_new[,k,j])==max(abs(Lambda_new[,k,j])))
			ll_mx[k] = max(abs(Lambda_new[,k,j]))
	}
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		for(k in 1:d){
			Lambda_new[,k,ind] = sign(Lambda_new[ll[k],k,ind]) * Lambda_new[,k,ind]
		}
		Theta_new[,,ind] = solve(S_new[,,ind])%*%W_new[,,ind]%*%t(Lambda_new[,,ind]%*%solve(S_new[,,ind])%*%W_new[,,ind])%*%
			Lambda_old[,,ind]%*%solve(t(Lambda_old[,,ind])%*%Lambda_old[,,ind])%*%Theta_old[,,ind]
	}
}

## Varimax & Promax Rotations
Lambda_vmx = Lambda_new
Theta_vmx = Theta_new
Rot_vmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
Lambda_pmx = Lambda_new
Theta_pmx = Theta_new
Rot_pmx = array(0, dim = c(d,d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		tmp = varimax(Lambda_new[,,ind])
		Rot_vmx[,,ind] = tmp$rot
		Lambda_vmx[,,ind] = Lambda_new[,,ind]%*%Rot_vmx[,,ind]
		Theta_vmx[,,ind] = t(Rot_vmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
		tmp = psych::Promax(Lambda_new[,,ind])
		Rot_pmx[,,ind] = tmp$rot
		Lambda_pmx[,,ind] = Lambda_new[,,ind]%*%Rot_pmx[,,ind]
		Theta_pmx[,,ind] = t(Rot_pmx[,,ind])%*%t(scale(t(Theta_new[,,ind])))
	}
}
## Permute Varimax factors
tr=principal(Y,nfactors=d,rotate="varimax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_vmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_vmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_vmx[,,ind] = Lambda_vmx[,ll,ind]
		Theta_vmx[,,ind] = Theta_vmx[ll,,ind]
		Rot_vmx[,,ind] = Rot_vmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_vmx[d1,,ind] = Theta_vmx[d1,,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Rot_vmx[,d1,ind] = Rot_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
			Lambda_vmx[,d1,ind] = Lambda_vmx[,d1,ind]*sign(Lambda_vmx[abs(Lambda_vmx[,d1,ind])==max(abs(Lambda_vmx[,d1,ind])),d1,ind])
		}
	}
}

## Permute Promax factors
tr=principal(Y,nfactors=d,rotate="promax")$loadings[,1:d]
for(d1 in 1:d){
	tr[,d1] = sign(tr[abs(tr[,d1])==max(abs(tr[,d1])),d1])*tr[,d1]
}
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		ll=rep(0,d)
		for(d1 in 1:d){
			tmp1=apply((Lambda_pmx[,d1,ind] - tr[,])^2,2,sum)
			tmp2=apply((Lambda_pmx[,d1,ind] + tr[,])^2,2,sum)
			if(min(tmp1) < min(tmp2)) ll[d1]=which(tmp1==min(tmp1))
			if(min(tmp1) > min(tmp2)) ll[d1]=which(tmp2==min(tmp2))
		}
		Lambda_pmx[,,ind] = Lambda_pmx[,ll,ind]
		Theta_pmx[,,ind] = Theta_pmx[ll,,ind]
		Rot_pmx[,,ind] = Rot_pmx[ll,ll,ind]
		for(d1 in 1:d){
			Theta_pmx[d1,,ind] = Theta_pmx[d1,,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
			Rot_pmx[,d1,ind] = Rot_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
			Lambda_pmx[,d1,ind] = Lambda_pmx[,d1,ind]*sign(Lambda_pmx[abs(Lambda_pmx[,d1,ind])==max(abs(Lambda_pmx[,d1,ind])),d1,ind])
		}
	}
}
## Variance explained by retained factors
var_new = array(0, dim=c(d,Nchains*Nsamples/2))
var_vmx = array(0, dim=c(d,Nchains*Nsamples/2))
var_pmx = array(0, dim=c(d,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		var_new[,ind] = diag(t(Lambda_new[,,ind])%*%(Lambda_new[,,ind]))
		var_vmx[,ind] = diag(t(Lambda_vmx[,,ind])%*%(Lambda_vmx[,,ind]))
		var_pmx[,ind] = diag(t(Lambda_pmx[,,ind])%*%(Lambda_pmx[,,ind]))
	}
}
## Communalities & Uniquenesses
comm_uniq_new = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_vmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
comm_uniq_pmx = array(0, dim=c(P,2,Nchains*Nsamples/2))
ind = 0
for(i in 1:dim(Q)[1]){
	for(j in 1:dim(Q)[2]){
		ind = ind + 1
		comm_uniq_new[,,ind] = cbind(diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])),1-diag((Lambda_new[,,ind])%*%t(Lambda_new[,,ind])))
		comm_uniq_vmx[,,ind] = cbind(diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])),1-diag((Lambda_vmx[,,ind])%*%t(Lambda_vmx[,,ind])))
		comm_uniq_pmx[,,ind] = cbind(diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])),1-diag((Lambda_pmx[,,ind])%*%t(Lambda_pmx[,,ind])))
	}
}
## Compute posterior median estimates and 95% posterior credible intervals
Q_m = array(0,dim=c(P,P,3))
for(p1 in 1:P){
	for(p2 in 1:P){
		Q_m[p1,p2,] = quantile(Q_new[p1,p2,],probs=c(.025,.5,.975))
	}
}
comm_uniq_new_m = array(0,dim=c(P,2,3))
comm_uniq_vmx_m = array(0,dim=c(P,2,3))
comm_uniq_pmx_m = array(0,dim=c(P,2,3))
Lambda_new_m = array(0,dim=c(P,d,3))
Lambda_vmx_m = array(0,dim=c(P,d,3))
Lambda_pmx_m = array(0,dim=c(P,d,3))
for(p in 1:P){
	comm_uniq_new_m[p,1,] = quantile(comm_uniq_new[p,1,],probs=c(.025,.5,.975))
	comm_uniq_vmx_m[p,1,] = quantile(comm_uniq_vmx[p,1,],probs=c(.025,.5,.975))
	comm_uniq_pmx_m[p,1,] = quantile(comm_uniq_pmx[p,1,],probs=c(.025,.5,.975))
	comm_uniq_new_m[p,2,] = quantile(comm_uniq_new[p,2,],probs=c(.025,.5,.975))
	comm_uniq_vmx_m[p,2,] = quantile(comm_uniq_vmx[p,2,],probs=c(.025,.5,.975))
	comm_uniq_pmx_m[p,2,] = quantile(comm_uniq_pmx[p,2,],probs=c(.025,.5,.975))
	for(k in 1:d){
		Lambda_new_m[p,k,] = quantile(Lambda_new[p,k,],probs=c(.025,.5,.975))
		Lambda_vmx_m[p,k,] = quantile(Lambda_vmx[p,k,],probs=c(.025,.5,.975))
		Lambda_pmx_m[p,k,] = quantile(Lambda_pmx[p,k,],probs=c(.025,.5,.975))
	}
}
Theta_new_m = array(0,dim=c(d,N,3))
Theta_vmx_m = array(0,dim=c(d,N,3))
Theta_pmx_m = array(0,dim=c(d,N,3))
for(k in 1:d){
	for(i in 1:N){
		Theta_new_m[k,i,] = quantile(Theta_new[k,i,],probs=c(.025,.5,.975))
		Theta_vmx_m[k,i,] = quantile(Theta_vmx[k,i,],probs=c(.025,.5,.975))
		Theta_pmx_m[k,i,] = quantile(Theta_pmx[k,i,],probs=c(.025,.5,.975))
	}
}		
var_new_m = array(0,dim=c(d,3))
var_vmx_m = array(0,dim=c(d,3))
var_pmx_m = array(0,dim=c(d,3))
for(k in 1:d){
	var_new_m[k,] = quantile(var_new[k,],c(.025,.5,.975))
	var_vmx_m[k,] = quantile(var_vmx[k,],c(.025,.5,.975))
	var_pmx_m[k,] = quantile(var_pmx[k,],c(.025,.5,.975))
}
var_expl_new = apply(var_new_m,2,sum)/P
var_expl_vmx = apply(var_vmx_m,2,sum)/P
var_expl_pmx = apply(var_pmx_m,2,sum)/P

## Create SM Table 1A
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_new_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_new_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_new_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm1a = knitr::kable(tabAsStringMatrix)

## Create SM Table 1B
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_pmx_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_pmx_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_pmx_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm1b = knitr::kable(tabAsStringMatrix)

## Create SM Table 3
comm_uniq.tab = as.data.frame(array(0,dim=c(P+1,6)))
names(comm_uniq.tab) = c("","Communality","","","Uniqueness","")
rownames(comm_uniq.tab) = c("Quantiles:",names(data)[ind_Y])
comm_uniq.tab[1,] = rep(c(".025","0.50",".975"),2)
comm_uniq.tab[2:(P+1),c(1,4)]=sprintf("%.3f", round(comm_uniq_vmx_m[,,1],3))
comm_uniq.tab[2:(P+1),c(2,5)]=sprintf("%.3f", round(comm_uniq_vmx_m[,,2],3))
comm_uniq.tab[2:(P+1),c(3,6)]=sprintf("%.3f", round(comm_uniq_vmx_m[,,3],3))
tabAsStringMatrix = print(comm_uniq.tab, printToggle = FALSE, noSpaces = TRUE)
tab_sm3 = knitr::kable(tabAsStringMatrix)

## Create Table 3
lambda.tab = as.data.frame(array(0,dim=c(P+1,d*3)))
names(lambda.tab) = c("","PC1","","","PC2","","","PC3","")
rownames(lambda.tab) = c("",names(data)[ind_Y])
lambda.tab[1,] = rep(c(".025","0.50",".975"),d)
for(d1 in 1:d){
	lambda.tab[2:(P+1),3*(d1-1)+1]=round(Lambda_vmx_m[,d1,1],3)
	lambda.tab[2:(P+1),3*(d1-1)+2]=round(Lambda_vmx_m[,d1,2],3)
	lambda.tab[2:(P+1),3*(d1-1)+3]=round(Lambda_vmx_m[,d1,3],3)
}
tabAsStringMatrix = print(lambda.tab, printToggle = FALSE, noSpaces = TRUE)
tab3 = knitr::kable(tabAsStringMatrix)

## Add PC scores into dataset
data1$pc1 = Theta_vmx_m[1,,2]
data1$pc2 = Theta_vmx_m[2,,2]
data1$pc3 = Theta_vmx_m[3,,2]