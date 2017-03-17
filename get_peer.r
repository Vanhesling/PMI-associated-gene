#2016-8-25
#Authered by ZYZ

library(peer)
library("R.matlab")
setwd("/media/data/ZYZ/GTEx2.0/2.results/peer_exclude_noneuropean")
rm(list=ls());
expr=dir(pattern = "expr_*");
smpl=dir(pattern = "smpl_*");
cov=dir(pattern = "cov_*");
for(i in 1:length(expr)){
    print(paste("now comes to",i));
    a=expr[i];
    b=unlist(strsplit(a,"_"));b=b[length(b)];
    n=as.numeric(unlist(strsplit(b,".m"))[1]);
    tissuename=sub('^.{5}',"",a);
    data <- readMat(expr[i]);
    expr_martrix=matrix(unlist(data), ncol = n, byrow = FALSE)
    covs = read.table(cov[i],header=FALSE)

    model = PEER()
    PEER_setPhenoMean(model,as.matrix(t(expr_martrix)))
    PEER_setAdd_mean(model, TRUE)
    PEER_setNk(model,15)
    PEER_setCovariates(model, as.matrix(covs))

    PEER_update(model)
    residuals = PEER_getResiduals(model)
    factors = PEER_getX(model)
    covsx = PEER_getCovariates(model)
    weigthx = PEER_getW(model)
    writeMat(paste("peer_",tissuename),residuals=residuals,factors=factors,covsx=covsx,weigthx=weigthx)
    print(paste("Finished",i));
}
