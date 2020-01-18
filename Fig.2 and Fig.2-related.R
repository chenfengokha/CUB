library(stats4)
library(emg)
library(magrittr)
library(plyr)
library(dplyr)
library(parallel)
library(readr)
library(bbmle)
library(stats4)
library(emg)
library(reshape2)
library(tidyr)
##load the human-flu ribo density data
a02 <- read.table("/mnt/data/home/phil/acting/virusTranslation/03.flu/05.res/SRR3623932.trim.rRNA_unmap.tophat.cov")
a22 <- read.table("/mnt/data/home/phil/acting/virusTranslation/03.flu/05.res/SRR3623937.trim.rRNA_unmap.tophat.cov")
rawdata1 <- rbind(a02 %>% cbind(time=0,type="36"),a22 %>% cbind(time=2,type="37"))
names(rawdata1)[1:7] <- c("cds","chr","sitechr","codon","pos","strand","ribo")
rawdata1$chr <- as.vector(rawdata1$chr)
rawdata1$cds <- as.vector(rawdata1$cds)
rawdata1$type <- as.vector(rawdata1$type)
##virus ribo
virusribo <- rawdata1 %>% dplyr::filter(chr %in% c("EF467817.1", "EF467818.1", "EF467819.1", "EF467820.1", "EF467821.1", "EF467822.1", "EF467823.1", "EF467824.1")) %>% separate(cds,c("cds","2"))
##human ribo
human <- paste("chr",c(1:22,"X","Y","M"),sep = "")
mydata1_2 <- rawdata1 %>% dplyr::filter(chr %in% human) %>% separate(cds,c("cds","2"))
#select the longest transcript of each gene
load("~/codonpaper/review/allhumangenome.Rdata")
allhumangenome %>% as.data.frame() %>% group_by(spe,name) %>% dplyr::filter(spe=="Homo_sapiens" & leng==max(leng)) -> ttt
idtrans <- ttt[,3:4] %>% group_by(name) %>% dplyr::summarize(namet=name1[1])
aa <- mydata1_2 %>% dplyr::filter(cds %in% idtrans$namet)
aa -> mydata1_2
##TA ratio
mydatah <- mydata1_2 %>%
  dplyr::filter(time==2) %>%
  group_by(time,cds) %>%
  dplyr::filter((length(codon) > 70) & (pos > 20) & (pos < (length(codon)-20))) %>%
  dplyr::summarize(txActi=mean(ribo,trim=0.25,na.rm=T))
mydatav <- virusribo %>%
  dplyr::filter(time==2) %>%
  group_by(time,cds) %>%
  dplyr::filter((length(codon) > 70) & (pos > 20) & (pos < (length(codon)-20))) %>%
  dplyr::summarize(txActi=mean(ribo,trim=0.25,na.rm=T))
tt <- (mydatav$txActi %>% mean())/(mydatah$txActi %>% mean())
##function to calculated typical decoding time
emgnb.mle <- function(x, expr, lower = NA, upper = NA) {
  if(length(x) != length(expr)) {
    stop("'x' and 'expr' must have the same length");
  }
  normX <- x / expr;
  rmIdx <- which(is.na(x));
  if(length(rmIdx) > 0) {
    x <- x[-rmIdx];
    expr <- expr[-rmIdx];
    normX <- normX[-rmIdx];
  }
  ## Limiting the space to search is very important here.
  ## Since the waiting time cannot be negative, the real sigma is expected to be < 1/4 of mu, which approximated by median(normX)
  if(is.na(lower)) {
    lower <- list(mu = min(normX[normX != 0])/2, sigma = sd(normX)/1000, lambda = 0.1/mean(normX),size=1e-3);
  }
  if(is.na(upper)) {
    upper <- list(mu = quantile(normX[normX != 0],p=.5), sigma = sd(normX),
                  lambda = 10/mean(normX),size=10);
  }
  startL <- list(mu = quantile(normX,p=.75)/2,
                 sigma = median(normX)/10, lambda = 1/mean(normX),size=1)
  mle2(function(size,mu,sigma,lambda) {
    emgnb.nllik(x,expr,size, mu, sigma, lambda)
  },
  method = "L-BFGS-B",  lower = lower,  upper = upper,  start = startL,default.start=F  );
}
emgnb.nllik <- function(x,expr,size,mu,sigma,lambda) {
  sigma <- max(c(sigma,1e-3));
  lambda <- max(c(lambda,1e-3));
  decodeTime <- seq(0,5,by=.01);
  vecRiboCnt <- rep(x,each=length(decodeTime));
  vecTxActi <- rep(expr,each=length(decodeTime));
  vecDecTime <- rep(decodeTime,length(x));
  vecIdx <- rep(seq(x),each=length(decodeTime));
  vecLogDemg <- demg(decodeTime,mu,sigma,lambda,log=T);
  sumExpVecLogDemg <- sum(exp(vecLogDemg));
  allP <- exp(dnbinom(vecRiboCnt,mu=vecDecTime * vecTxActi,size=size,log=T) + rep(vecLogDemg,length(x)) );
  if(length(x) == 1) {
    return(unname(-log(
      sum(allP) / sumExpVecLogDemg
    )));
  } else {
    tapply(allP,vecIdx,sum) -> tmpP;
    tmpP <- ifelse(tmpP < 1e-300,1e-300,tmpP);
    return(sum(unname(-log(
      tmpP))) - sumExpVecLogDemg
    );
  }
}
##calculate TDT and virus consumption
dfHuman.all <- mydata1_2 %>%
  group_by(time,cds) %>%
  dplyr::mutate(txActi=mean(ribo,trim=0.25,na.rm=T));  ## calculate translational activity
load("/mnt/data/home/chenfeng/codonpaper/tAIcodhum.Rdata")
tAIofcodonhuman$tai <- tAIofcodonhuman$Wi / max(tAIofcodonhuman$Wi);
dfVirusConsump.h <- virusribo %>%
  group_by(time,cds) %>%
  dplyr::filter(pos > 20 & pos < (length(codon)-20)) %>%
  dplyr::mutate(txActi=mean(ribo,trim=0.25,na.rm=T))%>%
  filter(substr(chr,1,2) == "EF") %>%
  group_by(time,codon) %>%
  dplyr::summarize(consump=sum(txActi)) %>%
  mutate(tai = tAIofcodonhuman$tai[match(codon,as.character(tAIofcodonhuman$codon))]) %>%
  mutate(RTSv = consump/tai);
dfAllRes <- dfHuman.all %>%
  filter(txActi > 1) %>%
  filter(substr(chr,1,2) != "EF") %>%
  group_by(cds) %>%
  dplyr::filter((length(codon) > 70) & (pos > 20) & (pos < (length(codon)-20))) %>% ## remove sites near ramp and termination
  split(paste0(.$time,.$codon)) %>%
  mclapply(mc.cores=60,function(x){
    if(dim(x)[1] < 10) {
      data.frame(NULL)
    } else {
      data.frame(t(attr(emgnb.mle(x$ribo,x$txActi),"coef"))) %>% cbind(codon=x$codon[1],time=x$time[1])
    }
  }) %>%
  rbind.fill();
dfVirusConsump.h$codon <- as.vector(dfVirusConsump.h$codon)
dfHuman.all$codon <- as.vector(dfHuman.all$codon)
dfAllRes$codon <- as.vector(dfAllRes$codon)
save(mydata1_2,virusribo,mydatah,mydatav,tt,dfHuman.all,dfAllRes,dfVirusConsump.h,file="~/codonpaper/code and data/fig2.finaldata.Rdata")
load("~/codonpaper/code and data/fig2.finaldata.Rdata")
mydf <- dfHuman.all %>% dplyr::filter(time==2)
mydf$txActi <- mydatah$txActi[match(mydf$cds,mydatah$cds)]
mydf <- mydf %>% dplyr::filter(!is.na(txActi)) %>% group_by(codon) %>% dplyr::summarize(conofhuman=sum(txActi)) %>% as.data.frame()
dfVirusConsump.h <- dfVirusConsump.h %>% dplyr::filter(!(codon %in% c("TGA","TAG","TAA")) & time==2) %>% as.data.frame()
dfVirusConsump.h$conofhuman <- mydf$conofhuman[match(dfVirusConsump.h$codon,mydf$codon)]
VtRNA0 <- unique(dfVirusConsump.h[,c(2,4)])
VtRNA0$time <- 0
VtRNA0$consump <- 0
VtRNA0$RTSv <- 0
VtRNA0$conofhuman <- 1
virustRNAconsumption <- rbind(VtRNA0,dfVirusConsump.h)
virustRNAconsumption$TDT <- dfAllRes$mu[match(paste(virustRNAconsumption$codon,virustRNAconsumption$time),paste(dfAllRes$codon,dfAllRes$time))]
virustRNAconsumption$rtsvofhuman <- virustRNAconsumption$consump/virustRNAconsumption$conofhuman
##calculated the sensitivity of codons after the virus infection, here slope is sensitivity
tmp2 <- mclapply(mc.cores=20,1:length(unique(virustRNAconsumption$codon)),function(x){
  mdf <- virustRNAconsumption %>% dplyr::filter(codon==unique(virustRNAconsumption$codon)[x])
  slop1 <- (mdf$TDT[2]/mdf$TDT[1])/mdf$consump[2]
  data.frame(codon=unique(mdf$codon),tAI=unique(mdf$tai),t2t0=(mdf$TDT[2]/mdf$TDT[1]),c2=mdf$consump[2],conofhuman=mdf$conofhuman[2],Sensitivity=slop1)
}) %>% rbind.fill()
##save the final data of fig2 and fig2-related
save(tmp2,virustRNAconsumption,file ="~/codonpaper/code and data/finaldata_forfig2.correlations.Rdata")
p <- format(cor.test(tmp2$Sensitivity,tmp2$tAI,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(tmp2$Sensitivity,tmp2$tAI,method = "spearman")$estimate[[1]], digits = 3)
source("~/Rfunction/style.print.R")
##draw Fig.2b
ggplot(data=tmp2, aes(x=Sensitivity,y=tAI))+
  geom_point() + 
  labs(title=paste("Rho=",rho,", P=",p),x="Sensitivity",y="tRNA supply of codons") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)
##draw Fig.S5b
p <- format(cor.test(tmp2$Sensitivity,tmp2$conofhuman,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(tmp2$Sensitivity,tmp2$conofhuman,method = "spearman")$estimate[[1]], digits = 3)
ggplot(data=tmp2, aes(x=Sensitivity,y=conofhuman))+
  geom_point() +
  labs(title=paste("Rho=",rho,", P=",p),x="Sensitivity",y="tRNA supply of codons") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)
new2a <- virustRNAconsumption %>% dplyr::filter(time==2)
p <- format(cor.test(new2a$RTSv,new2a$TDT,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(new2a$RTSv,new2a$TDT,method = "spearman")$estimate[[1]], digits = 3)
##Fig.2a
ggplot(new2a,aes(RTSv,TDT)) + geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  labs(title=paste("Rho=",rho,", P=",p),x="Relative tRNA shortages \n caused by viral translation",y="Typical decoding time") +
  style.print()+
  scale_x_continuous(trans = "log10")
##Fig.S5a
newS2a <- virustRNAconsumption %>% dplyr::filter(time==2)
p <- format(cor.test(newS2a$rtsvofhuman,newS2a$TDT,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(newS2a$rtsvofhuman,newS2a$TDT,method = "spearman")$estimate[[1]], digits = 3)
ggplot(newS2a,aes(rtsvofhuman,TDT)) + geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  labs(title=paste("Rho=",rho,", P=",p),x="Relative tRNA shortages \n caused by viral translation",y="Typical decoding time") +
  style.print()

##eight human gene,Fig.S5c
load("~/codonpaper/code and data/fig2.finaldata.Rdata")
human <- paste("chr",c(1:22,"X","Y","M"),sep = "")
mydatah <- mydata1_2 %>% 
  group_by(cds,time) %>%
  dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20) & chr %in% human) %>%
  dplyr::summarize(mfj=mean(ribo,trim=0.25,na.rm=T)) %>%
  as.data.frame() %>% arrange(desc(mfj))
highTAgenehuman <- (mydatah %>% dplyr::filter(time==0))[1:500,]$cds
##con of human gene
mydf <- mydata1_2 %>% dplyr::filter(chr %in% human)
mydf$TA <- mydatah$mfj[match(paste(mydf$cds,mydf$time),paste(mydatah$cds,mydatah$time))]
load("/mnt/data/home/chenfeng/codonpaper/tAIcodhum.Rdata")
tAIofcodonhuman$codon <- as.vector(tAIofcodonhuman$codon)
tAIofcodonhuman$wi <- tAIofcodonhuman$Wi/max(tAIofcodonhuman$Wi)
mydf$tai <- tAIofcodonhuman$wi[match(mydf$codon,tAIofcodonhuman$codon)]
##randomly select eight gene for 1000 times
aa <- mclapply(mc.cores = 60,1:1000,function(x){
  a <- mydf %>% dplyr::filter(!(is.na(TA)) & cds %in% sample(highTAgenehuman,8,replace = F)) %>% 
    group_by(codon,tai,time) %>% 
    #dplyr::filter(!(codon %in% c("TGA","TAG","TAA","ATG","TGG")))%>%
    dplyr::summarize(consump=sum(TA)) %>% as.data.frame() %>%
    dplyr::mutate(RTSv = consump/tai) 
  a$TDT <- dfAllRes$mu[match(paste(a$codon %>% as.vector(),a$time),paste(dfAllRes$codon %>% as.vector(),dfAllRes$time))]
  rho2a <- a %>% group_by(time) %>% dplyr::summarize(rho=format(cor.test(TDT,RTSv,method = "s")$estimate[[1]], digits = 3) %>% as.numeric())
  b0 <- a %>% dplyr::filter(time==0)  
  b2 <- a %>% dplyr::filter(time==2)
  b0$tdt2 <- b2$TDT[match(b0$codon,b2$codon)]
  b0$con2 <- b2$consump[match(b0$codon,b2$codon)]
  b0$sen <- b0$tdt2/b0$TDT/(b0$con2/b0$consump)
  rho2b <- format(cor.test(b0$sen,b0$tai,method = "s")$estimate[[1]], digits = 3) %>% as.numeric()
  p2b <- format(cor.test(b0$sen,b0$tai,method = "s")$p.value, digits = 3) %>% as.numeric()
  data.frame(r2atime0=rho2a$rho[which(rho2a$time==0)],r2atime2=rho2a$rho[which(rho2a$time==2)],rho2b,p2b)
}) %>% rbind.fill()
source("~/Rfunction/style.print.R")
data.frame(mean=c(mean(aa$r2atime2),mean(aa$rho2b)),sd=c(sd(aa$r2atime2),sd(aa$rho2b)),type=c("2","1")) %>%
  ggplot(aes(type,mean)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),width=0.3)+
  labs(y=expression(paste("ρ (",RTS[v], " ~ decoding time)")))+
  style.print() + 
  #theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5,angle = 25))+ 
  theme(legend.position = "none",axis.title.y=element_blank()) + 
  #scale_y_continuous(limits = c(0,0.4),breaks=c(0,0.1,0.2,0.3,0.4))+
  coord_flip()
