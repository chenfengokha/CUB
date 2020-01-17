##virus TA in yeast fig1a
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/datJustin.YPD1.RData");
datJustin.YPD1 <- datJustin; ## rename to avoid overwriting

load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/datJustin.YPD2.RData");
datJustin.YPD2 <- datJustin;

load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/datJustin.Sc-His.RData");
datJustin.ScHis <- datJustin;

load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/datJustin.Sc-Lys.RData");
datJustin.ScLys <- datJustin;

load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Ingolia2009_remap/12.2.datIngolia.RData");

load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Hunter2014/datHunter.RData");
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_McManus2014/datMcManus.RData");

datPop <- read.delim("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Pop2014/05.allData.tsv",header=F,as.is=T,sep="\t");
colnames(datPop) <- c("orf","chrom","chromPos","codon","pos","strand","ribo");

datZyn <- read.delim("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Zinshteyn2013/05.prepareDat.includeOffFrame.jitter2frame.tsv",header=F,as.is=T,sep="\t");
colnames(datZyn) <- c("orf","chrom","chromPos","codon","pos","strand","ribo")
mylist <- list(datJustin.ScHis,datJustin.ScLys,datJustin.YPD1,datJustin.YPD2,datIngolia,datHunter,datMcManus,datPop,datZyn)
names(mylist) <- c("ScHis","ScLys","YPD1","YPD2","Ingolia2009","Artieri2014","McManus2014","Pop2014","Zinshteyn2013")

yeastvirussite <- mclapply(mc.cores = 4,1:9,function(x){
  #fold change of TA in virus and yeast
  mydata <- mylist[[x]] %>%
      group_by(orf,chrom) %>%
      dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20)) %>%
      dplyr::summarize(TA=mean(ribo,trim=0.25,na.rm=T)) %>%
      as.data.frame()
  fc <- ((mydata %>% dplyr::filter(substr(chrom,1,2)!="ch"))$TA %>% mean()) / ((mydata %>% dplyr::filter(substr(chrom,1,2)=="ch"))$TA %>% mean())
  data.frame(fc,name=names(mylist[x]))
  # mydata <- mylist[[x]] %>%
  #   group_by(orf,chrom) %>%
  #   dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20)) %>%
  #   dplyr::summarize(TA=mean(ribo,trim=0.25,na.rm=T)) %>%
  #   as.data.frame()
  # virus <- mydata %>% dplyr::filter(chrom %in% c("SCU01060","NC_003745")) %>% dplyr::summarize(sumTA=sum(TA)) %>% as.vector()
  # yeastsite <- mydata %>% dplyr::filter(!(chrom %in% c("SCU01060","NC_003745"))) %>% arrange(desc(TA)) %>%
  #   dplyr::filter(TA>as.numeric(virus)) %>% nrow()
  # yeast <- mydata %>% dplyr::filter(!(chrom %in% c("SCU01060","NC_003745"))) %>% dplyr::summarize(sumTA=sum(TA)) %>% as.vector()
  # data.frame(name=names(mylist[x]),vTA=as.numeric(virus),yTA=as.numeric(yeast),yeastsite=yeastsite)
}) %>% rbind.fill()
yeastvirussite$per <- yeastvirussite$vTA/yeastvirussite$yTA
###########################################
##coverage. figs1
yeastvirus <- datJustin.ScHis %>% cbind(type="Sc-His") %>%
  rbind(datJustin.ScLys %>% cbind(type="Sc-Lys")) %>%
  rbind(datJustin.YPD1 %>% cbind(type="YPD1")) %>%
  rbind(datJustin.YPD2 %>% cbind(type="YPD2")) %>%
  rbind(datIngolia[,c("orf","chrom","chromPos","codon","pos","strand","ribo")] %>% cbind(type="Ingolia2009")) %>%
  rbind(datHunter %>% cbind(type="Artieri2014")) %>%
  rbind(datMcManus %>% cbind(type="McManus2014")) %>%
  rbind(datPop %>% cbind(type="Pop2014")) %>%
  rbind(datZyn %>% cbind(type="Zinshteyn2013")) 
seqdepth <- yeastvirus %>% 
  group_by(orf,type) %>%
  dplyr::filter(substr(chrom,1,2)=="ch") %>%
  dplyr::filter((length(codon) > 70) & (pos > 20) & (pos < (length(codon)-20))) %>%
  dplyr::mutate(mfj=mean(ribo,trim=0.25,na.rm=T)) %>%
  as.data.frame() %>%
  group_by(type) %>%
  dplyr::summarize(all=sum(mfj)) %>%
  as.data.frame()
virus <- yeastvirus %>% dplyr::filter(substr(chrom,1,2)!="ch") %>%
  dplyr::mutate(depth=seqdepth$all[match(as.vector(type),as.vector(seqdepth$type))]) %>%
  dplyr::mutate(ribonor=ribo*1000000/depth) %>% as.data.frame()

virus$type <- factor(virus$type,levels=c("Sc-Lys","Sc-His","YPD1","YPD2","Ingolia2009","Zinshteyn2013","Pop2014","Artieri2014","McManus2014"))

source("~/Rfunction/style.print.R")
##picture
a <- virus %>% group_by(type,chrom) %>% 
  dplyr::mutate(low=sort(ribonor)[ceiling(length(ribonor)*0.25)],high=sort(ribonor)[ceiling(length(ribonor)*0.75)]) %>%
  dplyr::filter(ribonor<=high & ribonor>=low) %>%
  dplyr::mutate(typesite="in")
a %>% group_by(type,typesite) %>% dplyr::summarize(total=sum(ribonor)) %>% as.data.frame() %>% arrange(typesite)
tmp <- a %>% dplyr::filter(typesite=="in")
win <- data.frame(stat=seq(from=1,to=4500-60+1,by=15),
                  end=seq(from=60,to=4500,by=15))
windowa <- mclapply(mc.cores = 20,1:nrow(win),function(x){
  tmp %>% dplyr::filter(chromPos>win[x,1] & chromPos<win[x,2]) %>%
    group_by(chrom,type) %>%
    dplyr::summarize(sum=sum(ribonor))  %>% 
    dplyr::mutate(site=win[x,1])
}) %>% rbind.fill()
save(windowa,file = "~/codonpaper/review/2a.coverage.windowa.Rdata")
load("~/codonpaper/review/2a.coverage.windowa.Rdata")
windowa %>%
  ggplot(aes(x=site,y=sum+1)) +
  geom_area() +
  labs(x="Position in ORF (nt)",y="Ribo-Seq coverage")+
  facet_grid(type~chrom)+
  style.print()+
  scale_x_continuous(breaks = c(0,2000,4000))+
  scale_y_log10()

###########################################
#figs2a and s4a 
load("~/codonpaper/review/1new.datjustin.sclys.Rdata")
load("/mnt/data/home/chenfeng/codonpaper/RTSvyeast.Rdata")
load("/mnt/data/home/phil/acting/virusTranslation/01.scratch/04.result.yeast.RData")
load("~/codonpaper/review/1BCDE.highTAgene.Rdata")

rhoce <- mclapply(mc.cores = 20,1:1000,function(x){
  consump <- datJustin.ScLys %>%
    dplyr::filter(orf %in% sample(highTAgene$orf,3,replace = F)) %>%
    group_by(codon,type) %>% 
    dplyr::filter(!is.na(TA) & !(codon %in% c("TAG","TGA","TAA"))) %>%
    dplyr::summarize(consump=sum(TA),tai=unique(wi)) %>%
    as.data.frame()
  consump$wi <- datJustin.ScLys$wi[match(paste(consump$codon,consump$type),paste(datJustin.ScLys$codon,datJustin.ScLys$type))]
  consump$rtsv <- consump$consump/consump$tai
  consump$rrt <- RTSv$rrt[match(consump$codon,RTSv$codon)]
  consump$tdt <- dfCodon.yeast$tdt[match(consump$codon,dfCodon.yeast$codon)]
  rhoC <- format(cor.test(consump$rrt,consump$rtsv,method = "s")$estimate[[1]], digits = 3)
  rhoE <- format(cor.test(consump$rtsv,consump$tdt,method = "spearman")$estimate[[1]], digits = 3)
  data.frame(rhoc=rhoC %>% as.numeric(),rhoe=rhoE %>% as.numeric(),x)
  
}) %>% rbind.fill() %>%
  dplyr::summarize(c=mean(rhoc),sdc=sd(rhoc),
                   e=mean(rhoe),sde=sd(rhoe)) %>% as.data.frame()

data.frame(mean=c(rhoce$e,rhoce$c),sd=c(rhoce$sde,rhoce$sdc),type=c("a","c")) %>%
  ggplot(aes(type,mean)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),width=0.3)+
  labs(y=expression(paste("ρ (",RTS[v], " ~ decoding time)")))+
  style.print() + 
  #theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5,angle = 25))+ 
  theme(legend.position = "none",axis.title.y=element_blank()) + 
  #scale_y_continuous(limits = c(0,0.4),breaks=c(0,0.1,0.2,0.3,0.4))+
  coord_flip()


#figs2c
load("/mnt/data/home/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/18.2.RData")
my1D <-  (dat.allRRT %>% dplyr::filter(names %in% c("McManus2014", "Hunter2014", "Pop2014", "Zinshteyn2013", "Sc-Lys", "Ingolia2009")) %>% as.data.frame())[,c(2,4,10)]
my1D$names[which(my1D$names=="Hunter2014")] <- "Artieri2014"

load("~/codonpaper/review/1Drawdata.datjustin.Rdata")
load("/mnt/data/home/chenfeng/codonpaper/tAIcodyeast.Rdata")
tAIofcodonyeast$wi <- tAIofcodonyeast$Wi/max(tAIofcodonyeast$Wi)
tAIofcodonyeast$codon <- as.vector(tAIofcodonyeast$codon)
datJustin$wi <- tAIofcodonyeast$wi[match(datJustin$codon,tAIofcodonyeast$codon)]
datJustin$rrt <- my1D$Pos6[match(paste(datJustin$exp,datJustin$codon),paste(my1D$names,my1D$codon))]
datJustin$type <- "1"
datJustin$type[which(datJustin$exp!="Sc-Lys")] <- "2"

data1D <- mclapply(mc.cores =60,1:1000,function(x){
  consump <- datJustin %>%
    dplyr::filter(orf %in% sample(highTAgene$orf,3,replace = F)) %>%
    group_by(orf,exp) %>% 
    dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20)) %>%
    dplyr::mutate(TA=mean(ribo,trim=0.25,na.rm=T)) %>%
    as.data.frame() %>%
    group_by(codon,exp,wi,rrt) %>%
    dplyr::filter(!is.na(TA) & !(codon %in% c("TAG","TGA","TAA"))) %>%
    dplyr::summarize(consump=sum(TA),type=unique(type)) %>%
    as.data.frame() %>%
    dplyr::mutate(RTSv = consump/wi)
  mclapply(mc.cores = 1,1:(consump$codon %>% unique() %>% length()),function(y){
    a <- unique(consump$codon)[y]
    #b <- mydf$tdt[which(mydf$exp=="Sc-Lys")]
    c <- consump %>% 
      dplyr::filter(codon==a) 
    c %>% group_by(codon) %>%
      dplyr::mutate(fcRTSv=c$RTSv[which(c$exp=="Sc-Lys")]/RTSv,fcrrt=c$rrt[which(c$exp=="Sc-Lys")]/rrt) %>%
      as.data.frame() %>% dplyr::filter(exp!="Sc-Lys")
  }) %>% rbind.fill() %>% group_by(exp) %>% 
    dplyr::summarize(rho=format(cor.test(fcRTSv,fcrrt,method = "s")$estimate[[1]], digits = 3) %>% as.numeric()) %>%
    as.data.frame() %>% cbind(x) 
}) %>% rbind.fill() %>% group_by(exp) %>% dplyr::summarize(mean=mean(rho),sd=sd(rho)) %>% as.data.frame()

data1D$exp <- factor(data1D$exp,levels = c("McManus2014", "Artieri2014","Pop2014", "Zinshteyn2013","Ingolia2009" ))
ggplot(data1D,aes(exp,mean)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),width=0.3)+
  labs(y=expression(paste("ρ (",RTS[v], " ~ decoding time)")))+
  style.print() + 
  #theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5,angle = 25))+ 
  #theme(legend.position = "none",axis.title.y=element_blank()) + 
  #scale_y_continuous(limits = c(0,0.4),breaks=c(0,0.1,0.2,0.3,0.4))+
  coord_flip()

##correlation of tdt~rtsv,rrt~rtsv of lower viral datasets
##rrt~rtsv figs2b
load("/mnt/data/home/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/18.2.RData")
tmp1 <- dat.allRRT %>% dplyr::filter(names %in% c("McManus2014", "Hunter2014", "Pop2014", "Zinshteyn2013", "Ingolia2009")) %>% as.data.frame()
allrrt <- tmp1[,c(2,4,10)]
allrrt$names[which(allrrt$names=="Hunter2014")] <- "Artieri2014"
load("/mnt/data/home/chenfeng/codonpaper/tAIcodyeast.Rdata")
tAIofcodonyeast$wi <- tAIofcodonyeast$Wi/max(tAIofcodonyeast$Wi)
tAIofcodonyeast$codon <- as.vector(tAIofcodonyeast$codon)
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Ingolia2009_remap/12.2.datIngolia.RData")
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Hunter2014/datHunter.RData");
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_McManus2014/datMcManus.RData");

datPop <- read.delim("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Pop2014/05.allData.tsv",header=F,as.is=T,sep="\t");
colnames(datPop) <- c("orf","chrom","chromPos","codon","pos","strand","ribo");

datZyn <- read.delim("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Zinshteyn2013/05.prepareDat.includeOffFrame.jitter2frame.tsv",header=F,as.is=T,sep="\t");
colnames(datZyn) <- c("orf","chrom","chromPos","codon","pos","strand","ribo")
lowyeast <- datIngolia[,c("orf","chrom","chromPos","codon","pos","strand","ribo")] %>% cbind(type="Ingolia2009") %>%
  rbind(datHunter %>% cbind(type="Artieri2014")) %>%
  rbind(datMcManus %>% cbind(type="McManus2014")) %>%
  rbind(datPop %>% cbind(type="Pop2014")) %>%
  rbind(datZyn %>% cbind(type="Zinshteyn2013")) 

lowyeast$wi <- tAIofcodonyeast$wi[match(lowyeast$codon,tAIofcodonyeast$codon)]

lowyeastRTSv <- lowyeast %>% dplyr::filter(substr(chrom,1,2)!="ch") %>%
  group_by(type,orf,chrom) %>%
  dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20)) %>%
  dplyr::mutate(TA=mean(ribo,trim=0.25,na.rm=T)) %>%
  as.data.frame() %>%
  group_by(type,codon,wi) %>%
  dplyr::summarize(consump=sum(TA)) %>%
  as.data.frame() %>%
  dplyr::mutate(RTSv=consump/wi)
lowyeastRTSv$type <- as.vector(lowyeastRTSv$type)
#save(lowyeastRTSv,file = "~/codonpaper/code and data/lowvirus.rrt.RTSv.Rdata")
load("~/codonpaper/code and data/lowvirus.rrt.RTSv.Rdata")
lowyeastRTSv$rrt <- allrrt$Pos6[match(paste(lowyeastRTSv$codon,lowyeastRTSv$type),paste(allrrt$codon,allrrt$names))]
source("~/Rfunction/style.print.R")

rrtrtsv <- lowyeastRTSv %>% group_by(type) %>% 
  dplyr::summarize(rho=cor.test(rrt,RTSv,method = "s")$estimate[1] %>% as.numeric(),
                   p=cor.test(rrt,RTSv,method = "s")$p.value %>% as.numeric())
rrtrtsv$type <- factor(rrtrtsv$type,levels=c("McManus2014", "Artieri2014", "Pop2014", "Zinshteyn2013", "Ingolia2009"))
ggplot(rrtrtsv,aes(type,rho)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  labs(x="",y="ρ(relative decoding time ~ relative tRNA \n shortages caused by viral translation)")+
  style.print() + 
  #theme(legend.position = "none",axis.title.y=element_blank()) + 
  coord_flip()

##rtsv-tdt figs4b
load("/mnt/data/home/phil/acting/virusTranslation/riboTranslocate/06.plot/03.RData")

tdt1 <- mclapply(mc.cores = 20,1:nrow(allRes),function(x){
  a <- allRes[x,2:4]
  a$type <-  sub("\\).*", "", sub(".*\\(", "", strsplit(a$exp %>% as.vector(),"~")[[1]][2]))
  a[,-1]
  }) %>% rbind.fill() 
tdt1 %>% dplyr::filter(type %in% c("McManus2014", "Artieri2014", "Pop2014", "Zinshteyn2013", "Ingolia2009")) -> tdt
lowyeastRTSv$tdt <- tdt$elongTime[match(paste(lowyeastRTSv$type,lowyeastRTSv$codon),paste(tdt$type,tdt$codon))]
tdtrtsv <- lowyeastRTSv %>% group_by(type) %>% 
  dplyr::summarize(rho=cor.test(tdt,RTSv,method = "s")$estimate[1] %>% as.numeric(),
                   p=cor.test(tdt,RTSv,method = "s")$p.value %>% as.numeric())
tdtrtsv$type <- factor(tdtrtsv$type,levels=c("McManus2014", "Artieri2014", "Pop2014", "Zinshteyn2013", "Ingolia2009"))
ggplot(tdtrtsv,aes(type,rho)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  labs(x="",y="ρ(typical decoding time ~ relative tRNA \n shortages caused by viral translation)")+
  style.print() + 
  #theme(legend.position = "none",axis.title.y=element_blank()) + 
  coord_flip()

## fig1e
load("/mnt/data/home/phil/acting/virusTranslation/01.scratch/04.result.yeast.RData")
p <- format(cor.test(dfCodon.yeast$RTSv,dfCodon.yeast$tdt,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(dfCodon.yeast$RTSv,dfCodon.yeast$tdt,method = "spearman")$estimate[[1]], digits = 3)
# fig3adata1 <- dfCodon.yeast %>% dplyr::filter(codon != "CGA")
# p1 <- format(cor.test(fig3adata1$RTSv,fig3adata1$tdt,method = "spearman")$p.value, digits = 3)
# rho1 <- format(cor.test(fig3adata1$RTSv,fig3adata1$tdt,method = "spearman")$estimate[[1]], digits = 3)
source("/mnt/data/home/chenfeng/Rfunction/style.print.R")
ggplot(dfCodon.yeast,aes(log(RTSv,10),tdt)) +
  labs(title=paste("Rho=",rho,", P=",p),x=expression(RTS[v]),y="TDT") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)+
  geom_point() +
  scale_x_continuous(limits = c(2,4),breaks=c(2,3,4))


##fig1C or s2a
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/Graph/13.RData");
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/Graph/13.3.RData");
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/Graph/18.2.slow_by_shortage.RRT.RData");

##rtsv=conofvirus/conofyeast
mydata <- plotCc2$data %>% as.data.frame()
load("/mnt/data/home/phil/acting/virusTranslation/01.scratch/datJustin.Sc-Lys.RData")
TA <- datJustin %>%
  group_by(orf) %>%
  dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20) & !(chrom %in% c("NC_003745","SCU01060"))) %>%
  dplyr::summarize(TA=mean(ribo,trim=0.25,na.rm=T)) %>%
  as.data.frame() 
datJustin$TA <- TA$TA[match(datJustin$orf,TA$orf)]
conofyeast <- datJustin %>% dplyr::filter(!(is.na(TA))) %>% group_by(codon) %>%
  dplyr::summarize(conofyeast=sum(TA)) %>% as.data.frame()
mydata$conofyeast <- conofyeast$conofyeast[match(mydata$codon,conofyeast$codon)]
mydata$rtsvofyeast <- mydata$cc/mydata$conofyeast

p <- format(cor.test(mydata$rtsvofyeast,mydata$rrt,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(mydata$rtsvofyeast,mydata$rrt,method = "spearman")$estimate[[1]], digits = 3)
ggplot(mydata,aes(log10(rtsvofyeast),rrt)) +
  labs(title=paste("Rho=",rho,", P=",p),x=expression(RTS[v]),y="Relative decoding time") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)+
  geom_point() +
  scale_x_continuous(limits = c(-3.3,-0.7),breaks=c(-3,-2,-1))

##rtsv=conofvirus/conofyeast figS4a
load("/mnt/data/home/phil/acting/virusTranslation/01.scratch/04.result.yeast.RData")
dfCodon.yeast$rtsvofyeast <- mydata$rtsvofyeast[match(as.vector(dfCodon.yeast$codon),mydata$codon)]
p <- format(cor.test(dfCodon.yeast$rtsvofyeast,dfCodon.yeast$tdt,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(dfCodon.yeast$rtsvofyeast,dfCodon.yeast$tdt,method = "spearman")$estimate[[1]], digits = 3)
# fig3adata1 <- dfCodon.yeast %>% dplyr::filter(codon != "CGA")
# p1 <- format(cor.test(fig3adata1$RTSv,fig3adata1$tdt,method = "spearman")$p.value, digits = 3)
# rho1 <- format(cor.test(fig3adata1$RTSv,fig3adata1$tdt,method = "spearman")$estimate[[1]], digits = 3)
source("/mnt/data/home/chenfeng/Rfunction/style.print.R")
ggplot(dfCodon.yeast,aes(log(rtsvofyeast,10),tdt)) +
  labs(title=paste("Rho=",rho,", P=",p),x=expression(RTS[v]),y="TDT") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)+
  geom_point() +
  scale_x_continuous(limits = c(-3.3,-0.7),breaks=c(-3,-2,-1))




