##virus TA in yeast, Fig.1a
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/mylist.RData");
##calculated the fold change of TA  
yeastvirussite <- mclapply(mc.cores = 4,1:9,function(x){
  mydata <- mylist[[x]] %>%
      group_by(orf,chrom) %>%
      dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20)) %>%
      dplyr::summarize(TA=mean(ribo,trim=0.25,na.rm=T)) %>%
      as.data.frame()
  fc <- ((mydata %>% dplyr::filter(substr(chrom,1,2)!="ch"))$TA %>% mean()) / ((mydata %>% dplyr::filter(substr(chrom,1,2)=="ch"))$TA %>% mean())
  data.frame(fc,name=names(mylist[x]))
}) %>% rbind.fill()
yeastvirussite$per <- yeastvirussite$vTA/yeastvirussite$yTA

##coverage. Fig.S1
yeastvirus <- mylist[[1]] %>% cbind(type="Sc-His") %>%
  rbind(mylist[[2]] %>% cbind(type="Sc-Lys")) %>%
  rbind(mylist[[3]] %>% cbind(type="YPD1")) %>%
  rbind(mylist[[4]] %>% cbind(type="YPD2")) %>%
  rbind(mylist[[5]] %>% cbind(type="Ingolia2009")) %>%
  rbind(mylist[[6]] %>% cbind(type="Artieri2014")) %>%
  rbind(mylist[[7]] %>% cbind(type="McManus2014")) %>%
  rbind(mylist[[8]] %>% cbind(type="Pop2014")) %>%
  rbind(mylist[[9]] %>% cbind(type="Zinshteyn2013")) 
seqdepth <- yeastvirus %>% 
  group_by(orf,type) %>%
  dplyr::filter(substr(chrom,1,2)=="ch") %>%
  dplyr::filter((length(codon) > 70) & (pos > 20) & (pos < (length(codon)-20))) %>%
  dplyr::mutate(mfj=mean(ribo,trim=0.25,na.rm=T)) %>%
  as.data.frame() %>%
  group_by(type) %>%
  dplyr::summarize(all=sum(mfj)) %>%
  as.data.frame()
##nromalized by sequence depth
virus <- yeastvirus %>% dplyr::filter(substr(chrom,1,2)!="ch") %>%
  dplyr::mutate(depth=seqdepth$all[match(as.vector(type),as.vector(seqdepth$type))]) %>%
  dplyr::mutate(ribonor=ribo*1000000/depth) %>% as.data.frame()
virus$type <- factor(virus$type,levels=c("Sc-Lys","Sc-His","YPD1","YPD2","Ingolia2009","Zinshteyn2013","Pop2014","Artieri2014","McManus2014"))
source("~/Rfunction/style.print.R")
##do the picture
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

#Fig.S2a and S4a 
load("~/codonpaper/review/1new.combined.rtsv.04.rr.Rdata")
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
##do the picture
data.frame(mean=c(rhoce$e,rhoce$c),sd=c(rhoce$sde,rhoce$sdc),type=c("a","c")) %>%
  ggplot(aes(type,mean)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),width=0.3)+
  labs(y=expression(paste("ρ (",RTS[v], " ~ decoding time)")))+
  style.print() +
  theme(legend.position = "none",axis.title.y=element_blank()) + 
  coord_flip()

##Fig.S2c
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
##draw the picture
data1D$exp <- factor(data1D$exp,levels = c("McManus2014", "Artieri2014","Pop2014", "Zinshteyn2013","Ingolia2009" ))
ggplot(data1D,aes(exp,mean)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),width=0.3)+
  labs(y=expression(paste("ρ (",RTS[v], " ~ decoding time)")))+
  style.print() + 
  coord_flip()

##correlation of tdt~rtsv,rrt~rtsv of lower viral datasets
##rrt~rtsv Fig.S2b
load("/mnt/data/home/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/18.2.RData")
tmp1 <- dat.allRRT %>% dplyr::filter(names %in% c("McManus2014", "Hunter2014", "Pop2014", "Zinshteyn2013", "Ingolia2009")) %>% as.data.frame()
allrrt <- tmp1[,c(2,4,10)]
allrrt$names[which(allrrt$names=="Hunter2014")] <- "Artieri2014"
load("/mnt/data/home/chenfeng/codonpaper/tAIcodyeast.Rdata")
tAIofcodonyeast$wi <- tAIofcodonyeast$Wi/max(tAIofcodonyeast$Wi)
tAIofcodonyeast$codon <- as.vector(tAIofcodonyeast$codon)
lowyeast <- mylist[[5]] %>% cbind(type="Ingolia2009") %>%
  rbind(mylist[[6]] %>% cbind(type="Artieri2014")) %>%
  rbind(mylist[[7]] %>% cbind(type="McManus2014")) %>%
  rbind(mylist[[8]] %>% cbind(type="Pop2014")) %>%
  rbind(mylist[[9]] %>% cbind(type="Zinshteyn2013")) 
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
load("~/codonpaper/code and data/lowvirus.rrt.RTSv.Rdata")
lowyeastRTSv$rrt <- allrrt$Pos6[match(paste(lowyeastRTSv$codon,lowyeastRTSv$type),paste(allrrt$codon,allrrt$names))]
rrtrtsv <- lowyeastRTSv %>% group_by(type) %>% 
  dplyr::summarize(rho=cor.test(rrt,RTSv,method = "s")$estimate[1] %>% as.numeric(),
                   p=cor.test(rrt,RTSv,method = "s")$p.value %>% as.numeric())
##draw the picture
rrtrtsv$type <- factor(rrtrtsv$type,levels=c("McManus2014", "Artieri2014", "Pop2014", "Zinshteyn2013", "Ingolia2009"))
ggplot(rrtrtsv,aes(type,rho)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  labs(x="",y="ρ(relative decoding time ~ relative tRNA \n shortages caused by viral translation)")+
  style.print() +  
  coord_flip()

##rtsv-tdt Fig.S4b
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
##draw the picture
tdtrtsv$type <- factor(tdtrtsv$type,levels=c("McManus2014", "Artieri2014", "Pop2014", "Zinshteyn2013", "Ingolia2009"))
ggplot(tdtrtsv,aes(type,rho)) + 
  geom_bar(stat = "identity",fill=gray.colors(11)[7]) +
  labs(x="",y="ρ(typical decoding time ~ relative tRNA \n shortages caused by viral translation)")+
  style.print() + 
  #theme(legend.position = "none",axis.title.y=element_blank()) + 
  coord_flip()

## Fig.1e
load("/mnt/data/home/phil/acting/virusTranslation/01.scratch/04.result.yeast.RData")
p <- format(cor.test(dfCodon.yeast$RTSv,dfCodon.yeast$tdt,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(dfCodon.yeast$RTSv,dfCodon.yeast$tdt,method = "spearman")$estimate[[1]], digits = 3)
source("/mnt/data/home/chenfeng/Rfunction/style.print.R")
##draw the picture
ggplot(dfCodon.yeast,aes(log(RTSv,10),tdt)) +
  labs(title=paste("Rho=",rho,", P=",p),x=expression(RTS[v]),y="TDT") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)+
  geom_point() +
  scale_x_continuous(limits = c(2,4),breaks=c(2,3,4))

##Fig.1C
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/Graph/13.RData");
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/Graph/13.3.RData");
load("/mnt/data3/disk/phil/acting/virusTranslation/riboTranslocate/data/riboSeq_Justin2014/Graph/18.2.slow_by_shortage.RRT.RData");
mydata <- plotCc2$data %>% as.data.frame()
load("/mnt/data/home/phil/acting/virusTranslation/01.scratch/datJustin.Sc-Lys.RData")
fig1c <- datJustin %>%
  group_by(orf) %>%
  dplyr::filter(length(codon)>70 & pos > 20 & pos < (length(codon)-20) & (chrom %in% c("NC_003745","SCU01060"))) %>%
  dplyr::mutate(TA=mean(ribo,trim=0.25,na.rm=T)) %>% 
  group_by(codon) %>% 
  dplyr::summarize(consump=sum(TA)) %>%
  dplyr::mutate(tai = tAIofcodonyeast$wi[match(codon,tAIofcodonyeast$codon)],
                rrt = mydata$rrt[match(codon,mydata$codon)]) %>%
  dplyr::mutate(RTSv=consump/tai) %>%
  as.data.frame()
##draw the picture
p <- format(cor.test(fig1c$RTSv,fig1c$rrt,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(fig1c$RTSv,fig1c$rrt,method = "spearman")$estimate[[1]], digits = 3)
ggplot(fig1c,aes(RTSv,rrt)) +
  labs(title=paste("Rho=",rho,", P=",p),x=expression(RTS[v]),x="Relative tRNA shortages \n caused by viral translation",y="Relative decoding time") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)+
  geom_point() +
  scale_x_continuous(trans = "log10",limits = c(100,10000),breaks=c(100,1000,10000))

##rtsv=conofvirus/conofyeast Fig.S4a
load("/mnt/data/home/phil/acting/virusTranslation/01.scratch/04.result.yeast.RData")
dfCodon.yeast$rtsvofyeast <- mydata$rtsvofyeast[match(as.vector(dfCodon.yeast$codon),mydata$codon)]
p <- format(cor.test(dfCodon.yeast$rtsvofyeast,dfCodon.yeast$tdt,method = "spearman")$p.value, digits = 3)
rho <- format(cor.test(dfCodon.yeast$rtsvofyeast,dfCodon.yeast$tdt,method = "spearman")$estimate[[1]], digits = 3)
source("/mnt/data/home/chenfeng/Rfunction/style.print.R")
ggplot(dfCodon.yeast,aes(log(rtsvofyeast,10),tdt)) +
  labs(title=paste("Rho=",rho,", P=",p),x=expression(RTS[v]),y="TDT") +
  style.print()+
  geom_smooth(method="lm",se=FALSE)+
  geom_point() +
  scale_x_continuous(limits = c(-3.3,-0.7),breaks=c(-3,-2,-1))
