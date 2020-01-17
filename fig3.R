library(Biostrings)

load("~/codonpaper/newdncu.Rdata")
load("/mnt/data/home/chenfeng/fcs/mcherry_yfp_dncu.Rdata")
mcherry_yfp_dncu <- mcherry_yfp_dncu[,c(1,2,4,6,8,10)]
mcherry_yfp_dncu$newdncutrans <- newdncu$newdncutrans[match(mcherry_yfp_dncu$dncu,newdncu$dncu)]
mcherry_yfp_dncu$newdncutai <- newdncu$newdncutai[match(mcherry_yfp_dncu$dncu,newdncu$dncu)]

load(file = "/mnt/data/home/chenfeng/codonpaper/fcsbeforegating.Rdata")
fcsbeforegating$mCherry <- fcsbeforegating$mCherryA/fcsbeforegating$sizeA
fcsbeforegating$YFP <- fcsbeforegating$YFPA/fcsbeforegating$sizeA
fcsbeforegating$newdncutrans <- newdncu$newdncutrans[match(fcsbeforegating$dncu,newdncu$dncu)]
fcsbeforegating$newdncutai <- newdncu$newdncutai[match(fcsbeforegating$dncu,newdncu$dncu)]
fcs <-  fcsbeforegating
fcsbeforegating <- mclapply(mc.cores = 40,1:length(unique(fcs$dncu)),function(x){
  mydf <- fcs %>% dplyr::filter(dncu==unique(fcs$dncu)[x])
  tmp <- mydf$mCherry
  tmp <- tmp[order(tmp,decreasing = F)]
  a <- tmp[as.integer(0.005*length(tmp))]
  b <- tmp[as.integer(0.995*length(tmp))]
  mydf %>% dplyr::filter(mCherry > a & mCherry < b)
}) %>% rbind.fill()

tmp <- fcsbeforegating %>% arrange(mCherry)
mi <- floor(nrow(tmp)*0.8)
ma <- nrow(tmp)
mydf <- tmp[mi:ma,]
bar1 <- mydf %>% group_by(newdncutrans,newdncutai) %>% dplyr::summarize(YFP=mean(YFP)) %>% as.data.frame() %>%
  dplyr::summarize(pdncumean1=format(cor.test(newdncutrans,max(YFP)-YFP,method = "spearman")$p.value, digits = 3),
                   rhodncumean1=format(cor.test(newdncutrans,max(YFP)-YFP,method = "spearman")$estimate[[1]], digits = 3),
                   pdncumean2=format(cor.test(newdncutai,max(YFP)-YFP,method = "spearman")$p.value, digits = 3),
                   rhodncumean2=format(cor.test(newdncutai,max(YFP)-YFP,method = "spearman")$estimate[[1]], digits = 3)) %>%
  as.data.frame()
bar2_4 <- fcsbeforegating %>% group_by(newdncutrans,newdncutai) %>%
  dplyr::summarize(YFP=mean(YFP),mCherry=mean(mCherry)) %>% as.data.frame() %>%
  dplyr::summarize(rho4=format(cor.test(newdncutrans,(max(YFP)-YFP)/mCherry,method = "spearman")$estimate[[1]], digits = 3),
                   p4=format(cor.test(newdncutrans,(max(YFP)-YFP)/mCherry,method = "spearman")$p.value, digits = 3),
                   rho5=format(cor.test(newdncutrans,(max(YFP)-YFP),method = "spearman")$estimate[[1]], digits = 3),
                   p5=format(cor.test(newdncutrans,(max(YFP)-YFP),method = "spearman")$p.value, digits = 3),
                   rho6=format(cor.test(newdncutai,(max(YFP)-YFP)/mCherry,method = "spearman")$estimate[[1]], digits = 3),
                   p6=format(cor.test(newdncutai,(max(YFP)-YFP)/mCherry,method = "spearman")$p.value, digits = 3),
                   rho7=format(cor.test(newdncutai,(max(YFP)-YFP),method = "spearman")$estimate[[1]], digits = 3),
                   p7=format(cor.test(newdncutai,(max(YFP)-YFP),method = "spearman")$p.value, digits = 3)
  )
bar3 <- c(format(cor.test(mcherry_yfp_dncu$newdncutrans,mcherry_yfp_dncu$meanmCherrydensityA,method = "spearman")$estimate[[1]], digits = 3),
          format(cor.test(mcherry_yfp_dncu$newdncutrans,mcherry_yfp_dncu$meanmCherrydensityA,method = "spearman")$p.value, digits = 3),
          format(cor.test(mcherry_yfp_dncu$newdncutai,mcherry_yfp_dncu$meanmCherrydensityA,method = "spearman")$estimate[[1]], digits = 3),
          format(cor.test(mcherry_yfp_dncu$newdncutai,mcherry_yfp_dncu$meanmCherrydensityA,method = "spearman")$p.value, digits = 3)
)
##figs8a
tmp <- fcsbeforegating %>% arrange(mCherry)
mi <- 0
ma <- nrow(tmp)
n=100
step <- (ma/n) %>% floor()
col3 <- mclapply(mc.cores = 40,1:n-1,function(x){
  a <- mi+step*(x-1)
  b <- a+step
  mydf <- tmp[a:b,]
  mydf %>% group_by(newdncutrans,newdncutai) %>% dplyr::summarize(YFP=mean(YFP)) %>% as.data.frame() %>%
    dplyr::summarize(pdncumean1=format(cor.test(newdncutrans,max(YFP)-YFP,method = "spearman")$p.value, digits = 3),
                     rhodncumean1=format(cor.test(newdncutrans,max(YFP)-YFP,method = "spearman")$estimate[[1]], digits = 3),
                     pdncumean2=format(cor.test(newdncutai,max(YFP)-YFP,method = "spearman")$p.value, digits = 3),
                     rhodncumean2=format(cor.test(newdncutai,max(YFP)-YFP,method = "spearman")$estimate[[1]], digits = 3)) %>%
    as.data.frame() %>% cbind(x,strainnumber=length(unique(mydf$newdncutai)))
}) %>% rbind.fill()
col3$rhodncumean2 <- as.numeric(col3$rhodncumean2)
col3 %>% ggplot(aes(x=x,y=rhodncumean2))+
  geom_bar(stat = "identity")+style.print()+
  labs(x="100 groups divided by mCherry expression from low to high",y=expression(paste("ρ (translation load ~ ",italic("D")[P],")")))

save(col3,bar1,bar3,bar2_4,file = "~/codonpaper/code and data/allcombinedfig3.Rdata")

##figs6
library(tidyverse);
library(reshape2);
library(scales);
library(RColorBrewer)
load("/mnt/data/home/chenfeng/codonpaper/review/fig3c_f.Rdata")

dfFacs <- figs3c_f %>%
  mutate(txLoad = max(YFP) - YFP,loadpermCherry=(max(YFP) - YFP)/mCherry) %>%
  rename(dpTai = newdncutai, dpTx = newdncutrans) %>%
  mutate(dpTaiRnk = ordered(dpTai,levels = sort(unique(dpTai))),
         dpTxRnk = ordered(dpTx,levels = sort(unique(dpTx))))

dfTry <- dfFacs 
# %>%
#   sample_n(50000);

discrete_gradient_pal <- function(colours, bins = 5) {
  ramp <- scales::colour_ramp(colours)
  
  function(x) {
    if (length(x) == 0) return(character())
    
    i <- floor(x * bins)
    i <- ifelse(i > bins-1, bins-1, i)
    ramp(i/(bins-1))
  }
}

scale_fill_discrete_gradient <- function(..., colours, bins = 5, na.value = "grey50", guide = "colourbar", aesthetics = "color", colors)  {
  colours <- if (missing(colours)) 
    colors
  else colours
  continuous_scale(
    aesthetics,
    "discrete_gradient",
    discrete_gradient_pal(colours, bins),
    na.value = na.value,
    guide = guide,
    ...
  )
}
useBreaks <- seq(0.1,0.9,by=0.1);
#C
source("~/Rfunction/style.print.R")
figC <- (dfTry %>% arrange(desc(mCherry)))[1:ceiling(0.2*nrow(dfTry)),] 
figC %>% group_by(dpTai) %>% dplyr::mutate(n=length(dpTai)) %>% dplyr::filter(n>50000) -> tt

#tt %>% group_by(dpTai) %>% dplyr::summarize(meanmcherry=mean(mCherry),semcherry=sd(mCherry)/(length(mCherry)^0.5),meanload=mean(txLoad),seload=sd(txLoad)/(length(txLoad)^0.5)) -> aaa


tt %>%
  mutate(dpTaiLvl = dpTai) %>%
  ggplot(aes(x=dpTaiRnk,y=txLoad,fill=dpTaiLvl)) +
  geom_violin(color="NA",scale="width",adjust=0.5,trim=F) +
  scale_y_continuous("Translation load \n (cells with top 20% mCherry expression)",limits = c(0.15,0.3),trans = "log10",breaks = c(0.15,0.2,0.3)) +
  scale_fill_discrete_gradient(expression(paste(italic(D)[P]," value")),
                               colours = brewer.pal(8,"Spectral"),
                               bins=8,
                               guide=guide_colorbar(raster=F,nbin=100,ticks.colour = NA),
                               breaks=c(0.1,0.5,0.9),
                               limits=c(0.1,0.9)) +
  scale_x_discrete(expression(paste("Rank of ",italic(D)[P])),labels=NULL) +
  geom_boxplot(width = 0.3, outlier.colour = NA, aes(fill = dpTaiLvl)) + 
  stat_summary(fun.y = 'mean', geom = 'point', shape = 18, colour = 'black')+
  style.print()+
  theme(legend.position="none")+
  #theme_classic()+
  theme(axis.ticks.x = element_blank())
#D
dfTry %>%
  mutate(dpTaiLvl = dpTai) %>%
  ggplot(aes(x=dpTaiRnk,y=loadpermCherry,fill=dpTaiLvl)) +
  geom_violin(color="NA",scale="width",adjust=0.5,trim=F) +
  scale_y_continuous("Translational load \n per unit of mCherry expression",trans = "log10") +
  scale_fill_discrete_gradient(expression(paste(italic(D)[P]," value")),
                               colours = brewer.pal(8,"Spectral"),
                               bins=8,
                               guide=guide_colorbar(raster=F,nbin=100,ticks.colour = NA),
                               breaks=c(0.1,0.5,0.9),
                               limits=c(0.1,0.9)) +
  scale_x_discrete(expression(paste("Rank of ",italic(D)[P])),labels=NULL) +
  geom_boxplot(width = 0.3, outlier.colour = NA, aes(fill = dpTaiLvl)) + 
  stat_summary(fun.y = 'mean', geom = 'point', shape = 18, colour = 'black')+
  style.print()+
  theme(legend.position="none")+
  theme(axis.ticks.x = element_blank())
#E
dfTry %>%
  mutate(dpTaiLvl = dpTai) %>%
  ggplot(aes(x=dpTaiRnk,y=mCherry,fill=dpTaiLvl)) +
  geom_violin(color="NA",scale="width",adjust=0.5,trim=F) +
  scale_y_continuous("mCherry expression",trans = "log10") +
  scale_fill_discrete_gradient(expression(paste(italic(D)[P]," value")),
                               colours = brewer.pal(8,"Spectral"),
                               bins=8,
                               guide=guide_colorbar(raster=F,nbin=100,ticks.colour = NA),
                               breaks=c(0.1,0.5,0.9),
                               limits=c(0.1,0.9)) +
  scale_x_discrete(expression(paste("Rank of ",italic(D)[P])),labels=NULL) +
  geom_boxplot(width = 0.3, outlier.colour = NA, aes(fill = dpTaiLvl)) + 
  stat_summary(fun.y = 'mean', geom = 'point', shape = 18, colour = 'black')+
  style.print()+
  theme(legend.position="none")+
  theme(axis.ticks.x = element_blank())
#F
dfTry %>%
  mutate(dpTaiLvl = dpTai) %>%
  ggplot(aes(x=dpTaiRnk,y=txLoad,fill=dpTaiLvl)) +
  geom_violin(color="NA",scale="width",adjust=0.5,trim=F) +
  scale_y_continuous("Translational load",limit=c(0.16,0.285),breaks = c(0.16,0.20,0.24,0.28)) +
  scale_fill_discrete_gradient(expression(paste(italic(D)[P]," value")),
                               colours = brewer.pal(8,"Spectral"),
                               bins=8,
                               guide=guide_colorbar(raster=F,nbin=100,ticks.colour = NA),
                               breaks=c(0.1,0.5,0.9),
                               limits=c(0.1,0.9)) +
  scale_x_discrete(expression(paste("Rank of ",italic(D)[P])),labels=NULL) +
  geom_boxplot(width = 0.3, outlier.colour = NA, aes(fill = dpTaiLvl)) + 
  stat_summary(fun.y = 'mean', geom = 'point', shape = 18, colour = 'black')+
  style.print()+
  theme(legend.position="none")+
  theme(axis.ticks.x = element_blank())

##Rtpcr plot
load("~/codonpaper/review/rtdata.Rdata")
rtdata$newdptai <- as.numeric(rtdata$newdptai)
rtdata %>% dplyr::filter(type=="mCherry/ACTB ratio") %>% arrange(newdptai) %>% dplyr::mutate(dpTaiRnk=1:36) %>%
  mutate(dpTaiLvl = newdptai) -> tttt 
tttt %>%
  ggplot(aes(x=dpTaiRnk,y=mean,color=dpTaiLvl)) + 
  geom_point(size=3)+
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),color="grey")+
  labs(y="mCherry/actin ratio") +
  scale_fill_discrete_gradient(expression(paste(italic(D)[P]," value")),
                               colours = brewer.pal(8,"Spectral"),
                               bins=8,
                               guide=guide_colorbar(raster=F,nbin=100,ticks.colour = NA),
                               breaks=c(0.1,0.5,0.9),
                               limits=c(0.1,0.9)) +
  scale_x_discrete(expression(paste("Rank of ",italic(D)[P])),labels=NULL) +
  style.print()+
  theme(legend.position="none")+
  theme(axis.ticks.x = element_blank())

rtdata %>% dplyr::filter(type=="YFP/ACTB ratio") %>% arrange(newdptai) %>% dplyr::mutate(dpTaiRnk=1:37) %>%
  mutate(dpTaiLvl = newdptai) -> tttt2 
tttt2 %>%
  ggplot(aes(x=dpTaiRnk,y=mean,color=dpTaiLvl)) + 
  geom_point(size=3)+
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se),color="grey")+
  labs(y="YFP/actin ratio") +
  scale_fill_discrete_gradient(expression(paste(italic(D)[P]," value")),
                               colours = brewer.pal(8,"Spectral"),
                               bins=8,
                               guide=guide_colorbar(raster=F,nbin=100,ticks.colour = NA),
                               breaks=c(0.1,0.5,0.9),
                               limits=c(0.1,0.9)) +
  scale_x_discrete(expression(paste("Rank of ",italic(D)[P])),labels=NULL) +
  ylim(0,0.6)+
  style.print()+
  theme(legend.position="none")+
  theme(axis.ticks.x = element_blank())


save(figs3c_f,rtdata,file = "~/codonpaper/code and data/allcombined.figs6.Rdata")















