library(plyr)
library(dplyr)
library(parallel)
library(readr)
library(LncFinder)
load("~/part3/new.Mh.test_plot.Rdata")
associate <- Mh.test_plot[,1:3]
names(associate) <- c("virus","symptomatic","natural")
host <- unique(c(associate$symptomatic %>% as.vector(),associate$natural %>% as.vector()))
load("~/codonpaper/review/allhumangenome.Rdata")
##xij of host gene
load("~/codonpaper/review/codon.Rdata")
xijofhost <- mclapply(mc.cores = 17,1:17,function(i){
  special <- unique(allhumangenome$spe)[i]
  myseq <- allhumangenome %>% dplyr::filter(spe==special)
  express <- read.csv(paste("~/codonpaper/review/transcriptome/express/",special,".abundance",sep = ""),sep = "\t",stringsAsFactors = F)
  myseq$express <- express$TPM[match(myseq$name1,express$Reference)]
  myseq <- (myseq %>% dplyr::filter(!is.na(express)) %>% arrange(desc(express)))[1:100,]
  xij <- mclapply(mc.cores = 4,1:nrow(myseq),function(yy){
    mydf <- myseq[yy,]
    b <- strsplit(mydf$seq,"_")[[1]] %>% table() %>% as.data.frame()
    names(b)[1] <- "codon"
    b$codon <- as.vector(b$codon)
    b$aa <- codon$aa[match(as.vector(b$codon),codon$codon)]
    b %>% dplyr::filter(!is.na(aa) & !(aa %in% c("W","M"))) %>% cbind(express=mydf$express)
  }) %>% rbind.fill() %>%
    group_by(codon,aa) %>%
    dplyr::summarize(fre=sum(Freq)) %>%
    as.data.frame() %>%
    group_by(aa) %>%
    dplyr::mutate(sumaa=sum(fre)) %>%
    dplyr::mutate(xij=fre/sumaa) %>%
    as.data.frame() %>% cbind(special)
  xij$special <- as.vector(xij$special)
  xij
}) %>% rbind.fill()
save(xijofhost,file = "~/codonpaper/code and data/1_2.xijofhost.finall.Rdata")
load("~/codonpaper/code and data/1_2.xijofhost.finall.Rdata")
## dp and or
VNS_DP_OR <- mclapply(mc.cores = 52,1:nrow(associate),function(x){
  asso <- associate[x,] 
  virusseq <- readDNAStringSet(paste("~/codonpaper/allcds/",as.vector(asso$virus),".txt",sep=""))
  yij <- mclapply(mc.cores=1,1:length(virusseq),function(y){
    seq <- substring(as.character(virusseq[[y]]),seq(1,(length(virusseq[[y]])-2),by=3),seq(3,length(virusseq[[y]]),by=3));
    a <- as.data.frame(table(seq));
    names(a)[1] <- c("codon");
    a$codon <- as.vector(a$codon)
    a$amino <- codon$aa[match(as.vector(a$codon),codon$codon)]
    a %>% dplyr::filter(!is.na(amino) & amino != "W")
  }) %>% rbind.fill() %>%
    group_by(codon,amino) %>%
    dplyr::summarize(fre=sum(Freq)) %>%
    as.data.frame() %>%
    group_by(amino) %>%
    dplyr::mutate(sumaa = sum(fre)) %>%
    dplyr::mutate(yij = fre/sumaa) %>%
    as.data.frame()
  syn <- xijofhost %>% dplyr::filter(special==as.vector(asso$symptomatic)) 
  nat <- xijofhost %>% dplyr::filter(special==as.vector(asso$natural))
  yij$xijsyn <- syn$xij[match(yij$codon,syn$codon)]  
  yij$xijnat <- nat$xij[match(yij$codon,nat$codon)]
  #DP
  Di <- yij %>% group_by(amino,sumaa) %>%
    dplyr::summarize(Divn=(sum((xijnat-yij)^2))^0.5,
                     Divs=(sum((xijsyn-yij)^2))^0.5) %>%
    as.data.frame()
  Dpvn <- prod(Di$Divn)^(1/nrow(Di)) 
  Dpvs <- prod(Di$Divs)^(1/nrow(Di)) 
  # #OR
  tmp <- mclapply(mc.cores=4,1:length(unique(yij$amino)),function(k){
    tmp <- unique(yij$amino)[k]
    myDf <- yij %>% dplyr::filter(amino == tmp)
    rowTrueOrFalse <- abs(myDf[,5]-myDf[,6]) < abs(myDf[,7]-myDf[,5]);
    colTrueOfFalse <- myDf[,5] > (1/length(myDf$codon));
    if(length(myDf$codon) <=1){return(data.frame(NULL));}
    data.frame(amino = myDf$amino[1],
               result = c(
                 length(which(rowTrueOrFalse & colTrueOfFalse)),
                 length(which(!rowTrueOrFalse & colTrueOfFalse)),
                 length(which(rowTrueOrFalse & !colTrueOfFalse)),
                 length(which(!rowTrueOrFalse & !colTrueOfFalse))))
  }) %>% rbind.fill()
  tmp2 <- array(tmp$result,dim = c(2,2,18))
  tmp3 <- mantelhaen.test(tmp2)       
  data.frame(Dpvn,Dpvs,P=tmp3$p.value,OR=as.numeric(tmp3$estimate)) %>% cbind(asso)
}) %>% rbind.fill()
save(VNS_DP_OR,file = "~/codonpaper/code and data/fig4.vns_realtranscriptome.Rdata")
load("~/codonpaper/code and data/fig4.vns_realtranscriptome.Rdata")
##group by translation selection
##Dp
max(VNS_DP_OR$Dpvn+VNS_DP_OR$Dpvs)
min(VNS_DP_OR$Dpvn+VNS_DP_OR$Dpvs)
tmp = seq(0.9,0.3,by=-.2)
alldata <- mclapply(mc.cores = 20,1:length(tmp),function(x){
  VNS_DP_OR %>% dplyr::filter(Dpvn+Dpvs <= tmp[x]) %>% 
    dplyr::summarize(fractionofcha=length(which(Dpvs-Dpvn<0))/length(Dpvn),
                     meanofdpvnjiandpvs=mean(Dpvn/Dpvs),
                     se=sd(Dpvn/Dpvs)/(length(Dpvn)^0.5),
                     n=length(Dpvn)) %>% as.data.frame() %>% cbind(type=paste("<",tmp[x]))
}) %>% rbind.fill()
source("~/Rfunction/style.print.R")
ggplot(data = alldata,aes(x=type,y=meanofdpvnjiandpvs))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymax=meanofdpvnjiandpvs+se,ymin=meanofdpvnjiandpvs-se),width=0.4)+
  labs(x=expression(italic("D")[P](V,S)+italic("D")[P](V,N)),
       y=expression(paste("Average ",italic("D")[P](V,N)/italic("D")[P](V,S)))) +
  #scale_y_continuous(limits = c(0,.55)) +
  style.print()
##or
allor <- mclapply(mc.cores = 4,1:length(tmp),function(x){
  associate <- VNS_DP_OR %>% dplyr::filter(Dpvn+Dpvs < tmp[x])
  #names(associate) <- c("virus","symptomatic","natural")
  tmp1 <- nrow(associate)
  mclapply(mc.cores = 2,1:1000, function(i){
    tmp2 <- sample(1:tmp1,tmp1,replace = T)
    mydf <- associate[tmp2,]
    #mydf$OR <- orA$OR[match(as.vector(mydf$virus),as.vector(orA$virus))]  
    mydf %>% dplyr::summarize(n=length(which(OR>1))/tmp1) %>% 
      as.data.frame() %>% cbind(type=paste("<",tmp[x]),type2=paste(tmp2,collapse = ""))
  }) %>% rbind.fill() %>% group_by(type) %>%
    dplyr::summarize(mean=mean(n),sd=sd(n)) %>% 
    as.data.frame()
}) %>% rbind.fill()
ggplot(data = allor,aes(x=type,y=mean))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),width=0.4)+
  labs(x=expression(italic("D")[P](V,S)+italic("D")[P](V,N)),
       y="Fraction of VNS-trio with OR > 1") +
  style.print()

##fig4c and figs9
library(readr);
library(plyr);
library(dplyr);
library(parallel);
library(Biostrings)
virusseq <- readDNAStringSet(filepath = "~/part4/virusseq.fas") 
alldncu1 <- mclapply(mc.cores=20,1:length(virusseq),function(x){
  seq <- substring(as.character(virusseq[[x]]),seq(1,(length(virusseq[[x]])-2),by=3),seq(3,length(virusseq[[x]]),by=3));
  dn <- data.frame(human1=dncu(seq)[1],wenzi1=dncu(seq)[2],human2=dncu(seq)[3],wenzi2=dncu(seq)[4])
  dn$gene <- strsplit(names(virusseq[x]),":")[[1]][1];
  dn
}) %>% rbind.fill()
dengue <- read.csv("/mnt/data/home/chenfeng/part4/newdeng2/dengue.txt")[,c(1,3,6,7,8)]
zika <- read.csv("/mnt/data/home/chenfeng/part4/newdeng2/zika.txt")[,c(1,4,5,6)]
dengue$virusname <- paste("Dengue virus type",dengue$type)
zika$virusname <- "Zika virus"
virus <- dengue[,-2] %>% rbind(zika)
mydncu <- alldncu1 %>% dplyr::filter(gene %in% as.vector(virus$accession))
mydncu$virusname <- virus$virusname[match(mydncu$gene,as.vector(virus$accession))]
mydncu$host <- virus$host[match(mydncu$gene,as.vector(virus$accession))]
hostype <- data.frame(host=mydncu$host %>% as.vector() %>% unique(),
                      type=c("Symptomatic host","Unknown","Unknown","Natural host","Symptomatic host","Natural host","Natural host","Symptomatic host","Symptomatic host",
                             rep("Natural host",5),"Symptomatic host", rep("Natural host",4)))
mydncu$hostype <- hostype$type[match(mydncu$host,hostype$host)]
mydncu$hostype <- factor(mydncu$hostype, levels=c("Symptomatic host","Unknown","Natural host"))
##Dp calculated by transcriptome
Dncu1 <- mydncu %>%  ggplot() +
  geom_point(aes(x=human1,y=wenzi1,color=hostype),size=0.5) +
  geom_point(data = base :: subset(mydncu, hostype == 'Natural host'),aes(x=human1,y=wenzi1,color=hostype),size=0.5) +
  labs(x=expression(paste(italic("D")[P]," (virus to human)")) , y=expression(paste(italic("D")[P]," (virus to mosquito)")))+
  #geom_vline(aes(xintercept=symnatdscu,colour=natname),linetype="dashed")+
  geom_abline(intercept=0,slope=1, linetype="dotted")+
  expand_limits(y=0.06) +
  scale_x_continuous(limits = c(0.06,0.3),breaks=c(0.1,0.15,0.2,0.25))+
  facet_wrap(~virusname,nrow = 2)+
  style.print()+
  scale_colour_manual(values=c("#6495ED","#5B88A0","red"))+
  theme(legend.position = "none")
mydncu %>% group_by(virusname) %>% dplyr::summarize(percent=length(which(wenzi1/human1>1))/length(wenzi1))
##Dp calculated by tAI of codons
Dncu2 <- mydncu %>%  ggplot() +
  geom_point(aes(x=human2,y=wenzi2,color=hostype)) +
  geom_point(data = base :: subset(mydncu, hostype == 'Natural host'),aes(x=human2,y=wenzi2,color=hostype)) +
  labs(x=expression(paste(italic("D")[P]," (virus to human)")) , y=expression(paste(italic("D")[P]," (virus to mosquito)")))+
  #geom_vline(aes(xintercept=symnatdscu,colour=natname),linetype="dashed")+
  geom_abline(intercept=0,slope=1, linetype="dotted")+
  expand_limits(y=0.1) +
  scale_x_continuous(limits = c(0.1,0.4),breaks=c(0.1,0.15,0.2,0.25))+
  facet_wrap(~virusname,nrow = 2)+
  style.print()+
  scale_colour_manual(values=c("#6495ED","#5B88A0","red"))+
  theme(legend.position = "none")
##print function
style.print <- function() {
  theme_classic() +
    theme(legend.text=element_text(size=unit(12,"bigpts")),
          legend.key.size=unit(12,"bigpts"),
          legend.key=element_blank(),
          legend.title=element_text(size=unit(12,"bigpts"),face="bold"),
          axis.title=element_text(size=unit(14,"bigpts")),
          axis.text=element_text(size=unit(12,"bigpts")),
          plot.title=element_text(hjust=-0.1,size=unit(16,"bigpts"),face="bold"),
          strip.text=element_text(size=unit(14,"bigpts")),
          panel.background=element_rect(fill="white"));
}
### dncu function, it means Dp in the NEE paper
library(dplyr);
library(plyr);
library(Biostrings);
load("/mnt/data/home/chenfeng/part4/codonfreq.Rdata")
##read and splite virus sequence to three nucletide element
##origin codon frequency in human, wenzi, virus
codonfre <- as.data.frame(codonfreq$A549)
codonfre$frewenzi <- codonfreq$Aedes$fraction[match(as.vector(codonfre$codon),as.vector(codonfreq$Aedes$codon))]
codonfre$frevirus <- codonfreq$Dengue_virus_2$fraction[match(as.vector(codonfre$codon),as.vector(codonfreq$Dengue_virus_2$codon))]
codonfreorigin <- codonfre[,-(3:4)]
names(codonfreorigin)[3] <- "frehuman"
##codon list table
codontable0 <- codonfreorigin[,(1:2)]
codontable1 <- data.frame(aa=c("M","W"),
                          codon=c("ATG","TGG"))
codontable <- codontable0 %>% rbind(codontable1)
load("/mnt/data/home/chenfeng/codonpaper/tAIcodhum.Rdata")
load("/mnt/data/home/chenfeng/codonpaper/tAIofcodonAno.Rdata")
codon <- read.table("~/fcs/codon.txt",header = TRUE,stringsAsFactors=F ) %>% rbind(data.frame(aa="M",codon="ATG"))
tAIofcodonhuman <- tAIofcodonhuman %>% dplyr::filter(Wi>0)
tAIofcodonhuman$aa <- codon$aa[match(as.vector(tAIofcodonhuman$codon),codon$codon)]
tAIofcodonhuman$wi <- tAIofcodonhuman$Wi/max(tAIofcodonhuman$Wi)
meanhuman <- tAIofcodonhuman %>% group_by(aa) %>% dplyr::summarize(all=sum(wi))
tAIofcodonhuman$all <- meanhuman$all[match(tAIofcodonhuman$aa,meanhuman$aa)]
tAIofcodonhuman$tai <- tAIofcodonhuman$wi/tAIofcodonhuman$all
tAIofcodonAno$aa <- codon$aa[match(as.vector(tAIofcodonAno$codon),codon$codon)]
tAIofcodonAno <- tAIofcodonAno %>% dplyr::filter(!(is.na(aa)) & !(aa %in% c("M","W")))
tAIofcodonAno$wi <- tAIofcodonAno$Wi/max(tAIofcodonAno$Wi)
meanAno <- tAIofcodonAno %>% group_by(aa) %>% dplyr::summarize(all=sum(wi))
tAIofcodonAno$all <- meanAno$all[match(tAIofcodonAno$aa,meanAno$aa)]
tAIofcodonAno$tai <- tAIofcodonAno$wi/tAIofcodonAno$all
dncu <- function(geneseq){
  seqnew <- geneseq[which(geneseq %in% as.vector(codontable0$codon))]
  virusfre <- as.data.frame(table(seqnew))
  virusfre$aa <- codontable0$aa[match(as.vector(virusfre$seqnew),as.vector(codontable0$codon))]
  virusaa <- as.data.frame(virusfre %>% group_by(aa) %>% dplyr::summarise(sum=sum(Freq)))
  virusfre$frevirus <- virusfre$Freq/virusaa$sum[match(virusfre$aa,virusaa$aa)]
  codonfrenew <- codonfreorigin[-5]
  codonfrenew$frevirus <- virusfre$frevirus[match(as.vector(codonfrenew$codon),as.vector(virusfre$seqnew))]
  codonfrenew$taihum <- tAIofcodonhuman$tai[match(as.vector(codonfrenew$codon),as.vector(tAIofcodonhuman$codon))]
  codonfrenew$taiAno <- tAIofcodonAno$tai[match(as.vector(codonfrenew$codon),as.vector(tAIofcodonAno$codon))]
  Di_aa <- codonfrenew %>% group_by(aa) %>%
    dplyr::summarize(aa_human1 = (sum((frehuman-frevirus)^2))^0.5,
                     aa_wenzi1 = (sum((frewenzi-frevirus)^2))^0.5,
                     aa_human2 = (sum((taihum-frevirus)^2))^0.5,
                     aa_wenzi2 = (sum((taiAno-frevirus)^2))^0.5) %>%
    as.data.frame()
  humandncu1 <- (Di_aa$aa_human1 %>% prod())^(1/nrow(Di_aa))
  wenzidncu1 <- (Di_aa$aa_wenzi1 %>% prod())^(1/nrow(Di_aa))
  humandncu2 <- (Di_aa$aa_human2 %>% prod())^(1/nrow(Di_aa))
  wenzidncu2 <- (Di_aa$aa_wenzi2 %>% prod())^(1/nrow(Di_aa))
  c(humandncu1,wenzidncu1,humandncu2,wenzidncu2)
}
##save all data
save(alldata,allor,mydncu,file = "~/codonpaper/code and data/allcombined.fig4.Rdata")

