# CUB
This repository contains core codes for our paper on "Dissimilation of synonymous codon usage bias in virus-host coevolution due to translational selection", as well as the processed data to facilitate direct reproduction of main figures and some supplementary figures.  Please refer to our paper for the detailed explanation of the variables appeared in the codes.

1)Fig.1, Fig.S2, and Fig.S4

"Fig.1 and Fig.1-related.R" is the code for Fig.1, Fig.S2 and Fig.S4.  Using this code, we calculated some mathematical parameters, such as "relative tRNA shortage caused by virus infection", "relative decoding time", "typical decoding time" using the raw yeast ribo-seq density data.

The final Rdataset of Fig.1, Fig.S2, and Fig.S4 are combinded in "combined_allfig1.Rdata".

2)Fig.2 and Fig.S5

"Fig.2 and Fig.2-related.R" is the code for Fig.2 and Fig.S5.  Using this code, we calculated some mathematical parameters, such as "relative tRNA shortage caused by virus infection", "typical decoding time", "Sensitivity (sensitivity of codons after the infection of Flu virus)" using the raw human ribo-seq density data.

The final datasets for Fig.2 and Fig.S5A-B were combined as "finaldata_forfig2.correlations.Rdata". "aa.eight.humangeneforfig2.Rdata" is the finally data for Fig.5c.  

3)Fig.3 and Fig.S6-8

"Fig.3 and Fig.3-related.R" is the code for Fig.3 and Fig.S6-8.  Using this code, we analyzed the trans regulation of CUB per unit expression of exogenous genes, cis-regulatory effect of CUB, net trans-regulatory effect of CUB by total expression of exogenous genes, and the point of "only top 208 high expression of exogenous gene can make significantly translation load to the host of yeast".

All final data used to draw the plot was combined in "allcombinedfig3.Rdata".

4)Fig.4 and Fig.S9

"Fig.4 and Fig.4-related.R" is the code for Fig.4 and Fig.S9.  From this code, we calculated the Dpvn and Dpvs values of the 52 VNS-trio and all mutants of Dengue and Zika virus.  

All final data used to draw the plot was combined in "allcombinedfig4.Rdata".
