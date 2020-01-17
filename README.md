# CUB
Here are all the codes of the paper "Dissimilation of synonymous codon usage bias in virus-host coevolution due to translational selection".  Details of explaination method for these mathematical parameters, please see our Nature Ecology & evolution paper.  

1)fig1, figs2, and figs4

"fig1 and fig1related.R" is the code for fig1, figs2 and figs4.  Using this code, we calculated some mathematical parameters, such as "RTSV", "rrt (relative decoding time)", "TDT (typical decoding time)" using the raw yeast ribo-seq density data.

The final Rdataset of fig1, figs2, and figs4 are combinded in "combined_allfig1.Rdata".

2)fig2 and figs5

"fig2A-B.final.R" is the code for fig2 and figs5.  Using this code, we calculated some mathematical parameters, such as "RTSV", "TDT (typical decoding time)", "Sensitivity (sensitivity of codons after the infection of Flu virus)" using the raw human ribo-seq density data.

The final datasets for fig2 and figs5A-B were combined as "finaldata_forfig2.correlations.Rdata". "aa.eight.humangeneforfig2.Rdata" is the finally data for fig.5c.  

3)fig3 and figs6-8

"fig3.R" is the code for fig3 and figs6-8.  Using this code, we analyzed the trans regulation of CUB per unit expression of exogenous genes, cis-regulatory effect of CUB, net trans-regulatory effect of CUB by total expression of exogenous genes, and the point of "only top 208 high expression of exogenous gene can make significantly translation load to the host of yeast".

All final data used to draw the plot was combined in "allcombinedfig3.Rdata".

4)fig4 and figs9

"fig4.R" is the code for fig4 and figs9.  From this code, we calculated the Dpvn and Dpvs values of the 52 VNS-trio and all mutants of Dengue and Zika virus.  

All final data used to draw the plot was combined in "allcombinedfig4.Rdata".








