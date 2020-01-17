# CUB
Here are all the codes for the paper "Dissimilation of synonymous codon usage bias in virus-host coevolution due to translational selection".  Details of explaination method for these mathematical parameters, please see our Nature Ecology & evolution paper.  

1)fig1, figs2, and figs4

"fig1 and fig1related.R" is the code for fig1, figs2 and figs4.  Using this code, you can get how we procede the raw yeast ribo-seq density data and the method of how to calculated these mathematical parameters, such as "RTSV", "rrt (relative decoding time)", "TDT (typical decoding time)".

The final Rdatasets of fig1, figs2, and figs4 are combinded in "combined_allfig1.Rdata".


2)fig2 and figs5

"fig2A-B.final.R" is the code for fig2 and figs5.  Using this code, you can get how we procede the raw human ribo-seq density data and the method of how to calculated these mathematical parameters, such as "RTSV", "TDT (typical decoding time)", "Sensitivity (sensitivity of codons after the infection of Flu virus)".

The final datasets for fig2 and figs5A-B were saved as "finaldata_forfig2.correlations.Rdata".  Using this data, you can get the
translation seletion landscape of 61 codons caused by virus infection."aa.eight.humangeneforfig2.Rdata" is the finally data for fig.5c.  This was used to show that the pattern observed in fig2 and figs5 was not caused by the small viral genome. 

3)fig3 and figs6-8

"fig3.R" is the code for fig3 and figs6-8.  Using this code, you can deeply understand the trans regulation of CUB per unit expression of exogenous genes, cis-regulatory effect of CUB, net trans-regulatory effect of CUB by total expression of exogenous genes, and the point of "only high engough expression of exogenous gene can make significantly translation load".

All final data used to draw the plot was combined in "allcombinedfig3.Rdata".

4)fig4 and figs9








