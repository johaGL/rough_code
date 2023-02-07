# playing fake gsea:  
# for myself:
# do a code to do a "fake gsea" nice visualization:
#     (fgsea in R) --> no because gsea was run with DEGs. its enough
#    - use a gmt file with the metabolic pathways (Thomas, Ben)
#     - use entire DESeq table (even with no significant, just to make the dots)
#     - make a list of dataframes, one dataframe by pathway
#     make, by each metabolic pathway a plot:
#           - build a df : 'gene' 'rank' 'log2FC' '
#           - rank the genes
#            - a dotplot with x the ranked , and y the log2FC 
#             - only label those having padj <= 0.05
#             - jitter ? how to deal with overlapping dots? think about it. 
#             -or do lollypop with invisible vertical bars ? 


fake_gmt = data.frame('pathway' =c('CELL MIGRATION', 'LIPID METABOLISM'),
  'leadingEdge'= c('AGT RARB APOE M N O P Q', 'SOX PAX RAS FAS MAD M P')
)  # leadingEdge?? confirm!!! or maybe the dataframe is absolutely otherwise

fake_DE_result = data.frame(
  'symbol'= c('AGT', 'RARB', 'APOE', 'M','N','O','P','Q', 'SOX', '', ''), # to complete
  'log2FoldChange'= c(1.9, 1.1, -1.1, 0.5, -0.6, 0, 0.04, ) ,  # to complete
  'padj'= c(0.004, 0.003)  # to complete
)