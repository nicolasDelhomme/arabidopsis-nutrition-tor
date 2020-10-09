library(pathview)
## Plot the genes and compounbd fold change in the pathway view
#gvals: log ratio in a given comparison, names are the EC code (enzyme code, you need to do tha mapping athId <-> EC)
# cpd.data- same for compounds
pathview(gene.data = gvals,
         cpd.data = cpd.data,
         species="ko",
         pathway.id = "00940",
         out.suffix="kegg00940-ligin-vs-nolignin")

# ---------------------------------------------------------------------------------------
#' #### KEGG
#' 
#' Save the KEGG enrichment to the files (for every clusters)
#' 
ltest <- ld.enrichments[2,] 
alpha=.1

dev.null <- sapply(1:5, function(k){
    ltestnew <- ltest[[k]]
    write.table(ltestnew,file=paste0("/mnt/picea/projects/spruce/meriksson/28_Circadian_Time_Series/analysis/KEGG_",k,".txt"),
                sep="\t",quote = FALSE, row.names=FALSE)
})

#' Find the pathway of the significantly enriched enzymes
pathways <- sapply(1:5, function(k){
    keggLink("pathway", ltest[[k]]$id[ltest[[k]]$padj<=alpha])
})

#' Pathways are duplicated, clean up
pathways <- sapply(1:5, function(k){
    pathways[[k]][grep("map",pathways[[k]])]
})

#' Proportion of the pathways
dev.null <- sapply(1:5, function(k){
    barplot(sort(table(sub("path:","",pathways[[k]])),decreasing = TRUE),las=2)
})

#' Get the pathway info
pinfo <- lapply(split(unique(unlist(pathways)), 
                      ceiling(seq_along(unique(unlist(pathways)))/10)),keggGet)

pathway.df <- as.data.frame(do.call(rbind,lapply(pinfo,function(x){data.frame(ENTRY=sapply(x,"[[","ENTRY"),
                                                                              NAME=sapply(x,"[[","NAME"))})))

tab <- sapply(1:5, function(k){
    table(sub("path:","",pathways[[k]]))
})

occurrence.df <- data.frame(sapply(1:5, function(k){
    tab[[k]][match(pathway.df$ENTRY,names(tab[[k]]))]
}))
occurrence.df[is.na(occurrence.df)] <- 0
colnames(occurrence.df) <- paste0("Cluster",5:1)

KEGG.results <- cbind(pathway.df,occurrence.df)

pander(KEGG.results)

write.table(KEGG.results,file="../KEGG_results.tsv",
            sep="\t",quote = FALSE, row.names=FALSE)

#' Create a wordcloud of the pathway names
dev.null <- sapply(1:5,function(i){
    occ <- KEGG.results[,paste0("Cluster",i)]
    wordcloud(KEGG.results$NAME,occ/sum(occ),
              colors = hpal,scale = c(2,.3),rot.per = 0)  
})

# ---------------------------------------------------------------------------------------
#' #### PFAM

ltest <- ld.enrichments[3,] 
alpha=.01

dev.null <- sapply(1:5, function(k){
    ltestnew <- ltest[[k]]
    write.table(ltestnew,file=paste0("/mnt/picea/projects/spruce/meriksson/28_Circadian_Time_Series/analysis/PFAM_",k,".txt"),
                sep="\t",quote = FALSE, row.names=FALSE)
    sel <- ltestnew$padj <= alpha
    wordcloud(ltestnew[sel,"name",drop=TRUE],
              unlist(ltestnew[sel,"nt"] / sum(ltestnew[sel,"nt"])),
              colors = hpal,scale = c(2,.3),rot.per = 0)
})

#' #### GO
#' Save the GO enrichment to the files (for every clusters)
ltest <- ld.enrichments[1,]
dev.null <- sapply(5:1, function(k){
    ltestnew <- ltest[[k]]
    write.table(ltestnew,file=paste0("/mnt/picea/projects/spruce/meriksson/28_Circadian_Time_Series/analysis/ld_enrichments_",k,".txt"),row.names=FALSE)
    write.table(data.frame(ltestnew$id,ltestnew$padj),file=paste0("/mnt/picea/projects/spruce/meriksson/28_Circadian_Time_Series/analysis/ld_",k,".txt"),row.names=FALSE, quote = FALSE)
})

#' ### In what cluster are clock genes
#' 
goi.presence <- sapply(5:1, function(k){
    gl <- ld$Probes[ld.sel][ctld == k]
    gene_list %in% gl
})
colnames(goi.presence) <- paste0("Cluster",1:5)
rownames(goi.presence) <- gene_list
barplot(colSums(goi.presence))

pander(goi.presence[rowSums(goi.presence) >0,])

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
