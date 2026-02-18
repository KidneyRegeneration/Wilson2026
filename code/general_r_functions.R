ggplotColors <- function(g){
  d <- 360/g
  h <- cumsum(c(15, rep(d,g - 1)))
  hcl(h = h, c = 100, l = 65)
}

## add genes in set to results table
addGenes <- function(results, cp_index, geneSymbols, topTab=NULL, pval=NULL){
  
  results = as.data.frame(results)
  results$genes = NA_character_
  
  for(i in 1:nrow(results)){
    symbols = geneSymbols[unlist(cp_index[rownames(results)[i]], use.names=FALSE)]
    symbols[symbols == ""] <- NA
    results$genes[i] = paste(symbols, collapse=",")
    if(!is.null(topTab)){
      tmp = topTab[, grepl("symbol",colnames(topTab),ignore.case = TRUE)]
      
      results$up[i] = paste(tmp[topTab$adj.P.Val < pval & topTab$logFC > 0 &
                                  tmp %in% symbols],
                            collapse = ",")
      results$down[i] = paste(tmp[topTab$adj.P.Val < pval & topTab$logFC < 0 &
                                    tmp %in% symbols],
                              collapse = ",")
    }
  }
  
  results
}

`%!in%` = Negate(`%in%`)


