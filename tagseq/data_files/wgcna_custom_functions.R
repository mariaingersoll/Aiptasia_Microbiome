### Custom functions for WGCNA pipeline

## Function to identify and select significant modules from WGCNA 
sigMods <- function(matrix, alpha = 0.0500001) {
  sig_mat <- (matrix < alpha)
  ind <- which(sig_mat == 'TRUE', arr.ind = TRUE)
  moduleCols <- gsub("ME", "", rownames(sig_mat)[ind[,"row"]])
  whichTrait <- colnames(sig_mat)[ind[,"col"]]
  out <- data.frame(moduleCols = as.data.frame(moduleCols),
                    whichTrait = as.data.frame(whichTrait))
  return(out)
}


## Function to do all the module/trait plots and GO MWU file creations
ModWGCNA <- function(wgcna_data, traitMods, traits, geneModuleMembership, MEs, vsdWG, corrplotpath, heatplotpath, goMWUpath) {
  
  ModWGCNA_out <- list()
  
  ########## Gene Significance
  
  # Create blank dataframes to populate with below forloop
  GSPvalue_df <- data.frame(matrix(nrow = ncol(datt), ncol = 0))
  GTSvalue_df <- data.frame(matrix(nrow = ncol(datt), ncol = 0))
  
  # Forloop to calculate gene trait significance and gene significance p values for each significant module/trait pairing
  for (t in 1:nrow(traitMods)) { #
    #t = 1
    ## trait selection/specification
    Trait <- traitMods$whichTrait[t] # pulls trait name
    selTrait <- as.data.frame(traits[,Trait]) # selects the trait data that matches the trait name
    names(selTrait) <- Trait # renames the selected trait data with trait name
    
    ## calculate trait correlation p values
    geneTraitSignificance <- as.data.frame(cor(datt, selTrait, use = "p")) # calculates correlation between GE data and trait data
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)) # Calculates Student asymptotic p-value
    names(geneTraitSignificance) <- names(selTrait) # rename geneTraitSignificance column with trait name
    names(GSPvalue) = paste("p.GS.", names(selTrait), sep = "") # rename GSPvalue colmn with trait name
    
    ## add the new trait columns to the existing dataframes
    GTSvalue_df <- cbind(GTSvalue_df, geneTraitSignificance)
    GSPvalue_df <- cbind(GSPvalue_df, GSPvalue)
    
  }
  
  ModWGCNA_out[["GTSig"]] <- GTSvalue_df
  ModWGCNA_out[["GSPvalue"]] <- GSPvalue_df
  
  
  ########## Correlation Plots
  
  # Create blank list to populate with correlation plots below
  plots <- list() # create an empty list to populate with plots
  
  # Forloop to create correlation plots per trait/module pairing (also saves as png)
  for (r in 1:nrow(traitMods)) {
    
    module <- traitMods$moduleCols[r]
    selectTrait <- traitMods$whichTrait[r]
    
    column <- match(module, modNames)
    moduleGenes <- (moduleColors == module)
    combo <- paste(module, selectTrait, sep = "; ")
    
    x <- abs(geneModuleMembership[moduleGenes, column])
    y <- abs(GTSvalue_df[moduleGenes, selectTrait])
    
    # plot the correlations
    plot <- ggplot() +
      theme_classic() +
      geom_point(aes_string(x = x, y = y), shape = 1, colour = module, size = 2) +
      xlab(paste(module, "module membership")) +
      ylab(paste("GS for", traitMods$whichTrait[r])) +
      stat_cor(aes_string(x = abs(geneModuleMembership[moduleGenes, column]), y = abs(GTSvalue_df[moduleGenes, selectTrait])), method = "pearson")
    ggsave(paste(corrplotpath, combo, ".png", sep = ""), width = 5, height = 4.3)
    
    plots[[combo]] <- plot
    
  }
  
  ModWGCNA_out[["CorrPlots"]] <- plots  
  
  
  ########## Module Heatmaps
  
  # Create blank list to populate with correlation plots below
  plots2 <- list() # create an empty list to populate with plots
  
  # Forloop to cycle through trait/module pairings to produce expression heatmaps and eigengene bar plots (also saves png)
  for (r in 1:nrow(traitMods)) {
    
    module <- traitMods$moduleCols[r] # pulls module colour
    selectTrait <- traitMods$whichTrait[r] # identified trait name
    combo2 <- paste(module, selectTrait, sep = "; ") # create name of module/trait pair
    module2 <- paste("ME", module, sep="") # name of module with 'ME' added for pulling correct data
    
    MEs$sample <- rownames(MEs) # name column for sample ID (from rownames) 
    ME <- MEs[, c("sample", module2)] # subset MEs for selected module only
    names <- rownames(arrange(ME, ME[,2])) #rearrange ME in ascending order and extract names
    datExpr <- cbind(wgcna_data, ME[,2]) # add mean expression to expression dataframe
    ME <- arrange(ME, ME[,2]) # order expression dataframe in ascending order
    
    ## Heatmap of isogroup expression per sample
    heatmap <- t(scale(datExpr[, moduleColors == module])) %>% 
      as.data.frame() %>%
      rownames_to_column("isogroup") %>%
      pivot_longer(-c(isogroup), names_to = "samples", values_to = "expr") %>%
      mutate(samples = fct_relevel(samples, names)) %>%
      ggplot(aes(x = samples, y = isogroup, fill = expr)) + 
      ggtitle(paste(selectTrait, "; ", module, " module", sep = "")) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position = "left", legend.title = element_blank(), axis.title = element_blank(), legend.text=element_text(size = 6)) +
      scale_fill_gradient2(low = "#5e3c99", high = "#e66101", mid = "gray95") +
      geom_tile()
    
    ## Barplot of eigengene expression per sample
    barplot <- ggplot(data = ME) +
      theme_cowplot() +
      geom_hline(aes_string(yintercept = 0)) +
      theme(text = element_text(size = 10), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_bar(aes_string(x = factor(ME$sample, levels = c(names)), y = ME[,2]), stat = "identity", fill = module, colour = "black", size = 0.2) +
      labs(x = "", y = "eigengene expression")
    
    ## Combine and save heatmap and bar plot
    heatbar <- cowplot::plot_grid(heatmap, barplot, nrow = 2, rel_heights = c(1, 1))
    ggsave(paste(heatplotpath, combo2, ".png", sep = ""), width = 8, height = 6)
    
    plots2[[combo2]] <- heatbar
    
  }
  
  ModWGCNA_out[["Heatmaps"]] <- plots2 
  
  
  ########## GO and KOG Files
  
  
  
  MEs <- MEs[,-length(MEs)] # removed the last column of the dataframe that should correspond with the sample IDs
  
  # calculating module memberships for all genes for all modules
  allkME <- as.data.frame(signedKME(datt, MEs))
  names(allkME) <- gsub("kME", "", names(allkME))
  
  # Forloop to cycle through trait/module pairings to produce GO ouput CSVs
  for (r in 1:nrow(traitMods)) {
    r = 1
    module <- traitMods$moduleCols[r] # pulls module colour
    selectTrait <- traitMods$whichTrait[r] # identified trait name
    
    # Saving data for Fisher-MWU combo test (GO_MWU)
    inModuleBinary <- as.numeric(moduleColors == module)
    combo3 <- data.frame("gene" = row.names(vsdWG), "Fish_kME" = allkME[,module] * inModuleBinary)
    write.csv(combo3, file = paste(goMWUpath, selectTrait, "_",module, ".csv", sep = ""), row.names = FALSE, quote = FALSE)
    
    # saving it also as a binary version
    combo4 <- combo3 %>% 
      mutate(binary = ifelse(Fish_kME > 0, 1, 0)) %>% 
      dplyr::select(-Fish_kME)
    write.csv(combo4, file = paste(goMWUpath, module, "_binary", ".csv", sep = ""), row.names = FALSE, quote = FALSE)
  }
  
  return(ModWGCNA_out)
  
}

