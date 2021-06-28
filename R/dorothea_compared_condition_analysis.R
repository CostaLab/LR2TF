if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("dorothea")
# BiocManager::install("viper")
# BiocManager::install('limma')

library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(clustree)
library(viridis)
library(tidyverse)
library(viper)
library(Dict)

#---- Differential gene expression function for integrated data

get_markers = function(sc,
                       condition,
                       control,
                       cell.type,
                       only.sig = TRUE,
                       adj.pval.cutoff = 0.05,
                       min.pct = 0.25,
                       logfc.threshold = 0.25,
                       out.dir = 'diff_genes') {
  # Take care if no cells
  markers = data.frame()
  markers = FindMarkers(
    sc,
    ident.1 = paste(cell.type, condition),
    ident.2 = paste(cell.type, control),
    min.pct = min.pct,
    logfc.threshold = logfc.threshold
  )
  
  if (only.sig) {
    markers = markers[markers$p_val_adj < adj.pval.cutoff, ]
  }
  
  markers$gene = rownames(markers)
  cols = colnames(markers)
  cell.type.name = gsub('/', '_', cell.type)
  cell.type.name = gsub(' ', '_', cell.type.name)
  
  if (nrow(markers) > 0) {
    write.table(
      markers[, c('gene', cols[-6])],
      file = paste0(
        out.dir,
        '/integrated.diff.genes.',
        cell.type.name,
        '.',
        condition,
        '.vs.',
        control,
        '.txt'
      ),
      quote = FALSE,
      sep = '\t',
      row.names = FALSE
    )
  }
  return(markers)
}




#---- TF activity inference with DoRothEA (regulon has to be built prior to this)

run_dorothea = function(d.diff_genes,
                        case,
                        control,
                        out.dir,
                        dorothea_regulon_human,
                        regulon,
                        cell.types) {
  suppressPackageStartupMessages(library(viper))
  
  TF_activities_df = data.frame(
    tf = unique(dorothea_regulon_human$tf),
    row.names = unique(dorothea_regulon_human$tf)
  )
  
  significant_TF_activities = list()
  
  for (cell.type in cell.types) {
    diff.genes <- d.diff_genes[[cell.type]]
    # Estimate z-score values for the gene expression signature (GES)
    myStatistics = matrix(diff.genes$avg_log2FC,
                          dimnames = list(diff.genes$gene, 'avg_log2FC'))
    myPvalue = matrix(diff.genes$p_val, dimnames = list(diff.genes$gene, 'p_val'))
    mySignature = (qnorm(myPvalue / 2, lower.tail = FALSE) * sign(myStatistics))[, 1]
    mySignature = mySignature[order(mySignature, decreasing = TRUE)]
    
    # Estimate TF activity
    mrs = msviper(
      ges = mySignature,
      regulon = regulon,
      minsize = 4,
      ges.filter = FALSE,
      verbose = FALSE
    )
    
    TF_activities = data.frame(
      Regulon = names(mrs$es$nes),
      Size = mrs$es$size[names(mrs$es$nes)],
      NES = mrs$es$nes,
      p.value = mrs$es$p.value,
      FDR = p.adjust(mrs$es$p.value, method = 'fdr')
    )
    TF_activities = TF_activities[order(TF_activities$FDR), ]
    
    tf = rownames(filter(TF_activities, FDR < 0.05))
    
    if (length(tf) != 0) {
      df_end = as.data.frame(tf)
      
      if (length(significant_TF_activities) == 0) {
        df_end['first'] = cell.type
      } else {
        df_end[[cell.type]] = cell.type
      }
      
      significant_TF_activities[[cell.type]] <- df_end
    }
    
    # Write to file per cell type
    write.table(
      TF_activities,
      file = paste0(
        out.dir,
        case,
        '_vs_',
        control,
        '_tf_activity_',
        gsub(' ', '_', cell.type),
        '.txt'
      ),
      sep = '\t',
      row.names = FALSE,
      quote = FALSE
    )
    
    TF_activities_df[[cell.type]] = TF_activities[rownames(TF_activities_df), 'NES']
  }
  
  results = list(TF_activities_df, significant_TF_activities)
  
  return(results)
}




#---- Plot results from run_dorothea

plot_dorothea = function(df, case, control) {
  suppressPackageStartupMessages(library(pheatmap))
  
  
  
  # Order heatmap by TF activity
  df$tf = NULL
  df$diff = apply(
    df,
    1,
    FUN = function(x)
      diff(range(x))
  )
  df = df[order(df$diff), ]
  df.plot = df[complete.cases(df), ]
  tf_cluster = hclust(dist(t(df.plot)))
  tfs_ordered = tf_cluster$labels[tf_cluster$order]
  celltype_cluster = hclust(dist(df.plot))
  celltype_ordered = celltype_cluster$labels[celltype_cluster$order]
  
  paletteLength = 100
  myColor = colorRampPalette(c('Darkblue', 'white', 'red'))(paletteLength)
  viperBreaks = c(
    seq(min(df.plot), 0,
        length.out = ceiling(paletteLength / 2) + 1),
    seq(
      max(df.plot) / paletteLength,
      max(df.plot),
      length.out = floor(paletteLength / 2)
    )
  )
  
  p = pheatmap(
    t(df.plot)[tfs_ordered, celltype_ordered],
    fontsize_col = 8,
    fontsize_row = 8,
    color = myColor,
    breaks = viperBreaks,
    main = paste0(case, ' vs ', control),
    angle_col = 90,
    border_color = NA
  )
  
  return(p)
}

rds_path <- '~/dorothea_execution/20200325_step5.cca.anno.rds'
out_path <- '~/R_results/results_dorothea_combined_analysis_ABC'
confidence_level <- 'ABC'

dir.create(out_path)
dir.create(paste0(out_path, '/diff_genes'))

load(rds_path)
seuratobject <- UpdateSeuratObject(step5.cca)

new.cluster.ids <- c("Neural", "MSC", "Fibroblast", "Megakaryocyte", "Myeloid")
names(new.cluster.ids) <- levels(seuratobject)
seuratobject <- RenameIdents(seuratobject,  new.cluster.ids)
seuratobject[["new_annotation"]] <- Idents(object = seuratobject)
seuratobject <- StashIdent(object = seuratobject, save.name = "new_annotation")

seuratobject$celltype.condition = paste(seuratobject$new_annotation, seuratobject$protocol)
Idents(seuratobject) = 'celltype.condition'


markers = data.frame()
condition = "PMF,MF2"
control = "control"

#neural_cells <- subset(x = seuratobject, idents = cell)

#markers_cells = FindMarkers(seuratobject, ident.1 = paste(cell, condition),
#                      ident.2 = paste(cell, control), min.pct = 0.25, logfc.threshold = 0.25)

#markers = FindConservedMarkers(seuratobject, ident.1 = cell, grouping.var = "protocol", verbose = FALSE)

list.marker_dfs <- list()

for (cell in new.cluster.ids) {
  df.markers <-
    get_markers(
      seuratobject,
      condition,
      control,
      cell,
      only.sig = FALSE,
      adj.pval.cutoff = 0.05,
      min.pct = 0.0,
      logfc.threshold = 0.0,
      out.dir = paste0(out_path, '/diff_genes')
    )
  
  list.marker_dfs[[cell]] = df.markers
  
}

dorothea.path = 'https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv'
dorothea_regulon_human = read_csv(dorothea.path)

# Group regulons
regulon = dorothea_regulon_human %>%
  dplyr::filter(confidence %in% as.list(confidence_level)) %>%
  split(.$tf) %>%
  map(function(dat) {
    tf = dat %>% distinct(tf) %>% pull()
    targets = setNames(dat$mor, dat$target)
    likelihood = dat$likelihood
    list(tfmode = targets, likelihood = likelihood)
  })

dorothea_results = run_dorothea(
  list.marker_dfs,
  condition,
  control,
  paste0(out_path, '/diff_genes'),
  dorothea_regulon_human,
  regulon,
  new.cluster.ids
)

#dorothea_results$tf = NULL
#dorothea_results$diff = apply(dorothea_results, 1, FUN = function(x) diff(range(x)))
#dorothea_results = dorothea_results[order(dorothea_results$diff),]

svg(filename = paste0(out_path, '/heatmap_dorothea_combined.svg'))
plot = plot_dorothea(dorothea_results[[1]], condition, control)
print(plot)
dev.off()

if (length(dorothea_results[[2]]) > 1) {
  combined_significants = Reduce(function(dtf1, dtf2)
    merge(dtf1, dtf2, by = "tf", all = TRUE),
    dorothea_results[[2]])
  combined_significants = combined_significants %>%
    unite(sig_in,
          first:ncol(combined_significants),
          na.rm = TRUE,
          sep = " ")
  
  end_result = merge(dorothea_results[[1]],
                     combined_significants,
                     by = "tf",
                     all = TRUE)
  write.csv(end_result, paste0(out_path, '/all_tfs.csv'))
  end_result = end_result[complete.cases(end_result), ]
  write.csv(end_result, paste0(out_path, '/significant_tfs.csv'))
} else if (length(dorothea_results[[2]]) == 0) {
  write.csv(dorothea_results[[2]], paste0(out_path, '/all_tfs.csv'))
} else {
  colnames(dorothea_results[[2]][[1]])[colnames(dorothea_results[[2]][[1]]) != 'tf'] <-
    'sig_in'
  end_result = merge(dorothea_results[[1]],
                     dorothea_results[[2]][[1]],
                     by = "tf",
                     all = TRUE)
  write.csv(end_result, paste0(out_path, '/all_tfs.csv'))
  end_result = end_result[complete.cases(end_result), ]
  write.csv(end_result, paste0(out_path, '/significant_tfs.csv'))
}
