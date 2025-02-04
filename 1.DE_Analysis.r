#Run Rv4.1.0
########## -- Libraries --##########
options(mc.cores = 8)
library(SetRank)
library(Polychrome)
library(GSEABase)
library(pheatmap)
library(UpSetR)
library(GeneSets.Homo.sapiens)
library(limma)
library(edgeR)
library(stats)
library(factoextra)
library(umap)
library(Rtsne)
library(ggplot2)
library(ExpressionNormalizationWorkflow)
library(stringr)
library(pvca)
library(mclust)
library(fossil)
library(amap)
library(dendextend)
library(data.table)
library(rrcov)
library("Hmisc")
library(ggrepel)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# TITLE: Rachel Rutishauser HIV study
# AUTHOR: Leanne Whitmore
########## -- Files & Objects --##########
ENSEMBL2human <- read.csv(
    file = "./Esenmbl2Symbol.csv",
    header = TRUE,
    stringsAsFactors = FALSE
)

# When set to true code will run enrichment analysis (takes a bit of time so have added this option to easily turn off)
SetRankRun <- FALSE
if (isTRUE(SetRankRun)) {
    # converters
    symbol2EntrezID <- createIDConverter(
        "Homo.sapiens", "SYMBOL",
        "ENTREZID"
    )
    IDConverter <- createIDConverter(
        "Homo.sapiens", "ENTREZID",
        "SYMBOL"
    )
}

########## -- FUNCTIONS --##########
theme_Publication <- function(base_size = 14) {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

vizualize_DE_genes_bp <- function(results, plot_file) {
    print("STATUS: Generating bar plot of number of DE genes...")
    results_t <- t(summary(results))
    results_t <- results_t[, -2]

    for (i in 1:(length(row.names(results_t)))) {
        results_t[i, 1] <- results_t[i, 1] * -1
    }

    DE <- as.data.frame(results_t)
    DE <- setnames(DE,
        old = c("Var1", "Var2", "Freq"),
        new = c("Time_Point", "grodwn", "DE_genes")
    )

    # Create plot
    ggplot(DE, aes(
        x = Time_Point, y = DE_genes, fill = grodwn,
        label = DE$DE_genes
    )) +
        geom_bar(stat = "identity", position = "identity") +
        # geom_text(size = 5, position = position_stack(vjust = 0) )+
        # theme_light() +
        theme_minimal() +
        scale_fill_manual(values = c("#0808c4", "#da9618")) +
        # xlab("Time point")
        ylab("# of DE Genes") + labs(fill="") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 10)
        )
    ggsave(plot_file, width = 6, height = 4, units="in", bg="white", dpi = 300)
}

generate_boxplots_voom <- function(data, labels, filename, figres, maintitle, ylabtitle) {
    png(filename, width = 10, height = 8, units = "in", res = figres)
    # par(mar=c(1,1,1,1))
    minvalue <- min(data)
    maxvalue <- max(data)
    boxplot(data,
        labels = labels, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = ylabtitle, main = maintitle, cex.axis = .6, las = 2,
        frame = FALSE
    )
    dev.off()
}

generate_density_plot <- function(data, labels, filename, figres) {
    png(filename, res = figres)
    par(xpd = TRUE)
    if (length(labels) > 10) {
        plotDensities(data, legend = FALSE)
    } else {
        plotDensities(data,
            legend = "topright",
            inset = c(-0.2, 0), levels(labels)
        )
    }
    dev.off()
}

normalize_data <- function(CM2, targetfile) {

    # order target and count matrix so they are the same (THIS IS IMPORTANT)
    CM2 <- CM2[, rownames(targetfile)]

    # CHECK IF ORDER IS THE SAME
    if (all.equal(colnames(CM2), rownames(targetfile)) != TRUE) {
        print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
        print(rownames(targetfile))
        print(colnames(CM2))
    }
    # biological reps
    targetfile$bioreps <- targetfile$RR.Name
    bioreps <- factor(targetfile$bioreps)
    biorepsdesign <- model.matrix(~ 0 + bioreps)

    # normalize
    CM2 <- DGEList(counts = CM2)
    CM2 <- calcNormFactors(CM2, method = "TMM") # TMM normalization
    message("STATUS: normalizing")
    png(file.path(norm_results, "mean_variance_norm.png"))
    Pi.CPM <- voom(counts = CM2, normalize.method = "none",
        design = biorepsdesign, plot = T, span = 0.1, save.plot=T)
    dev.off()
    write.csv(Pi.CPM$E, file.path(norm_results, "1.norm_matrix.csv"))
    message("STATUS: getting the corfit")
    corfit <- duplicateCorrelation(CM2$counts,
        block = factor(targetfile$Scope_ID
      )
    ) # account for repeated sampling of individuals

    message("STATUS: renormalizing with corfit")
    png(file.path(norm_results, "mean_variance_norm_corfit.png"))
    Pi.CPM <- voom(counts = CM2, normalize.method = "none",
                design = biorepsdesign,
                correlation=corfit$consensus.correlation,
                plot = T, span = 0.1, save.plot=T)
    dev.off()

    message("STATUS: recalculating corfit")
    corfit <- duplicateCorrelation(Pi.CPM,
        block = factor(targetfile$Scope_ID)
    )
    write.csv(Pi.CPM$E, file.path(norm_results, "1.norm_matrix_corfit.csv"))

    sig_HGNC <- merge(ENSEMBL2human, Pi.CPM$E,
        by.x = "Gene.stable.ID",
        by.y = "row.names",
        all.X = T, all.Y = T
    )

    sig_HGNC <- sig_HGNC[, !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
    sig_HGNC <- avereps(sig_HGNC,
        ID = sig_HGNC$HGNC.symbol
    )
    rownames(sig_HGNC) <- sig_HGNC[, "HGNC.symbol"]
    sig_HGNC <- sig_HGNC[, !(colnames(sig_HGNC) %in% c("HGNC.symbol"))]
    sig_HGNC <- as.matrix(data.frame(sig_HGNC))
    write.csv(sig_HGNC, file.path(norm_results, "1.norm_matrix_HGNC.csv"), quote = FALSE)
    return(list("norm"=Pi.CPM, "corfit"=corfit, "TMM"=CM2))
}

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca=FALSE, legend="right",
                    morePCs=FALSE) {

    # Run PCA/SVD reduction
    if (isFALSE(pca)) {
        pca <- prcomp(t(exprs))
    }
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)


    vizualize_pca(
        file.path(results_path, paste0("svd_", base_file_name)),
        sv$u, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, E, size, legend
    )
    vizualize_pca(
        file.path(results_path, paste0("pca_", base_file_name)),
        pca$x, labels[, target_columns[1]],
        labels[, target_columns[2]],
        figres, E, size, legend
    )
    vizualize_scree_plot(
        file.path(
            results_path,
            paste0("scree_", base_file_name)
        ), pca, figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
        loadingscores["PC1"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
        loadingscores["PC2"], figres
    )
    if (isTRUE(morePCs)) {
        is_pc3_0 <- loadingscores$PC3 > 0
        is_pc4_0 <- loadingscores$PC4 > 0
        is_pc5_0 <- loadingscores$PC5 > 0
        is_pc5_0 <- loadingscores$PC5 > 0
        is_pc6_0 <- loadingscores$PC6 > 0
        is_pc7_0 <- loadingscores$PC6 > 0

        loadingscores <- as.data.frame(pca$rotation)
        loadingscores <- loadingscores[is_pc3_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC3)), ]
        save_loading_scores(
            file.path(results_path, paste0("loadingscores_pc3", base_file_name, ".txt")),
              loadingscores["PC3"], figres
        )

        loadingscores <- as.data.frame(pca$rotation)
        loadingscores <- loadingscores[is_pc4_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC4)), ]
        save_loading_scores(
            file.path(results_path, paste0("loadingscores_pc4", base_file_name, ".txt")),
              loadingscores["PC4"], figres
        )

        loadingscores <- as.data.frame(pca$rotation)
        loadingscores <- loadingscores[is_pc5_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC5)), ]
        save_loading_scores(
            file.path(results_path, paste0("loadingscores_pc5", base_file_name, ".txt")),
            loadingscores["PC5"], figres
        )

        loadingscores <- as.data.frame(pca$rotation)
        loadingscores <- loadingscores[is_pc6_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC6)), ]
        save_loading_scores(
            file.path(results_path, paste0("loadingscores_pc6", base_file_name, ".txt")),
              loadingscores["PC6"], figres
        )

        loadingscores <- as.data.frame(pca$rotation)
        loadingscores <- loadingscores[is_pc7_0, ]
        loadingscores <- loadingscores[with(loadingscores, order(-PC7)), ]
        save_loading_scores(
            file.path(results_path, paste0("loadingscores_pc7", base_file_name, ".txt")),
            loadingscores["PC7"], figres
        )

    }
    return(pca)
}

umap_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns,
                     figres = 100, size = 1, UMAP=FALSE, legend="right") {
    # Runs default paramaters of umap
    if (isFALSE(UMAP)) {
        UMAP <- umap(t(exprs))
    }
    vizualize_umap(
        file.path(results_path, paste0("umap_", base_file_name)),
        UMAP$layout, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size, legend
    )

    return(UMAP)
}

vizualize_umap <- function(plot_file, U, class1, class2, figres, size, legend) {
    # Vizualize umap reduction
    library(Polychrome)
    minx <- min(U[, 1])
    maxx <- max(U[, 1])
    miny <- min(U[, 2])
    maxy <- max(U[, 2])
    P36 <- createPalette(length(levels(factor(class2))), c("#6f0c0c", "#127012", "#090971"))

        qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
            theme_Publication() + theme(legend.title = element_blank()) +
            xlab("UMAP 1") +
            ylab("UMAP 2") +
            scale_color_manual(values = as.character(P36)) +
            scale_fill_manual(values = as.character(P36)) +
            xlim(minx, maxx) + ylim(miny, maxy) +
            theme(legend.position = legend, legend.key.size = unit(0.4, "cm")) +
            scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
    
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size, legend) {
    # Vizualize PCA  results
    library(Polychrome)
    minx <- min(PCA[, 1])
    maxx <- max(PCA[, 1])
    miny <- min(PCA[, 2])
    maxy <- max(PCA[, 2])

    P36 <- createPalette(length(levels(factor(class2))), c("#6f0c0c", "#127012", "#090971"))
    qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
        theme_Publication() +
        theme(legend.title = element_blank()) +
        xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
        ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
        theme(legend.position = legend) +
        scale_color_manual(values = as.character(P36)) +
        scale_fill_manual(values = as.character(P36)) +
        scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
       
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 300)
}

tsne_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns, figres = 300, size = 1, T = FALSE, legend="right") {
    # Runs default paramaters of umap
    if (isFALSE(T)) {
        T <- Rtsne(t(exprs), perplexity = 1)
    }
    vizualize_tSNE(
        file.path(results_path, paste0("tsne_", base_file_name)),
        T$Y, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size, legend
    )
    return (T)
}

vizualize_tSNE <- function(plot_file, T, class1, class2, figres, size, legend) {
    # Vizualize tsne reduction
    library(Polychrome)
    minx <- min(T[, 1])
    maxx <- max(T[, 1])
    miny <- min(T[, 2])
    maxy <- max(T[, 2])
    P36 <- createPalette(length(levels(factor(class2))), c("#6f0c0c", "#127012", "#090971"))
    qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
        theme_Publication() + theme(legend.title = element_blank()) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        scale_color_manual(values = as.character(P36)) +
        scale_fill_manual(values = as.character(P36)) +
        xlim(minx, maxx) + ylim(miny, maxy) +
        theme(legend.position = legend) +
        scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_scree_plot <- function(plot_file, PCA, figres) {
    # Vizualize principle component variation results
    scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
    png(plot_file, width = 7, height = 6, units = "in", res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores <- function(write_file, df, figres) {
    # Save list of genes that have a positive effect on variation of principle
    # component 1 and 2 sorted from most influential
    write.table(df, file = write_file)
}

filter_read_counts_mean <- function(cm, filter_cutoff) {
    # Filter value was calculated by:
    #### Filters by row means usually set at 10 reads per gene across all samples

    A <- rowMeans(cm)
    isexpr <- A >= filter_cutoff
    cmfl <- cm[isexpr, ]
    return(cmfl)
}

rename_samples <- function(samples) {
    newsampleIDs <- c()
    for (i in samples) {
        i <- str_remove(i, "_DARE\\S+$")
        newsampleIDs <- c(newsampleIDs, i)
    }
    return(newsampleIDs)
}

generate_design_matrix_2 <- function(normmatrix, target) {
    target$bioreps <- paste(target$Controllers, target$RR.Name,  sep = "_")
    com <- factor(target$bioreps)
    batch <- paste(target$Date.Sample.Rec, sep = "_")
    batch <- factor(batch)

    mm <- model.matrix(~ 0 + com + batch)

    rownames(mm) <- colnames(normmatrix)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
    mm <- mm[, colSums(mm) > 0]

    excludeAll <- nonEstimable(mm)
    if (length(excludeAll) > 0) {
        message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
    }

    if ("ti" %in% excludeAll) {
        return("interactions term non estimable")
    }
    mm <- mm[, !colnames(mm) %in% excludeAll]
    if (!is.fullrank(mm)) {
        return("not full rank")
    }
    return(mm)
}

generate_clusterdendograms <- function(hc, plotfilename1, adjvalue, labelsize = 0.7) {
    counter <- 0
    labelsf <- c()
    colors <- c()
    dend <- as.dendrogram(hc)
    dend_labels <- labels(dend)
    P36 <- createPalette(length(levels(factor(dend_labels))),
         c("#6f0c0c", "#127012", "#090971"))
    names(P36) <- unique(dend_labels)
    for (i in dend_labels) {
        labelsf <- c(labelsf, i)
        colors <- c(colors, P36[[i]])
    }
    labels_colors(dend) <- colors
    labels_cex(dend) <- labelsize
    png(plotfilename1,
        units = "in", # bg = "transparent",
        width = 14.5, height = 5, res = 300
    )
    par(mar = c(6, 3, 2, 0.5), xpd = TRUE)
    plot(dend, xlab = "", main = "")
    mtext(paste0("Adj Rand index ", round(adjvalue, 3)))
    dev.off()
}

gene_enrichment <- function(genes, results_folder, cluster, geneTranslation=TRUE) {
    if (isTRUE(geneTranslation)) {
        inputGenesTrans <- ENSEMBL2human[ENSEMBL2human$Gene.stable.ID %in% genes, ]
        inputGenesHGNC <- unique(unlist(inputGenesTrans$HGNC.symbol))
    } else {
       inputGenesHGNC <- genes
    }
    inputGenes <- symbol2EntrezID(inputGenesHGNC)
    network <- setRankAnalysis(inputGenes, collection,
        use.ranks = FALSE,
        setPCutoff = 0.01,
        fdrCutoff = 0.05
    )

    generate_folder(results_folder)
    #### IMPORTANT OUTPUT INFORMATION###
    # SetRank value -  value reflects the prominence of a gene set in a gene set network (based on PageRank algorithm developed by google)
    # pSetRank value - significance of the SetRank value
    exportSingleResult(network, inputGenes,
        collection, paste0(results_folder, "/de_unranked_", cluster),
        IDConverter = IDConverter
    )
    # png(file.path(results_folder, paste0("de_unranked_", cluster,".png")), res=100)
    # plot(network, layout = layout.spring)
    # dev.off()
    return(network)
}

run_gsea_and_ora <- function(finalrankedgenes, gmt.file, universe, region, results_folder, GSEA=TRUE) {
    library(fgsea)
    pathways <- gmtPathways(gmt.file)
    if (isTRUE(GSEA)) {
        fgseaRes <- fgsea(
            pathways = pathways,
            stats = finalrankedgenes,
            minSize = 15,
            maxSize = 1000
        )
        fgseaRes <- fgseaRes[order(pval), ]
        fgseaResSig <- fgseaRes[fgseaRes$padj <= 0.05, ]
        fgseaResdf <- as.data.frame(fgseaRes)
        fwrite(as.data.frame(fgseaResdf), file.path(results_folder, paste0(region, "_GSEA_allResults.csv")))
        fwrite(as.data.frame(fgseaResSig), file.path(results_folder, paste0(region, "_GSEA_allResultsSig.csv")))
    }
    foraRes <- fora(pathways, finalrankedgenes, universe, minSize = 5, maxSize = Inf)
    foraRes <- foraRes[order(pval), ]
    foraResSig <- foraRes[foraRes$padj <= 0.05, ]
    fwrite(as.data.frame(foraRes), file.path(results_folder, paste0(region, "_ORA_allResults.csv")))
    fwrite(as.data.frame(foraResSig), file.path(results_folder,paste0(region, "_ORA_allResultsSig.csv")))
    # return(list("ORA" = foraRes, "sigORA" = foraResSig, "sigGSEA" =fgseaResSig, "GSEA" = fgseaRes))
}

volcano_de_nompval <- function(results, filename="", lfccutoff=0.25,
     pointsize=2.3, maxoverlaps=35) {
    p <- ggplot(
        results,
        aes(
            x = LFC, y = -log10(pval),
            color = Sign, label = genes
        )
    ) +
        geom_vline(xintercept = c(lfccutoff, -1*lfccutoff), lty = "dashed") +
        geom_hline(yintercept = -log10(0.05), lty = "dashed") +
        geom_point(size = pointsize, alpha=0.5) +
        labs(
            x = "log2(FC)",
            y = "Significance, -log10(P)",
            color = "Sig"
        ) +
        scale_color_manual(
            values = c(
                "Pval < 0.01" = "orange",
                "Pval < 0.05" = "dodgerblue",
                "NS" = "gray"
            ),
            guide = guide_legend(override.aes = list(size = 2))
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        geom_text_repel(
            data =subset(results, pval  < 0.05), aes(labels=genes),
           #data = results, #subset(results, pval.adj < 0.05),
            size = 2,  point.padding = 0.15, color = "black",
            min.segment.length = .1, box.padding = .2, lwd = 1,
            # max.overlaps = 25
            #point.padding = 0.15, color = "black",
            #min.segment.length = .1, box.padding = .1, lwd = 1,
            max.overlaps = maxoverlaps
        ) +
        theme_classic(base_size = 12) +
        theme(legend.position = "bottom",legend.key.size = unit(0.4, "cm")) +
        theme(strip.background = element_rect(fill = "white")) # c("DWM"='red', "GM"='darkblue',"SWM"="purple")
    ggsave(filename, width = 3.5, height = 3, units = "in", dpi = 300, bg = "white")
    ggsave(str_replace(filename,".png", ".pdf"), width = 3.5, height = 3, units = "in", dpi = 300, bg = "white")
}

convertensembl2hgnc <- function(df) {
    df_HGNC <- merge(ENSEMBL2human, df,
        by.x = "Gene.stable.ID",
        by.y = "row.names",
        all.X = T, all.Y = T
    )

    df_HGNC <- df_HGNC[, !(names(df_HGNC) %in% c("Gene.stable.ID"))]
    df_HGNC <- avereps(df_HGNC,
        ID = df_HGNC$HGNC.symbol
    )
    rownames(df_HGNC) <- df_HGNC[, "HGNC.symbol"]
    df_HGNC <- df_HGNC[, !(colnames(df_HGNC) %in% c("HGNC.symbol"))]
    df_HGNC <- as.matrix(data.frame(df_HGNC, check.names = FALSE))
    class(df_HGNC) <- "numeric"
    return(df_HGNC)
}

########## -- Load Data --##########

# --Read in target files
message("STATUS: Load tables")
cm <- read.table("count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("Target-file.csv", row.names = 1,sep = ",", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
PINC <- c("1292" ,"1408", "8106")
PIC <- c("2251", "2594", "2887", "2908", "8012", "8092")
x <- c()
for (i in target$Scope_ID) {
    if (i %in% PINC) {
        x <- c(x, "PINC")
    } else if (i %in% PIC) { 
        x <- c(x, "PIC")
    } else { 
        x <- c(x, "NA")

    }
}
target$Controllers <- x
target$RR.Name <- str_remove_all(target$RR.Name, " ")
target$RR.Name <- str_replace_all(target$RR.Name, "R3", "R2")
x <- c()
for (i in target$VL) {
    if (i == "<30") {
        x <- c(x, 15)
    } else if (i == "<40") { 
        x <- c(x, 15)
    } else if (i == "0") { 
        x <- c(x, 15)
    } else if (is.na(i)) { 
        x <- c(x, 0)
    } else { 
        x <- c(x, i)

    }
}
target$VLadj <- as.numeric(x)

# --Rename sample names in count matrix
newsampleIDs <- rename_samples(colnames(cm))
colnames(cm) <- newsampleIDs
newsampleIDstarg <- rename_samples(rownames(target))
rownames(target) <- newsampleIDstarg
cm <- cm[,colnames(cm) %in% rownames(target)]
# remove the samples because of contamination or unclear patient assignmnt (G092)
removetest <- c("G092")
cmtest <- cm[, !(colnames(cm) %in% removetest)]
targettest <- target[!(rownames(target) %in% removetest), ]
if (all.equal(colnames(cmtest), rownames(targettest)) != TRUE) {
    print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
    targettest <- targettest[colnames(cmtest), ]
    if (all.equal(colnames(cmtest), rownames(targettest)) == TRUE) {
        print("ISSUE HAS BEEN FIXED")
    } else {
        print("ISSUE IS NOT FIXED, PLEASE LOOK AT MANUALLY")
    }
}

remove1 <- c("G100", "G130","G018", "G092", "G133")
cm <- cm[, !(colnames(cm) %in% remove1)]
target <- target[!(rownames(target) %in% remove1), ]
if (all.equal(colnames(cm), rownames(target)) != TRUE) {
    print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
    target <- target[colnames(cm), ]
    if (all.equal(colnames(cm), rownames(target)) == TRUE) {
        print("ISSUE HAS BEEN FIXED")
    } else {
        print("ISSUE IS NOT FIXED, PLEASE LOOK AT MANUALLY")
    }
}

# --Write new renamed files to file
count_results <- "1.count_data/"
generate_folder(count_results)
write.csv(cm, file.path(count_results, "count_matrix_renamed.csv"))
write.csv(target, file.path(count_results, "target_renamed.csv"))
# --generate figure of all counts
generate_density_plot(
    cm, rownames(target), file.path(count_results, "de_intensities_raw_counts.png"), 100
)

########## -- Normalize Data --##########

# for (cutoff in c(0, 1, 2,3,4,6,7,8,9,10, 12,14,15,16,18,20, 25, 30, 35, 40, 45, 50)) {
#     cmfl_counts <- filter_read_counts_mean(cm, cutoff)
#     write.csv(cmfl_counts, file.path(count_results, paste0("count_matrix_renamed_fl_", cutoff, ".csv")))

#     norm_results <- paste0("1.norm_data_fl_",cutoff, "/")
#     generate_folder(norm_results)
#     nd2 <- normalize_data(cmfl_counts, target)
#     Pi.CPM = nd2$norm
#     corfit = nd2$corfit
#     generate_boxplots_voom(Pi.CPM$E, target$Scope_ID,
#         file.path(norm_results, "boxplot_vnorm_all.png"),
#         100,
#         maintitle = "Normalized count matrix",
#         ylabtitle = "voom normalized expression"
#     )
#     generate_density_plot(
#         Pi.CPM$E, rownames(target), file.path(norm_results, "de_intensities_norm_counts.png"), 100
#     )
#     saveRDS(Pi.CPM, file.path(norm_results, paste0("normobject_", cutoff, ".rds")))
#     saveRDS(corfit, file.path(norm_results, paste0("corfitobject_", cutoff, ".rds")))
# }

# I chose a cutoff value of 50 because density and variance plot 
# Doing normalization including 3 repeat samples to see if there was difference due to batch
for (cutoff in c(50)) {
    cmfl_counts <- filter_read_counts_mean(cmtest, cutoff)
    write.csv(cmfl_counts, file.path(count_results, paste0("count_matrix_renamed_fl_test_", cutoff, ".csv")))

    norm_results <- paste0("1.norm_data_fl_test_",cutoff, "/")
    generate_folder(norm_results)
    nd2 <- normalize_data(cmfl_counts, targettest)
    Pi.CPMtest = nd2$norm
    corfittest = nd2$corfit
    generate_boxplots_voom(Pi.CPMtest$E, targettest$Scope_ID,
        file.path(norm_results, "boxplot_vnorm_all.png"),
        100,
        maintitle = "Normalized count matrix",
        ylabtitle = "voom normalized expression"
    )
    generate_density_plot(
        Pi.CPMtest$E, rownames(targettest), file.path(norm_results, "de_intensities_norm_counts.png"), 100
    )
    saveRDS(Pi.CPMtest, file.path(norm_results, paste0("normobject_", cutoff, ".rds")))
    saveRDS(corfittest, file.path(norm_results, paste0("corfitobject_", cutoff, ".rds")))
}

# I chose a cutoff value of 50 because density and variance plot 
norm_results <- "1.norm_data_fl_50/"
Pi.CPM <- readRDS(file.path(norm_results, "normobject_50.rds"))
corfit <- readRDS(file.path(norm_results, "corfitobject_50.rds"))

cibersort <- Pi.CPM$E
cibersort <- merge(ENSEMBL2human, cibersort,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)
cibersort <- cibersort[, !(names(cibersort) %in% c("Gene.stable.ID"))]
cibersort <- avereps(cibersort,
  ID = cibersort$HGNC.symbol
)
rownames(cibersort) <- cibersort[, "HGNC.symbol"]
cibersort <- cibersort[, !(colnames(cibersort) %in% c("HGNC.symbol"))]
cibersort <- as.matrix(data.frame(cibersort))
class(cibersort) <- "numeric"
write.table(cibersort, file.path(norm_results, "1.norm_HGNC_cibersort.txt"), sep="\t", quote=F)

########## -- Feature Reduction --########## 

message("STATUS: Run Feature reduction")
feature_results <- "1.feature_red"
generate_folder(feature_results)

pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_norm.png",
    c("RR.Name", "Scope_ID"), 300, 4
)
pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_normControllers.png",
    c("Scope_ID", "Controllers"), 300, 4, pca=pca
)
pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_normControllers2.png",
    c("Controllers", "Controllers"), 300, 4, pca=pca
)
pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_normControllersdata.png",
    c("Date.Sample.Rec", "Controllers"), 300, 4, pca=pca
)
pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_normControllersdate2.png",
    c("Date.Sample.Rec", "Date.Sample.Rec"), 300, 4, pca=pca
)
pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_normControllersdate2RRName.png",
    c("Date.Sample.Rec", "RR.Name"), 300, 4, pca=pca
)

umap <- umap_fun(
    Pi.CPM$E, target, 
    feature_results, "_norm.png",
    c("Scope_ID", "RR.Name"), 300, 4)
umap <- umap_fun(
    Pi.CPM$E, target, 
    feature_results, "_normControllers.png",
    c("Scope_ID", "Controllers"), 300, 4, UMAP=umap)
umap <- umap_fun(
    Pi.CPM$E, target, 
    feature_results, "_normControllersandDate.png",
    c("Date.Sample.Rec", "Controllers"), 300, 4, UMAP=umap)

########## -- Hierarchal clustering --##########

message("STATUS: Running hierarchal clustering...")
hi_cluster_results <- "1.hi_cluster_results"
generate_folder(hi_cluster_results)
d <- Pi.CPM$E
colnames(d) <- target$Scope_ID
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(target$Scope_ID)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(target$Scope_ID))))
generate_clusterdendograms(hc,
    paste0(hi_cluster_results, "/dendogram_Scope_ID.png"), m, labelsize = 1)

########## -- de controllers individual time points --##########
message("STATUS: Running de controllers individual time points")

deresults_path <- "1.de_controllers_timeComparison"
generate_folder(deresults_path)
Pi.lmfit <- lmFit(Pi.CPM, design = mm_all, block=target$Scope_ID,
             correlation=corfit$consensus)
saveRDS(Pi.lmfit, file.path("1.de_controllers", "lmfitobject.RDS"))
# Pi.lmfit <- readRDS( file.path("1.de_controllers", "lmfitobject.RDS"))

mm_all <- generate_design_matrix_2(Pi.CPM$E, target)

contrastsmatrix <- c(
    "comPINC_BL -comPIC_BL",
    "comPINC_preATI -comPIC_preATI",
    "comPINC_preR -comPIC_preR",
    "comPINC_postR1 -comPIC_postR1",
    "comPINC_postR2 -comPIC_postR2"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm_all)
contrast <- contrasts.fit(Pi.lmfit, contrasts = contr)
# Pi.contrasts <- eBayes(contrast, robust = TRUE, trend = TRUE)
Pi.contrasts <- treat(contrast, lfc = 0.58, robust = TRUE, trend = TRUE)
png(file.path(deresults_path, "SAplot.png"))
plotSA(Pi.contrasts, main = "Residual standard deviation vs. average log expression")
dev.off()

results <- decideTests(Pi.contrasts, lfc = 0.58, method = "separate", adjust.method = "none", p.value = 0.05)
write.csv(Pi.contrasts$coefficients, file = file.path(deresults_path, "coefficients.csv"), quote = F)
write.csv(Pi.contrasts$t, file = file.path(deresults_path, "t_stats.csv"), quote = F)
write.csv(Pi.contrasts$p.value, file = file.path(deresults_path, "p_value.csv"), quote = F)
Pi.contrasts$p.value.adj <- Pi.contrasts$p.value
for (i in 1:ncol(Pi.contrasts$p.value.adj)) Pi.contrasts$p.value.adj[, i] <- p.adjust(Pi.contrasts$p.value.adj[, i], method = "BH")
write.csv(Pi.contrasts$p.value.adj, file = file.path(deresults_path, "p_value_adj.csv"), quote = F)

dataMatrixde <- Pi.contrasts$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

Pi.contrasts$genes <- data.frame(ID_REF = rownames(Pi.contrasts))

write.csv(ExpressMatrixde, file = file.path(deresults_path, "expression_matrix_de.csv"), quote = F)
write.csv(results, file = file.path(deresults_path, "results_de.csv"), quote = F)
write.csv(dataMatrixde, file = file.path(deresults_path, "full_expression_matrix_de.csv"), quote = F)
message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))
ExpressMatrixde_HGNC <- convertensembl2hgnc(ExpressMatrixde)
write.csv(ExpressMatrixde_HGNC, file.path(deresults_path, "expression_matrix_lm_hgnc.csv"), quote = F)
dataMatrixde_HGNC <- convertensembl2hgnc(dataMatrixde)
write.csv(dataMatrixde_HGNC, file.path(deresults_path, "full_expression_matrix_de_lm_hgnc.csv"), quote = F)
pvalhgnc <- convertensembl2hgnc(Pi.contrasts$p.value)

vizualize_DE_genes_bp(results, file.path(deresults_path, "barplot.pdf"))
vizualize_DE_genes_bp(results, file.path(deresults_path, "barplot.svg"))
vizualize_DE_genes_bp(results, file.path(deresults_path, "barplot.png"))

### COMPLEX HEATMAP OF MULTIPLE HEATMAPS
row_ha <- rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
col_fun <- colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), c("darkblue", "mediumblue", "dodgerblue", "white", "orange", "red", "darkred"))
# Stage <- factor(rep(c("PIC", "PINC"), c(4, 4)), 
#     levels = c(c("PIC", "PINC")))

hmC <- Heatmap(ExpressMatrixde,
    name = "LFC",
    col = col_fun,
    show_row_names = FALSE,
    # row_split = PCs,
    # row_km = length(unique(c(clustergenesdf))),
    cluster_row_slices = T,
    cluster_columns = FALSE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    # column_split = Stage,
    cluster_column_slices = FALSE,
    # column_labels = collabels,
    column_names_rot = 90,
    border = F, width = unit(5, "cm"),
    column_names_gp = grid::gpar(fontsize = 8),
    column_title_gp = grid::gpar(fontsize = 8.5, fontface ="bold"),
    # row_names_gp = grid::gpar(fontsize = 20),
    heatmap_legend_param = list(
        legend_height = unit(3, "cm"),
        labels_gp = gpar(fontsize = 14),
        title_gp = gpar(fontsize = 14, fontface = "bold")
    ), use_raster = TRUE,  raster_quality = 5
)

dsdf <- as.dist(1 - cor(t(ExpressMatrixde), method = "pearson"))
# dsdf <- dist(ExpressMatrixde, method = "euclidean")
hcdf <- hclust(dsdf, method = "ward.D2")
xdf <- cutree(hcdf, h = 12)
geneorderdf <- hcdf$labels[hcdf$order]
clustergenesdf <- xdf[geneorderdf]

P36 <- createPalette(length(unique(c(clustergenesdf))), c("#ff0000", "#00ff00", "#0000ff"))
names(P36) <- unique(c(clustergenesdf))

clustergenesdf <- clustergenesdf[hmC@row_names_param$labels]
ht3 <- Heatmap(clustergenesdf,
    name = "clusters",
    col = P36, width = unit(0.3, "cm"), show_row_names = FALSE
)

ht_list <- hmC + ht3
png(file.path(deresults_path, "ComplexHeatmap.png"), width = 6, height = 6, units = "in", res = 700)
# pdf(file.path(deresults_path, "ComplexHeatmap.pdf")) # , width = 6, height = 6, units = "in", res = 700)

heatmapCOV <- draw(ht_list,
    heatmap_legend_side = "left"
)
dev.off()
ExpressMatrixdecluster <- ExpressMatrixde[row_order(heatmapCOV), ]
clustergenesdfcluster <- clustergenesdf[row_order(heatmapCOV)]

dotplot4bulkheatmap(clustergenesdfcluster, ExpressMatrixdecluster,
    Pi.contrasts$p.value,
    breaks = NULL, width = 3.5, labelsize = 8, height=6, figurename=file.path(deresults_path, "dotplot4pearsonComplexHM.png"),
    scalesize=8, sizebreaks=6
)

line_plots_for_clusters(clustergenesdfcluster, Pi.CPM$E, orderclusters = unique(clustergenesdfcluster))

##############--volcano plot --##############
target$RR.Name <- factor(target$RR.Name, levels=c("BL", "vacc", "preATI", "preR", "postR1", "postR2"))
for (i in 1:length(levels(target$RR.Name)[-2])) {
    time <- levels(target$RR.Name)[-2][i]
    
    total <- merge(dataMatrixde_HGNC[, i], pvalhgnc[, i], by = "row.names")
    colnames(total) <- c("genes", "LFC", "pval")
    total$Sign <- "NS"
    total$Sign[total$pval < 0.05] <- "Pval < 0.05"
    total$Sign[total$pval < 0.01] <- "Pval < 0.01"
    total$Sign <- factor(total$Sign,
        levels = c(
            "NS", "Pval < 0.05", "Pval < 0.01"
        )
    )
    volcano_de_nompval(total, file.path(deresults_path, paste0("volcano",time, ".png")))
    write.csv(total, file.path(deresults_path, paste0(time, "totalinfo.csv")))
    if (isTRUE(SetRankRun)) {
        print("STATUS: gene enrichments for up regulated")
        genesup <- total[(total$pval < 0.05 & total$LFC > 0.26), "genes"]
        network <- gene_enrichment(genesup, 
            file.path(deresults_path, paste0("SetRank_results_", time)), "up",  geneTranslation=FALSE)
        print("STATUS: gene enrichments for module down regulated")
        genesdwn <- total[(total$pval < 0.05 & total$LFC < 0.26), "genes"]
        network <- gene_enrichment(genesdwn, 
            file.path(deresults_path, paste0("SetRank_results_", time)), "dwn",  geneTranslation=FALSE)
   }
}
# figures for set rank pathways
dfdwn <- read.table("/share/lwhitmo/projects/DARE_01/2nd_analysis/1.de_controllers_timeComparison/SetRank_results_BL/de_unranked_dwn_pathways.txt", sep="\t", header=TRUE)
BSLNdwn = c("R-HSA-2172127", "GO:0045087", "GO:0044548", "META_PWY-4081", "hsa04510")
dfdwn4fig <- dfdwn[dfdwn$name %in% BSLNdwn, ]
dfdwn4fig$type <- rep("dwn", nrow(dfdwn4fig))
dfdwn4fig$logadjustedPValue <- log10(dfdwn4fig$adjustedPValue) 
dfup <- read.table("/share/lwhitmo/projects/DARE_01/2nd_analysis/1.de_controllers_timeComparison/SetRank_results_BL/de_unranked_up_pathways.txt", sep="\t", header=TRUE)
BSLNup =c("GO:0007205", "GO:0046339")
dfup4fig <- dfup[dfup$name %in% BSLNup, ]
dfup4fig$type <- rep("up", nrow(dfup4fig))
dfup4fig$logadjustedPValue <- -log10(dfup4fig$adjustedPValue) 

df4fig <- rbind(dfdwn4fig,dfup4fig)
df4fig$description <- factor(df4fig$description, levels=c(df4fig[df4fig$name %in% BSLNup, "description"], df4fig[df4fig$name %in% BSLNdwn,"description"])) 

ggplot(df4fig, aes(x=logadjustedPValue, y=description, color=type)) +
    # geom_bar(stat='identity') + 
    geom_segment( aes(yxend=name, xend=0)) +
    geom_point( size=4) +
    theme_Publication() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(deresults_path, "BSLNenrich.png"), width=8, height=4, units="in", dpi=300)


