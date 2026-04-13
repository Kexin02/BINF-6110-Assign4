# Temporal Dynamics of Cell-Type Composition and Epithelial Remodeling in the Respiratory Mucosa Revealed by Single-Cell RNA Sequencing

## 1. Introduction
The respiratory mucosa serves as a critical interface between the external environment and the host immune system, providing both a physical barrier and an active participant in immune defense (Şenel, 2021). Upon infection, epithelial cells not only maintain barrier integrity but also produce cytokines, chemokines, and interferon-stimulated genes that orchestrate innate and adaptive immune responses (Holtzman et al., 2002). These responses are highly dynamic and involve complex interactions between epithelial and immune cell populations, which are difficult to resolve using bulk transcriptomic approaches.

Single-cell RNA sequencing (scRNA-seq) has emerged as a powerful tool for dissecting cellular heterogeneity within complex tissues (Potter, 2018). Unlike bulk RNA sequencing, which averages gene expression across many cells, scRNA-seq enables the identification of rare cell types, transient cellular states, and dynamic transcriptional programs (Khozyainova et al., 2023). In the context of infection, scRNA-seq allows simultaneous analysis of immune activation, epithelial remodeling, and shifts in cell type composition over time (J. Li et al., 2025).

Several computational frameworks have been developed for analyzing scRNA-seq data, with Seurat being one of the most widely used tools for normalization, clustering, and visualization (Gribov et al., 2010). Graph-based clustering approaches, such as those implemented in Seurat, provide robust identification of transcriptionally distinct populations, though they require careful parameter tuning and biological interpretation. Alternative approaches, such as trajectory inference or probabilistic models, can capture continuous cellular transitions but may be more sensitive to noise and computational complexity (Jin et al., 2018; Y. Li et al., 2015).

In this study, we analyzed scRNA-seq data from mouse respiratory mucosa across multiple timepoints (D02, D05, D08, D14, and naive) to characterize infection-associated transcriptional changes. The objectives were to (1) identify and annotate cell populations, (2) examine infection-induced transcriptional programs, and (3) investigate temporal changes in cell composition and the emergence of specific cell states. Particular attention was given to a remodeling epithelial population (cluster 24), which exhibited dynamic abundance and distinct gene expression patterns.


## 2. Methods
Single-cell RNA sequencing data were analyzed using R (version 4.5.3) with the Seurat framework (Seurat v5.4.0) for preprocessing, clustering, and visualization (Hao et al., 2024). Additional packages included dplyr (v2.5.2) for data manipulation (Wickham et al., 2014), ggplot2 (v4.0.2) for visualization (Wickham, 2016), clusterProfiler (v4.16.0) for functional enrichment analysis (Yu, 2024), and org.Mm.eg.db (v3.21.0) for mouse gene annotation (Marc Carlson, 2025). Plot composition and auxiliary visualizations were generated using patchwork (v1.3.2) (Pedersen, 2019) and pheatmap (v1.0.13) (Kolde, 2010). All analyses were performed under a sequential computation plan to avoid parallelization-related memory issues.

Quality control was assessed using standard single-cell metrics, including the number of detected genes per cell (nFeature_RNA) and the total transcript counts per cell (nCount_RNA). These metrics were visualized across experimental timepoints (D02, D05, D08, D14, and naive) using violin plots and feature scatter plots to evaluate data consistency and potential outliers. No additional filtering thresholds were applied, as the distributions were comparable across conditions and did not indicate substantial technical bias.

For downstream analysis, cells from the respiratory mucosa (RM) were subsetted to focus on the tissue most relevant to infection-induced remodeling. Data normalization was performed using Seurat’s NormalizeData function with default parameters, followed by identification of highly variable genes using the variance-stabilizing transformation (VST) method, selecting the top 2000 features. Scaled expression values were obtained using ScaleData, and principal component analysis (PCA) was conducted to reduce dimensionality. Uniform Manifold Approximation and Projection (UMAP) was then applied using the first 20 principal components to visualize the global transcriptional structure of the dataset.

Graph-based clustering was performed using Seurat’s nearest-neighbor approach (FindNeighbors, dims = 1:20) followed by modularity optimization (FindClusters, resolution = 0.5), resulting in discrete transcriptional clusters. Cluster identities were manually annotated based on canonical marker genes, including epithelial (Epcam), T cell (Cd3e), myeloid (Lyz2), and olfactory neuron (Cnga4) markers. While this approach provides biologically interpretable labels, it relies on prior knowledge and may not fully capture intermediate or transitional cellular states.

Differential expression analysis was performed using the FindMarkers function to identify genes enriched in specific clusters relative to all other cells, using a log2 fold-change threshold of 0.25 and a minimum expression fraction of 25%. Functional enrichment analysis was conducted using the clusterProfiler package. Over-representation analysis (ORA) was applied to significantly upregulated genes using Gene Ontology biological process (GO-BP) terms, while gene set enrichment analysis (GSEA) was performed on a ranked list of genes based on log2 fold-change values. Both approaches rely on curated pathway annotations and may be biased toward well-characterized biological processes.

To assess temporal dynamics, the relative abundance of specific clusters and annotated cell types was calculated across timepoints by normalizing cell counts within each condition. These proportions were visualized to identify changes in cellular composition associated with infection and recovery. All visualizations were generated using ggplot2 with consistent formatting for clarity and reproducibility.


## 3. Results

### 3.1 Quality control of single-cell RNA-seq data

To ensure the reliability of downstream analyses, standard quality control metrics were first evaluated. The number of detected genes per cell (nFeature_RNA) exhibited comparable distributions across all timepoints, with similar medians and interquartile ranges (Figure 1). This indicates consistent library complexity and suggests that no subset of cells suffered from reduced gene detection.

Similarly, total transcript counts per cell (nCount_RNA) showed consistent distributions across conditions (Figure 2), supporting uniform sequencing depth and minimal technical variability. The absence of pronounced outliers further indicates that low-quality or multiplet cells were not a major confounding factor in the dataset.

The relationship between sequencing depth and gene detection was further assessed by examining the correlation between nCount_RNA and nFeature_RNA. A strong positive correlation (r = 0.827) was observed (Figure 3), confirming that increased sequencing depth resulted in improved gene detection efficiency. Collectively, these results demonstrate that the dataset is of high quality and suitable for subsequent analyses.

![QC_nFeature](figures/1.QC_nFeature.png)

**Figure 1. Distribution of detected genes per cell (nFeature_RNA).**  
Violin plots show the number of genes detected per cell across all timepoints. Comparable distributions indicate consistent library complexity and high-quality data across conditions.

![QC_nCount](figures/2.QC_nCount.png)

**Figure 2. Distribution of total RNA counts per cell (nCount_RNA).**  
Violin plots illustrate sequencing depth across samples. Similar distributions across timepoints suggest minimal technical variation.

![QC_scatter](figures/3.QC_scatter.png)

**Figure 3. Correlation between sequencing depth and gene detection.**  
Scatter plot of nCount_RNA versus nFeature_RNA shows a strong positive correlation (r = 0.827), indicating efficient gene capture and robust data quality.

### 3.2 Clustering reveals transcriptional heterogeneity

Dimensionality reduction using UMAP revealed a complex transcriptional landscape composed of multiple distinct clusters (Figure 4). These clusters represent transcriptionally heterogeneous cell populations within the respiratory mucosa.

The spatial separation between clusters suggests substantial biological diversity, while some degree of continuity between clusters indicates transitional or intermediate states. This organization reflects both discrete cell identities and continuous transcriptional variation.

![UMAP_clusters](figures/4.UMAP_clusters.png)

**Figure 4. UMAP visualization of unsupervised clustering.**  
Cells are grouped into transcriptionally distinct clusters based on gene expression profiles, revealing cellular heterogeneity within the dataset.

### 3.3 Cell type annotation based on canonical markers

To assign biological identities to clusters, canonical marker genes were examined. Violin plots demonstrated distinct expression patterns of key markers, including Epcam (epithelial), Cd3e (T cells), Lyz2 (myeloid), Krt13 (epithelial subset), and Cnga4 (olfactory neurons) (Figure 5). These markers showed strong specificity and minimal overlap across clusters.

Feature plots further confirmed that these markers localized to distinct regions of the UMAP embedding (Figure 7), reinforcing the spatial coherence of the identified populations. Based on these patterns, clusters were annotated into major cell types, including epithelial, IFN-responsive epithelial, myeloid, olfactory neuron, T cells, and remodeling epithelial populations (Figure 6).

![Marker_violin](figures/5.Marker_violin.png)

**Figure 5. Expression of canonical marker genes across clusters.**  
Violin plots show distinct expression patterns of key markers used for cell type annotation, supporting accurate classification of cell populations.

![UMAP_annotation](figures/6.UMAP_annotation.png)

**Figure 6. UMAP visualization of annotated cell types.**  
Clusters are assigned to major biological cell types based on canonical marker expression.

![Feature_annotation](figures/7.Feature_annotation.png)

**Figure 7. Spatial distribution of marker gene expression.**  
Feature plots illustrate localization of canonical markers within the UMAP, confirming cluster identities.

### 3.4 Infection-related transcriptional programs

To investigate infection-associated responses, the expression of interferon-stimulated and immune-related genes was examined. Genes such as Isg15, Ifit1, and Rsad2 exhibited elevated expression in specific epithelial and immune clusters (Figure 8), indicating activation of antiviral signaling pathways.

Additional genes, including Cxcl16 and Cd274, displayed more localized expression patterns, suggesting roles in immune modulation and cell–cell communication. These findings highlight functional heterogeneity within epithelial populations and demonstrate that subsets of epithelial cells adopt immune-responsive transcriptional states.

![Feature_story](figures/8.Feature_story.png)

**Figure 8. Expression of infection-related genes.**  
Feature plots show expression of interferon-stimulated and immune-related genes, indicating activation of antiviral and inflammatory pathways.

### 3.5 Temporal dynamics of cluster 24

Cluster 24 was identified as a distinct epithelial subpopulation exhibiting dynamic changes over time. Its relative abundance increased from 0.38% at D02 to 1.25% at D14, followed by a decrease in the naive condition (Figure 9).

Although cluster 24 represents a small fraction of the total cell population, its consistent expansion at later timepoints suggests a biologically meaningful role. This temporal pattern is consistent with involvement in tissue remodeling or recovery processes rather than early immune activation.

![Cluster24_abundance](figures/9.Cluster24_abundance.png)

**Figure 9. Temporal dynamics of cluster 24 abundance.**  
The proportion of cluster 24 cells increases over time and peaks at D14, suggesting involvement in late-stage biological processes.

### 3.6 Functional enrichment of cluster 24

To further characterize cluster 24, Gene Ontology enrichment analysis was performed. Over-representation analysis (ORA) revealed significant enrichment of biological processes related to epidermal development, cell–cell junction organization, and wound healing (Figure 10).

Complementary GSEA analysis identified enrichment of processes such as keratinization, RNA processing, translation, and ribonucleoprotein complex biogenesis (Figure 11). The convergence of these results suggests that cluster 24 is a metabolically active epithelial population undergoing differentiation and structural remodeling.

![ORA_dotplot](figures/10.ORA_dotplot.png)

**Figure 10. GO enrichment analysis (ORA) of cluster 24 marker genes.**  
Dot plot showing significantly enriched biological processes. Dot size represents gene count, and color indicates adjusted p-values.

![GSEA_dotplot](figures/11.Cluster24_GSEA_dotplot.png)

**Figure 11. GSEA of cluster 24 marker genes.**  
Gene set enrichment analysis reveals activation of pathways related to keratinization, RNA processing, and translation, indicating active epithelial remodeling.

### 3.7 Cell type composition changes over time

Analysis of overall cell type composition revealed dynamic but coordinated changes across timepoints (Figure 12). Epithelial cells remained the dominant population in all conditions, while IFN-responsive epithelial cells increased during intermediate stages (D05–D08) and declined thereafter.

Immune populations, including myeloid cells and T cells, also exhibited temporal variation, suggesting coordinated immune activation. These results indicate that infection induces both compositional and transcriptional changes within the tissue.

![Celltype_composition](figures/12.Celltype_composition.png)

**Figure 12. Cell type composition across timepoints.**  
Stacked bar plots show relative proportions of major cell types. Temporal shifts reflect coordinated immune and epithelial responses.

### 3.8 Marker gene expression defines cluster 24 identity

Finally, the expression of representative marker genes was examined across clusters. Cluster 24 showed strong and specific expression of genes such as Krt13, Krt6b, Plet1, Csta1, Prss27, Calml3, and Pglyrp4 (Figure 13).

These genes are associated with epithelial differentiation and barrier function, further supporting the classification of cluster 24 as a specialized remodeling epithelial population.

![Marker_dotplot](figures/13.Marker_dotplot.png)

**Figure 13. Expression of representative cluster 24 marker genes.**  
Dot plot shows expression level and proportion of expressing cells across clusters. Cluster 24 displays strong enrichment of epithelial remodeling markers.


## 4. Discussion
This study provides a detailed characterization of cellular heterogeneity and transcriptional dynamics in the respiratory mucosa during infection. The results highlight both stable cell identities and dynamic transcriptional responses across time.

A key finding is the identification of a remodeling epithelial population (cluster 24), which exhibits a clear temporal increase and peaks at D14. This population shows strong enrichment of genes associated with epithelial differentiation, including Krt13 and Krt6b, as well as proteases such as Prss27. These features are consistent with epithelial remodeling and barrier repair processes observed in inflamed tissues (Garcia-Hernandez et al., 2023; Liao et al., 2025). The enrichment of keratinization and epidermal development pathways further supports this interpretation, suggesting that cluster 24 represents a differentiation-associated epithelial state.

In addition to structural remodeling, cluster 24 also shows enrichment of RNA processing and translation-related pathways. This indicates high metabolic and transcriptional activity, which is often associated with rapidly differentiating or regenerating cells (Ayyaz et al., 2019). Similar transcriptional signatures have been observed in epithelial regeneration following injury, where cells undergo coordinated changes in gene expression to restore tissue integrity.

Interferon-stimulated genes, including Isg15 and Ifit1, were upregulated in subsets of epithelial and immune cells, indicating activation of antiviral responses. These findings are consistent with previous studies demonstrating that epithelial cells play an active role in innate immunity through interferon signaling (Aarreberg et al., 2019; Liao et al., 2025). The presence of IFN-responsive epithelial cells highlights the functional plasticity of epithelial populations during infection.

Changes in cell composition further support a coordinated response. While epithelial cells remained dominant, IFN-responsive epithelial and immune populations exhibited temporal variation. This is consistent with studies showing that infection induces both immune infiltration and epithelial state transitions (Lindeboom et al., 2024; Melms et al., 2021). The increase in cluster 24 at later timepoints suggests a shift from immune activation to tissue repair.

Despite these insights, several limitations should be considered. First, cell type annotation relied on canonical markers, which may not fully capture intermediate or transitional states. Second, clustering results depend on parameters such as resolution and dimensionality, which can influence interpretation. Third, enrichment analyses rely on existing annotations and may not identify novel biological processes.

Furthermore, this study analyzes discrete timepoints rather than continuous trajectories. Methods such as pseudotime analysis or RNA velocity could provide additional insight into lineage relationships and state transitions (Hou et al., 2023; La Manno et al., 2018). Future work integrating these approaches, along with experimental validation, would strengthen the conclusions.

Overall, this study demonstrates the power of single-cell transcriptomics in uncovering dynamic cellular responses to infection. The identification of a remodeling epithelial population provides new insight into tissue repair mechanisms and highlights potential targets for further investigation.


## 5. References
Aarreberg, L. D., Esser-Nobis, K., Driscoll, C., Shuvarikov, A., Roby, J. A., & Gale, M. (2019). Interleukin-1β Induces mtDNA Release to Activate Innate Immune Signaling via cGAS-STING. Molecular Cell, 74(4), 801-815.e6. https://doi.org/10.1016/j.molcel.2019.02.038
Ayyaz, A., Kumar, S., Sangiorgi, B., Ghoshal, B., Gosio, J., Ouladan, S., Fink, M., Barutcu, S., Trcka, D., Shen, J., Chan, K., Wrana, J. L., & Gregorieff, A. (2019). Single-cell transcriptomes of the regenerating intestine reveal a revival stem cell. Nature, 569(7754), 121–125. https://doi.org/10.1038/s41586-019-1154-y
Garcia-Hernandez, V., Raya-Sandino, A., Azcutia, V., Miranda, J., Kelm, M., Flemming, S., Birkl, D., Quiros, M., Brazil, J. C., Parkos, C. A., & Nusrat, A. (2023). Inhibition of Soluble Stem Cell Factor Promotes Intestinal Mucosal Repair. Inflammatory Bowel Diseases, 29(7), 1133–1144. https://doi.org/10.1093/ibd/izad003
Gribov, A., Sill, M., Lück, S., Rücker, F., Döhner, K., Bullinger, L., Benner, A., & Unwin, A. (2010). SEURAT: Visual analytics for the integrated analysis of microarray data. BMC Medical Genomics, 3(1), 21. https://doi.org/10.1186/1755-8794-3-21
Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293–304. https://doi.org/10.1038/s41587-023-01767-y
Holtzman, M. J., Morton, J. D., Shornick, L. P., Tyner, J. W., O’Sullivan, M. P., Antao, A., Lo, M., Castro, M., & Walter, M. J. (2002). Immunity, Inflammation, and Remodeling in the Airway Epithelial Barrier: Epithelial-Viral-Allergic Paradigm. Physiological Reviews, 82(1), 19–46. https://doi.org/10.1152/physrev.00020.2001
Hou, W., Ji, Z., Chen, Z., Wherry, E. J., Hicks, S. C., & Ji, H. (2023). A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. Nature Communications, 14(1), 7286. https://doi.org/10.1038/s41467-023-42841-y
Jin, S., MacLean, A. L., Peng, T., & Nie, Q. (2018). scEpath: Energy landscape-based inference of transition probabilities and cellular trajectories from single-cell transcriptomic data. Bioinformatics, 34(12), 2077–2086. https://doi.org/10.1093/bioinformatics/bty058
Khozyainova, A. A., Valyaeva, A. A., Arbatsky, M. S., Isaev, S. V., Iamshchikov, P. S., Volchkov, E. V., Sabirov, M. S., Zainullina, V. R., Chechekhin, V. I., Vorobev, R. S., Menyailo, M. E., Tyurin-Kuzmin, P. A., & Denisov, E. V. (2023). Complex Analysis of Single-Cell RNA Sequencing Data. Biochemistry (Moscow), 88(2), 231–252. https://doi.org/10.1134/S0006297923020074
Kolde, R. (2010). pheatmap: Pretty Heatmaps (p. 1.0.13) [Dataset]. https://doi.org/10.32614/CRAN.package.pheatmap
La Manno, G., Soldatov, R., Zeisel, A., Braun, E., Hochgerner, H., Petukhov, V., Lidschreiber, K., Kastriti, M. E., Lönnerberg, P., Furlan, A., Fan, J., Borm, L. E., Liu, Z., van Bruggen, D., Guo, J., He, X., Barker, R., Sundström, E., Castelo-Branco, G., … Kharchenko, P. V. (2018). RNA velocity of single cells. Nature, 560(7719), 494–498. https://doi.org/10.1038/s41586-018-0414-6
Li, J., Jiang, Y., Ma, M., Wang, L., Jing, M., Yang, Z., Zhang, M., Chen, K., & Fan, J. (2025). Epithelial cell diversity and immune remodeling in bladder cancer progression: Insights from single-cell transcriptomics. Journal of Translational Medicine, 23(1), 135. https://doi.org/10.1186/s12967-025-06138-6
Li, Y., Baccelli, F., Dhillon, H. S., & Andrews, J. G. (2015). Statistical Modeling and Probabilistic Analysis of Cellular Networks With Determinantal Point Processes. IEEE Transactions on Communications, 63(9), 3405–3422. https://doi.org/10.1109/TCOMM.2015.2456016
Liao, G., Nakayama, T., Zhu, B., Lee, I. T., Yeung, J., Yeo, Y. Y., Chang, Y., Wang, C., Liao, S. C.-K., Nkosi, D., Renteria, A., Bravo, D. T., Overdevest, J. B., Yan, C. H., Zarabanda, D., Gall, P. A., Dholakia, S. S., Borchard, N. A., Yang, A., … Jiang, S. (2025). Multi-scaled transcriptomics of chronically inflamed nasal epithelium reveals immune-epithelial dynamics and tissue remodeling in nasal polyp formation. Immunity, 58(10), 2593-2608.e6. https://doi.org/10.1016/j.immuni.2025.08.009
Lindeboom, R. G. H., Worlock, K. B., Dratva, L. M., Yoshida, M., Scobie, D., Wagstaffe, H. R., Richardson, L., Wilbrey-Clark, A., Barnes, J. L., Kretschmer, L., Polanski, K., Allen-Hyttinen, J., Mehta, P., Sumanaweera, D., Boccacino, J. M., Sungnak, W., Elmentaite, R., Huang, N., Mamanova, L., … Teichmann, S. A. (2024). Human SARS-CoV-2 challenge uncovers local and systemic response dynamics. Nature, 631(8019), 189–198. https://doi.org/10.1038/s41586-024-07575-x
Marc Carlson. (2025). org.Mm.eg.db: Genome wide annotation for Mouse (Version v3.21.0) [Computer software].
Melms, J. C., Biermann, J., Huang, H., Wang, Y., Nair, A., Tagore, S., Katsyv, I., Rendeiro, A. F., Amin, A. D., Schapiro, D., Frangieh, C. J., Luoma, A. M., Filliol, A., Fang, Y., Ravichandran, H., Clausi, M. G., Alba, G. A., Rogava, M., Chen, S. W., … Izar, B. (2021). A molecular single-cell lung atlas of lethal COVID-19. Nature, 595(7865), 114–119. https://doi.org/10.1038/s41586-021-03569-1
Pedersen, T. L. (2019). patchwork: The Composer of Plots (p. 1.3.2) [Dataset]. https://doi.org/10.32614/CRAN.package.patchwork
Potter, S. S. (2018). Single-cell RNA sequencing for the study of development, physiology and disease. Nature Reviews Nephrology, 14(8), 479–492. https://doi.org/10.1038/s41581-018-0021-7
Şenel, S. (2021). An Overview of Physical, Microbiological and Immune Barriers of Oral Mucosa. International Journal of Molecular Sciences, 22(15), 7821. https://doi.org/10.3390/ijms22157821
Wickham, H. (2016). Data Analysis. In H. Wickham, Ggplot2 (pp. 189–201). Springer International Publishing. https://doi.org/10.1007/978-3-319-24277-4_9
Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2014). dplyr: A Grammar of Data Manipulation (p. 1.2.1) [Dataset]. https://doi.org/10.32614/CRAN.package.dplyr
Yu, G. (2024). Thirteen years of clusterProfiler. The Innovation, 5(6), 100722. https://doi.org/10.1016/j.xinn.2024.100722





