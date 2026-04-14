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

### 3.1 Quality control indicates overall comparable data quality across conditions

Quality control metrics showed that the dataset was broadly comparable across all timepoints. The distribution of detected genes per cell (nFeature_RNA) was similar among conditions, although D05 and D08 tended to show slightly higher medians, whereas D14 displayed a lower median and somewhat reduced upper range (Figure 1). A similar pattern was observed for total transcript counts per cell (nCount_RNA), with D05 and D08 exhibiting relatively higher central values and D14 shifted slightly lower (Figure 2). Despite these modest differences, the overall distributions remained within a similar range, suggesting that the samples were technically comparable and suitable for downstream analysis.

The relationship between sequencing depth and gene detection was further evaluated using a feature scatter plot. A strong positive correlation was observed between nCount_RNA and nFeature_RNA (r = 0.827), indicating that cells with higher transcript counts generally had more detected genes (Figure 3). The overall smooth relationship and lack of a dominant low-quality population support the conclusion that the dataset has acceptable quality for clustering and annotation.

![QC_nFeature](figures/1.QC_nFeature.png)

**Figure 1. Distribution of detected genes per cell across timepoints.** Violin plots show the distribution of nFeature_RNA in D02, D05, D08, D14, and Naive samples, with embedded boxplots indicating the median and interquartile range. D05 and D08 display slightly higher central values, whereas D14 shows a lower median, suggesting modest time-dependent variation in transcriptional complexity while overall data quality remains comparable across conditions. 

![QC_nCount](figures/2.QC_nCount.png)

**Figure 2. Distribution of total transcript counts per cell across timepoints.** Violin plots show nCount_RNA for each condition. D05 and D08 exhibit relatively higher count distributions, while D14 is shifted modestly lower. Black points indicate cells with unusually high count values. Overall, the comparable distributions support consistent sequencing depth across samples.

![QC_scatter](figures/3.QC_scatter.png)

**Figure 3. Relationship between transcript counts and detected genes.** Scatter plot of nCount_RNA versus nFeature_RNA across all cells. Each point represents one cell and is colored by condition. The strong positive Pearson correlation (r = 0.827) indicates that gene detection scales appropriately with sequencing depth, supporting overall dataset quality.

### 3.2 Unsupervised clustering identifies multiple transcriptionally distinct populations

Unsupervised clustering of respiratory mucosa cells revealed extensive transcriptional heterogeneity, with 30 Seurat clusters resolved in UMAP space (Figure 4). Several large populations were clearly separated, including a major left-sided population, multiple right-sided clusters, and several smaller isolated groups. Notably, cluster 24 formed a small and highly separated population at the bottom of the embedding, suggesting a distinct transcriptional identity rather than a transitional state between major clusters.

These clusters were subsequently annotated into broader biological categories using canonical marker genes. The annotated UMAP resolved six major populations: epithelial cells, IFN-responsive epithelial cells, myeloid cells, olfactory neurons, T cells, and a small remodeling epithelial population (Figure 6). The annotated populations were spatially coherent and broadly consistent with the original unsupervised clustering structure, supporting the biological validity of the annotation.

![UMAP_clusters](figures/4.UMAP_clusters.png)

**Figure 4. UMAP visualization of transcriptionally defined Seurat clusters.** Cells were clustered using a graph-based Seurat workflow and visualized by UMAP. Each color and label corresponds to one Seurat cluster. The clear separation between major groups indicates substantial transcriptional heterogeneity within the respiratory mucosa, while the presence of small isolated clusters suggests rare but biologically distinct cell states.

![UMAP_annotation](figures/6.UMAP_annotation.png)
**Figure 6. UMAP visualization of annotated major cell types.** Clusters were grouped into major biological categories based on canonical marker expression. Annotated populations include epithelial, IFN-responsive epithelial, myeloid, olfactory neuron, T-cell, and remodeling epithelial populations. The strong spatial coherence of annotated groups supports the robustness of the manual annotation strategy.

### 3.3 Canonical marker genes support cell type annotation

The annotation was supported by canonical marker expression patterns across clusters. Broad *Epcam* expression marked many epithelial-associated clusters, whereas *Cd3e* was restricted to a small subset corresponding to T cells (Figure 5). *Lyz2* was enriched in several right-sided clusters consistent with myeloid identity, while Cnga4 was concentrated in the large left-sided population annotated as olfactory neurons. In contrast, *Krt13* expression was highly restricted to a small isolated cluster, supporting its interpretation as a specialized epithelial subset rather than a general epithelial program.

Feature plots confirmed these spatial relationships on the UMAP embedding (Figure 7). Epcam was widely distributed across epithelial regions, Cd3e localized to the T-cell compartment, Lyz2 was concentrated in myeloid populations, and Cnga4 marked the olfactory neuron territory. Most importantly, Krt13 was sharply localized to the isolated cluster 24, further supporting its classification as a remodeling epithelial population.

![Marker_violin](figures/5.Marker_violin.png)

**Figure 5. Canonical marker expression across Seurat clusters.** Stacked violin plots show log-normalized expression of *Epcam*, *Cd3e*, *Lyz2*, *Krt13*, and *Cnga4* across clusters. *Epcam* broadly marks epithelial-associated clusters, *Cd3e* is restricted to T-cell clusters, *Lyz2* is enriched in myeloid populations, *Cnga4* marks olfactory neurons, and *Krt13* is highly specific to a small epithelial subset, supporting cluster-level annotation.


![Feature_annotation](figures/7.Feature_annotation.png)

**Figure 7. Feature plots of canonical marker genes used for annotation.** UMAP feature plots show the spatial distribution of *Epcam*, *Cd3e*, *Lyz2*, *Krt13*, and *Cnga4*. Marker expression localizes to distinct UMAP regions, confirming the annotation of epithelial, T-cell, myeloid, olfactory neuron, and remodeling epithelial populations.

### 3.4 Infection-related genes define activated epithelial and immune states

Genes associated with interferon signaling and immune activation showed clear cell type–restricted expression patterns. *Isg15*, *Ifit1*, and *Rsad2* were enriched primarily in upper-right and right-central regions of the UMAP, corresponding to IFN-responsive epithelial and immune-associated populations (Figure 8). *Cxcl16* showed a similar but somewhat more localized pattern, whereas *Cd274* was more sparsely expressed but still concentrated in activated regions.

These results indicate that the response to infection is not uniform across all epithelial cells. Instead, a subset of epithelial and immune cells adopts a distinct antiviral or inflammatory transcriptional program, which supports separating IFN-responsive epithelial cells from baseline epithelial populations in the annotation.

![Feature_story](figures/8.Feature_story.png)

**Figure 8. Feature plots of infection- and immune-associated genes.** UMAP feature plots show the expression of *Isg15*, *Ifit1*, *Rsad2*, *Cxcl16*, and *Cd274*. These genes are enriched in specific upper-right and right-central regions of the UMAP, indicating activation of interferon-related and immune-associated transcriptional programs in a subset of epithelial and immune cells.

### 3.5 A rare remodeling epithelial population expands at later timepoints

Cluster 24 remained rare across all conditions but displayed a clear temporal trend. Its relative abundance was low at D02 (0.00376) and lowest at D05 (0.00196), increased modestly at D08 (0.00499), and then rose sharply at D14 (0.0125) before returning to a lower level in the Naive condition (0.00350) (Figure 9). Although the absolute fraction of cells in this cluster was small, the substantial enrichment at D14 suggests that this population emerges preferentially during a later phase of the response.

The delayed increase in cluster 24 is consistent with a role in epithelial remodeling or tissue recovery rather than immediate antiviral activation. Its temporal behavior also complements the IFN-related expression patterns shown in Figure 8, supporting a model in which early immune activation is followed by later epithelial restructuring.

![Cluster24_abundance](figures/9.Cluster24_abundance.png)

**Figure 9. Relative abundance of cluster 24 across timepoints.** Line plot showing the proportion of cells assigned to cluster 24 in each condition. Cluster 24 remains rare overall but increases markedly at D14, consistent with a late-stage remodeling-associated epithelial response.

### 3.6 Over-representation analysis reveals wound healing and epithelial organization programs in cluster 24

Functional enrichment of cluster 24 marker genes by over-representation analysis (ORA) revealed strong enrichment for biological processes related to epithelial differentiation, tissue repair, and structural organization (Figure 10). The most prominent terms included actin filament organization, regulation of proteolysis, epidermis development, and wound healing, together with cell-cell junction organization, cell-cell junction assembly, apical junction assembly, and keratinocyte differentiation.

These pathways strongly support the interpretation that cluster 24 represents an epithelial population engaged in structural remodeling and restoration of barrier architecture. Enrichment of junction- and membrane-localization terms further suggests that these cells are not merely stressed epithelial cells but are actively participating in rebuilding organized epithelial tissue.

![ORA_dotplot](figures/10.ORA_dotplot.png)

**Figure 10. Gene Ontology enrichment analysis of cluster 24 marker genes by over-representation analysis.** Dot plot showing the top enriched biological processes identified by ORA. The x-axis represents gene ratio, dot size reflects the number of genes contributing to each term, and color indicates adjusted p-value significance. Enrichment is dominated by epithelial differentiation, cytoskeletal organization, junction assembly, and wound-healing processes, supporting a remodeling-associated epithelial phenotype.

### 3.7 GSEA highlights biosynthetic and differentiation-associated transcriptional programs

Gene set enrichment analysis (GSEA) provided complementary evidence that cluster 24 is transcriptionally active and functionally specialized (Figure 11). The most enriched pathways included keratinization, mRNA processing, translation, ribonucleoprotein complex biogenesis, RNA splicing, and cytoplasmic translation. Additional enrichment for protein-containing complex assembly and organization indicates active intracellular restructuring.

Compared with the ORA results, which emphasized extracellular and tissue-level remodeling processes, GSEA highlights the intracellular biosynthetic programs that likely support this remodeling phenotype. Together, these results suggest that cluster 24 is not only structurally specialized but also metabolically active, with coordinated activation of differentiation and protein synthesis programs.

![GSEA_dotplot](figures/11.Cluster24_GSEA_dotplot.png)

**Figure 11. GSEA of cluster 24 marker genes.**  
Gene set enrichment analysis reveals activation of pathways related to keratinization, RNA processing, and translation, indicating active epithelial remodeling.

### 3.8 Cell type composition shifts over time indicate coordinated immune activation and epithelial remodeling

Changes in overall cell type composition revealed coordinated temporal dynamics across the tissue (Figure 12). Epithelial cells remained the dominant population at all timepoints, but their relative abundance decreased at intermediate stages, particularly D05 and D08. In contrast, IFN-responsive epithelial cells increased noticeably at D05 and remained elevated at D08, consistent with activation of antiviral signaling during the early-to-mid response.

Myeloid cells also increased during the time course and were particularly prominent at D14, suggesting sustained or secondary immune involvement. T cells showed their strongest relative representation at D08, indicating an intermediate-stage adaptive immune component. Olfactory neurons remained comparatively stable, whereas remodeling epithelial cells remained a minor population overall, consistent with their rarity in Figure 9. Taken together, these shifts suggest a progression from early epithelial/immune activation toward later tissue remodeling.

![Celltype_composition](figures/12.Celltype_composition.png)

**Figure 12. Cell type composition across timepoints.** Stacked bar plot showing the relative proportions of annotated cell types in each condition. Epithelial cells dominate across all timepoints, while IFN-responsive epithelial, myeloid, and T-cell populations show dynamic changes over time. The overall pattern suggests coordinated progression from immune activation at intermediate timepoints to later epithelial remodeling.

### 3.9 Marker gene expression confirms the distinct identity of cluster 24

The dot plot of representative marker genes confirmed that cluster 24 has a highly specific transcriptional signature (Figure 13). *Krt13*, *Krt6b*, *Plet1*, *Csta1*, *Prss27*, *Calml3*, and *Pglyrp4* all showed strong enrichment in cluster 24, with both high average expression and high percentages of expressing cells. Most of these markers showed minimal expression in other clusters, although Plet1 displayed weaker expression in a small number of non-cluster 24 cells.

The combination of keratin-associated genes, epithelial differentiation markers, and barrier-associated genes strongly supports the interpretation that cluster 24 represents a distinct remodeling epithelial population rather than a generic epithelial or inflammatory state.

![Marker_dotplot](figures/13.Marker_dotplot.png)

**Figure 13. Expression of representative marker genes for cluster 24 across all clusters.** Dot plot showing average expression (color intensity) and percentage of expressing cells (dot size) for selected marker genes across clusters. Cluster 24 exhibits strong and highly specific enrichment of *Krt13*, *Krt6b*, *Plet1*, *Csta1*, *Prss27*, *Calml3*, and *Pglyrp4*, confirming its identity as a distinct remodeling epithelial population.


## 4. Discussion
This study provides a detailed characterization of cellular heterogeneity and transcriptional dynamics in the respiratory mucosa during infection. By combining clustering, marker-based annotation, temporal composition analysis, and functional enrichment, the results reveal both relatively stable major cell identities and dynamic state changes associated with immune activation and epithelial remodeling. In particular, the data support a model in which early and intermediate responses are dominated by interferon-associated epithelial and immune programs, whereas later stages are marked by the emergence of a distinct remodeling epithelial population.

A key finding of this study is the identification of a rare but transcriptionally distinct remodeling epithelial population, cluster 24, which shows a clear temporal increase and peaks at D14. This cluster is characterized by strong expression of genes associated with epithelial differentiation and remodeling, including *Krt13*, *Krt6b*, *Plet1*, *Csta1*, and *Prss27*. The high specificity of these markers across clusters argues that cluster 24 is not simply a stressed epithelial subset, but rather a biologically distinct cell state. These features are consistent with epithelial remodeling and barrier repair processes observed in inflamed or regenerating mucosal tissues and align with recent work showing that epithelial injury responses are often accompanied by keratin-associated transcriptional reprogramming (Garcia-Hernandez et al., 2023; Liao et al., 2025). The ORA results further strengthen this interpretation by showing enrichment for wound healing, epidermis development, keratinocyte differentiation, and cell-cell junction organization, all of which are hallmarks of structural epithelial repair.

In addition to structural remodeling signatures, cluster 24 also showed enrichment for RNA-processing- and translation-related pathways in GSEA, including mRNA processing, RNA splicing, translation, and ribonucleoprotein complex biogenesis. This suggests that these cells are not only differentiated in identity but also highly biosynthetically active. Such a profile is often associated with epithelial populations undergoing active regeneration, where cellular restructuring is coupled with elevated transcriptional and translational demand (Ho et al., 2021). The coexistence of keratinization-related pathways and biosynthetic pathways is particularly notable, because it implies that cluster 24 may represent an actively differentiating repair state rather than a terminally quiescent epithelial population. Similar transcriptional signatures have been reported in regenerating epithelial compartments following injury, where restoration of tissue integrity requires both lineage-specific differentiation and coordinated macromolecular synthesis.

Another important result is the presence of interferon-responsive epithelial and immune-associated populations marked by genes such as *Isg15*, *Ifit1*, and *Rsad2*. These genes were enriched in restricted subsets of cells rather than uniformly across the tissue, indicating that infection-induced antiviral responses are spatially and cell-type specific. This selective activation is biologically plausible, as epithelial cells are increasingly recognized not only as passive barrier cells but also as active participants in innate immune sensing and signaling. Previous studies have shown that epithelial cells can initiate and amplify antiviral responses through interferon-mediated programs, thereby shaping the broader tissue immune environment (Dalskov et al., 2023; Liu et al., 2024; Major et al., 2020). In this context, the IFN-responsive epithelial compartment identified here likely represents an activated epithelial state that precedes or accompanies later tissue repair programs.

Changes in cell composition across time further support the interpretation of a coordinated and staged response. While epithelial cells remained the dominant population throughout the dataset, IFN-responsive epithelial cells were expanded at intermediate stages, and myeloid and T-cell populations also showed temporal variation. These patterns are consistent with a transition from early inflammatory and antiviral activation toward later epithelial restructuring. Similar shifts in cell abundance and state have been reported in infection and inflammation models, where epithelial injury is accompanied by immune recruitment, transient activation states, and subsequent restoration programs (Kim et al., 2023; Oikonomou et al., 2021; Wynn & Vannella, 2016). The increase of cluster 24 at D14 is therefore especially informative, because it suggests that remodeling-associated epithelial programs intensify after the peak of early immune signaling rather than simultaneously with it.

Despite these insights, several limitations should be considered. First, cell type annotation relied primarily on canonical markers, which may not fully resolve intermediate or transitional states. This is particularly relevant for epithelial populations, where activation, injury, and differentiation states can partially overlap transcriptionally. Second, clustering outcomes depend on analytical choices such as dimensionality and resolution, and alternative parameter settings could merge or split some of the observed populations. Third, enrichment analyses depend on existing gene annotations and curated pathway databases, which may bias interpretation toward well-characterized biological processes while underrepresenting novel or poorly annotated functions.

A further limitation is that the present study examines discrete timepoints rather than continuous cellular trajectories. As a result, the data support temporal association but do not directly establish lineage relationships between IFN-responsive epithelial states and remodeling epithelial states. Methods such as pseudotime inference, RNA velocity, or lineage tracing would be valuable for testing whether cluster 24 emerges from pre-existing epithelial populations that undergo progressive state transitions during recovery (Hou et al., 2023; La Manno et al., 2018). In addition, experimental validation of marker genes such as *Krt13*, *Krt6b*, and *Plet1* would strengthen the interpretation that cluster 24 represents a true repair-associated epithelial population in vivo.

Overall, this study demonstrates the power of single-cell transcriptomics in uncovering dynamic cellular responses to infection in the respiratory mucosa. Rather than revealing only static cell identities, the analysis captures a temporal sequence of epithelial and immune state changes, from interferon-associated activation to late remodeling-associated differentiation. The identification of a remodeling epithelial population that peaks at D14 and is enriched for wound-healing, keratinization, and RNA-processing programs provides new insight into how mucosal tissues may transition from immune defense to structural repair. These findings highlight cluster 24 as a biologically interesting candidate for future investigation into epithelial recovery, barrier restoration, and post-infectious tissue remodeling.


## 5. References
Dalskov, L., Gad, H. H., & Hartmann, R. (2023). Viral recognition and the antiviral interferon response. *The EMBO Journal*, *42*(14), e112907. https://doi.org/10.15252/embj.2022112907

Garcia-Hernandez, V., Raya-Sandino, A., Azcutia, V., Miranda, J., Kelm, M., Flemming, S., Birkl, D., Quiros, M., Brazil, J. C., Parkos, C. A., & Nusrat, A. (2023). Inhibition of Soluble Stem Cell Factor Promotes Intestinal Mucosal Repair. *Inflammatory Bowel Diseases*, *29*(7), 1133–1144. https://doi.org/10.1093/ibd/izad003

Gribov, A., Sill, M., Lück, S., Rücker, F., Döhner, K., Bullinger, L., Benner, A., & Unwin, A. (2010). SEURAT: Visual analytics for the integrated analysis of microarray data. *BMC Medical Genomics*, *3*(1), 21. https://doi.org/10.1186/1755-8794-3-21

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, *42*(2), 293–304. https://doi.org/10.1038/s41587-023-01767-y

Ho, J. J. D., Man, J. H. S., Schatz, J. H., & Marsden, P. A. (2021). Translational remodeling by RNA-binding proteins and noncoding RNAs. *WIREs RNA*, *12*(5), e1647. https://doi.org/10.1002/wrna.1647

Holtzman, M. J., Morton, J. D., Shornick, L. P., Tyner, J. W., O’Sullivan, M. P., Antao, A., Lo, M., Castro, M., & Walter, M. J. (2002). Immunity, Inflammation, and Remodeling in the Airway Epithelial Barrier: Epithelial-Viral-Allergic Paradigm. *Physiological Reviews*, *82*(1), 19–46. https://doi.org/10.1152/physrev.00020.2001

Hou, W., Ji, Z., Chen, Z., Wherry, E. J., Hicks, S. C., & Ji, H. (2023). A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. *Nature Communications*, *14*(1), 7286. https://doi.org/10.1038/s41467-023-42841-y

Jin, S., MacLean, A. L., Peng, T., & Nie, Q. (2018). scEpath: Energy landscape-based inference of transition probabilities and cellular trajectories from single-cell transcriptomic data. *Bioinformatics*, *34*(12), 2077–2086. https://doi.org/10.1093/bioinformatics/bty058

Khozyainova, A. A., Valyaeva, A. A., Arbatsky, M. S., Isaev, S. V., Iamshchikov, P. S., Volchkov, E. V., Sabirov, M. S., Zainullina, V. R., Chechekhin, V. I., Vorobev, R. S., Menyailo, M. E., Tyurin-Kuzmin, P. A., & Denisov, E. V. (2023). Complex Analysis of Single-Cell RNA Sequencing Data. *Biochemistry (Moscow)*, *88*(2), 231–252. https://doi.org/10.1134/S0006297923020074

Kim, J., Kim, S., Lee, S.-Y., Jo, B.-K., Oh, J.-Y., Kwon, E.-J., Kim, K.-T., Adpaikar, A. A., Kim, E.-J., Jung, H.-S., Kim, H.-R., Roe, J.-S., Hong, C. P., Kim, J. K., Koo, B.-K., & Cha, H.-J. (2023). Partial in vivo reprogramming enables injury-free intestinal regeneration via autonomous Ptgs1 induction. *Science Advances*, *9*(47), eadi8454. https://doi.org/10.1126/sciadv.adi8454

Kolde, R. (2010). pheatmap: Pretty Heatmaps (p. 1.0.13) [Dataset]. https://doi.org/10.32614/CRAN.package.pheatmap

La Manno, G., Soldatov, R., Zeisel, A., Braun, E., Hochgerner, H., Petukhov, V., Lidschreiber, K., Kastriti, M. E., Lönnerberg, P., Furlan, A., Fan, J., Borm, L. E., Liu, Z., van Bruggen, D., Guo, J., He, X., Barker, R., Sundström, E., Castelo-Branco, G., … Kharchenko, P. V. (2018). RNA velocity of single cells. *Nature*, *560*(7719), 494–498. https://doi.org/10.1038/s41586-018-0414-6

Li, J., Jiang, Y., Ma, M., Wang, L., Jing, M., Yang, Z., Zhang, M., Chen, K., & Fan, J. (2025). Epithelial cell diversity and immune remodeling in bladder cancer progression: Insights from single-cell transcriptomics. *Journal of Translational Medicine*, *23*(1), 135. https://doi.org/10.1186/s12967-025-06138-6

Li, Y., Baccelli, F., Dhillon, H. S., & Andrews, J. G. (2015). Statistical Modeling and Probabilistic Analysis of Cellular Networks With Determinantal Point Processes. *IEEE Transactions on Communications*, *63*(9), 3405–3422. https://doi.org/10.1109/TCOMM.2015.2456016

Liao, G., Nakayama, T., Zhu, B., Lee, I. T., Yeung, J., Yeo, Y. Y., Chang, Y., Wang, C., Liao, S. C.-K., Nkosi, D., Renteria, A., Bravo, D. T., Overdevest, J. B., Yan, C. H., Zarabanda, D., Gall, P. A., Dholakia, S. S., Borchard, N. A., Yang, A., … Jiang, S. (2025). Multi-scaled transcriptomics of chronically inflamed nasal epithelium reveals immune-epithelial dynamics and tissue remodeling in nasal polyp formation. *Immunity*, *58*(10), 2593-2608.e6. https://doi.org/10.1016/j.immuni.2025.08.009

Liu, Y.-G., Jin, S.-W., Zhang, S.-S., Xia, T.-J., Liao, Y.-H., Pan, R.-L., Yan, M.-Z., & Chang, Q. (2024). Interferon lambda in respiratory viral infection: Immunomodulatory functions and antiviral effects in epithelium. *Frontiers in Immunology*, 15, 1338096. https://doi.org/10.3389/fimmu.2024.1338096

Major, J., Crotta, S., Llorian, M., McCabe, T. M., Gad, H. H., Priestnall, S. L., Hartmann, R., & Wack, A. (2020). Type I and III interferons disrupt lung epithelial repair during recovery from viral infection. *Science (New York, N.y.)*, *369*(6504), 712–717. https://doi.org/10.1126/science.abc2061

Marc Carlson. (2025). org.Mm.eg.db: Genome wide annotation for Mouse (Version v3.21.0) [Computer software].

Oikonomou, N., Schuijs, M. J., Chatzigiagkos, A., Androulidaki, A., Aidinis, V., Hammad, H., Lambrecht, B. N., & Pasparakis, M. (2021). Airway epithelial cell necroptosis contributes to asthma exacerbation in a mouse model of house dust mite-induced allergic inflammation. *Mucosal Immunology*, *14*(5), 1160–1171. https://doi.org/10.1038/s41385-021-00415-5

Pedersen, T. L. (2019). patchwork: The Composer of Plots (p. 1.3.2) [Dataset]. https://doi.org/10.32614/CRAN.package.patchwork

Potter, S. S. (2018). Single-cell RNA sequencing for the study of development, physiology and disease. *Nature Reviews Nephrology*, *14*(8), 479–492. https://doi.org/10.1038/s41581-018-0021-7

Şenel, S. (2021). An Overview of Physical, Microbiological and Immune Barriers of Oral Mucosa. *International Journal of Molecular Sciences*, *22*(15), 7821. https://doi.org/10.3390/ijms22157821

Wickham, H. (2016). Data Analysis. In H. Wickham, Ggplot2 (pp. 189–201). Springer International Publishing. https://doi.org/10.1007/978-3-319-24277-4_9

Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2014). dplyr: A Grammar of Data Manipulation (p. 1.2.1) [Dataset]. https://doi.org/10.32614/CRAN.package.dplyr

Wynn, T. A., & Vannella, K. M. (2016). Macrophages in Tissue Repair, Regeneration, and Fibrosis. *Immunity*, *44*(3), 450–462. https://doi.org/10.1016/j.immuni.2016.02.015

Yu, G. (2024). Thirteen years of clusterProfiler. *The Innovation*, *5*(6), 100722. https://doi.org/10.1016/j.xinn.2024.100722

