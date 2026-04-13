# BINF-6110-Assign4

# 3. Results
The overall quality of the single-cell RNA sequencing dataset was first assessed using standard quality control metrics. The number of detected genes per cell (nFeature_RNA) showed broadly consistent distributions across all timepoints, with similar median values and variability (Figure 1). Likewise, total transcript counts per cell (nCount_RNA) were comparable across conditions, suggesting minimal technical bias introduced during sequencing (Figure 2). A strong positive correlation between nCount_RNA and nFeature_RNA was observed (r = 0.827), indicating that increased sequencing depth was associated with higher gene detection efficiency (Figure 3). Together, these results confirm the overall robustness and consistency of the dataset.

![Distribution of detected features per cell](figures/1.nFeature_RNA.png)

**Figure 1. Distribution of detected genes (nFeature_RNA) across timepoints.** 
Violin plots show comparable gene detection across conditions.

![Distribution of total counts per cell](figures/2.nCount_RNA.png)

**Figure 2. Distribution of total RNA counts (nCount_RNA) across timepoints.** 
Similar distributions suggest consistent sequencing depth.

![Correlation between nCount and nFeature](figures/3.FeatureScatter.png)

**Figure 3. Correlation between nCount_RNA and nFeature_RNA.** 
A strong positive correlation (r = 0.827) indicates good data quality.

Dimensionality reduction and clustering revealed the global transcriptional structure of the dataset. UMAP visualization showed that cells formed multiple distinct clusters, reflecting transcriptionally heterogeneous populations (Figure 5). When colored by timepoint, cells exhibited partial overlap but also displayed shifts in distribution, suggesting dynamic transcriptional changes over time (Figure 4). These results indicate that both cell identity and temporal effects contribute to the overall structure of the data.

![UMAP colored by time](figures/4.UMAP_time.png)

**Figure 4. UMAP colored by timepoint.** 
Cells show partial temporal separation.

![UMAP clusters](figures/5.UMAP_clusters.png)

**Figure 5. UMAP of Seurat clusters.** 
Distinct clusters represent transcriptionally defined cell populations.

To assign biological meaning to these clusters, cell types were annotated based on canonical marker gene expression. Distinct expression patterns of key markers—including Epcam (epithelial), Cd3e (T cells), Lyz2 (myeloid), Krt13 (epithelial subset), and Cnga4 (olfactory neurons)—were observed across clusters (Figure 7). These markers exhibited largely non-overlapping distributions and were spatially localized to specific regions on the UMAP (Figure 8), supporting the accuracy of the annotation. As a result, clusters were classified into major cell types, including epithelial, IFN-responsive epithelial, myeloid, olfactory neuron, T cells, and a remodeling epithelial population (Figure 6).

![Annotated cell types](figures/6.UMAP_celltype.png)

**Figure 6. UMAP of annotated cell types.** 
Clusters are grouped into major biological cell types.

![Marker gene violin plots](figures/7.VlnPlot_markers.png)

**Figure 7. Expression of canonical marker genes.** 
Distinct marker expression supports cell type annotation.

![Marker gene feature plots](figures/8.FeaturePlot_annotation.png)

**Figure 8. Spatial expression of marker genes on UMAP.** 
Markers localize to specific clusters.

In addition to canonical markers, genes associated with infection and immune responses showed distinct expression patterns. Interferon-stimulated genes such as Isg15, Ifit1, and Rsad2 were enriched in specific epithelial and immune clusters (Figure 9), indicating activation of antiviral pathways. These expression patterns suggest that subsets of epithelial cells actively respond to infection through immune-related signaling mechanisms.

![Infection-related genes](figures/9.FeaturePlot_IFN.png)

**Figure 9. Expression of infection-related genes.** 
Interferon-related genes are enriched in specific clusters.

A key focus of this study was cluster 24, which was identified as a remodeling epithelial population. The abundance of this cluster varied across timepoints, increasing from 0.38% at D02 to a peak of 1.25% at D14, followed by a decrease in the naive condition (Figure 10). This temporal pattern suggests that cluster 24 is associated with later stages of infection and may play a role in tissue remodeling or recovery.

![Cluster 24 abundance](figures/10.Cluster24_abundance.png)

**Figure 10. Relative abundance of cluster 24 across timepoints.** 
Cluster 24 peaks at D14.

To further characterize the biological role of cluster 24, functional enrichment analysis was performed on its marker genes. Gene Ontology analysis revealed significant enrichment for processes related to keratinization, RNA processing, translation, and epidermal development (Figure 11). These functions are consistent with active epithelial differentiation and structural remodeling, supporting the classification of cluster 24 as a specialized epithelial subpopulation involved in tissue repair.

![GO enrichment](figures/11.Cluster24_GSEA_dotplot.png)

**Figure 11. GO enrichment of cluster 24 markers.** 
Cluster 24 is enriched for keratinization and RNA-related processes.

Changes in overall cell type composition were also examined across timepoints. Epithelial cells remained the dominant population under all conditions, while IFN-responsive epithelial cells increased during intermediate stages and declined at later timepoints (Figure 12). Immune populations, including myeloid cells and T cells, showed moderate temporal variation, suggesting coordinated immune and epithelial responses during infection progression.

![Cell composition](figures/12.Celltype_composition.png)

**Figure 12. Cell type composition across time.** 
Dynamic changes reflect immune activation and epithelial remodeling.

Finally, marker gene expression analysis confirmed the distinct identity of cluster 24. Genes such as Krt13, Krt6b, Plet1, Csta1, Prss27, and Pglyrp4 were highly enriched in this cluster compared to others (Figure 13). The strong and specific expression of these genes further supports the role of cluster 24 as a remodeling epithelial population with specialized functional properties.

![Dot plot markers](figures/13.DotPlot_cluster24.png)

**Figure 13. Expression of cluster 24 marker genes.** 
Markers are highly enriched in cluster 24.
