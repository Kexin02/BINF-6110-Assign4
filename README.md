# BINF-6110-Assign4

## 3. Results

# 3.1 Quality control of single-cell RNA-seq data

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
