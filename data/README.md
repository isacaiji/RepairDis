# RepairDis Data Files

This directory contains curated tables and compact source-data files used by
the RepairDis manuscript and static web release.

- `tables/TableS1_236_curated_DDR_genes.csv`: curated DDR gene list.
- `tables/TableS2_Final_DPRS_model_genes_and_coefficients.csv`: genes and
  coefficients retained in the final LUAD DPRS model.
- `tables/TableS3_Published_signature_benchmark_comparison.csv`: published
  DDR- or DNA repair-related prognostic signatures used for benchmarking.
- `source_data/Figure2`: plot data and statistical summaries for Figure 2.
- `source_data/Figure3`: plot data and statistical summaries for Figure 3.
- `source_data/Figure4`: DPRS model, validation and benchmark source data.
- `source_data/SupplementaryFigures`: compact source data for Supplementary
  Figures S1-S5 where available.

Large raw and processed omics datasets, including full TCGA/GEO/GDSC matrices
and intermediate analysis workspaces, are not included in this repository.
They should be obtained from the original public resources or regenerated with
the provided R scripts.
