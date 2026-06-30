# RepairDis

RepairDis is a web-accessible DNA damage repair (DDR) resource for pan-cancer
DDR gene annotation, multi-omics exploration, DDR-state scoring, immune
microenvironment analysis, drug association browsing, interaction networks,
synthetic-lethality information and evolutionary analysis.

This repository contains the RepairDis front-end code, selected R analysis
scripts used for the manuscript figures, small curated summary tables and
compact source-data files for the final figures. Large public datasets and
generated intermediate matrices are not included in the repository.

## Repository Structure

```text
.
|-- src/                         # Vue 3 front-end source code
|-- public/                      # Static assets and precomputed display files
|-- Rcode/
|   |-- manuscript_figures/       # R scripts for manuscript figure generation
|   `-- dprs_luad/                # R scripts for LUAD DPRS modeling/benchmarking
`-- data/
    |-- tables/                   # Curated supplementary tables
    `-- source_data/              # Compact source data for main/supplementary figures
```

## Front-End

The web interface is implemented with Vue 3 and Vite.

```bash
npm install
npm run dev
npm run build
```

The current repository is arranged as a front-end/static-resource release.
Database deployment files and private server-side configuration are not
included.

## R Analysis Code

R scripts are provided as analysis code associated with the manuscript:

- `Rcode/manuscript_figures`: pan-cancer DDR-state, immune remodeling,
  immunotherapy-response and gene-drug-pathway figure scripts.
- `Rcode/dprs_luad`: LUAD DPRS construction, validation, model comparison and
  related analysis scripts.

The scripts assume that TCGA, GEO, GDSC and processed RepairDis analysis files
are available locally. Paths may need to be adjusted before rerunning.

## Data Availability

Small curated summary tables are included in `data/tables`, and compact source
data for the final figures are included in `data/source_data`. Large-scale
transcriptomic, mutation, copy-number, clinical and drug-response datasets
should be obtained from their original public resources, including TCGA/GDC,
UCSC Xena, GEO and GDSC.

## Reference Implementation

The repository organization was cleaned with reference to the public adiDB
project structure while retaining RepairDis-specific Vue 3 code and DDR
analysis scripts.
