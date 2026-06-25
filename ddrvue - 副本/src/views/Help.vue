<template>
  <div class="help">
    <div class="help-container">
      <div class="sidebar">
        <el-menu
            :default-active="activeIndex"
            background-color="#f4f4f5"
            text-color="#333"
            active-text-color="#409eff"
            @select="handleMenuSelect"
        >
          <!-- 1. Introduction -->
          <el-menu-item :index="0" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            1. Introduction
          </el-menu-item>
          <!-- 2. Home -->
          <el-menu-item :index="1" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            2. Home
          </el-menu-item>
          <!-- 3. Browser -->
          <el-menu-item :index="2" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            3. Browser
          </el-menu-item>
          <!-- 4. Molecular（含4-1、4-2子项） -->
          <el-menu-item :index="3" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            4. Molecular
          </el-menu-item>
          <el-menu-item :index="3.1" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            4-1 Base Analysis
          </el-menu-item>
          <!-- 4-1的子项（4.1.1等） -->
          <el-menu-item :index="'3.1.1'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.1.1 Copy Number Variation
          </el-menu-item>
          <el-menu-item :index="'3.1.2'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.1.2 mRNA Expression
          </el-menu-item>
          <el-menu-item :index="'3.1.3'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.1.3 DNA Methylation
          </el-menu-item>
          <el-menu-item :index="'3.1.4'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.1.4 Mutation
          </el-menu-item>
          <el-menu-item :index="'3.1.5'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.1.5 Survival analysis
          </el-menu-item>
          <el-menu-item :index="'3.1.6'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.1.6 Experssion
          </el-menu-item>
          <el-menu-item :index="'3.1.7'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.1.6 CLUSTER
          </el-menu-item>

          <el-menu-item :index="3.2" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            4-2 Immune Analysis
          </el-menu-item>
          <!-- 4-2的子项 -->
          <el-menu-item :index="'3.2.1'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.2.1 Checkpoint
          </el-menu-item>
          <el-menu-item :index="'3.2.2'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.2.2 Chemokine
          </el-menu-item>
<!--          <el-menu-item :index="'3.2.3'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">-->
<!--            4.2.3 Immu cell-->
<!--          </el-menu-item>-->
<!--          <el-menu-item :index="'3.2.4'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">-->
<!--            4.2.4 Immune score-->
<!--          </el-menu-item>-->
          <el-menu-item :index="'3.2.5'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.2.3 Immune inhibitor
          </el-menu-item>
          <el-menu-item :index="'3.2.6'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.2.4 Immune stimulator
          </el-menu-item>
          <el-menu-item :index="'3.2.7'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.2.5 Receptor
          </el-menu-item>
          <el-menu-item :index="'3.2.8'" style="margin-left: 30px;font-size: 14px;height: 35px;display: flex; align-items: center;">
            4.2.6 ssGSEA analysis
          </el-menu-item>

          <!-- 5. Drug（含5-1/5-2/5-3子项） -->
          <el-menu-item :index="4" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            5. Drug
          </el-menu-item>
          <el-menu-item :index="4.1" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            5-1 Drug Overview
          </el-menu-item>
          <el-menu-item :index="4.2" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            5-2 Drug Sensitivity
          </el-menu-item>
          <el-menu-item :index="4.3" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            5-3 Drug Pathway
          </el-menu-item>

          <!-- 6. Network及以下（保持原有） -->
          <el-menu-item :index="5" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            6. Network
          </el-menu-item>
          <el-menu-item :index="5.1" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            6-1 PPI Network
          </el-menu-item>
          <el-menu-item :index="5.2" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            6-2 Cross-talk Network
          </el-menu-item>
          <el-menu-item :index="5.3" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            6-3 TF Network
          </el-menu-item>
          <el-menu-item :index="5.4" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            6-4 ncRNA Network
          </el-menu-item>
          <el-menu-item :index="6" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            7. SL
          </el-menu-item>
          <el-menu-item :index="6.1" style="margin-left: 10px;font-size: 15px;height: 40px;display: flex; align-items: center;">
            7-1 Interaction
          </el-menu-item>
          <el-menu-item :index="7" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            8. Evolution Analysis
          </el-menu-item>
          <el-menu-item :index="8" style="font-weight: bold;font-size: 18px;color: #051d31;height: 50px; display: flex; align-items: center;">
            9. Contact us
          </el-menu-item>
        </el-menu>
      </div>
      <div class="content-wrapper" ref="contentWrapperRef">
        <div class="content">
          <!-- 1. Introduction（保持原有） -->
          <el-card id="introduction" class="content-card">
            <span style="font-weight: bold;font-size: 24px">1. Introduction</span>
            <el-divider></el-divider>
            <p style="font-size: 17px">
              <span style="color: #3a8ee6;font-weight: bold">RepairDis</span>
              RepairDis is a comprehensive online analytical platform dedicated to dissecting DNA damage repair
              mechanisms in cancer. It systematically integrates multi-dimensional authoritative biomedical datasets,
              including TCGA multi-omics data (genomics, transcriptomics, epigenomics, proteomics) for profiling
              repair-related genes, DrugBank drug information for exploring repair gene-drug interactions, biological
              network data from KEGG, Reactome and STRING for constructing regulatory models, and multi-species
              sequence data from NCBI for cross-species evolutionary analysis. Based on standardized processing of
              these multi-source heterogeneous data, RepairDis enables users to perform diverse analyses
              (survival analysis, ssGSEA, drug-target association analysis, pathway visualization) across cancer types,
              providing a one-stop solution to elucidate the roles of DNA damage repair in cancer initiation,
              progression and therapeutic resistance.
            </p>
          </el-card>

          <!-- 2. Home（保持原有） -->
          <el-card id="home" class="content-card">
            <span style="font-weight: bold;font-size: 24px">2. Home</span>
            <el-divider></el-divider>
            <p>
              The <span style="font-weight: bold;color: #3a8ee6">Home</span> page is the entry point of the system, providing a quick search function and an overview of the platform's key features.
              Users can directly search for genes of interest in the search box on the homepage to quickly access detailed information.
            </p>
            <br>
            <div>
              <el-card class="image-card">
                <img src="@/assets/help/home/home.png"
                     alt="Home Page"
                     style="width: 100%">
              </el-card>
              <p style="margin-top: 10px; color: #666; font-size: 14px;">
                <i class="el-icon-info"></i> You can enter gene names (e.g., TP53, BRCA1) in the search box to quickly access gene details.
              </p>
            </div>
          </el-card>

          <!-- 3. Browser（保持原有） -->
          <el-card id="browser" class="content-card">
            <span style="font-weight: bold;font-size: 24px">3. Browser</span>
            <el-divider></el-divider>
            <p>
              The <span style="font-weight: bold;color: #3a8ee6">Browser</span> module allows users to browse and query comprehensive gene information in the database.
              It provides a list of all genes with basic information and supports fuzzy search functionality to quickly locate genes of interest.
            </p>
            <br>
            <div>
              <el-card class="image-card">
                <img src="@/assets/help/browser/Browser.png"
                     alt="Browser Page"
                     style="width: 100%">
              </el-card>
            </div>
          </el-card>

          <!-- 4. Molecular（含4-1、4-2子内容） -->
          <el-card id="molecular" class="content-card">
            <span style="font-weight: bold;font-size: 24px">4. Molecular Landscape</span>
            <el-divider></el-divider>
            <p>
              The <span style="font-weight: bold;color: #3a8ee6">Molecular landscape</span> module provides multi-omics & immune correlation analysis results of DNA damage repair genes,
              helping users understand the molecular mechanisms and immune regulatory roles of genes in cancer.
            </p>
            <br>

            <!-- 4-1 Base Analysis -->
            <el-card id="base-analysis" style="width: 98%;cursor: pointer" class="sub-content-card">
              <span style="font-weight: bold;font-size: 20px;color: #55475b" @click="baseAnalysis = !baseAnalysis">4-1 Base Analysis</span>
              <div v-show="baseAnalysis">
                <el-divider></el-divider>
                <p style="font-size: 17px">
                  This section includes 6 types of basic multi-omics analysis for DNA damage repair genes, covering core molecular characteristics across different cancer types.
                </p>
                <br>
                <div>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/menus.png"
                         alt="Molecular landscape menu Page"
                         style="width: 100%">
                  </el-card>
                </div>

                <!-- 4.1.1 CNV Analysis -->
                <el-card id="cnv-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.1.1 CNV Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The Copy Number Variation results are displayed in this segment.
                    You can either enter the gene you wish to study or use the fuzzy search function for genes
                    in our library to obtain relevant findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/cnv.png" alt="CNV Analysis" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.1.2 mRNA Expression Analysis -->
                <el-card id="mRNA-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.1.2 mRNA Expression Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The Log2FC results of mRNA differential expression analysis are demonstrated in this segment.
                    You can either type in a target gene or conduct a fuzzy search on the genes stored in our database
                    to get the desired findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/analysis.png" alt="CNV Analysis" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.1.3 Methylation Analysis -->
                <el-card id="methylation-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.1.3 Methylation Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The DNA methylation analysis results are displayed in this segment.
                    You can either enter a target gene or use the fuzzy search function for genes in our library to obtain
                    relevant findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/Methy.png" alt="Methylation Analysis" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.1.4 Mutation Analysis -->
                <el-card id="mutation-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.1.4 Mutation Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The gene mutation analysis results are displayed in this segment.
                    You can either enter a target gene or use the fuzzy search function for genes in our library
                    to obtain relevant findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/mutations.png" alt="Mutation Analysis" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.1.5 Survival Analysis -->
                <el-card id="survival-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.1.5 Survival Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The survival analysis results are displayed in this segment. You can either enter the cancer type
                    and a target gene, or use the fuzzy search function for genes in our library to obtain relevant
                    findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/analysis-sur.png" alt="Survival Analysis" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.1.6 Expression Analysis -->
                <el-card id="expression-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.1.6 Expression Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The gene expression profiles across 33 cancer types are displayed in this segment, presented in
                    the form of violin plots. You can either enter a target gene and click the "Click" button,
                    or use the fuzzy search function for genes in our library to obtain relevant findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/Expression.png" alt="Expression Analysis" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.1.7 CLUSTER Analysis -->
                <el-card id="CLUSTER-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.1.7 CLUSTER Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The cancer typing results across 33 cancer types based on the collected genes are displayed in
                    this segment, including the following content: "Sample Consensus Clustering Heatmap (k=2)",
                    "Cluster Stability Analysis Chart", and "ICL Optimal Cluster Number Determination Chart".
                    You can either select a cancer type of interest to view the corresponding results,
                    or use the fuzzy search function for cancer types in our library to obtain relevant findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/CLUSTER.png" alt="Expression Analysis" style="width: 100%">
                  </el-card>
                </el-card>
              </div>
            </el-card>
            <br>

            <!-- 4-2 Immune Analysis -->
            <el-card id="immune-analysis" style="width: 98%;cursor: pointer" class="sub-content-card">
              <span style="font-weight: bold;font-size: 20px;color: #55475b" @click="immuneAnalysis = !immuneAnalysis">4-2 Immune Analysis</span>
              <div v-show="immuneAnalysis">
                <el-divider></el-divider>
                <p style="font-size: 17px">
                  This section focuses on the correlation between DNA damage repair genes and tumor immune microenvironment, including 6 types of immune-related analysis.
                </p>
                <br>
                <div>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/ImmuneMenu.png"
                         alt="Molecular landscape menu Page"
                         style="width: 100%">
                  </el-card>
                </div>

                <!-- 4.2.1 Immune Checkpoint -->
                <el-card id="immune-checkpoint" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.2.1 Immune Checkpoint</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    This section allows you to enter a gene of interest, with the resulting output demonstrating
                    the correlation between the expression levels of this gene and those of immune checkpoint genes,
                    along with a corresponding heatmap graph.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/CheckPoint.png" alt="Immune Checkpoint" style="width: 100%">
                  </el-card>
                </el-card>


                <!-- 4.2.2 Chemokine -->
                <el-card id="chemokine" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.2.2 Chemokine</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    This section allows you to enter a gene of interest, with the resulting output
                    demonstrating the correlation between the expression levels of this gene and those of
                    chemokine genes, along with a corresponding heatmap.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/Chemokine.png" alt="Chemokine" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.2.3 Immune Cell Enrichment Analysis -->
<!--                <el-card id="immune-cell-analysis" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">-->
<!--                  <span style="font-weight: bold;font-size: 18px">4.2.3 Immune Cell</span>-->
<!--                  <el-divider></el-divider>-->
<!--                  <p style="font-size: 16px">-->
<!--                    This section allows you to enter a gene of interest, with the resulting output-->
<!--                    demonstrating the correlation between the expression levels of this gene and immune-->
<!--                    cell percentages, along with a corresponding heatmap graph.-->
<!--                  </p>-->
<!--                  <br>-->
<!--                  <el-card class="image-card">-->
<!--                    <img src="@/assets/help/analysis/imm7.png" alt="Immune Cell " style="width: 100%">-->
<!--                  </el-card>-->
<!--                </el-card>-->

<!--                &lt;!&ndash; 4.2.4 Immune Score &ndash;&gt;-->
<!--                <el-card id="immune-score" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">-->
<!--                  <span style="font-weight: bold;font-size: 18px">4.2.4 Immune Score</span>-->
<!--                  <el-divider></el-divider>-->
<!--                  <p style="font-size: 16px">-->
<!--                    This section allows you to enter a gene of interest, with the resulting output demonstrating-->
<!--                    the correlation between the expression levels of this gene and the three immune scores,-->
<!--                    along with a corresponding heatmap graph.-->
<!--                  </p>-->
<!--                  <br>-->
<!--                  <el-card class="image-card">-->
<!--                    <img src="@/assets/help/analysis/imm6.png" alt="Immune Stimulator" style="width: 100%">-->
<!--                  </el-card>-->
<!--                </el-card>-->

                <!-- 4.2.5 Immune Inhibitor -->
                <el-card id="immune-inhibitor" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.2.3 Immune Inhibitor</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    This section allows you to enter a gene of interest, with the resulting output demonstrating
                    the correlation between the expression levels of this gene and those of immunosuppressive factor
                    genes, along with a corresponding heatmap.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/ImmuneInhibitor.png" alt="Immune Inhibitor" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.2.6 Immune Stimulator -->
                <el-card id="immune-stimulator" style="width: 98%;margin-bottom: 10px" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.2.4 Immune Stimulator</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    This section allows you to enter a gene of interest, with the resulting output
                    demonstrating the correlation between the expression levels of this gene and those
                    of immunostimulatory factor genes, along with a corresponding heatmap.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/ImmuneSt.png" alt="Immune Stimulator" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.2.7 Receptor -->
                <el-card id="receptor" style="width: 98%" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.2.7 Receptor</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    This section allows you to enter a gene of interest, with the resulting output demonstrating
                    the correlation between the expression levels of this gene and those of chemokine receptors,
                    along with a corresponding heatmap.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/ImmuneRec.png" alt="Receptor" style="width: 100%">
                  </el-card>
                </el-card>

                <!-- 4.2.7 ssGSEA Analysis -->
                <el-card id="ssgsea-analysis" style="width: 98%" class="sub-sub-content-card">
                  <span style="font-weight: bold;font-size: 18px">4.2.7 ssGSEA Analysis</span>
                  <el-divider></el-divider>
                  <p style="font-size: 16px">
                    The ssGSEA analysis results are displayed in this segment. You can either enter the cancer
                    type and a target gene, or use the fuzzy search function for genes in our library to obtain
                    relevant findings.
                  </p>
                  <br>
                  <el-card class="image-card">
                    <img src="@/assets/help/analysis/analysis-gsea.png" alt="ssGSEA Analysis" style="width: 100%">
                  </el-card>
                </el-card>
              </div>
            </el-card>
          </el-card>
          <br>
          <!-- 5. Drug -->
          <el-card id="drug" class="content-card">
            <span style="font-weight: bold;font-size: 24px">5. Drug</span>
            <el-divider></el-divider>
            <p>
              The <span style="font-weight: bold;color: #3a8ee6">Drug</span> module provides 3 types of drug-related analysis for DNA damage repair genes, supporting exploration of potential therapeutic strategies.
            </p>
            <br>
            <div>
              <el-card class="image-card">
                <img src="@/assets/help/analysis/DrugMenu.png"
                     alt="Drug menu Page"
                     style="width: 100%">
              </el-card>
            </div>

            <!-- 5-1 Drug Overview -->
            <el-card id="drug-overview" style="width: 98%;cursor: pointer" class="sub-content-card">
              <span style="font-weight: bold;font-size: 20px;color: #55475b" @click="drugOverview = !drugOverview">5-1 Drug Overview</span>
              <div v-show="drugOverview">
                <el-divider></el-divider>
                <p style="font-size: 17px">
                  Displays basic information of drugs targeting DNA damage repair genes, including drug name,
                  target, indication and clinical trial phase.
                </p>
                <br>
                <el-card class="image-card">
                  <img src="@/assets/help/drug/drug.png" alt="Drug Overview" style="width: 100%">
                </el-card>
              </div>
            </el-card>
            <br>
            <!-- 5-2 Drug Sensitivity -->
            <el-card id="drug-sensitivity" style="width: 98%;cursor: pointer" class="sub-content-card">
              <span style="font-weight: bold;font-size: 20px;color: #55475b" @click="drugSensitivity = !drugSensitivity">5-2 Drug Sensitivity</span>
              <div v-show="drugSensitivity">
                <el-divider></el-divider>
                <p style="font-size: 17px">
                  This section allows you to enter a gene of interest to retrieve drug correlation data,
                  where the top six correlation results are presented in the visualization.
                </p>
                <br>
                <el-card class="image-card">
                  <img src="@/assets/help/drug/Re.png" alt="Drug Sensitivity" style="width: 100%">
                </el-card>
              </div>
            </el-card>
            <br>
            <!-- 5-3 Drug Pathway -->
            <el-card id="drug-pathway" style="width: 98%;cursor: pointer" class="sub-content-card">
              <span style="font-weight: bold;font-size: 20px;color: #55475b" @click="drugPathway = !drugPathway">5-3 Drug Pathway</span>
              <div v-show="drugPathway">
                <el-divider></el-divider>
                <p style="font-size: 17px">
                  This section allows you to enter a gene list to retrieve drug pathway data
                  related to the submitted genes.
                </p>
                <br>
                <el-card class="image-card">
                  <img src="@/assets/help/drug/pathway.png" alt="Drug Pathway" style="width: 100%">
                </el-card>
              </div>
            </el-card>
          </el-card>
          <br>

          <!-- 6. Network及以下（保持原有） -->
          <el-card id="network" class="content-card">
            <span style="font-weight: bold;font-size: 24px">6. Network</span>
            <el-divider></el-divider>
            <p>
              The <span style="font-weight: bold;color: #3a8ee6">Network</span> module includes multiple sub-networks to visualize interactions between genes:
            </p>
            <ul style="margin-left: 20px; line-height: 1.8;">
              <li>PPI Network: Protein-Protein Interaction network showing physical interactions between proteins</li>
              <li>Cross-talk Network: Cross-talk interactions between different signaling pathways</li>
              <li>TF Network: Transcription Factor regulatory network demonstrating gene regulation by transcription factors</li>
              <li>ncRNA Network: Non-coding RNA regulatory network showing regulatory relationships involving ncRNAs</li>
            </ul>
            <br>
            <div>
              <el-card id="ppi-network" class="sub-content-card">
                <span style="font-size: 18px; font-weight: bold;">6-1 PPI Network</span>
                <br><br>
                <img src="@/assets/help/network/ppi.png"
                     alt="PPI Network"
                     style="width: 100%">
              </el-card>

              <el-card id="cross-network" class="sub-content-card">
                <span style="font-size: 18px; font-weight: bold;">6-2 Cross-talk Network</span>
                <br><br>
                <img src="@/assets/help/network/cross.png"
                     alt="Cross-talk Network"
                     style="width: 100%">
              </el-card>

              <el-card id="tf-network" class="sub-content-card">
                <span style="font-size: 18px; font-weight: bold;">6-3 TF Network</span>
                <br><br>
                <img src="@/assets/help/network/tf.png"
                     alt="TF Network"
                     style="width: 100%">
              </el-card>

              <el-card id="ncrna-network" class="sub-content-card">
                <span style="font-size: 18px; font-weight: bold;">6-4 ncRNA Network</span>
                <br><br>
                <img src="@/assets/help/network/nc.png"
                     alt="ncRNA Network"
                     style="width: 100%">
              </el-card>
            </div>
          </el-card>

          <el-card id="sl" class="content-card">
            <span style="font-weight: bold;font-size: 24px">7. SL</span>
            <el-divider></el-divider>
            <p>
              The <span style="font-weight: bold;color: #3a8ee6">SL (Synthetic Lethality)</span> module provides information on synthetic lethal interactions between genes.
              Synthetic lethality refers to a genetic interaction where the simultaneous loss of two genes leads to cell death,
              while inactivation of either gene alone does not. This has important implications for cancer therapy.
            </p>
            <br>
            <div>
              <el-card id="interaction" class="sub-content-card">
                <span style="font-size: 18px; font-weight: bold;">7-1 Interaction</span>
                <br><br>
                <img src="@/assets/help/network/sl.png"
                     alt="Synthetic Lethality Interaction"
                     style="width: 100%">
              </el-card>
            </div>
          </el-card>

          <el-card id="evolution" class="content-card">
            <span style="font-weight: bold;font-size: 24px">8. Evolution Analysis</span>
            <el-divider></el-divider>
            <p>
              In the <span style="font-weight: bold;color: #3a8ee6">Evolution Analysis</span> section, you can query cross-species DDR gene conservation using phylogenetic trees.
              The evolution tree supports circular and rectangular layouts for comparing conserved DDR genes across representative species.
            </p>
            <br>
            <el-card class="image-card" style="width: 98%">
              <img src="@/assets/help/evolution/help-evo.png" style="width: 100%">
            </el-card>
          </el-card>

          <el-card id="contact" class="content-card">
            <span style="font-weight: bold;font-size: 24px">9. Contact us</span>
            <el-divider></el-divider>
            <p>
              If you have any questions, suggestions, or need technical support, please contact us through the following ways:
            </p>
            <br>
            <div style="line-height: 1.8; font-size: 16px;">
              <p><i class="el-icon-message"></i> Email: lguo@njupt.edu.cn</p>
              <p><i class="el-icon-info"></i> Please include "RepairDis Inquiry" in the subject line for faster response.</p>
            </div>
            <br>
            <div>
              <el-card class="image-card">
                <img src="@/assets/help/other/other.png"
                     alt="Contact Page"
                     style="width: 100%">
              </el-card>
            </div>
          </el-card>
          <!-- 已删除多余的50px空白div -->
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted, onUnmounted } from 'vue';

// 菜单激活状态
const activeIndex = ref('0');
// 折叠状态变量
const baseAnalysis = ref(false);
const immuneAnalysis = ref(false);
const drugOverview = ref(false);
const drugSensitivity = ref(false);
const drugPathway = ref(false);
// 右侧内容容器引用
const contentWrapperRef = ref(null);
// 防抖定时器
let scrollTimer = null;

// 所有可滚动到的区域ID映射
const sectionIds = {
  '0': 'introduction',
  '1': 'home',
  '2': 'browser',
  '3': 'molecular',
  '3.1': 'base-analysis',
  '3.1.1': 'cnv-analysis',
  '3.1.2': 'mRNA-analysis',
  '3.1.3': 'methylation-analysis',
  '3.1.4': 'mutation-analysis',
  '3.1.5': 'survival-analysis',
  '3.1.6': 'expression-analysis',
  '3.1.7': 'CLUSTER-analysis',
  '3.2': 'immune-analysis',
  '3.2.1': 'immune-checkpoint',
  '3.2.2': 'chemokine',
  // '3.2.3': 'immune-cell-analysis',
  // '3.2.4': 'immune-score',
  '3.2.5': 'immune-inhibitor',
  '3.2.6': 'immune-stimulator',
  '3.2.7': 'receptor',
  '3.2.8': 'ssgsea-analysis',
  '4': 'drug',
  '4.1': 'drug-overview',
  '4.2': 'drug-sensitivity',
  '4.3': 'drug-pathway',
  '5': 'network',
  '5.1': 'ppi-network',
  '5.2': 'cross-network',
  '5.3': 'tf-network',
  '5.4': 'ncrna-network',
  '6': 'sl',
  '6.1': 'interaction',
  '7': 'evolution',
  '8': 'contact'
};

// 滚动到指定区域
const scrollTo = (id) => {
  const element = document.getElementById(id);
  if (element && contentWrapperRef.value) {
    // 滚动到右侧容器内的指定元素
    element.scrollIntoView({
      behavior: 'smooth',
      block: 'start'
    });
    // 点击菜单时自动展开对应内容
    if (id === 'base-analysis') baseAnalysis.value = true;
    if (id === 'immune-analysis') immuneAnalysis.value = true;
    if (id === 'drug-overview') drugOverview.value = true;
    if (id === 'drug-sensitivity') drugSensitivity.value = true;
    if (id === 'drug-pathway') drugPathway.value = true;
    // 4.1.x/4.2.x子项自动展开上级
    if (id.includes('analysis') && !id.includes('immune-') && !id.includes('ssgsea')) {
      baseAnalysis.value = true;
    }
    if (id.includes('immune-') || id === 'chemokine' || id === 'receptor' || id === 'ssgsea-analysis') {
      immuneAnalysis.value = true;
    }
  }
};

// 处理菜单选择
const handleMenuSelect = (index) => {
  activeIndex.value = index.toString();
  const targetId = sectionIds[index];
  if (targetId) {
    scrollTo(targetId);
  }
};

// 处理右侧内容滚动，同步侧边栏高亮
const handleContentScroll = () => {
  // 防抖处理，避免高频触发
  if (scrollTimer) clearTimeout(scrollTimer);
  scrollTimer = setTimeout(() => {
    if (!contentWrapperRef.value) return;

    const scrollTop = contentWrapperRef.value.scrollTop + 50; // 偏移量，优化高亮触发时机
    let currentActiveIndex = '0';
    let maxMatchOffset = Infinity;

    // 遍历所有区域，找到当前可视区域对应的菜单索引
    Object.entries(sectionIds).forEach(([index, id]) => {
      const element = document.getElementById(id);
      if (!element) return;

      const elementTop = element.offsetTop;
      const elementBottom = elementTop + element.offsetHeight;

      // 判断元素是否在可视区域内
      if (scrollTop >= elementTop && scrollTop < elementBottom) {
        const offset = Math.abs(scrollTop - elementTop);
        // 优先选择距离滚动顶部最近的元素
        if (offset < maxMatchOffset) {
          maxMatchOffset = offset;
          currentActiveIndex = index;
        }
      }
    });

    // 更新激活的菜单索引
    if (activeIndex.value !== currentActiveIndex) {
      activeIndex.value = currentActiveIndex;
    }
  }, 100); // 100ms防抖延迟
};

onMounted(() => {
  // 初始滚动到introduction
  scrollTo('introduction');
  // 监听右侧内容容器的滚动事件
  if (contentWrapperRef.value) {
    contentWrapperRef.value.addEventListener('scroll', handleContentScroll);
  }
});

onUnmounted(() => {
  // 移除滚动监听，清除定时器
  if (contentWrapperRef.value) {
    contentWrapperRef.value.removeEventListener('scroll', handleContentScroll);
  }
  clearTimeout(scrollTimer);
});
</script>

<style scoped>
.help {
  width: 100%;
  max-width: 1358px;
  margin: 0 auto;
  background-color: #ffffff;
  padding: 10px 0;
  height: 100vh;
}

.help-container {
  display: flex;
  flex-direction: row;
  gap: 20px;
  padding: 0 20px;
  height: 100%;
}
.sidebar {
  width: 280px;
  padding: 10px;
  border-right: 1px solid #ccc;
  top: 75px;
  display: block;
  text-align: left;
  overflow-y: auto;
  max-height: calc(100% - 20px);
}

.content-wrapper {
  width: 1120px;
  overflow-y: auto;
  max-height: calc(100% - 20px);
  top: 75px;
}

.content {
  width: 100%;
  padding: 10px;
  text-align: left;
  scroll-behavior: smooth;
  height: 100%;
  overflow-y: auto;
}

/* 滚动条美化 */
.content::-webkit-scrollbar, .sidebar::-webkit-scrollbar {
  width: 8px;
  height: 8px;
}
.content::-webkit-scrollbar-track, .sidebar::-webkit-scrollbar-track {
  background: #f1f1f1;
  border-radius: 4px;
}
.content::-webkit-scrollbar-thumb, .sidebar::-webkit-scrollbar-thumb {
  background: #c1c1c1;
  border-radius: 4px;
}
.content::-webkit-scrollbar-thumb:hover, .sidebar::-webkit-scrollbar-thumb:hover {
  background: #a8a8a8;
}

.content p{
  line-height: 1.5;
  font-size: 18px;
  display: block;
  text-align: left;
}

.content ul {
  line-height: 1.8;
  font-size: 17px;
}



/* 激活菜单样式 */
::v-deep .el-menu-item.is-active {
  background-color: #e6f7ff !important;
  border-left: 3px solid #409eff !important;
}

/* 子内容卡片样式 */
.sub-content-card {
  border: 1px solid #e6e6e6;
  border-radius: 8px;
  transition: all 0.3s ease;
  margin-bottom: 15px;
}
.sub-content-card:hover {
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

/* 子子内容卡片样式 */
.sub-sub-content-card {
  border: 1px solid #f0f0f0;
  border-radius: 6px;
  padding: 10px;
}
</style>
