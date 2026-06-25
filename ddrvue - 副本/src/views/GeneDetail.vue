<template>
  <div class="gene-detail-wrapper" v-if="gene">
    <div class="detail-header">
      <div class="header-inner">
        <el-button class="btn-back" @click="handleBack" plain>
          <el-icon><ArrowLeft /></el-icon>
          <span>Return</span>
        </el-button>

        <div class="header-title-group">
          <h1 class="gene-main-title">{{ gene.geneName }}</h1>
        </div>

        <div class="header-divider"></div>
      </div>
    </div>

    <main class="detail-content">
      <div class="info-card">
        <el-descriptions :column="2" border class="custom-descriptions">
          <el-descriptions-item label="Gene Symbol" label-class-name="label-bg">
            <span class="gene-highlight">{{ gene.geneName }}</span>
          </el-descriptions-item>

          <el-descriptions-item label="Mean MO-DDRweight" label-class-name="label-bg" label-align="center">
            <el-tooltip content="Mean gene-level MO-DDRweight across available TCGA cancer types. Click to view cancer-specific values." placement="top" effect="light">
              <div class="score-interactive-zone" @click="isDialogVisible = true">
                <el-tag :type="getScoreType(getMeanWeight(gene))" effect="dark" class="score-pill">
                  {{ formatScore(getMeanWeight(gene)) }}
                </el-tag>
                <div class="score-meta">
                  <span class="score-status" :class="getScoreType(getMeanWeight(gene))">
                    {{ getScoreStatusText(getMeanWeight(gene)) }}
                  </span>
                  <el-icon class="icon-expand"><Histogram /></el-icon>
                </div>
              </div>
            </el-tooltip>
          </el-descriptions-item>

          <el-descriptions-item :rowspan="4" label="Protein Conformation" label-align="center" label-class-name="label-bg">
            <div class="structure-viewport">
              <Protein :name="gene.geneName" class="viewer-instance" />
              <div class="viewer-label">Source: AlphaFold / PDB</div>
            </div>
          </el-descriptions-item>

          <el-descriptions-item label="Ensembl ID" label-class-name="label-bg">
            <div @click="goToEns(gene.ensembl)" class="data-link">
              {{ gene.ensembl }} <el-icon><TopRight /></el-icon>
            </div>
          </el-descriptions-item>

          <el-descriptions-item label="Functional Role" label-class-name="label-bg">
            <span class="text-secondary">{{ gene.function }}</span>
          </el-descriptions-item>

          <el-descriptions-item label="Primary Literature" label-class-name="label-bg">
            <div @click="goToPubMed(gene.pmid)" class="data-link">
              PMID: {{ gene.pmid }} <el-icon><Link /></el-icon>
            </div>
          </el-descriptions-item>

          <el-descriptions-item label="Reference Title" label-class-name="label-bg">
            <span class="text-secondary italic">{{ gene.title }}</span>
          </el-descriptions-item>

          <el-descriptions-item label="Scientific Abstract" :colspan="2" label-class-name="label-bg">
            <div v-html="gene.abstract" class="abstract-body"></div>
          </el-descriptions-item>
        </el-descriptions>
      </div>

      <div class="analysis-card">
        <div class="card-header-simple">
          <el-icon><DataLine /></el-icon>
          <span>Mutation Heatmap Analysis</span>
        </div>

        <div class="heatmap-container">
          <el-loading v-if="isLoadingHeatmap" text="Rendering Heatmap..." />
          <div v-else-if="heatmapUrl" class="heatmap-img-box">
            <img :src="heatmapUrl" class="heatmap-img" />
            <div class="heatmap-caption">Pan-cancer genomic alteration frequency for {{ gene.geneName }}</div>
          </div>
          <el-empty v-else description="No mutation data available" :image-size="100" />
        </div>
      </div>
    </main>

    <CancerScoreDialog
        v-model="isDialogVisible"
        :gene-name="gene.geneName"
        :score-data="scoreData"
    />
  </div>
</template>

<script setup>
import { ref, onMounted, onUnmounted, watch } from 'vue';
import { useRoute, useRouter } from 'vue-router';
import axios from 'axios';
import Protein from "@/Evolution/structure/Protein.vue";
import CancerScoreDialog from './CancerScoreDialog.vue';
import { ElMessage } from 'element-plus';
import {
  ArrowLeft, Histogram, TopRight, Link, DataLine
} from '@element-plus/icons-vue';

const route = useRoute();
const router = useRouter();
const gene = ref(null);
const isDialogVisible = ref(false);
const isLoadingHeatmap = ref(false);
const heatmapUrl = ref('');
const scoreData = ref([]);
const apiBase = ref('http://121.37.88.191:9016');

const getMeanWeight = (item) => {
  const value = item?.meanMoDdrWeight ?? item?.drfsScore;
  const numericValue = Number(value);
  return Number.isFinite(numericValue) ? numericValue : null;
};
const formatScore = (s) => {
  if (s == null || Number.isNaN(Number(s))) return 'NA';
  return Number(s).toFixed(2);
};
const getScoreType = (s) => {
  if (s == null || Number.isNaN(Number(s))) return 'info';
  return 'primary';
};
const getScoreStatusText = (s) => {
  if (s == null || Number.isNaN(Number(s))) return 'Unavailable';
  return 'Mean across cancers';
};

const fetchGeneDetail = async () => {
  try {
    const response = await axios.get(`${apiBase.value}/api/genes/${route.query.id}`);
    const data = response.data;
    gene.value = data;
    if (data.geneName) {
      fetchHeatmapData(data.geneName);
      fetchCancerScoreData(data.geneName);
    }
  } catch (e) { ElMessage.error('Failed to load gene data'); }
};

const fetchCancerScoreData = async (name) => {
  try {
    const res = await axios.get(`${apiBase.value}/api/genes/${name}/score`);
    const cancerOrder = [
      'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM',
      'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC',
      'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD',
      'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM'
    ];
    scoreData.value = cancerOrder
        .map(name => ({
          name,
          value: Number(res.data?.[name.toLowerCase()])
        }))
        .filter(item => Number.isFinite(item.value));
  } catch (e) {
    scoreData.value = [];
  }
};

const fetchHeatmapData = async (name) => {
  isLoadingHeatmap.value = true;
  try {
    const res = await axios.get(`${apiBase.value}/api/genes/${name}/heatmutation-image`, { responseType: 'blob' });
    heatmapUrl.value = URL.createObjectURL(new Blob([res.data], { type: 'image/png' }));
  } catch (e) {} finally { isLoadingHeatmap.value = false; }
};

const handleBack = () => router.back();
const goToEns = (id) => window.open(`https://www.ensembl.org/id/${id}`, '_blank');
const goToPubMed = (id) => window.open(`https://pubmed.ncbi.nlm.nih.gov/${id}`, '_blank');

onMounted(fetchGeneDetail);
onUnmounted(() => heatmapUrl.value && URL.revokeObjectURL(heatmapUrl.value));
</script>

<style scoped>
/* 核心容器：使用中性灰蓝，建立专业感 */
.gene-detail-wrapper {
  background-color: #f5f8fa;
  min-height: 100vh;
  padding-bottom: 50px;
  font-family: 'Inter', -apple-system, system-ui, sans-serif;
}

/* 页头：保持经典渐变，但微调细节 */
.detail-header {
  background: linear-gradient(135deg, #002766 0%, #004494 100%);
  padding: 40px 0;
  position: relative;
  box-shadow: 0 4px 12px rgba(0,0,0,0.1);
}

.header-inner {
  max-width: 1300px;
  margin: 0 auto;
  padding: 0 40px;
  text-align: center;
}

.gene-main-title {
  color: #fff;
  font-size: 32px;
  font-weight: 800;
  margin: 10px 0;
  letter-spacing: -0.5px;
}

.header-subtitle {
  color: rgba(255,255,255,0.7);
  font-size: 14px;
  text-transform: uppercase;
  letter-spacing: 1px;
}

.btn-back {
  position: absolute;
  left: 40px;
  top: 50%;
  transform: translateY(-50%);
  background: rgba(255,255,255,0.1);
  border: 1px solid rgba(255,255,255,0.2);
  color: #fff;
}

.btn-back:hover { background: rgba(255,255,255,0.2); color: #fff; }

/* 内容卡片 */
.detail-content {
  max-width: 1300px;
  margin: -30px auto 0;
  padding: 0 40px;
}

.info-card, .analysis-card {
  background: #fff;
  border-radius: 12px;
  padding: 30px;
  box-shadow: 0 4px 20px rgba(0, 40, 85, 0.05);
  margin-bottom: 25px;
  border: 1px solid #eef2f7;
}

/* 表格定制 */
:deep(.label-bg) {
  background-color: #f9fbfe !important;
  color: #475a80 !important;
  font-weight: 700 !important;
  width: 160px;
}

.gene-highlight {
  font-size: 20px;
  font-weight: 800;
  color: #003a8c;
}

/* 评分区域交互 */
.score-interactive-zone {
  display: inline-flex;
  align-items: center;
  gap: 15px;
  cursor: pointer;
  padding: 8px 12px;
  border-radius: 8px;
  transition: all 0.2s;
  border: 1px solid transparent;
}

.score-interactive-zone:hover {
  background: #f0f7ff;
  border-color: #d0e3ff;
}

.score-pill { font-size: 18px; font-weight: 800; height: 32px; padding: 0 15px; }

.score-meta { display: flex; flex-direction: column; }
.score-status { font-size: 12px; font-weight: 700; text-transform: uppercase; }
.score-status.primary { color: #1d5fa7; }
.score-status.danger { color: #f56c6c; }
.score-status.warning { color: #e6a23c; }
.score-status.success { color: #67c23a; }
.icon-expand { font-size: 14px; color: #909399; margin-top: 2px; }

/* 蛋白质视图 */
.structure-viewport {
  height: 400px;
  background: #fafbfc;
  border: 1px solid #f0f2f5;
  border-radius: 8px;
  position: relative;
}

.viewer-instance { height: 100%; width: 100%; }
.viewer-label {
  position: absolute; bottom: 10px; right: 10px;
  font-size: 10px; color: #99a9bf; font-weight: 600;
}

/* 链接与文字 */
.data-link {
  color: #1890ff;
  font-weight: 700;
  cursor: pointer;
  display: inline-flex;
  align-items: center;
  gap: 5px;
}
.data-link:hover { text-decoration: underline; color: #003a8c; }

.text-secondary { color: #4b5a75; font-size: 14.5px; line-height: 1.6; }
.italic { font-style: italic; }

.abstract-body {
  line-height: 1.8;
  font-size: 15px;
  color: #334a75;
  text-align: justify;
  background: #fcfdfe;
  padding: 20px;
  border-radius: 8px;
  border: 1px solid #f0f4f8;
}

/* 分析区 */
.card-header-simple {
  display: flex;
  align-items: center;
  gap: 10px;
  margin-bottom: 25px;
  font-size: 18px;
  font-weight: 800;
  color: #003a8c;
}

.heatmap-container {
  min-height: 400px;
  background: #fdfdfe;
  border-radius: 8px;
  border: 1px dashed #d1d9e6;
  display: flex;
  justify-content: center;
  align-items: center;
}

.heatmap-img {
  max-width: 100%;
  border-radius: 4px;
  box-shadow: 0 10px 30px rgba(0,0,0,0.08);
}

.heatmap-caption {
  text-align: center;
  font-size: 12px;
  color: #94a3b8;
  margin-top: 15px;
  font-style: italic;
}
</style>
