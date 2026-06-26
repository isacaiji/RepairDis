<template>
  <div class="expression-analysis-container">
    <!-- Title Section - 简化标题描述 -->
    <div class="analysis-header">
      <h2>Expression Analysis</h2>
    </div>

    <!-- Search Section - 搜索框和按钮在同一行 -->
    <div class="search-container">
      <div class="search-row">
        <!-- 自定义样式的 el-input 搜索框 -->
        <el-input
            v-model="searchInput"
            placeholder="Search gene symbol..."
            class="custom-search-input"
            @input="handleInputDebounce"
            @focus="showDropdown = true"
            @blur="handleBlur"
        >
          <template #prefix>
            <el-icon class="search-icon">
              <Search />
            </el-icon>
          </template>
          <template #suffix>
            <el-icon
                class="clear-icon"
                v-if="searchInput"
                @click="clearSearch"
            >
              <CircleClose />
            </el-icon>
          </template>
        </el-input>

        <!-- 分析按钮 - 与搜索框同行 -->
        <el-button
            type="primary"
            @click="fetchExpressionData"
            :disabled="!selectedGene || isLoading"
            class="analysis-button"
        >
          <el-icon v-if="isLoading" class="loading-icon">
            <Loading />
          </el-icon>
          <span v-else>Click</span>
        </el-button>
      </div>
      <!--      下拉框-->
      <div
          class="dropdown-container"
          v-show="showDropdown && filteredGenes.length > 0"
      >
        <div
            class="dropdown-item"
            v-for="(gene, index) in filteredGenes"
            :key="index"
            @click="selectGene(gene)"
            :class="{ active: gene === highlightedGene }"
            @mouseenter="highlightedGene = gene"
        >
          {{ gene }}
        </div>
      </div>
    </div>

    <!-- Result Display Section -->
    <div class="result-container">
      <!-- Empty State -->
      <el-empty
          v-if="!showResult && !isLoading && !errorMsg"
          description="No data to display"
          class="empty-state"
      >
        <template #image>
          <el-icon class="empty-icon">
            <Picture/>
          </el-icon>
        </template>
        <template #description>
          <p>Please select a gene and click "Analyze"</p>
        </template>
      </el-empty>

      <!-- Error State -->
      <div v-if="errorMsg && !isLoading" class="error-state">
        <el-icon class="error-icon"><Warning /></el-icon>
        <p class="error-text">{{ errorMsg }}</p>
      </div>

      <!-- Loading State -->
      <div v-if="isLoading" class="loading-state">
        <el-loading
            v-loading="isLoading"
            text="Retrieving data..."
            spinner="el-icon-loading"
            background="rgba(255, 255, 255, 0.9)"
            class="loading-overlay"
        >
        </el-loading>
      </div>

      <!-- Result Content - 下载按钮在图片下方 -->
      <div v-if="showResult && !isLoading && !errorMsg" class="result-content">
        <h3 class="result-title">{{ selectedGene }} Expression Analysis</h3>

        <!-- 图片容器（直接用el-image加载后端PNG接口） -->
        <div class="image-container">
          <el-image
              :src="`${baseUrl}/expression/img?gene=${selectedGene}`"
              fit="contain"
              class="expression-image"
              :preview-src-list="['' + baseUrl + '/expression/img?gene=' + selectedGene]"
              lazy
          >
            <template #error>
              <div class="image-error">
                <el-icon class="error-icon"><Warning /></el-icon>
                <p>Failed to load image. Please try again.</p>
              </div>
            </template>
          </el-image>
        </div>

        <!-- 下载按钮 - 调用后端PDF下载接口 -->
        <div class="download-container">
          <el-button
              type="success"
              icon="Download"
              @click="downloadPDF"
              class="download-button"
              :loading="isDownloading"
          >
            <el-icon v-if="isDownloading" class="loading-icon">
              <Loading />
            </el-icon>
            <span v-else>Download PDF</span>
          </el-button>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import {ref} from 'vue';
import axios from 'axios';
import {Search, CircleClose, Loading, Picture, Warning, Download} from '@element-plus/icons-vue';
import {ElMessage} from 'element-plus';

// 基础URL配置（与后端对应）
const baseUrl = ref('http://121.37.88.191:9016');

// DDR基因列表
const allGenes = ref([
  "ALKBH2", "ALKBH3", "APEX1", "APEX2", "APLF", "APTX", "ATM", "ATR", "ATRIP", "ATRX",
  "BARD1", "BLM", "BRCA1", "BRCA2", "BRIP1", "CCNH", "CDK7", "CETN2", "CHAF1A", "CHEK1",
  "CHEK2", "CLK2", "DCLRE1A", "DCLRE1B", "DCLRE1C", "DDB1", "DDB2", "DMC1", "DNA2", "DNPH1",
  "DNTT", "DUT", "EME1", "EME2", "ENDOV", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5",
  "ERCC6", "ERCC8", "EXO1", "EXO5", "FAAP100", "FAAP20", "FAAP24", "FAN1", "FANCA", "FANCB",
  "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "FEN1", "GEN1",
  "GTF2E2", "GTF2H1", "GTF2H2", "GTF2H3", "GTF2H4", "GTF2H5", "HELQ", "HERC2", "HFM1", "HLTF",
  "HMCES", "HUS1", "LIG1", "LIG3", "LIG4", "MAD2L2", "MBD4", "MDC1", "MGMT", "MLH1",
  "MLH3", "MMS19", "MNAT1", "MPG", "MPLKIP", "MSH2", "MSH3", "MSH4", "MSH5", "MSH6",
  "MUS81", "MUTYH", "NABP2", "NBN", "NEIL1", "NEIL2", "NEIL3", "NHEJ1", "NTHL1", "NUDT1",
  "NUDT15", "NUDT18", "OGG1", "PALB2", "PARG", "PARK7", "PARP1", "PARP2", "PARP3", "PARPBP",
  "PAXIP1", "PCNA", "PDS5B", "PER1", "PMS1", "PMS2", "PNKP", "POLA1", "POLB", "POLD1",
  "POLD2", "POLD3", "POLD4", "POLE", "POLE2", "POLE3", "POLE4", "POLG", "POLH", "POLI",
  "POLK", "POLL", "POLM", "POLN", "POLQ", "PRIMPOL", "PRKDC", "PRPF19", "RAD1", "RAD17",
  "RAD18", "RAD23A", "RAD23B", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAD52", "RAD54B",
  "RAD54L", "RAD9A", "RBBP8", "RDM1", "RECQL", "RECQL4", "RECQL5", "REV1", "REV3L", "RIF1",
  "RMI1", "RNF168", "RNF4", "RNF8", "RPA1", "RPA2", "RPA3", "RPA4", "RRM2B", "SETMAR",
  "SHPRH", "SLX1A", "SLX1B", "SLX4", "SMC5", "SMC6", "SMUG1", "SPIDR", "SPO11", "SPRTN",
  "SWI5", "SWSAP1", "TDG", "TDP1", "TDP2", "TOP3A", "TOPBP1", "TP53", "TP53BP1", "TREX1",
  "TREX2", "UBE2A", "UBE2B", "UBE2N", "UBE2T", "UBE2V2", "UNG", "USP1", "UVSSA", "WDR48",
  "WRN", "XAB2", "XPA", "XPC", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "XRCC5", "XRCC6", "ZSWIM7"
]);

// 响应式状态
const searchInput = ref('');
const selectedGene = ref('');
const filteredGenes = ref([]);
const showDropdown = ref(false);
const highlightedGene = ref('');
const isLoading = ref(false); // 图片加载状态
const isDownloading = ref(false); // 下载状态
const showResult = ref(false);
const debounceTimer = ref(null);
const errorMsg = ref('');

// 防抖处理输入事件
const handleInputDebounce = (value) => {
  if (debounceTimer.value) clearTimeout(debounceTimer.value);

  if (!value.trim()) {
    filteredGenes.value = [];
    showDropdown.value = false;
    return;
  }

  debounceTimer.value = setTimeout(() => {
    const lowerCaseValue = value.toLowerCase();
    filteredGenes.value = allGenes.value.filter(gene =>
        gene.toLowerCase().includes(lowerCaseValue)
    );
    showDropdown.value = true;
  }, 300);
};

// 失去焦点延迟关闭下拉框
const handleBlur = () => {
  setTimeout(() => {
    showDropdown.value = false;
  }, 200);
};

// 清空搜索
const clearSearch = () => {
  searchInput.value = '';
  filteredGenes.value = [];
  showDropdown.value = false;
  selectedGene.value = '';
  errorMsg.value = '';
  showResult.value = false;
};

// 选择基因
const selectGene = (gene) => {
  selectedGene.value = gene;
  searchInput.value = gene;
  showDropdown.value = false;
  highlightedGene.value = '';
  errorMsg.value = '';
};

// 获取基因图片（仅验证接口可用性，图片直接通过el-image src加载）
const fetchExpressionData = async () => {
  if (!selectedGene.value) {
    ElMessage.warning('Please select a gene first');
    return;
  }

  isLoading.value = true;
  showResult.value = false;
  errorMsg.value = '';

  try {
    const response = await axios.head(`${baseUrl.value}/expression/img?gene=${selectedGene.value}`);
    showResult.value = true;
  } catch (error) {
    console.error('Error fetching expression data:', error);
    if (error.response) {
      if (error.response.status === 404) {
        // 解析后端返回的404提示
        errorMsg.value = error.response.data
            ? new TextDecoder().decode(error.response.data)
            : `Image not found for gene: ${selectedGene.value}`;
      } else {
        errorMsg.value = 'Failed to retrieve data. Server error.';
      }
    } else if (error.request) {
      errorMsg.value = 'Failed to connect to server. Please check network.';
    } else {
      errorMsg.value = 'An unexpected error occurred. Please try again.';
    }
    ElMessage.error(errorMsg.value);
  } finally {
    isLoading.value = false;
  }
};

// 下载PDF（调用后端专门的PDF接口）
const downloadPDF = () => {
  if (!selectedGene.value) return;

  isDownloading.value = true;
  errorMsg.value = '';

  try {
    axios.get(`${baseUrl.value}/expression/pdf?gene=${selectedGene.value}`, {
      responseType: 'blob' // 接收二进制文件
    }).then(response => {
      const blob = new Blob([response.data], {type: 'application/pdf'});
      const url = window.URL.createObjectURL(blob);
      const link = document.createElement('a');

      // 从响应头获取文件名（后端已设置）
      const disposition = response.headers['content-disposition'];
      let fileName = `${selectedGene.value}_expression.pdf`;
      if (disposition) {
        const match = disposition.match(/filename="?([^";]+)"?/);
        if (match && match[1]) fileName = match[1];
      }

      link.href = url;
      link.download = fileName;
      document.body.appendChild(link);
      link.click();

      // 释放资源
      document.body.removeChild(link);
      window.URL.revokeObjectURL(url);

      ElMessage.success('PDF downloaded successfully');
    }).catch(error => {
      console.error('Error downloading PDF:', error);
      let errMsg = 'Failed to download PDF. Please try again.';
      if (error.response) {
        if (error.response.status === 404) {
          errMsg = error.response.data
              ? new TextDecoder().decode(error.response.data)
              : `PDF not found for gene: ${selectedGene.value}`;
        } else if (error.response.status === 500) {
          errMsg = 'Server error. Please try again later.';
        }
      } else if (error.request) {
        errMsg = 'Failed to connect to server. Please check network.';
      }
      errorMsg.value = errMsg;
      ElMessage.error(errMsg);
    }).finally(() => {
      isDownloading.value = false;
    });
  } catch (error) {
    console.error('Download error:', error);
    errorMsg.value = 'Failed to download PDF. Please try again.';
    ElMessage.error(errorMsg.value);
    isDownloading.value = false;
  }
};
</script>

<style scoped>
/* 全局容器 */
.expression-analysis-container {
  max-width: 1200px;
  margin: 0 auto;
  padding: 30px 24px;
  background-color: #ffffff;
  border-radius: 16px;
  box-shadow: 0 10px 30px rgba(0, 0, 0, 0.08);
}

/* 标题区域 */
.analysis-header {
  text-align: center;
  margin-bottom: 30px;
  padding-bottom: 15px;
  border-bottom: 2px solid #e0e7ff;
}

.analysis-header h2 {
  color: #003f88;
  font-size: 26px;
  margin: 0;
  font-weight: 600;
  letter-spacing: 0.5px;
}

/* 搜索区域 */
.search-container {
  background-color: rgba(248, 250, 252, 0.8);
  padding: 24px;
  border-radius: 12px;
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.03);
  margin-bottom: 30px;
  position: relative;
}

/* 搜索行 - 搜索框和按钮同行 */
.search-row {
  display: flex;
  gap: 12px;
  max-width: 700px;
  margin: 0 auto;
  align-items: center;
}

/* 自定义搜索框样式 */
.custom-search-input {
  flex: 1;
  height: 52px !important;
  font-size: 16px;
  border-radius: 10px !important;
  border: 1px solid #e0e7ff !important;
  background-color: #ffffff !important;
  transition: all 0.3s ease;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.03);
}

.custom-search-input:focus-within {
  border-color: #003f88 !important;
  box-shadow: 0 0 0 3px rgba(0, 63, 136, 0.15), 0 2px 8px rgba(0, 0, 0, 0.05) !important;
  outline: none;
  border-width: 1.5px !important;
}

.custom-search-input .el-input__inner {
  height: 52px !important;
  line-height: 52px !important;
  padding: 0 20px !important;
  border-radius: 10px !important;
  border: none !important;
}

.custom-search-input .el-input__inner::placeholder {
  color: #a0aec0;
  font-style: italic;
  font-size: 15px;
}

/* 图标样式 */
.search-icon {
  color: #003f88;
  font-size: 19px;
  margin-left: 8px;
}

.clear-icon {
  color: #a0aec0;
  font-size: 17px;
  cursor: pointer;
  margin-right: 8px;
  transition: color 0.2s ease;
}

.clear-icon:hover {
  color: #003f88;
}

/* 分析按钮样式 */
.analysis-button {
  padding: 0 28px;
  height: 52px !important;
  font-size: 16px;
  font-weight: 500;
  background: linear-gradient(to right, #003f88, #00509d) !important;
  border: none !important;
  border-radius: 10px !important;
  transition: all 0.3s ease;
  box-shadow: 0 4px 12px rgba(0, 63, 136, 0.15);
  white-space: nowrap;
}

.analysis-button:hover {
  background: linear-gradient(to right, #002855, #003f88) !important;
  transform: translateY(-2px);
  box-shadow: 0 8px 20px rgba(0, 63, 136, 0.2);
}

.analysis-button:disabled {
  background: linear-gradient(to right, #c5d5e6, #d1dce8) !important;
  transform: none;
  box-shadow: none;
  cursor: not-allowed;
}

/* 加载图标动画 */
.loading-icon {
  margin-right: 8px;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  from {
    transform: rotate(0deg);
  }
  to {
    transform: rotate(360deg);
  }
}

/* 下拉框样式 - 关键修改：贴紧搜索框 */
.dropdown-container {
  position: absolute;
  top: 85%; /* 完全贴合搜索行底部 */
  left: 47%;
  transform: translateX(-50%);
  width: calc(100% - 48px);
  max-width: 652px;
  background-color: #ffffff;
  border-radius: 0 0 10px 10px !important; /* 上圆角取消，与搜索框呼应 */
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); /* 只显示下阴影 */
  max-height: 250px;
  overflow-y: auto;
  z-index: 1000;
  border: 1px solid #e0e7ff;
  border-top: none !important; /* 取消上边框，避免双重边框 */
  margin-top: -1px; /* 微调1px，完全贴合无缝隙 */
}

.dropdown-item {
  padding: 12px 20px;
  font-size: 15px;
  color: #334155;
  cursor: pointer;
  transition: all 0.2s ease;
}

.dropdown-item:hover, .dropdown-item.active {
  background-color: #f0f7ff;
  color: #003f88;
  font-weight: 500;
}

.dropdown-container::-webkit-scrollbar {
  width: 6px;
}

.dropdown-container::-webkit-scrollbar-track {
  background: #f8fafc;
  border-radius: 0 0 3px 3px;
}

.dropdown-container::-webkit-scrollbar-thumb {
  background: #cbd5e1;
  border-radius: 3px;
}

/* 结果容器 */
.result-container {
  min-height: 600px;
  background-color: #ffffff;
  border-radius: 12px;
  overflow: hidden;
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.03);
}

/* 空状态 */
.empty-state {
  height: 600px;
  display: flex;
  flex-direction: column;
  justify-content: center;
}

.empty-icon {
  font-size: 52px;
  color: #c5d5e6;
  margin-bottom: 18px;
}

.el-empty__description {
  color: #6c757d !important;
  font-size: 15px !important;
}

/* 错误状态样式 */
.error-state {
  height: 600px;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  text-align: center;
  padding: 20px;
}

.error-text {
  color: #dc2626;
  font-size: 16px;
  margin-top: 20px;
  max-width: 500px;
  line-height: 1.6;
}

/* 加载状态 */
.loading-state {
  height: 600px;
  position: relative;
}

.loading-overlay {
  height: 100%;
}

.el-loading__text {
  color: #003f88 !important;
  font-size: 16px !important;
}

/* 结果内容区域 */
.result-content {
  padding: 30px;
  text-align: center;
}

.result-title {
  color: #003f88;
  font-size: 18px;
  margin: 0 0 24px;
  font-weight: 600;
}

/* 图片容器 */
.image-container {
  height: 800px;
  display: flex;
  justify-content: center;
  align-items: center;
  background-color: #f8fafc;
  border-radius: 10px;
  overflow: hidden;
  padding: 25px;
  box-shadow: inset 0 2px 8px rgba(0, 0, 0, 0.05);
  margin-bottom: 24px;
}

.expression-image {
  max-width: 100%;
  max-height: 100%;
  transition: all 0.3s ease;
  border-radius: 8px;
  box-shadow: 0 8px 24px rgba(0, 0, 0, 0.08);
}

.expression-image:hover {
  transform: scale(1.01);
  box-shadow: 0 12px 32px rgba(0, 0, 0, 0.12);
}

/* 图片加载中状态 */
.image-loading {
  text-align: center;
  color: #003f88;
  padding: 50px 20px;
}

/* 图片错误状态 */
.image-error {
  text-align: center;
  padding: 50px 20px;
  color: #dc2626;
}

/* 下载按钮容器 - 图片正下方居中 */
.download-container {
  text-align: center;
  margin-top: 20px;
}

.download-button {
  padding: 8px 32px;
  font-size: 16px;
  background: linear-gradient(to right, #059669, #047857) !important;
  border: none !important;
  border-radius: 8px !important;
  transition: all 0.3s ease;
  box-shadow: 0 4px 12px rgba(5, 150, 105, 0.15);
}

.download-button:hover {
  background: linear-gradient(to right, #065f46, #047857) !important;
  transform: translateY(-2px);
  box-shadow: 0 8px 20px rgba(5, 150, 105, 0.2);
}
</style>