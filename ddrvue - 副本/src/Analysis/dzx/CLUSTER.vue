<template>
  <div class="cluster-container">
    <!-- 标题区域 -->
    <div class="cluster-header">
      <h2>Cluster</h2>
    </div>

    <!-- 选择区域 -->
    <div class="select-container">
      <div class="input-button-group">
        <!-- 输入框 + 下拉包裹容器 -->
        <div class="input-with-dropdown">
          <el-input
              v-model="selectedCancer"
              placeholder="Enter cancer type (e.g. KIRC)..."
              class="cancer-input"
              @focus="showDropdown = true"
              @blur="handleBlur"
              @input="handleInputDebounce"
          >
            <template #prefix>
              <el-icon class="search-icon"><Search /></el-icon>
            </template>
            <template #suffix>
              <el-icon
                  class="clear-icon"
                  v-if="selectedCancer"
                  @click="clearSelection"
              >
                <CircleClose />
              </el-icon>
            </template>
          </el-input>

          <!-- Cancer下拉建议框（紧贴输入框下方） -->
          <div
              class="dropdown-container"
              v-show="showDropdown && filteredCancers.length > 0"
          >
            <div
                class="dropdown-item"
                v-for="(cancer, index) in filteredCancers"
                :key="index"
                @click="selectCancer(cancer)"
                :class="{ active: cancer === highlightedCancer }"
                @mouseenter="highlightedCancer = cancer"
            >
              {{ cancer }}
            </div>
          </div>
        </div>

        <!-- View Cluster 按钮 -->
        <el-button
            type="primary"
            @click="loadClusterImages"
            :disabled="!selectedCancer || isLoading"
            class="confirm-button"
        >
          <el-icon v-if="isLoading" class="loading-icon"><Loading /></el-icon>
          <span v-else>View Cluster</span>
        </el-button>
      </div>
    </div>

    <!-- 结果展示区域 -->
    <div class="result-container">
      <!-- 空状态 -->
      <el-empty
          v-if="!isLoaded && !isLoading && !errorMsg"
          description="No cluster data to display"
          class="empty-state"
      >
        <template #image>
          <el-icon class="empty-icon"><Picture /></el-icon>
        </template>
        <template #description>
          <p>Please select a cancer type and click "View Cluster"</p>
        </template>
      </el-empty>

      <!-- 错误状态 -->
      <div v-if="errorMsg && !isLoading" class="error-state">
        <el-icon class="error-icon"><Warning /></el-icon>
        <p class="error-text">{{ errorMsg }}</p>
      </div>

      <!-- 加载状态 -->
      <div v-if="isLoading" class="loading-state">
        <el-loading
            v-loading="isLoading"
            text="Loading cluster images..."
            spinner="el-icon-loading"
            background="rgba(255, 255, 255, 0.9)"
            class="loading-overlay"
        >
        </el-loading>
      </div>

      <!-- 图片平铺展示区域（扩展为4张图） -->
      <div v-if="isLoaded && !isLoading && !errorMsg && clusterImages.length === 4" class="image-grid-container">
        <h3 class="images-title">
          {{ selectedCancer }} Cluster Analysis
        </h3>

        <!-- 4张图片网格布局：横向平铺，自适应响应 -->
        <div class="image-grid">
          <div
              class="image-card"
              v-for="(img, index) in clusterImages"
              :key="index"
          >
            <div class="card-image-container">
              <el-image
                  :src="img"
                  :preview-src-list="clusterImages"
                  fit="contain"
                  :alt="imageTitles[index]"
                  style="width: 100%; height: 100%;"
              >
                <template #error>
                  <div class="image-error">
                    <el-icon><Warning /></el-icon>
                    <p>{{ imageTitles[index] }} 加载失败</p>
                  </div>
                </template>
                <template #loading>
                  <div class="image-loading">
                    <el-icon class="loading-icon"><Loading /></el-icon>
                  </div>
                </template>
              </el-image>
            </div>
            <!-- 自定义标题展示 -->
            <div class="image-caption">{{ imageTitles[index] }}</div>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, watch } from 'vue';
import axios from 'axios';
import { Search, CircleClose, Loading, Picture, Warning } from '@element-plus/icons-vue';
import { ElMessage } from 'element-plus';

// 基础URL
const baseUrl = ref('http://121.37.88.191:9016');
// const baseUrl = ref('http://localhost:9016');
// 支持的癌症类型
const allCancers = ref([
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC",
  "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
  "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
  "UCEC", "UCS", "UVM"
]);

// 自定义图片标题（新增survival对应标题，扩展为4个）
const imageTitles = ref([
  "Sample Consensus Clustering Heatmap (k=2)",
  "Cluster Stability Analysis Chart",
  "ICL Optimal Cluster Number Determination Chart",
  "Prognostic Survival Analysis of Cluster Subtypes"
]);

// 响应式数据
const selectedCancer = ref('');
const filteredCancers = ref([]);
const showDropdown = ref(false);
const highlightedCancer = ref('');
const isLoading = ref(false);
const isLoaded = ref(false);
const errorMsg = ref('');
const clusterImages = ref([]); // 存储4张图片（2张consensus + 1张icl + 1张survival）
const debounceTimer = ref(null);

// 防抖输入
const handleInputDebounce = (value) => {
  if (debounceTimer.value) clearTimeout(debounceTimer.value);
  if (!value?.trim()) {
    filteredCancers.value = [];
    showDropdown.value = false;
    return;
  }
  debounceTimer.value = setTimeout(() => {
    const lower = value.toLowerCase();
    filteredCancers.value = allCancers.value.filter(c =>
        c.toLowerCase().includes(lower)
    );
    showDropdown.value = true;
  }, 300);
};

// 失焦延迟关闭下拉
const handleBlur = () => {
  setTimeout(() => {
    showDropdown.value = false;
  }, 200);
};

// 清空
const clearSelection = () => {
  selectedCancer.value = '';
  filteredCancers.value = [];
  showDropdown.value = false;
  clusterImages.value = [];
  isLoaded.value = false;
  errorMsg.value = '';
};

// 选择癌症类型
const selectCancer = (cancer) => {
  selectedCancer.value = cancer;
  filteredCancers.value = [];
  showDropdown.value = false;
  highlightedCancer.value = '';
  errorMsg.value = '';
};

// 监听图片数组变化，确认加载状态（改为4张图）
watch(
    () => clusterImages.value,
    (newVal) => {
      isLoaded.value = newVal.length === 4; // 改为4张图才标记加载完成
    }
);

// 加载图片（新增survival接口请求）
const loadClusterImages = async () => {
  if (!selectedCancer.value) {
    ElMessage.warning('Please select a cancer type first');
    return;
  }

  isLoading.value = true;
  errorMsg.value = '';
  clusterImages.value = [];

  try {
    // 并行请求三个接口（新增survival接口）
    const [consensusRes, iclRes, survivalRes] = await Promise.all([
      axios.get(`${baseUrl.value}/cluster/consensus/images`, {
        params: {cancer: selectedCancer.value},
        responseType: 'json'
      }),
      axios.get(`${baseUrl.value}/cluster/icl/images`, {
        params: {cancer: selectedCancer.value},
        responseType: 'json'
      }),
      axios.get(`${baseUrl.value}/cluster/survival/images`, { // 新增survival接口请求
        params: {cancer: selectedCancer.value},
        responseType: 'json'
      })
    ]);

    // 提取有效图片（过滤空字符串）
    const consensusImgs = Array.isArray(consensusRes.data) ? consensusRes.data.filter(img => img) : [];
    const iclImgs = Array.isArray(iclRes.data) ? iclRes.data.filter(img => img) : [];
    const survivalImgs = Array.isArray(survivalRes.data) ? survivalRes.data.filter(img => img) : []; // 新增survival图片提取

    // 合并为4张图（2张consensus + 1张icl + 1张survival）
    clusterImages.value = [...consensusImgs, ...iclImgs, ...survivalImgs];

    // 验证图片数量（确保是4张）
    if (clusterImages.value.length !== 4) {
      throw new Error(`Expected 4 images, but got ${clusterImages.value.length} images`);
    }

    ElMessage.success(`Successfully loaded 4 cluster images`);

  } catch (error) {
    console.error('load error:', error);
    errorMsg.value = error.response?.data || error.message || 'Failed to load cluster images. Please try again.';
    ElMessage.error(errorMsg.value);
  } finally {
    isLoading.value = false;
  }
};
</script>

<style scoped>
.cluster-container {
  max-width: 1800px;
  margin: 0 auto;
  padding: 30px 24px;
  background-color: #ffffff;
  border-radius: 16px;
  box-shadow: 0 10px 30px rgba(0, 0, 0, 0.08);
}

.cluster-header {
  text-align: center;
  margin-bottom: 30px;
  padding-bottom: 15px;
  border-bottom: 2px solid #e0e7ff;
}

.cluster-header h2 {
  color: #003f88;
  font-size: 26px;
  margin: 0;
  font-weight: 600;
}

.select-container {
  background-color: rgba(248, 250, 252, 0.8);
  padding: 24px;
  border-radius: 12px;
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.03);
  margin-bottom: 30px;
}

.input-button-group {
  display: flex;
  gap: 16px;
  justify-content: center;
  align-items: center;
  flex-wrap: wrap;
}

/* 包裹输入框和下拉的容器 */
.input-with-dropdown {
  position: relative;
  width: 420px;
}

.cancer-input {
  width: 100% !important;
  height: 52px !important;
  font-size: 16px;
  border-radius: 10px !important;
  border: 1px solid #e0e7ff !important;
  background-color: #ffffff !important;
}

.cancer-input .el-input__inner {
  height: 52px !important;
  line-height: 52px !important;
  padding: 0 20px !important;
}

.confirm-button {
  height: 52px !important;
  font-size: 16px;
  padding: 0 28px;
  background: linear-gradient(to right, #003f88, #00509d) !important;
  border: none !important;
  border-radius: 10px !important;
  white-space: nowrap;
}

.confirm-button:hover {
  background: linear-gradient(to right, #002855, #003f88) !important;
  transform: translateY(-2px);
  box-shadow: 0 8px 20px rgba(0, 63, 136, 0.2);
}

/* 下拉框样式 */
.dropdown-container {
  position: absolute;
  top: calc(100% + 4px);
  left: 0;
  right: 0;
  width: 100%;
  background-color: #ffffff;
  border-radius: 0 0 10px 10px;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
  max-height: 250px;
  overflow-y: auto;
  z-index: 1000;
  border: 1px solid #e0e7ff;
  border-top: none;
}

.dropdown-item {
  padding: 12px 20px;
  font-size: 15px;
  color: #334155;
  cursor: pointer;
}

.dropdown-item:hover,
.dropdown-item.active {
  background-color: #f0f7ff;
  color: #003f88;
  font-weight: 500;
}

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
}

.result-container {
  min-height: 800px; /* 调整最小高度适配4张图 */
  background-color: #ffffff;
  border-radius: 12px;
  overflow: hidden;
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.03);
}

.empty-state {
  height: 800px; /* 同步调整空状态高度 */
  display: flex;
  flex-direction: column;
  justify-content: center;
}

.empty-icon {
  font-size: 52px;
  color: #c5d5e6;
  margin-bottom: 18px;
}

.error-state {
  height: 800px; /* 同步调整错误状态高度 */
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
  margin: 20px 0 0;
  max-width: 500px;
}

.loading-state {
  height: 800px; /* 同步调整加载状态高度 */
  position: relative;
}

.loading-overlay {
  height: 100%;
}

/* 图片网格容器样式（适配4张图） */
.image-grid-container {
  padding: 30px;
}

.images-title {
  color: #003f88;
  font-size: 18px;
  text-align: center;
  margin: 0 0 24px;
  font-weight: 600;
}

/* 4张图片网格布局：优化响应式，2列/4列自适应 */
.image-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); /* 缩小最小宽度适配4张图 */
  gap: 20px;
  justify-content: center;
  align-items: center;
}

.image-card {
  display: flex;
  flex-direction: column;
  align-items: center;
  background-color: #f8fafc;
  border-radius: 12px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  overflow: hidden;
}

.card-image-container {
  width: 100%;
  height: 380px; /* 微调图片容器高度适配布局 */
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 20px;
}

.image-caption {
  padding: 12px 0;
  font-size: 15px;
  color: #003f88;
  font-weight: 500;
}

.image-loading,
.image-error {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  background-color: #f8fafc;
  color: #dc2626;
  font-size: 14px;
  gap: 8px;
}

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
</style>