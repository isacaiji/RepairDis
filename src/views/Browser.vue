<template>
  <div class="browser-container">
    <!-- 搜索表单区域 -->
    <div class="search-form-wrapper">
      <div class="search-form">
        <el-input
            v-model="searchQuery"
            placeholder="Enter gene name (e.g., TP53, BRCA1)"
            :loading="loading"
            clearable
            @input="handleInput"
            @clear="resetSearch"
            @keyup.enter="searchGenes"
            class="search-input"
        >
        </el-input>

        <!-- 实时搜索建议下拉框 - 宽度与input保持一致 -->
        <div
            v-if="showSuggestions && options.length > 0"
            class="suggestions-dropdown"
        >
          <div
              v-for="option in options"
              :key="option"
              class="suggestion-item"
              @click="selectSuggestion(option)"
          >
            <span
                v-html="highlightMatch(option)"
            ></span>
          </div>
        </div>

        <!-- 移除了搜索按钮的图标 -->
        <el-button type="primary" @click="searchGenes" class="search-btn">
          Search
        </el-button>
        <el-button
            type="danger"
            @click="resetSearch"
            v-if="isSearchActive"
            class="reset-btn"
        >
          <el-icon><Refresh /></el-icon> Reset
        </el-button>
      </div>
      <!-- 搜索示例 -->
      <div class="search-examples">
        <span>eg: </span>
        <span @click="fillSearchQuery('TP53')" class="example-item">TP53</span>
        <span> | </span>
        <span @click="fillSearchQuery('BRCA1')" class="example-item">BRCA1</span>
        <span> | </span>
        <span @click="fillSearchQuery('ATM')" class="example-item">ATM</span>
      </div>
      <el-divider class="custom-divider"></el-divider>
    </div>

    <!-- 内容区域保持不变 -->
    <div class="content-wrapper">
      <!-- 搜索无结果 -->
      <div class="empty-res" v-if="isSearchActive && searchResults.length === 0">
        <el-empty
            description="No matching DNA repair genes found"
            style="background-color: white; width: 1200px; height: 750px; margin: 0 auto;"
        />
      </div>

      <!-- 搜索结果展示 -->
      <div class="search-result-wrapper" v-if="isSearchActive && searchResults.length > 0">
        <div class="gene-table-container">
          <div class="table-wrapper">
            <el-table
                :data="searchResults"
                stripe
                border
                highlight-current-row
                :empty-text="loading ? 'Searching...' : 'No data'"
                class="gene-table"
                :header-cell-style="tableHeaderStyle"
                :row-class-name="tableRowClassName"
                :row-height="70"
            >
              <el-table-column prop="geneName" label="Gene" align="center" width="150"></el-table-column>
              <el-table-column prop="ensembl" label="Ensembl ID" align="center" width="200"></el-table-column>
              <el-table-column prop="pathway" label="DDR Pathway" align="center"></el-table-column>
              <el-table-column prop="pmid" label="PMID" align="center" width="150"></el-table-column>
              <el-table-column label="Operation" align="center" width="180">
                <template #default="{ row }">
                  <el-button
                      type="primary"
                      size="small"
                      @click="goToGenedetail(row.id)"
                      class="detail-btn"
                  >
                    <el-icon><ArrowRight /></el-icon> Details
                  </el-button>
                </template>
              </el-table-column>
            </el-table>
          </div>
        </div>
      </div>

      <!-- 基因列表（分页）展示 -->
      <div class="list-result-wrapper" v-else>
        <!-- 列表内容保持不变 -->
        <div class="table-container">
          <!-- 列表加载态 -->
          <el-loading
              v-if="isFirstPageLoading"
              target=".table-wrapper"
              text="Loading DNA repair genes..."
              background="rgba(255, 255, 255, 0.8)"
          ></el-loading>

          <!-- 基因列表表格 -->
          <div v-else class="gene-table-container">
            <div class="table-wrapper">
              <el-table
                  :data="currentPageGenes"
                  stripe
                  border
                  highlight-current-row
                  empty-text="No gene data available"
                  :loading="isFullSummaryLoading && currentPage > 1"
                  class="gene-table"
                  :header-cell-style="tableHeaderStyle"
                  :row-class-name="tableRowClassName"
                  :row-height="70"
              >
                <el-table-column prop="geneName" label="Gene" align="center" width="150"></el-table-column>
                <el-table-column prop="ensemblId" label="Ensembl ID" align="center" width="200"></el-table-column>
                <el-table-column prop="pathway" label="Repair Pathway" align="center"></el-table-column>
                <el-table-column prop="pmid" label="PMID" align="center" width="150"></el-table-column>
                <el-table-column label="Operation" align="center" width="180">
                  <template #default="{ row }">
                    <el-button
                        type="primary"
                        size="small"
                        @click="goToGenedetail(row.id)"
                        :disabled="!isFullSummaryLoaded"
                        class="detail-btn"
                    >
                      <el-icon><ArrowRight /></el-icon> Details
                    </el-button>
                  </template>
                </el-table-column>
              </el-table>
            </div>
          </div>
        </div>

        <!-- 列表分页控件 -->
        <el-pagination
            v-if="!isFirstPageLoading"
            @current-change="handleCurrentChange"
            :current-page="currentPage"
            :page-size="pageSize"
            layout="total, prev, pager, next, jumper"
            :total="total"
            :disabled="!isFullSummaryLoaded && currentPage > 1"
            class="pagination"
        ></el-pagination>
      </div>
    </div>
  </div>
</template>

<script setup>
// 脚本部分保持不变
import { ref, onMounted, computed, nextTick, watch } from 'vue';
import { useRouter } from 'vue-router';
import axios from 'axios';
import { Search, Refresh, ArrowRight } from '@element-plus/icons-vue';
import { ElMessage } from 'element-plus';

const router = useRouter();
// const apiBase = ref('http://localhost:9016/api');
const apiBase = ref('http://121.37.88.191:9016');
// 搜索相关状态
const searchQuery = ref('');
const searchResults = ref([]);
const allGenes = ref([]);
const options = ref([]);
const loading = ref(false);
const isSearchActive = ref(false);
const showSuggestions = ref(false);

// 分页相关状态
const pageSize = ref(15);
const currentPage = ref(1);
const total = ref(218); // 明确显示DNA损伤修复基因总数
const firstPageSummary = ref([]);
const fullSummaryData = ref([]);
const isFirstPageLoading = ref(true);
const isFullSummaryLoading = ref(false);
const isFullSummaryLoaded = ref(false);

// 表格样式相关
const tableHeaderStyle = {
  'background-color': '#003F88',
  'color': '#ffffff',
  'font-weight': '600',
  'font-size': '16px',
  'padding': '16px 0',
  'border-bottom': '2px solid #002855'
};

const tableRowClassName = ({ rowIndex }) => {
  return rowIndex % 2 === 0 ? 'table-row-even' : 'table-row-odd';
};

// 计算属性
const totalGeneCount = computed(() => {
  return fullSummaryData.value.length > 0 ? fullSummaryData.value.length : total.value;
});

const currentPageGenes = computed(() => {
  if (!isFullSummaryLoaded.value) return firstPageSummary.value;
  const start = (currentPage.value - 1) * pageSize.value;
  const end = Math.min(start + pageSize.value, fullSummaryData.value.length);
  return fullSummaryData.value.slice(start, end);
});

// 处理输入事件，实现实时搜索建议
const handleInput = (query) => {
  if (!query.trim()) {
    options.value = [];
    showSuggestions.value = false;
    return;
  }

  showSuggestions.value = true;
  remoteMethod(query);
};

// 选择建议项
const selectSuggestion = (option) => {
  searchQuery.value = option;
  showSuggestions.value = false;
  searchGenes();
};

// 高亮匹配的文本
const highlightMatch = (option) => {
  if (!searchQuery.value) return option;
  const lowerQuery = searchQuery.value.toLowerCase();
  const lowerOption = option.toLowerCase();
  const index = lowerOption.indexOf(lowerQuery);

  if (index === -1) return option;

  const before = option.substring(0, index);
  const match = option.substring(index, index + searchQuery.value.length);
  const after = option.substring(index + searchQuery.value.length);

  return `${before}<span class="highlight">${match}</span>${after}`;
};

// 生命周期与方法
onMounted(async () => {
  try {
    const [allGenesRes, firstPageRes] = await Promise.all([
      axios.get(`${apiBase.value}/api/genes/all`),
      axios.get(`${apiBase.value}/api/genes/summary`, {
        params: { pageNum: currentPage.value, pageSize: pageSize.value }
      })
    ]);
    allGenes.value = allGenesRes.data || [];
    if (firstPageRes.data?.list && firstPageRes.data.total !== undefined) {
      firstPageSummary.value = firstPageRes.data.list;
      total.value = firstPageRes.data.total;
    }
    isFirstPageLoading.value = false;

    await nextTick();
    isFullSummaryLoading.value = true;
    const fullSummaryRes = await axios.get(`${apiBase.value}/api/genes/summary/all`);
    fullSummaryData.value = fullSummaryRes.data || [];
    isFullSummaryLoaded.value = true;
    isFullSummaryLoading.value = false;

  } catch (error) {
    console.error('Initialization error:', error);
    ElMessage.error('Failed to load data. Please refresh the page.');
    isFirstPageLoading.value = false;
    isFullSummaryLoading.value = false;
  }
});

const searchGenes = async () => {
  if (!searchQuery.value.trim()) {
    resetSearch();
    return;
  }

  try {
    loading.value = true;
    const response = await axios.get(`${apiBase.value}/api/genes/search`, {
      params: { query: searchQuery.value.trim() }
    });
    searchResults.value = response.data || [];
    isSearchActive.value = true;
    showSuggestions.value = false; // 搜索后隐藏建议
    if (searchResults.value.length === 0) {
      ElMessage.info(`No results found for "${searchQuery.value}"`);
    }
  } catch (error) {
    console.error('Search error:', error);
    ElMessage.error('Search failed. Please try again.');
    searchResults.value = [];
    isSearchActive.value = true;
  } finally {
    loading.value = false;
  }
};

const resetSearch = () => {
  searchQuery.value = '';
  searchResults.value = [];
  isSearchActive.value = false;
  options.value = [];
  showSuggestions.value = false;
  currentPage.value = 1;
};

const remoteMethod = (query) => {
  if (query.trim() === "") {
    options.value = [];
    return;
  }

  loading.value = true;
  const lowerQuery = query.toLowerCase();
  const matched = allGenes.value
      .filter(item => item.toLowerCase().includes(lowerQuery))
      .sort((a, b) => {
        const aStarts = a.toLowerCase().startsWith(lowerQuery) ? -1 : 1;
        const bStarts = b.toLowerCase().startsWith(lowerQuery) ? -1 : 1;
        return aStarts - bStarts;
      })
      .slice(0, 10);

  options.value = matched;
  loading.value = false;
};

const fillSearchQuery = (example) => {
  searchQuery.value = example;
  searchGenes();
};

const handleCurrentChange = (newPage) => {
  if (!isFullSummaryLoaded.value && newPage !== 1) {
    ElMessage.warning('Please wait for full data loading');
    return;
  }
  currentPage.value = newPage;
  document.querySelector('.table-container')?.scrollIntoView({ behavior: 'smooth' });
};

const goToGenedetail = (id) => {
  if (!id) {
    ElMessage.error('Invalid gene ID');
    return;
  }
  router.push({ path: '/detail', query: { id } });
};
</script>

<style scoped>
/* 主要调整的样式 */
.search-form {
  display: flex;
  width: 70%;
  gap: 15px;
  align-items: center;
  margin: 0 auto;
  position: relative; /* 确保下拉建议相对于搜索框定位 */
}

/* 搜索建议下拉框 - 确保宽度与输入框一致 */
.suggestions-dropdown {
  position: absolute;
  top: 100%;
  left: 0;
  right: 0; /* 左右都设为0，与父容器等宽 */
  margin-top: 5px;
  background-color: white;
  border-radius: 8px;
  box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.1), 0 8px 10px -6px rgba(0, 0, 0, 0.1);
  z-index: 1000;
  max-height: 300px;
  overflow-y: auto;
  border: 1px solid #e2e8f0;
  width: calc(100% - 160px); /* 计算宽度，减去按钮占用的空间 */
}

/* 搜索按钮样式 - 移除图标后调整 */
.search-btn {
  background: linear-gradient(to right, #003F88, #00509D);
  border: none;
  color: white;
  height: 52px;
  padding: 0 30px;
  font-size: 16px;
  border-radius: 8px;
  transition: all 0.3s ease;
  font-weight: 500;
  min-width: 100px; /* 确保按钮宽度一致 */
}

/* 其他样式保持不变 */
.browser-container {
  max-width: 1400px;
  margin: 0 auto;
  padding: 40px;
  background-color: #f2f4f7;
  min-height: calc(100vh - 80px);
  box-shadow: 0 0 30px rgba(0, 0, 0, 0.08);
}

.search-form-wrapper {
  background-color: #ffffff;
  border-radius: 12px;
  box-shadow: 0 5px 20px rgba(0, 0, 0, 0.1);
  padding: 40px;
  margin-bottom: 40px;
  border-top: 5px solid #003F88;
  position: relative;
}

.search-input {
  flex: 1;
  height: 52px;
  border-radius: 8px;
  border: 1px solid #c9d1d9;
  background-color: #fefefe;
  box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.05);
  transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
}

.search-input .el-input__wrapper {
  height: 100%;
  border: none;
  box-shadow: none;
  background-color: transparent;
}

.search-input .el-input__input {
  height: 100%;
  font-size: 16px;
  padding: 0 18px 0 45px;
  font-family: 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  color: #2d3748;
  letter-spacing: 0.3px;
}

.search-input .el-input__input::placeholder {
  color: #a0aec0;
  font-style: italic;
}

.search-input:focus-within {
  border-color: #003F88;
  box-shadow: 0 0 0 3px rgba(0, 63, 136, 0.15), 0 1px 3px rgba(0, 0, 0, 0.1);
  outline: none;
  border-width: 1.5px;
}

.search-icon {
  color: #003F88;
  font-size: 18px;
  left: 15px;
}

.suggestion-item {
  padding: 12px 20px;
  font-size: 15px;
  color: #2d3748;
  cursor: pointer;
  transition: all 0.2s ease;
}

.suggestion-item:hover {
  background-color: #edf2f7;
  color: #003F88;
  padding-left: 24px;
}

.suggestion-item .highlight {
  color: #003F88;
  font-weight: 600;
  background-color: rgba(0, 63, 136, 0.1);
  padding: 2px 4px;
  border-radius: 3px;
}

.search-btn:hover {
  background: linear-gradient(to right, #002855, #003F88);
  transform: translateY(-2px);
  box-shadow: 0 5px 15px rgba(0, 0, 0, 0.15);
}

.reset-btn {
  height: 52px;
  padding: 0 30px;
  font-size: 16px;
  border-radius: 8px;
  transition: all 0.3s ease;
}

.reset-btn:hover {
  transform: translateY(-2px);
  box-shadow: 0 5px 15px rgba(0, 0, 0, 0.1);
}

/* 其他样式保持不变 */
.search-examples {
  margin-top: 22px;
  text-align: center;
  color: #4A5568;
  font-size: 16px;
}

.example-item {
  color: #003F88;
  font-weight: 500;
  cursor: pointer;
  padding: 0 10px;
  transition: all 0.2s;
}

.example-item:hover {
  color: #002855;
  text-decoration: underline;
}

.custom-divider {
  margin: 30px 0 0;
  background: linear-gradient(90deg, rgba(0,63,136,0) 0%, rgba(0,63,136,0.3) 50%, rgba(0,63,136,0) 100%);
  height: 1px;
}

.content-wrapper {
  background-color: #ffffff;
  border-radius: 12px;
  box-shadow: 0 5px 20px rgba(0, 0, 0, 0.1);
  overflow: hidden;
}

.result-header {
  display: flex;
  justify-content: space-between;
  padding: 22px 35px;
  background-color: #f7fafc;
  border-bottom: 1px solid #e2e8f0;
  font-size: 18px;
  color: #1A202C;
  font-weight: 600;
}

.search-content strong {
  color: #003F88;
}

.gene-table-container {
  padding: 30px 35px;
  max-width: 1200px;
  margin: 0 auto;
}

.table-wrapper {
  box-shadow: 0 3px 15px rgba(0, 0, 0, 0.07);
  border-radius: 8px;
  overflow: hidden;
  border: 1px solid #e2e8f0;
}

.gene-table {
  border-radius: 8px 8px 0 0;
  font-size: 15px;
  width: 100%;
}

.gene-table th {
  background-color: #003F88 !important;
  color: white !important;
  font-weight: 600;
  font-size: 16px;
  padding: 16px 0 !important;
  position: relative;
}

.gene-table td {
  color: #2D3748;
  padding: 20px 15px !important;
  border-bottom: 1px solid #e2e8f0;
  vertical-align: middle;
  font-size: 15px;
}

.gene-table tr {
  transition: all 0.3s ease;
}

.gene-table tr:hover > td {
  background-color: #edf2f7;
  transform: translateX(3px);
}

.table-row-even {
  background-color: #f8fafc;
}

.table-row-odd {
  background-color: #ffffff;
}

.pagination {
  margin: 35px auto;
  display: flex;
  justify-content: center;
  padding-bottom: 30px;
}

.el-pagination button, .el-pagination span:not([class*=suffix]) {
  color: #2D3748;
}

.el-pagination .btn-next, .el-pagination .btn-prev {
  background-color: #f8fafc;
  border: 1px solid #e2e8f0;
}

.el-pagination .el-pager li {
  border-radius: 4px;
  transition: all 0.2s ease;
}

.el-pagination .el-pager li:hover:not(.active) {
  color: #003F88;
  border-color: #003F88;
}

.el-pagination .el-pager li.active {
  background-color: #003F88;
  color: white;
}

.empty-res {
  padding: 100px 0;
}

.detail-btn {
  background-color: #003F88 !important;
  border-color: #003F88 !important;
  color: white !important;
  transition: all 0.3s ease;
  padding: 8px 16px !important;
  font-size: 14px !important;
}

.detail-btn:hover {
  background-color: #002855 !important;
  border-color: #002855 !important;
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
}
</style>