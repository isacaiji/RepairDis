<template>
  <div class="drug-target-container">
    <div class="content-wrapper">
      <!-- 搜索区域卡片 -->
      <el-card class="search-card">
        <!-- 搜索模式切换 -->
        <div class="search-mode-radio">
          <el-radio-group v-model="searchMode" class="radio-group" @change="handleSearchModeChange">
            <el-radio :label="'drug'" class="radio-item">
              <span class="radio-icon"><i class="el-icon-capsule"></i></span>
              <span class="radio-text">By Drug</span>
            </el-radio>
            <el-radio :label="'target'" class="radio-item">
              <span class="radio-icon"><i class="el-icon-dna"></i></span>
              <span class="radio-text">By Target Gene</span>
            </el-radio>
          </el-radio-group>
        </div>

        <!-- 带联想的搜索框 -->
        <div class="search-controls">
          <el-autocomplete
              v-model="searchQuery"
              :placeholder="searchMode === 'drug' ? 'Enter drug name' : 'Enter target gene name'"
              clearable
              size="large"
              @keyup.enter="performSearch"
              class="search-input"
              :fetch-suggestions="querySearchAsync"
              :loading="suggestionLoading"
              @select="handleSelectSuggestion"
              highlight-first-item
              popper-class="search-suggestion-popper"
          ></el-autocomplete>

          <div class="search-buttons">
            <el-button
                type="primary"
                size="default"
                @click="performSearch"
                class="search-btn"
                :loading="isLoading"
            >
              Search
            </el-button>
            <el-button
                type="default"
                size="default"
                @click="resetSearch"
                class="reset-btn"
            >
              Reset
            </el-button>
          </div>
        </div>
      </el-card>

      <!-- 结果统计 -->
      <div v-if="isSearchPerformed && searchQuery.trim() && totalItems > 0" class="result-stats">
        <p>
          <span class="stats-label">Search results:</span>
          <span class="stats-count">{{ totalItems }}</span>
          <span class="stats-desc">matches for "{{ searchQuery }}"</span>
        </p>
      </div>

      <!-- 结果表格 -->
      <el-card class="result-card">
        <div v-if="isLoading" class="loading-container">
          <el-loading-spinner size="large"></el-loading-spinner>
          <p>Loading data...</p>
        </div>

        <el-table
            :data="currentData"
            border
            stripe
            highlight-current-row
            :header-cell-style="{ 'background-color': '#f0f5ff', 'color': '#003a8c', 'font-size': '16px' }"
            :row-style="{ 'font-size': '15px' }"
            :empty-text="isSearchPerformed ? 'No matching data found' : 'Loading initial data...'"
            class="result-table"
        >
          <el-table-column
              prop="drugName"
              label="Drug Name"
              align="center"
              min-width="200"
              :cell-style="{ 'color': '#00509d', 'font-weight': '500', 'cursor': 'pointer', 'padding': '12px 0' }"
              @cell-click="handleCellClick"
          ></el-table-column>
          <el-table-column
              prop="geneName"
              label="Target Name"
              align="center"
              min-width="200"
              :cell-style="{ 'color': '#00509d', 'font-weight': '500', 'cursor': 'pointer', 'padding': '12px 0' }"
              @cell-click="handleCellClick"
          ></el-table-column>
          <el-table-column
              prop="uniprotId"
              label="Uniprot Id"
              align="center"
              min-width="200"
              :cell-style="{ 'padding': '12px 0' }"
          ></el-table-column>
          <el-table-column
              prop="action"
              label="Action"
              align="center"
              min-width="200"
              show-overflow-tooltip
              :cell-style="{ 'padding': '12px 0' }"
          ></el-table-column>
          <el-table-column
              prop="approved"
              label="Approved"
              align="center"
              min-width="200"
              :cell-style="{ 'padding': '12px 0' }"
          ></el-table-column>
        </el-table>

        <!-- 导出按钮 -->
        <div v-if="totalItems > 0" class="export-container">
          <el-button
              type="success"
              size="default"
              @click="exportCurrentData"
              class="export-btn"
          >
            Export Current Results
          </el-button>
          <el-button
              type="primary"
              size="default"
              @click="exportAllData"
              class="export-btn"
              style="margin-left: 15px"
          >
            Export All Data
          </el-button>
        </div>

        <!-- 分页 -->
        <div class="pagination-container" v-if="totalItems > 0">
          <el-pagination
              @size-change="handleSizeChange"
              @current-change="handleCurrentChange"
              :current-page="currentPage"
              :page-sizes="[10, 20, 50]"
              :page-size="pageSize"
              layout="total, sizes, prev, pager, next, jumper"
              :total="totalItems"
              style="font-size: 14px"
          ></el-pagination>
        </div>
      </el-card>

      <!-- 网络可视化 -->
      <el-card class="network-card" v-if="showNetwork">
        <div class="network-header">
          <h3>Drug-Target Interaction Network</h3>
          <el-button
              type="text"
              @click="toggleNetwork"
              class="toggle-network-btn"
          >
            {{ isNetworkVisible ? 'Hide Network' : 'Show Network' }}
          </el-button>
        </div>
        <div
            id="drug-target-network"
            class="network-container"
            :style="{ height: isNetworkVisible ? '500px' : '0' }"
        >
          <el-loading
              v-if="isNetworkLoading"
              target="#drug-target-network"
              text="Loading network data..."
              background="rgba(255, 255, 255, 0.8)"
          />
          <el-empty
              v-else-if="!networkData.nodes || networkData.nodes.length === 0"
              description="No network data available"
              class="network-empty"
          />
        </div>
      </el-card>
    </div>
  </div>
</template>

<script setup>
import { nextTick, onMounted, onUnmounted, ref } from 'vue';
import { ElMessage, ElNotification } from 'element-plus';
import axios from 'axios';
import * as echarts from 'echarts';

// 基础配置 - 替换为你的后端地址
const baseURL = 'http://121.37.88.191:9016';
// const baseURL = 'http://localhost:9016';
// 搜索相关变量
const searchMode = ref('drug');
const searchQuery = ref('');
const currentPage = ref(1);
const pageSize = ref(10);
const currentData = ref([]);
const totalItems = ref(0);
const isSearchPerformed = ref(false);
const isLoading = ref(false);
const suggestionLoading = ref(false);
const suggestionList = ref([]);

// 网络相关变量
const isNetworkVisible = ref(true);
const showNetwork = ref(true);
const isNetworkLoading = ref(false);
const networkData = ref({ nodes: [], links: [] });
let networkChart = null;

// 执行搜索
const performSearch = () => {
  isLoading.value = true;
  const params = {
    pageNum: currentPage.value,
    pageSize: pageSize.value
  };

  // 根据搜索模式构建条件
  if (searchMode.value === 'drug' && searchQuery.value.trim()) {
    params.drugName = searchQuery.value.trim();
  } else if (searchMode.value === 'target' && searchQuery.value.trim()) {
    params.geneName = searchQuery.value.trim();
  }

  axios.post(`${baseURL}/api/drug-target/query`, params)
      .then(response => {
        const result = response.data;
        currentData.value = result.list || [];
        totalItems.value = result.total || 0;
        isSearchPerformed.value = true;
        isLoading.value = false;

        // 加载网络数据
        loadNetworkData({
          name: searchMode.value === 'drug' ? params.drugName : params.geneName
        });
      })
      .catch(error => {
        console.error('Error fetching data:', error);
        isLoading.value = false;
        ElMessage.error('Failed to load data');
      });
};

// 重置搜索
const resetSearch = () => {
  searchQuery.value = '';
  currentPage.value = 1;
  performSearch();
};

// 分页处理
const handleSizeChange = (val) => {
  pageSize.value = val;
  currentPage.value = 1;
  performSearch();
};

const handleCurrentChange = (val) => {
  currentPage.value = val;
  performSearch();
  document.querySelector('.result-card')?.scrollIntoView({ behavior: 'smooth' });
};

// 单元格点击事件
const handleCellClick = (row, column) => {
  if (column.property === 'drugName') {
    searchMode.value = 'drug';
    searchQuery.value = row.drugName;
  } else if (column.property === 'geneName') {
    searchMode.value = 'target';
    searchQuery.value = row.geneName;
  }
  performSearch();
};

// 导出当前结果
const exportCurrentData = () => {
  const params = {};
  if (searchMode.value === 'drug' && searchQuery.value.trim()) {
    params.drugName = searchQuery.value.trim();
  } else if (searchMode.value === 'target' && searchQuery.value.trim()) {
    params.geneName = searchQuery.value.trim();
  }

  axios({
    url: `${baseURL}/api/drug-target/export/conditions`,
    method: 'post',
    data: params,
    responseType: 'blob'
  })
      .then(response => {
        const url = window.URL.createObjectURL(new Blob([response.data]));
        const a = document.createElement('a');
        a.href = url;
        a.download = `drug_target_interactions_${searchQuery.value.trim() || 'all'}.tsv`;
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        document.body.removeChild(a);

        ElNotification({
          title: 'Export Successful',
          message: `Exported ${totalItems.value} records`,
          type: 'success',
          duration: 3000
        });
      })
      .catch(error => {
        console.error('Error exporting data:', error);
        ElMessage.error('Failed to export data');
      });
};

// 导出所有数据
const exportAllData = () => {
  axios({
    url: `${baseURL}/api/drug-target/export/all`,
    method: 'get',
    responseType: 'blob'
  })
      .then(response => {
        const url = window.URL.createObjectURL(new Blob([response.data]));
        const a = document.createElement('a');
        a.href = url;
        a.download = 'all_drug_target_interactions.tsv';
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        document.body.removeChild(a);

        ElNotification({
          title: 'Export Successful',
          message: 'All drug-target interactions exported',
          type: 'success',
          duration: 3000
        });
      })
      .catch(error => {
        console.error('Error exporting all data:', error);
        ElMessage.error('Failed to export all data');
      });
};

// 加载网络数据
const loadNetworkData = (params) => {
  isNetworkLoading.value = true;
  axios.get(`${baseURL}/api/drug-target/network`, { params })
      .then(response => {
        networkData.value = response.data || { nodes: [], links: [] };
        isNetworkLoading.value = false;
        nextTick(() => {
          initNetworkChart();
        });
      })
      .catch(error => {
        console.error('Error loading network data:', error);
        isNetworkLoading.value = false;
        networkData.value = { nodes: [], links: [] };
        ElMessage.error('Failed to load network data');
      });
};

// 初始化网络图表
const initNetworkChart = () => {
  const container = document.getElementById('drug-target-network');
  if (!container) return;

  if (networkChart) {
    networkChart.dispose();
  }

  networkChart = echarts.init(container);

  const formattedNodes = (networkData.value.nodes || []).map(node => ({
    ...node,
    category: node.type === 'drug' ? 0 : 1
  }));

  const option = {
    tooltip: {
      trigger: 'item',
      formatter: params => {
        if (params.dataType === 'node') {
          return `${params.name} (${params.data.type === 'drug' ? 'Drug' : 'Target Gene'})`;
        } else {
          return `${params.data.source} → ${params.data.target}: ${params.data.action || 'Unknown'}`;
        }
      }
    },
    legend: {
      top: 10,
      left: 10,
      orient: 'vertical',
      backgroundColor: 'rgba(255, 255, 255, 0.8)',
      borderColor: '#ddd',
      borderWidth: 1,
      borderRadius: 4,
      padding: 10,
      itemGap: 15,
      textStyle: {
        fontSize: 12,
        color: '#333'
      },
      data: ['Drug', 'Target Gene']
    },
    series: [{
      type: 'graph',
      layout: 'force',
      force: {
        repulsion: 300,
        edgeLength: 150,
        gravity: 0.2
      },
      roam: true,
      draggable: true,
      label: {
        show: true,
        fontSize: 12,
        color: '#333'
      },
      edgeSymbol: ['none', 'arrow'],
      edgeSymbolSize: [4, 10],
      data: formattedNodes,
      links: networkData.value.links || [],
      categories: [
        {
          name: 'Drug',
          itemStyle: { color: '#1E90FF' }
        },
        {
          name: 'Target Gene',
          itemStyle: { color: '#32CD32' }
        }
      ],
      emphasis: {
        focus: 'adjacency',
        lineStyle: { width: 5 }
      }
    }]
  };

  networkChart.setOption(option);

  const handleResize = () => {
    networkChart && networkChart.resize();
  };
  window.addEventListener('resize', handleResize);

  onUnmounted(() => {
    window.removeEventListener('resize', handleResize);
    if (networkChart) {
      networkChart.dispose();
      networkChart = null;
    }
  });
};

// 切换网络显示状态
const toggleNetwork = () => {
  isNetworkVisible.value = !isNetworkVisible.value;
  nextTick(() => {
    networkChart && networkChart.resize();
  });
};

// 搜索模式切换
const handleSearchModeChange = () => {
  searchQuery.value = '';
  suggestionList.value = [];
};

// 异步联想搜索（核心：仅传递keyword，无limit）
const querySearchAsync = (queryString, callback) => {
  if (!queryString.trim()) {
    callback([]);
    return;
  }

  suggestionLoading.value = true;

  const apiUrl = searchMode.value === 'drug'
      ? `${baseURL}/api/drug-target/suggest/drug`
      : `${baseURL}/api/drug-target/suggest/gene`;

  axios.get(apiUrl, {
    params: {
      // 仅传递关键词，无limit参数（后端已去掉limit）
      keyword: queryString.trim()
    }
  })
      .then(response => {
        const results = response.data || [];
        const formattedResults = results.map(item => ({ value: item }));
        suggestionList.value = formattedResults;
        callback(formattedResults);
        suggestionLoading.value = false;
      })
      .catch(error => {
        console.error('Error fetching suggestions:', error);
        callback([]);
        suggestionLoading.value = false;
        ElMessage.warning('Failed to load');
      });
};

// 选择联想建议
const handleSelectSuggestion = (item) => {
  searchQuery.value = item.value;
  performSearch();
};

// 初始化
onMounted(() => {
  currentPage.value = 1;
  performSearch();
});
</script>

<style scoped>
.drug-target-container {
  margin: 0;
  padding: 20px;
  background-color: #f9f9f9;
  min-height: 100vh;
}

.content-wrapper {
  max-width: 1800px;
  margin: 0 auto;
}

.search-card {
  margin-bottom: 25px;
  border-radius: 10px;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
  border: none;
  padding: 25px;
  background-color: #fff;
}

.search-mode-radio {
  margin-bottom: 25px;
  padding: 18px;
  background-color: #f0f5ff;
  border-radius: 10px;
}

.radio-group {
  display: flex;
  gap: 35px;
}

.radio-item {
  display: flex;
  align-items: center;
  cursor: pointer;
  padding: 10px 18px;
  border-radius: 25px;
  transition: all 0.2s;
  font-size: 16px;
}

.radio-item:hover {
  background-color: rgba(255, 255, 255, 0.7);
}

.radio-icon {
  color: #00509d;
  font-size: 19px;
  margin-right: 8px;
}

.radio-text {
  font-size: 17px;
  color: #333;
}

.search-controls {
  display: flex;
  gap: 18px;
  align-items: center;
  flex-wrap: nowrap;
}

.search-input {
  flex: 1;
  min-width: 400px;
}

.search-suggestion-popper {
  --el-autocomplete-input-height: 40px !important;
  font-size: 15px;
}

.search-suggestion-popper .el-popper__inner {
  padding: 5px 0;
  border-radius: 8px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
}

.search-suggestion-popper .el-autocomplete-suggestion__list {
  max-height: 300px;
}

.search-suggestion-popper .el-autocomplete-suggestion__item {
  padding: 10px 15px;
  cursor: pointer;
}

.search-suggestion-popper .el-autocomplete-suggestion__item:hover {
  background-color: #f0f5ff;
}

.search-suggestion-popper .el-autocomplete-suggestion__item--highlighted {
  background-color: #e6f0ff;
  color: #00509d;
}

.search-buttons {
  display: flex;
  gap: 12px;
}

.search-btn, .reset-btn, .export-btn {
  padding: 0 18px;
  height: 40px;
  font-size: 14px;
  border-radius: 6px;
  transition: all 0.2s;
}

.search-btn {
  background: linear-gradient(to right, #003f88, #00509d);
  border: none;
  color: #fff;
}

.search-btn:hover {
  background: linear-gradient(to right, #001d3d, #002855);
  transform: translateY(-2px);
}

.reset-btn {
  background-color: #fff;
  color: #003f88;
  border: 1px solid #00509d;
}

.reset-btn:hover {
  background-color: #f0f5ff;
  transform: translateY(-2px);
}

.result-stats {
  margin: 0 0 15px 0;
  color: #666;
  font-size: 16px;
  text-align: left;
}

.stats-label {
  font-weight: 500;
  color: #333;
}

.stats-count {
  color: #00509d;
  font-weight: 600;
  font-size: 17px;
  margin: 0 5px;
}

.result-card {
  border-radius: 10px;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
  border: none;
  padding: 25px;
  margin-bottom: 25px;
  background-color: #fff;
}

.loading-container {
  text-align: center;
  padding: 70px 0;
  color: #666;
}

.loading-container p {
  margin-top: 18px;
  font-size: 16px;
}

.result-table {
  margin-top: 8px;
  margin-bottom: 18px;
}

.el-table th {
  height: 48px;
  font-weight: 600;
}

.el-table td {
  height: 58px;
}

.export-container {
  text-align: right;
  margin: 0 0 18px 0;
  padding-top: 5px;
  border-top: 1px solid #f0f0f0;
}

.export-btn:first-child {
  background: linear-gradient(to right, #0077b6, #0096c7);
  border: none;
  color: #fff;
}

.export-btn:first-child:hover {
  background: linear-gradient(to right, #005f8c, #007aa3);
  transform: translateY(-2px);
}

.export-btn:last-child {
  background: linear-gradient(to right, #003f88, #00509d);
  border: none;
  color: #fff;
}

.export-btn:last-child:hover {
  background: linear-gradient(to right, #001d3d, #002855);
  transform: translateY(-2px);
}

.pagination-container {
  display: flex;
  justify-content: center;
  margin-top: 15px;
  padding: 8px 0;
}

.network-card {
  margin-bottom: 25px;
  border-radius: 10px;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
  border: none;
  background-color: #fff;
  overflow: hidden;
  transition: all 0.3s ease;
}

.network-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 14px 20px;
  border-bottom: 1px solid #f0f0f0;
}

.network-header h3 {
  margin: 0;
  color: #003f88;
  font-size: 18px;
  font-weight: 600;
}

.toggle-network-btn {
  color: #00509d;
  padding: 5px 10px;
  font-size: 14px;
}

.network-container {
  width: 100%;
  transition: height 0.3s ease;
  overflow: hidden;
}

.network-empty {
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
}
</style>