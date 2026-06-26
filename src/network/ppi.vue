<template>
  <div class="bg-gray-50 min-h-screen flex flex-col" style="height: 100%; background-color: #f8f9fb;">
    <!-- 1. 顶部控制区 -->
    <div class="mb-4 flex-shrink-0" style="width: 90%; margin: 0 auto; padding-top: 30px;">
      <div class="header" style="margin-bottom: 10px; text-align: left;">
        <h1 style="color: #003f88; margin-bottom: 5px; font-size: 24px; font-weight: 600;">Protein-Protein Interaction Network</h1>
      </div>
      <el-divider></el-divider>
      <div class="sl-hunt">
        <p><strong>Protein-Protein Interaction Network:</strong> Protein-Protein Interaction Network is a graph-theoretic construct encoding direct physical associations between proteins, with nodes as molecular entities and edges as empirically validated binding events.</p>
        <p><strong>Data Source:</strong> The data is derived from literature collected in the <a href="https://cn.string-db.org/" target="_blank" style="color: #00509d;">STRING</a> database.</p>
      </div>

      <!-- 搜索和控制区域卡片 -->
      <el-card class="search-card" style="margin-top: 20px; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px;">
        <div class="flex flex-col md:flex-row items-start md:items-center gap-6 p-2 w-full" style="width: 100%; margin: 0 auto">
          <!-- 左侧：蛋白质输入框 + 高亮按钮 -->
          <div style="margin: 0 auto; display: flex; align-items: center; gap: 10px; width: 100%">
            <span style="color: #333; font-weight: 500;">Gene:</span>
            <el-input
                v-model="targetProtein"
                placeholder="Enter gene symbol"
                size="default"
                clearable
                @clear="handleInputClear"
                style="height: 40px; flex: 1;"
            />
            <el-button
                type="primary"
                @click="handleProteinHighlight"
                :disabled="!hasData || !targetProtein.trim()"
                size="default"
                style="height: 40px; padding: 0 20px; background: linear-gradient(to right, #003f88, #00509d); border: none; font-size: 14px;"
            >
              Input
            </el-button>
          </div>

          <!-- 中间：Slider + 输入框组 - 调整滑块长度和按钮大小 -->
          <div style="width: 100%; white-space: normal; padding: 15px; border: 1px solid #eee; margin-top: 10px; border-radius: 8px; background-color: #f9f9f9;">
            <!-- 滑块和输入框 - 优化布局 -->
            <div style="display: flex; flex-wrap: wrap; align-items: center; gap: 20px; width: 100%;">
              <!-- 左侧滑块区域 - 缩短滑块长度 -->
              <div style="flex: 1; min-width: 300px;">
                <div style="display: flex; align-items: center; gap: 15px;">
                  <span style="font-weight: 500; color: #333; white-space: nowrap;">Degree:</span>
                  <el-slider
                      v-model.number="currentSliderValue"
                      :min="ppiData.minCount"
                      :max="ppiData.maxCount"
                      :disabled="!hasData"
                      style="flex: 1; max-width: 500px; margin: 0;"
                  />
                  <el-input-number
                      v-model.number="currentSliderValue"
                      :min="ppiData.minCount"
                      :max="ppiData.maxCount"
                      :disabled="!hasData"
                      @change="handleGoClick"
                      controls-position="right"
                      style="width: 100px;"
                  />
                </div>
                <div class="flex justify-between w-full text-xs text-gray-500" style="width: 100%; max-width: 550px; margin: 5px 0 0 45px;">  <!-- 调整标签位置 -->
                  <span>min: {{ ppiData.minCount }}</span>
                  <span>max: {{ ppiData.maxCount }}</span>
                </div>
              </div>

              <!-- 右侧按钮组 - 增大按钮尺寸 -->
              <div style="display: flex; gap: 15px; margin-left: auto; white-space: nowrap;">
                <el-button
                    type="primary"
                    @click="handleGoClick"
                    :disabled="!hasData || currentSliderValue === appliedSliderValue"
                    size="default"
                    style="padding: 0 20px; height: 45px; background: linear-gradient(to right, #003f88, #00509d); border: none; font-size: 15px;"
                >
                  <el-icon style="margin-right: 5px;"><Position /></el-icon> Go
                </el-button>
                <el-button
                    type="default"
                    @click="handleResetClick"
                    :disabled="!hasData || currentSliderValue === ppiData.minCount"
                    size="default"
                    style="padding: 0 20px; height: 45px; color: #003f88; border-color: #00509d; font-size: 15px;"
                >
                  <el-icon style="margin-right: 5px;"><RefreshRight /></el-icon> Reset
                </el-button>
                <el-button
                    type="success"
                    @click="handleExportTxt"
                    :disabled="!hasData"
                    :loading="isExporting"
                    size="default"
                    style="padding: 0 20px; height: 45px; background: linear-gradient(to right, #0077b6, #0096c7); border-color: #0077b6; font-size: 15px;"
                >
                  <ElIcon style="margin-right: 5px;"><Download /></ElIcon> Download
                </el-button>
              </div>
            </div>
          </div>
        </div>
      </el-card>
    </div>

    <el-divider></el-divider>
    <!-- 2. 下方视图区 -->
    <div class="flex-1 flex flex-col gap-4 p-0 md:p-4 overflow-hidden" style="margin-bottom: 10px">
      <!-- 核心图表容器 - 包含图例 -->
      <el-card class="result-card" style="border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px; padding: 20px;">
        <div
            id="chart-container"
            class="w-full border rounded overflow-hidden flex-1"
            :style="{ minHeight: '500px' }"
        >
          <el-loading
              v-if="isLoadingData"
              target="#chart-container"
              text="Loading data..."
              background="rgba(255, 255, 255, 0.8)"
          />
          <el-alert
              v-else-if="loadError"
              type="error"
              :description="`${loadError}`"
              show-icon
              class="m-4"
              dangerously-use-html-string
          />
          <el-empty
              v-else-if="!hasData"
              description="No data available"
              class="h-full flex items-center justify-center"
          />
          <el-empty
              v-else-if="!filteredNodes.length"
              description="No data after filtering"
              class="h-full flex items-center justify-center"
          >
            <template #image><ElIcon class="text-5xl text-gray-300"><Filter /></ElIcon></template>
            <template #footer><el-button type="primary" @click="handleResetClick" size="small">Reset Filter</el-button></template>
          </el-empty>
        </div>
      </el-card>

      <!-- 表格区域（带分页） -->
      <el-card class="result-card" style="border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px; padding: 20px; background-color: #fff;">
        <div
            class="mb-2 font-bold text-[#00509d] flex items-center justify-between"
            v-if="selectedNode"
            style="margin-bottom: 10px"
        >
          <span>Selected: {{ selectedNode }}</span>
          <el-button
              type="default"
              size="small"
              @click="clearSelection"
              style="background-color: #f0f5ff; color: #00509d; border-color: #00509d; padding: 0 15px;"
          >
            Cancel
          </el-button>
        </div>

        <!-- 分页表格 -->
        <div>
          <el-table
              ref="unifiedTable"
              :data="tableData"
              v-if="hasData"
              :key="tableRefreshKey"
              border
              stripe
              max-height="400px"
              style="width: 100%;"
              :header-cell-style="{ 'background-color': '#f0f5ff', 'color': '#003a8c', 'font-size': '14px' }"
              :row-style="{ 'font-size': '14px' }"
              v-loading="isLoadingTable"
              element-loading-text="Loading table data..."
          >
            <el-table-column prop="node1" label="Node1" align="center" />
            <el-table-column prop="node2" label="Node2" align="center" />
            <el-table-column prop="combined_score" label="Combined Score" align="center" />
            <el-table-column prop="node1_string_id" label="Node1 STRING ID" align="center" />
            <el-table-column prop="node2_string_id" label="Node2 STRING ID" align="center" />
          </el-table>

          <!-- 分页控件 -->
          <div class="pagination-container" v-if="hasData && total > 0">
            <el-pagination
                @size-change="handlePageSizeChange"
                @current-change="handleCurrentPageChange"
                :current-page="currentPage"
                :page-sizes="[10, 20, 50]"
                :page-size="pageSize"
                layout="total, sizes, prev, pager, next, jumper"
                :total="total"
                style="font-size: 14px"
                background
            />
          </div>

          <el-empty
              v-if="!hasData ||
                    (!selectedNode && total === 0) ||
                    (selectedNode && tableData.length === 0)"
              class="py-4"
          />
        </div>
      </el-card>
    </div>
  </div>
</template>

<script setup>
import { ref, computed, onMounted, onUnmounted, nextTick } from 'vue';
import * as echarts from 'echarts';
import { ElMessage, ElIcon } from 'element-plus';
import { Filter, Download, Position, RefreshRight } from '@element-plus/icons-vue';

// 状态管理
const targetProtein = ref('');
const isLoadingData = ref(false);
const isLoadingTable = ref(false);
const loadError = ref('');
const isExporting = ref(false);
const selectedNode = ref('');
const currentSliderValue = ref(0);
const appliedSliderValue = ref(0);
const tableRefreshKey = ref(true);

// 分页相关状态
const currentPage = ref(1);    // 当前页码
const pageSize = ref(10);      // 每页条数
const total = ref(0);          // 总数据条数
const tableData = ref([]);     // 分页表格数据

// 连接数-颜色映射表（用于图表图例）
const colorScale = ref([
  { min: 1, max: 5, color: 'hsla(210, 80%, 70%, 0.9)', name: '1-5' },
  { min: 6, max: 15, color: 'hsla(210, 80%, 55%, 0.9)', name: '6-15' },
  { min: 16, max: 30, color: 'hsla(240, 80%, 60%, 0.9)', name: '16-30' },
  { min: 31, max: 50, color: 'hsla(300, 80%, 60%, 0.9)', name: '31-50' },
  { min: 51, max: 62, color: 'hsla(340, 80%, 55%, 0.9)', name: '51-62' }
]);

// 数据存储
const filteredNodes = ref([]);
const filteredLinks = ref([]);
const ppiData = ref({
  nodes: [],
  links: [],
  node1Counts: {},
  maxCount: 62,
  minCount: 0,
  originalRecords: []
});

// 计算属性：判断数据是否存在
const hasData = computed(() =>
    ppiData.value.nodes.length > 0 &&
    ppiData.value.links.length > 0 &&
    total.value > 0
);

let chartInstance = null;
let resizeDebounceTimer = null;
let highlightedNodeIndex = -1;

// 通用节点选择处理函数
const selectNode = (nodeName) => {
  if (!hasData.value || !chartInstance) return false;

  // 查找节点
  const matchedNode = filteredNodes.value.find(node => node.name === nodeName);
  if (!matchedNode) {
    clearSelection();
    ElMessage.warning(`${nodeName} Not Found`);
    return false;
  }

  // 取消之前的高亮
  if (highlightedNodeIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'downplay',
      dataIndex: highlightedNodeIndex
    });
  }

  // 高亮当前节点并记录索引
  highlightedNodeIndex = filteredNodes.value.findIndex(node => node.name === nodeName);
  chartInstance.dispatchAction({
    type: 'highlight',
    dataIndex: highlightedNodeIndex
  });

  // 更新选中状态并重新加载分页数据
  selectedNode.value = nodeName;
  currentPage.value = 1;  // 重置为第一页
  loadPPITable(appliedSliderValue.value);

  // 更新输入框
  targetProtein.value = nodeName;

  // 强制表格刷新
  tableRefreshKey.value = !tableRefreshKey.value;

  return true;
};

// 按度数计算颜色
const getNodeColorByDegree = (degree) => {
  if (degree === 0) return 'hsla(0, 0%, 80%, 0.9)';

  for (const scale of colorScale.value) {
    if (degree >= scale.min && degree <= scale.max) {
      return scale.color;
    }
  }
  return 'hsla(0, 0%, 60%, 0.9)';
};

// 清除节点选择
const clearSelection = () => {
  if (chartInstance && highlightedNodeIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'downplay',
      dataIndex: highlightedNodeIndex
    });
    highlightedNodeIndex = -1;
  }

  // 重置选中状态并重新加载分页数据
  selectedNode.value = '';
  currentPage.value = 1;  // 重置为第一页
  loadPPITable(appliedSliderValue.value);

  targetProtein.value = '';
  tableRefreshKey.value++;
};

// 输入框清空处理
const handleInputClear = () => {
  if (selectedNode.value) {
    clearSelection();
  }
};

// 蛋白质高亮处理
const handleProteinHighlight = async () => {
  if (!hasData.value || !targetProtein.value.trim() || !chartInstance) {
    return;
  }
  selectNode(targetProtein.value.trim());
};

// 计算节点度数
const calculateNodeDegrees = (links) => {
  const degrees = {};
  if (!links || links.length === 0) return degrees;

  links.forEach(link => {
    const sourceName = typeof link.source === 'string' ? link.source : link.source.name;
    const targetName = typeof link.target === 'string' ? link.target : link.target.name;
    degrees[sourceName] = 0;
    degrees[targetName] = 0;
  });

  links.forEach(link => {
    const sourceName = typeof link.source === 'string' ? link.source : link.source.name;
    const targetName = typeof link.target === 'string' ? link.target : link.target.name;
    degrees[sourceName]++;
    degrees[targetName]++;
  });

  return degrees;
};
const apiBase = 'http://121.37.88.191:9016';
// 加载图表数据
const loadPPIData = async (degreeThreshold) => {
  isLoadingData.value = true;
  loadError.value = '';

  try {
    const response = await fetch(`${apiBase}/api/network/ppi?degreeThreshold=${degreeThreshold}`);
    if (!response.ok) throw new Error(`Failed to fetch graph data (${response.status})`);

    const data = await response.json();
    if (data.error) throw new Error(data.error);

    const currentNodeCounts = calculateNodeDegrees(data.links || []);
    const formattedNodes = data.nodes.map(node => {
      const nodeDegree = currentNodeCounts[node.name] || 0;
      return {
        ...node,
        symbolSize: 55,
        itemStyle: { color: getNodeColorByDegree(nodeDegree), borderWidth: 2, borderColor: '#fff' },
        label: { color: '#fff', fontSize: 12, fontWeight: 'bold', overflow: 'truncate', width: 80 },
        tooltip: { formatter: `<b>${node.name}</b><br/>Degree: ${nodeDegree}` },
        __debug: { degree: nodeDegree, color: getNodeColorByDegree(nodeDegree) }
      };
    });

    ppiData.value = {
      ...data,
      nodes: formattedNodes,
      node1Counts: currentNodeCounts,
      minCount: data.minCount ?? 1,
      maxCount: data.maxCount ?? Math.max(...Object.values(currentNodeCounts), 62)
    };

    filteredNodes.value = formattedNodes;
    filteredLinks.value = data.links || [];

    await nextTick();
    return true;
  } catch (error) {
    loadError.value = `Fetch error: ${error.message}`;
    ElMessage.error(loadError.value);
    return false;
  } finally {
    isLoadingData.value = false;
  }
};

// 加载分页表格数据
const loadPPITable = async (degreeThreshold) => {
  isLoadingTable.value = true;
  loadError.value = '';

  try {
    // 构建分页请求参数
    const params = new URLSearchParams({
      degreeThreshold,
      pageNum: currentPage.value,
      pageSize: pageSize.value,
      ...(selectedNode.value && { selectedNode: selectedNode.value })
    });

    const response = await fetch(`${apiBase}/api/network/ppi/table?${params}`);
    if (!response.ok) throw new Error(`Failed to fetch table data (${response.status})`);

    const res = await response.json();
    if (res.error) throw new Error(res.error);

    // 更新分页数据
    total.value = res.total || 0;
    tableData.value = res.list || [];
    return true;
  } catch (error) {
    loadError.value = `Fetch error: ${error.message}`;
    ElMessage.error(loadError.value);
    return false;
  } finally {
    isLoadingTable.value = false;
  }
};

// 组合加载数据
const loadAllData = async (degreeThreshold) => {
  // 重置分页状态
  highlightedNodeIndex = -1;
  selectedNode.value = '';
  currentPage.value = 1;

  try {
    const chartLoaded = await loadPPIData(degreeThreshold);
    if (!chartLoaded) return;

    const tableLoaded = await loadPPITable(degreeThreshold);
    if (!tableLoaded) return;

    await nextTick();
    initEcharts();
    tableRefreshKey.value++;
    appliedSliderValue.value = degreeThreshold;
  } catch (error) {
    console.error('Data loading failed:', error);
  }
};

// 初始化图表 - 包含图例配置
const initEcharts = () => {
  if (chartInstance) {
    chartInstance.dispose();
  }
  chartInstance = echarts.init(document.getElementById('chart-container'));

  const chartOption = {
    animationDuration: 1800,
    animationEasingUpdate: 'quinticInOut',
    // 图例配置 - 整合到网络图表中
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
      // 使用颜色映射表作为图例数据
      data: colorScale.value.map(item => ({
        name: `Degree: ${item.name}`,
        icon: 'circle',
        textStyle: {color: '#333'}
      }))
    },
    tooltip: {
      trigger: 'item',
      formatter: (params) => {
        if (params.dataType === 'node') {
          const nodeDegree = params.data.__debug.degree || 1;
          return `<b>${params.name}</b><br/>Degree: ${nodeDegree}`;
        }
        return params.name;
      }
    },
    series: [{
      type: 'graph',
      layout: 'force',
      force: {repulsion: 800, edgeLength: 220, gravity: 0.1, friction: 0.6},
      roam: true,
      draggable: true,
      label: {show: true},
      edgeSymbol: ['none', 'none'],
      edgeLabel: {show: false},
      data: filteredNodes.value,
      links: filteredLinks.value,
      // 使用颜色映射表作为类别
      categories: colorScale.value.map(item => ({
        name: `Degree: ${item.name}`,
        itemStyle: {color: item.color}
      })),
      lineStyle: {
        width: 2,
        color: 'rgba(0, 0, 0, 0.2)',
        curveness: 0.1
      },
      emphasis: {
        focus: 'adjacency',
        lineStyle: {width: 4, color: '#165DFF'},
        itemStyle: {shadowBlur: 15, shadowColor: 'rgba(22, 93, 255, 0.8)'},
        label: {fontSize: 14, color: '#fff'}
      },
      inactive: {
        itemStyle: {opacity: 0.2},
        linkStyle: {opacity: 0.2},
        label: {opacity: 0.2}
      }
    }]
  };

  chartInstance.setOption(chartOption);

  // 节点点击事件
  chartInstance.on('click', (params) => {
    if (params.dataType === 'node') {
      selectNode(params.name);
    }
  });

  // 空白处点击取消选择
  chartInstance.getZr().on('click', (e) => {
    if (selectedNode.value && !e.target) {
      clearSelection();
    }
  });

  window.addEventListener('resize', () => {
    if (chartInstance) {
      chartInstance.resize();
    }
  });
};

// 分页交互 - 页码改变
const handleCurrentPageChange = (val) => {
  currentPage.value = val;
  loadPPITable(appliedSliderValue.value);
};

// 分页交互 - 每页条数改变
const handlePageSizeChange = (val) => {
  pageSize.value = val;
  currentPage.value = 1; // 条数改变时重置为第一页
  loadPPITable(appliedSliderValue.value);
};

// 筛选交互
const handleGoClick = () => {
  if (hasData.value && currentSliderValue.value !== appliedSliderValue.value) {
    loadAllData(currentSliderValue.value);
  }
};

const handleResetClick = () => {
  if (hasData.value && currentSliderValue.value !== ppiData.value.minCount) {
    currentSliderValue.value = ppiData.value.minCount;
    loadAllData(currentSliderValue.value);
  }
};

// 导出文件
const handleExportTxt = async () => {
  if (!hasData.value) return;

  isExporting.value = true;
  // 如果输入框有值，则export-current
  if (targetProtein.value) {
    await fetch(`${apiBase}/api/network/ppi/export-current?name=${targetProtein.value}`, {
      method: 'GET',
      headers: {'Accept': 'text/tab-separated-values'}
    }).then(response => {
      if (!response.ok) throw new Error(`Export failed (${response.status})`);
      return response.blob();
    }).then(blob => {
      const timeStr = new Date().toISOString().replace(/[-:\.T]/g, '').slice(0, 14);
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `ppi_current_${selectedNode.value}_${timeStr}.tsv`;
      a.click();
      URL.revokeObjectURL(url);
      ElMessage.success('File exported successfully');
    }).catch(error => {
      loadError.value = error.message;
    }).finally(() => {
      isExporting.value = false;
    })
  }
  // 否则全量导出
  else {
    await fetch(`${apiBase}/api/network/ppi/export`, {
      method: 'GET',
      headers: {'Accept': 'text/tab-separated-values'}
    }).then(response => {
      if (!response.ok) throw new Error(`Export failed (${response.status})`);

      return response.blob();
    }).then(blob => {
      const timeStr = new Date().toISOString().replace(/[-:\.T]/g, '').slice(0, 14);
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `ppi_all_${timeStr}.tsv`;
      a.click();
      URL.revokeObjectURL(url);
      ElMessage.success('File exported successfully');
    }).catch(error => {
      loadError.value = error.message;
    }).finally(() => {
      isExporting.value = false;
    })
  }
};
// 窗口 resize 防抖处理
const handleWindowResize = () => {
  clearTimeout(resizeDebounceTimer);
  resizeDebounceTimer = setTimeout(() => {
    if (chartInstance) {
      chartInstance.resize();
    }
  }, 100);
};

// 生命周期
onMounted(() => {
  loadAllData(currentSliderValue.value);
  window.addEventListener('resize', handleWindowResize);
});

onUnmounted(() => {
  window.removeEventListener('resize', handleWindowResize);
  if (chartInstance) {
    chartInstance.dispose();
  }
  clearSelection();
  clearTimeout(resizeDebounceTimer);
  Object.assign(ppiData.value, {nodes: [], links: [], node1Counts: {}, maxCount: 0, minCount: 0, originalRecords: []});
  [filteredNodes.value, filteredLinks.value, tableData.value] = [[], [], []];
  total.value = 0;
});
</script>

<style scoped>
#chart-container {
  width: 100%;
  height: 100%;
}

/* 提示区 */
.sl-hunt {
  font-size: 16px;
  text-align: left;
  font-weight: 500;
  line-height: 1.8;
  color: #2c3e50;
  background: #eaf3fc;
  padding: 18px 22px;
  border-left: 4px solid #2a5eaa;
  border-radius: 8px;
  margin: -10px 0 10px;
}

.sl-hunt a {
  color: #1150bd;
  text-decoration: none;
}

.sl-hunt a:hover {
  text-decoration: underline;
}

.el-table .cell {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

/* 分页控件样式优化 */
::v-deep .el-pagination {
  margin-top: 15px;
}

::v-deep .el-pagination__total {
  color: #666;
}

::v-deep .el-pager li {
  margin: 0 5px;
}

::v-deep .el-pager li.active {
  background-color: #1150bd;
  color: white;
}

/* 自定义样式类 */
.search-card {
  padding: 20px;
  background-color: #fff;
}

.result-card {
  margin-bottom: 30px;
}

.pagination-container {
  display: flex;
  justify-content: center;
  margin-top: 25px;
  text-align: center;
  padding: 10px 0;
}

/* 按钮悬停效果 */
::v-deep .el-button--primary {
  transition: all 0.2s;
}

::v-deep .el-button--primary:hover {
  transform: translateY(-2px);
}

::v-deep .el-button--default {
  transition: all 0.2s;
}

::v-deep .el-button--default:hover {
  background-color: #f0f5ff;
  transform: translateY(-2px);
}
</style>