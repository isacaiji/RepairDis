<template>
  <div class="bg-gray-50 min-h-screen flex flex-col" style="height: 100%; background-color: #f8f9fb;">
    <div class="mb-4 flex-shrink-0" style="width: 90%; margin: 0 auto; padding-top: 30px;">
      <div class="header" style="margin-bottom: 10px; text-align: left;">
        <h1 style="color: #003f88; margin-bottom: 5px; font-size: 24px; font-weight: 600;">Synthetic-Lethality Gene Interaction</h1>
      </div>
      <el-divider></el-divider>
      <div class="sl-hunt">
        <p><strong>Synthetic Lethality:</strong> Synthetic lethality describes a genetic interaction
          where concurrent loss of two genes induces cell death, whereas inactivation of either alone
          does not—an essential concept in cancer therapy research.</p>
        <p><strong>Data Source:</strong> The data is derived from literature collected in the
          <a href="https://pubmed.ncbi.nlm.nih.gov/" target="_blank" style="color: #00509d"> PubMed </a>
          database and the
          <a href="https://spidrweb.org/app/spidrweb/" target="_blank" style="color: #00509d"> SPIDRweb </a>
          database.</p>
      </div>

      <el-card class="search-card" style="margin-top: 20px; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px;">
        <div class="flex flex-col md:flex-row items-start md:items-center gap-6 p-2 w-full" style="width: 100%; margin: 0 auto">
          <div style="margin: 0 auto; display: flex; align-items: center; gap: 10px; width: 100%; flex-wrap: wrap;">
            <span style="color: #333; font-weight: 500; white-space: nowrap;"> Gene:</span>
            <el-input
                v-model="targetProtein"
                placeholder="Enter gene symbol (e.g., BRCA1, ATM)"
                style="height: 40px; flex: 1; min-width: 250px;"
                clearable
                @clear="handleInputClear"
            />
            <el-button
                type="primary"
                @click="handleProteinHighlight"
                :disabled="!hasOriginalData || !targetProtein.trim()"
                style="height: 40px; padding: 0 20px; background: linear-gradient(to right, #003f88, #00509d); border: none; font-size: 14px; white-space: nowrap;"
            >
              Locate Gene
            </el-button>
          </div>
        </div>

        <div style="width: 100%; white-space: nowrap; padding: 10px; border: 1px solid #eee; margin-top: 10px;">
          <div style="display: inline-block; vertical-align: middle; width: 30%; white-space: normal;">
            <div class="flex justify-start gap-4 flex-wrap">
              <div v-for="(item, index) in colorScale" :key="index" style="display: flex; align-items: center; gap: 10px;">
                <div :style="{width: '16px', height: '16px', borderRadius: '4px', backgroundColor: item.color}"></div>
                <span class="text-sm text-gray-700">
                  {{ item.min }}-{{ item.max }}
                </span>
              </div>
            </div>
          </div>

          <div style="display: inline-block; vertical-align: middle; width: 40%; text-align: center; white-space: normal;">
            <div style="display: inline-flex; align-items: center; gap: 15px;">
              <span>Sensitivity:</span>
              <el-slider
                  v-model.number="currentSliderValue"
                  :min="originalPpiData.minSensitivity"
                  :max="Math.min(originalPpiData.maxSensitivity, -1.0)"
                  :disabled="!hasData"
                  style="width: 200px; margin: 0;"
                  :format-tooltip="(value) => `${value.toFixed(1)}`"
                  step="0.1"
              />
              <el-input-number
                  v-model.number="currentSliderValue"
                  :min="originalPpiData.minSensitivity"
                  :max="Math.min(originalPpiData.maxSensitivity, -1.0)"
                  :disabled="!hasData"
                  @change="handleSliderValueChange"
                  controls-position="right"
                  style="width: 100px;"
                  step="0.1"
                  precision="1"
              />
            </div>
            <div class="flex justify-between w-full text-xs text-gray-500" style="width: 100%; max-width: 350px; margin: 5px auto 0;">
              <span>Range: {{ originalPpiData.minSensitivity.toFixed(1) }} —— -1.0</span>
            </div>
          </div>

          <div style="display: inline-block; vertical-align: middle; width: 30%; text-align: right; white-space: normal;">
            <el-space direction="horizontal" size="medium">
              <el-button
                  type="success"
                  @click="handleGoClick"
                  :disabled="!hasData || currentSliderValue === appliedSliderValue"
                  style="padding: 0 15px; background: linear-gradient(to right, #0077b6, #0096c7); border-color: #0077b6; font-size: 14px;"
              >
                <el-icon><Position /></el-icon> Filter
              </el-button>
              <el-button
                  type="warning"
                  @click="handleResetClick"
                  :disabled="!hasData || currentSliderValue === -1.0"
                  style="padding: 0 15px; font-size: 14px;"
              >
                <el-icon><RefreshRight /></el-icon> Reset
              </el-button>
              <el-button
                  type="primary"
                  @click="handleExportTxt"
                  :disabled="!hasData"
                  :loading="isExporting"
                  style="padding: 0 15px; background: linear-gradient(to right, #003f88, #00509d); border: none; font-size: 14px;"
              >
                <ElIcon><Download /></ElIcon> Export Data
              </el-button>
            </el-space>
          </div>
        </div>
      </el-card>
    </div>

    <el-divider></el-divider>
    <div class="flex-1 flex flex-col gap-4 p-0 md:p-4 overflow-hidden" style="padding-bottom: 10px">
      <div class="flex-1 flex flex-col md:flex-row gap-4" style="width: 98%; margin: 0 auto">
        <div
            id="chart-container"
            class="w-full md:w-3/4 border rounded overflow-hidden flex-1"
            :style="{ minHeight: '500px', backgroundColor: '#fff' }"
            v-loading="isLoadingData"
            element-loading-text="Loading Synthetic Lethality Network Data..."
            element-loading-background="rgba(255, 255, 255, 0.8)"
        >
          <el-alert
              v-if="1"
              type="error"
              :description="`${loadError}`"
              show-icon
              class="m-4"
              dangerously-use-html-string
          />
          <el-empty
              v-else-if="!hasData && !isLoadingData"
              description="No Synthetic Lethality Data Available"
              class="h-full flex items-center justify-center"
          />
          <el-empty
              v-else-if="ppiData.links.length === 0 && !isLoadingData"
              description="No Data After Filtering, Please Adjust Threshold"
              class="h-full flex items-center justify-center"
          >
            <template #image><ElIcon class="text-5xl text-gray-300"><Filter /></ElIcon></template>
            <template #footer><el-button type="primary" @click="handleResetClick" size="small">Reset to -1.0</el-button></template>
          </el-empty>
        </div>
      </div>

      <el-divider></el-divider>
      <el-card class="flex-shrink-0 result-card" style="min-height: 500px; background-color: #fff; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px; padding: 20px;">
        <div
            class="mb-2 font-bold text-[#00509d] flex items-center justify-between"
            v-if="selectedNode"
            style="margin-bottom: 10px"
        >
          <span>Synthetic Lethality Interactions Related to {{ selectedNode }}</span>
          <el-button
              type="info"
              size="small"
              @click="clearSelection"
              style="background-color: #1F77B4;margin-left: 20px"
          >
            Clear Selection
          </el-button>
        </div>
        <div v-else class="mb-2 font-bold text-[#00509d]">
          Synthetic Lethality Interaction List (Total {{ total }} pairs)
        </div>

        <div>
          <el-table
              ref="unifiedTable"
              :data="tableData"
              v-if="hasData"
              :key="tableRefreshKey"
              border
              max-height="600px"
              style="width: 100%;"
              :header-cell-style="{ 'background-color': '#f0f5ff', 'color': '#003a8c', 'font-size': '14px' }"
              :row-style="{ 'font-size': '14px' }"
              v-loading="isLoadingTable"
              element-loading-text="Loading table data..."
          >
            <el-table-column prop="geneA" label="Gene A" align="center" />
            <el-table-column prop="geneB" label="Gene B" align="center" />
            <el-table-column prop="geminiSensitive" label="GEMINI Sensitivity Score" align="center" />
            <el-table-column prop="cellLine" label="Evidence Source" align="center" />
            <el-table-column label="Action" align="center">
              <template #default="scope">
                <el-button
                    type="text"
                    @click="selectNodePair(scope.row.geneA, scope.row.geneB)"
                    style="color: #1150bd"
                    class="hover:text:#003f88 transition-colors"
                >
                  View in Network
                </el-button>
              </template>
            </el-table-column>
          </el-table>

          <div class="pagination-container flex justify-center mt-4" v-if="hasData && total > 0" style="padding-bottom: 20px">
            <el-pagination
                @size-change="handlePageSizeChange"
                @current-change="handleCurrentPageChange"
                :current-page="currentPage"
                :page-sizes="[10, 20, 50]"
                :page-size="pageSize"
                layout="total, sizes, prev, pager, next, jumper"
                :total="total"
                background
                style="justify-content: center; font-size: 14px"
            />
          </div>

          <el-empty
              v-if="!hasData && !isLoadingData && !loadError"
              description="No synthetic lethality interactions found"
              class="py-4"
          />
        </div>
      </el-card>
    </div>
  </div>
</template>

<script setup>
import {ref, computed, onMounted, onUnmounted, nextTick} from 'vue';
import * as echarts from 'echarts';
import {ElMessage, ElIcon} from 'element-plus';
import {Filter, Download, Position, RefreshRight} from '@element-plus/icons-vue';

// 后端API基础地址
const apiBase = 'http://121.37.88.191:9016';

// 状态管理
const targetProtein = ref('');
const isLoadingData = ref(false);
const isLoadingTable = ref(false);
const loadError = ref('');
const isExporting = ref(false);
const selectedNode = ref('');
const selectedNodeInfo = ref(null);
const currentSliderValue = ref(-1.0);
const appliedSliderValue = ref(-1.0);
const tableRefreshKey = ref(true);
const filteredGene = ref(''); // 记录当前过滤的基因

// 分页相关状态
const currentPage = ref(1);
const pageSize = ref(10);
const total = ref(0);
const tableData = ref([]);

// 颜色映射表
const colorScale = ref([
  {min: 1, max: 3, color: 'hsla(210, 80%, 70%, 0.9)'},
  {min: 4, max: 10, color: 'hsla(210, 80%, 55%, 0.9)'},
  {min: 11, max: 20, color: 'hsla(240, 80%, 60%, 0.9)'},
  {min: 21, max: 50, color: 'hsla(300, 80%, 60%, 0.9)'},
  {min: 51, max: 100, color: 'hsla(340, 80%, 55%, 0.9)'}
]);

// 原始数据存储（未过滤的完整网络）
const originalPpiData = ref({
  nodes: [],
  links: [],
  minSensitivity: -3,
  maxSensitivity: 1,
  nodeLethalCounts: {}
});

// 当前显示数据（过滤后的目标基因关联网络）
const ppiData = ref({
  nodes: [],
  links: [],
  minSensitivity: -3,
  maxSensitivity: 1,
  nodeLethalCounts: {}
});

// 计算属性：判断原始数据是否存在
const hasOriginalData = computed(() =>
    originalPpiData.value.nodes.length > 0 &&
    originalPpiData.value.links.length > 0 &&
    total.value > 0
);

// 计算属性：判断当前显示数据是否存在
const hasData = computed(() =>
    ppiData.value.nodes.length > 0 &&
    ppiData.value.links.length > 0 &&
    total.value > 0
);

// 使用普通变量存储 ECharts 实例，避免 Vue 响应式代理造成的严重性能损耗
let chartInstance = null;
let resizeDebounceTimer = null;
let highlightedNodeIndex = -1;
let highlightedLinkIndices = [];

// 工具函数：统一转换为大写（不区分大小写核心）
const normalizeGeneName = (name) => {
  return name ? name.trim().toUpperCase() : '';
};

// 节点选择和高亮
const selectNode = (nodeName) => {
  if (!hasData.value || !chartInstance) return false;

  const targetName = normalizeGeneName(nodeName);
  const nodeIndex = ppiData.value.nodes.findIndex(node =>
      normalizeGeneName(node.name) === targetName
  );

  if (nodeIndex === -1) {
    clearSelection();
    ElMessage.warning(`${nodeName} not found in network`);
    return false;
  }
  const matchedNode = ppiData.value.nodes[nodeIndex];
  const actualNodeName = matchedNode.name;

  clearHighlights();

  // 高亮当前节点
  highlightedNodeIndex = nodeIndex;
  chartInstance.dispatchAction({type: 'highlight', dataIndex: highlightedNodeIndex});

  // 高亮关联的边
  highlightedLinkIndices = ppiData.value.links
      .map((link, idx) => {
        const source = typeof link.source === 'object' ? link.source.name : link.source;
        const target = typeof link.target === 'object' ? link.target.name : link.target;
        return (normalizeGeneName(source) === targetName || normalizeGeneName(target) === targetName) ? idx : -1;
      })
      .filter(idx => idx !== -1);

  highlightedLinkIndices.forEach(idx => {
    chartInstance.dispatchAction({type: 'highlight', dataIndex: idx, seriesIndex: 0});
  });

  // 更新节点信息与表格
  selectedNode.value = actualNodeName;
  selectedNodeInfo.value = {
    pathway: matchedNode.pathway || 'Not provided',
    lethalCount: ppiData.value.nodeLethalCounts[actualNodeName] || 0,
    lethalGenes: getLethalGenesForNode(actualNodeName),
    function: matchedNode.function || 'No information available'
  };
  currentPage.value = 1;
  loadLethalityTable(appliedSliderValue.value);
  targetProtein.value = actualNodeName;
  tableRefreshKey.value = !tableRefreshKey.value;

  return true;
};

// 选择基因对并在网络中高亮
const selectNodePair = (node1, node2) => {
  if (!hasData.value || !chartInstance) return;

  clearHighlights();
  const target1 = normalizeGeneName(node1);
  const target2 = normalizeGeneName(node2);

  const indices = [node1, node2]
      .map(name => ppiData.value.nodes.findIndex(node =>
          normalizeGeneName(node.name) === normalizeGeneName(name)
      ))
      .filter(idx => idx !== -1);

  indices.forEach(idx => {
    chartInstance.dispatchAction({type: 'highlight', dataIndex: idx});
  });

  const linkIndex = ppiData.value.links.findIndex(link => {
    const source = typeof link.source === 'object' ? link.source.name : link.source;
    const target = typeof link.target === 'object' ? link.target.name : link.target;
    const s = normalizeGeneName(source);
    const t = normalizeGeneName(target);
    return (s === target1 && t === target2) || (s === target2 && t === target1);
  });

  if (linkIndex !== -1) {
    chartInstance.dispatchAction({type: 'highlight', dataIndex: linkIndex, seriesIndex: 0});
    highlightedLinkIndices = [linkIndex];
  }

  highlightedNodeIndex = indices[0];
  const matchedNode = ppiData.value.nodes[highlightedNodeIndex];
  selectedNode.value = matchedNode.name;
  selectedNodeInfo.value = {
    pathway: matchedNode.pathway || 'Not provided',
    lethalCount: ppiData.value.nodeLethalCounts[matchedNode.name] || 0,
    lethalGenes: getLethalGenesForNode(matchedNode.name),
    function: matchedNode.function || 'No information available'
  };

  document.getElementById('chart-container').scrollIntoView({behavior: 'smooth'});
};

// 获取节点的合成致死基因列表
const getLethalGenesForNode = (nodeName) => {
  const targetName = normalizeGeneName(nodeName);
  const lethalGenes = new Set();
  ppiData.value.links.forEach(link => {
    const source = typeof link.source === 'object' ? link.source.name : link.source;
    const target = typeof link.target === 'object' ? link.target.name : link.target;
    const s = normalizeGeneName(source);
    const t = normalizeGeneName(target);

    if (s === targetName) lethalGenes.add(target);
    if (t === targetName) lethalGenes.add(source);
  });
  return Array.from(lethalGenes);
};

// 根据互作数量计算节点颜色
const getNodeColorByLethalCount = (count) => {
  if (count === 0) return 'hsla(0, 0%, 80%, 0.9)';
  for (const scale of colorScale.value) {
    if (count >= scale.min && count <= scale.max) return scale.color;
  }
  return colorScale.value[colorScale.value.length - 1].color;
};

// 清除所有高亮
const clearHighlights = () => {
  if (chartInstance) {
    if (highlightedNodeIndex !== -1) {
      chartInstance.dispatchAction({type: 'downplay', dataIndex: highlightedNodeIndex});
      highlightedNodeIndex = -1;
    }
    highlightedLinkIndices.forEach(idx => {
      chartInstance.dispatchAction({type: 'downplay', dataIndex: idx, seriesIndex: 0});
    });
    highlightedLinkIndices = [];
  }
};

// 清除节点选择
const clearSelection = () => {
  clearHighlights();
  selectedNode.value = '';
  selectedNodeInfo.value = null;
  currentPage.value = 1;
  loadLethalityTable(appliedSliderValue.value);
  targetProtein.value = '';
  tableRefreshKey.value++;
};

// 输入框清空处理 - 恢复显示完整网络
const handleInputClear = () => {
  if (filteredGene.value) {
    filteredGene.value = '';
    ppiData.value = JSON.parse(JSON.stringify(originalPpiData.value));
    initEcharts();
    loadLethalityTable(appliedSliderValue.value);
  }
  clearSelection();
};

// 基因定位处理
const handleProteinHighlight = async () => {
  if (!hasOriginalData.value || !targetProtein.value.trim() || !chartInstance) return;

  const inputGene = normalizeGeneName(targetProtein.value);
  const geneExists = originalPpiData.value.nodes.some(node =>
      normalizeGeneName(node.name) === inputGene
  );

  if (!geneExists) {
    ElMessage.warning(`Gene ${targetProtein.value} not found in network`);
    return;
  }

  const relatedLinks = originalPpiData.value.links.filter(link => {
    const source = typeof link.source === 'object' ? link.source.name : link.source;
    const target = typeof link.target === 'object' ? link.target.name : link.target;
    const s = normalizeGeneName(source);
    const t = normalizeGeneName(target);
    return s === inputGene || t === inputGene;
  });

  if (relatedLinks.length === 0) {
    ElMessage.info(`No direct interactions found for ${targetProtein.value}`);
    return;
  }

  const relatedNodeNames = new Set();
  const originalGeneName = originalPpiData.value.nodes.find(node =>
      normalizeGeneName(node.name) === inputGene
  ).name;
  relatedNodeNames.add(originalGeneName);

  relatedLinks.forEach(link => {
    const source = typeof link.source === 'object' ? link.source.name : link.source;
    const target = typeof link.target === 'object' ? link.target.name : link.target;
    relatedNodeNames.add(source);
    relatedNodeNames.add(target);
  });

  const relatedNodes = originalPpiData.value.nodes.filter(node =>
      relatedNodeNames.has(node.name)
  );

  const nodeLethalCounts = calculateLethalCounts(relatedLinks);
  const formattedNodes = relatedNodes.map(node => {
    const count = nodeLethalCounts[node.name] || 0;
    return {
      ...node,
      symbolSize: Math.min(10 + count * 1.5, 40),
      itemStyle: {
        color: getNodeColorByLethalCount(count),
        borderWidth: 1,
        borderColor: '#ddd'
      },
      label: {
        color: '#333',
        fontSize: Math.min(12, 8 + count * 0.2),
        fontWeight: 'bold',
        overflow: 'truncate',
        width: 70
      },
      tooltip: {}
    };
  });

  ppiData.value = {
    ...originalPpiData.value,
    nodes: formattedNodes,
    links: relatedLinks.map(link => ({
      ...link,
      lineStyle: {
        width: 1.5,           // 优化点：变细
        color: '#b0b8c1',     // 优化点：高级灰蓝色
        curveness: 0.2,
        opacity: 0.6          // 局部过滤时透明度适中
      },
      tooltip: {}
    })),
    nodeLethalCounts
  };

  filteredGene.value = originalGeneName;
  await nextTick();
  initEcharts();
  selectNode(originalGeneName);
};

// 计算节点互作数量
const calculateLethalCounts = (links) => {
  const counts = {};
  if (!links || links.length === 0) return counts;

  links.forEach(link => {
    const source = typeof link.source === 'object' ? link.source.name : link.source;
    const target = typeof link.target === 'object' ? link.target.name : link.target;
    counts[source] = 0;
    counts[target] = 0;
  });

  links.forEach(link => {
    const source = typeof link.source === 'object' ? link.source.name : link.source;
    const target = typeof link.target === 'object' ? link.target.name : link.target;
    counts[source]++;
    counts[target]++;
  });

  return counts;
};

// 处理滑块值变更
const handleSliderValueChange = (value) => {
  if (value > -1.0) {
    currentSliderValue.value = -1.0;
    ElMessage.info('Maximum sensitivity value is -1.0');
  }
};

// 加载图表数据
const loadLethalityData = async (sensitivityThreshold) => {
  isLoadingData.value = true;
  loadError.value = '';

  try {
    const threshold = Math.min(sensitivityThreshold, -1.0);
    const response = await fetch(`${apiBase}/api/network/sl?sensitivityThreshold=${threshold}`);
    if (!response.ok) throw new Error(`Network request failed (Status code: ${response.status})`);
    const data = await response.json();

    if (!data || !Array.isArray(data.links) || !Array.isArray(data.nodes)) {
      throw new Error("Backend returned data format error, missing links or nodes");
    }

    if (data.links.length === 0) {
      const hint = data.minSensitivity >= threshold
          ? `All data sensitivity scores are ≥ ${threshold}, please increase threshold`
          : `No data with sensitivity score < ${threshold} found`;
      loadError.value = hint;
      ElMessage.warning(hint);
      return false;
    }

    const nodeLethalCounts = calculateLethalCounts(data.links);
    const formattedNodes = data.nodes
        .map(node => {
          const count = nodeLethalCounts[node.name] || 0;
          return {
            ...node,
            symbolSize: Math.min(10 + count * 1.5, 40),
            itemStyle: {
              color: getNodeColorByLethalCount(count),
              borderWidth: 1,
              borderColor: '#ddd'
            },
            label: {
              color: '#333',
              fontSize: Math.min(12, 8 + count * 0.2),
              fontWeight: 'bold',
              overflow: 'truncate',
              width: 70
            },
            tooltip: {}
          };
        })
        .filter(node => {
          const count = nodeLethalCounts[node.name] || 0;
          return count > 0;
        });

    originalPpiData.value = {
      ...data,
      nodes: formattedNodes,
      links: data.links.map(link => ({
        ...link,
        lineStyle: {
          width: 1.5,           // 优化点：全局连线变细
          color: '#b0b8c1',     // 优化点：浅灰蓝底色
          curveness: 0.2,
          opacity: 0.4          // 优化点：降低透明度，背景线不抢眼
        },
        tooltip: {}
      })),
      nodeLethalCounts,
      minSensitivity: data.minSensitivity || -3,
      maxSensitivity: data.maxSensitivity || 1
    };

    ppiData.value = JSON.parse(JSON.stringify(originalPpiData.value));
    await nextTick();
    return true;
  } catch (error) {
    loadError.value = `Data loading failed: ${error.message}`;
    ElMessage.error(loadError.value);
    return false;
  } finally {
    isLoadingData.value = false;
  }
};

// 加载表格数据
const loadLethalityTable = async (sensitivityThreshold) => {
  isLoadingTable.value = true;
  loadError.value = '';

  try {
    const threshold = Math.min(sensitivityThreshold, -1.0);
    const params = new URLSearchParams({
      threshold: threshold,
      pageNum: currentPage.value,
      pageSize: pageSize.value
    });

    const target = filteredGene.value || selectedNode.value;
    if (target && target.trim()) {
      params.append('gene', target.trim());
    }

    const requestUrl = `${apiBase}/api/network/sl/table?${params.toString()}`;
    const response = await fetch(requestUrl);
    if (!response.ok) throw new Error(`Table data request failed (Status code: ${response.status})`);
    const res = await response.json();

    if (res.error) throw new Error(`Backend error: ${res.error}`);

    tableData.value = (res.list || []).map(item => ({
      geneA: item.geneA,
      geneB: item.geneB,
      geminiSensitive: item.geminiSensitive,
      cellLine: item.cellLine || 'Not provided'
    }));
    total.value = res.total || 0;
    return true;
  } catch (error) {
    loadError.value = `Table data loading failed: ${error.message}`;
    ElMessage.error(loadError.value);
    return false;
  } finally {
    isLoadingTable.value = false;
  }
};

// 导出数据
const handleExportTxt = async () => {
  if (!hasData.value) return;

  isExporting.value = true;
  try {
    const threshold = Math.min(appliedSliderValue.value, -1.0);
    const params = new URLSearchParams();
    params.append('sensitivity', threshold);

    let exportUrl = '';
    const inputValue = targetProtein.value.trim();
    const normalizedInput = normalizeGeneName(inputValue);

    let geneValue = '';
    if (normalizedInput) {
      const matchedNode = originalPpiData.value.nodes.find(node =>
          normalizeGeneName(node.name) === normalizedInput
      );
      geneValue = matchedNode ? matchedNode.name : inputValue;
    }

    if (geneValue) {
      exportUrl = `${apiBase}/api/network/sl/export-current`;
      params.append('name', geneValue);
    } else {
      exportUrl = `${apiBase}/api/network/sl/export`;
    }

    const response = await fetch(`${exportUrl}?${params.toString()}`, {
      method: 'GET',
      headers: {'Accept': 'text/tab-separated-values'}
    });

    if (!response.ok) throw new Error(`Export request failed (Status code: ${response.status})`);

    const blob = await response.blob();
    const timeStr = new Date().toISOString().replace(/[-:\.T]/g, '').slice(0, 14);
    const fileName = geneValue
        ? `synthetic_lethality_${geneValue}_${timeStr}.tsv`
        : `synthetic_lethality_network_${timeStr}.tsv`;

    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = fileName;
    a.click();
    URL.revokeObjectURL(url);

    ElMessage.success('Synthetic lethality data exported successfully');
  } catch (error) {
    ElMessage.error('Export failed: ' + error.message);
  } finally {
    isExporting.value = false;
  }
};

// 组合数据加载
const loadAllData = async (sensitivityThreshold = -1.0) => {
  clearHighlights();
  selectedNode.value = '';
  selectedNodeInfo.value = '';
  filteredGene.value = '';
  currentPage.value = 1;

  try {
    const threshold = Math.min(sensitivityThreshold, -1.0);
    const chartLoaded = await loadLethalityData(threshold);
    if (!chartLoaded) return;

    currentSliderValue.value = threshold;
    appliedSliderValue.value = threshold;

    const tableLoaded = await loadLethalityTable(threshold);
    if (!tableLoaded) return;

    await nextTick();
    initEcharts();
    tableRefreshKey.value++;
  } catch (error) {
    console.error('Data loading failed:', error);
  }
};

// 初始化图表
const initEcharts = () => {
  if (chartInstance) chartInstance.dispose();
  const container = document.getElementById('chart-container');
  if (!container) return;

  chartInstance = echarts.init(container);

  const chartOption = {
    animationDuration: 1800,
    animationEasingUpdate: 'quinticInOut',
    tooltip: {trigger: 'none'},
    series: [{
      type: 'graph',
      layout: 'force',
      force: {
        repulsion: 300,
        edgeLength: 150,
        gravity: 0.1,
        friction: 0.6
      },
      roam: true,
      draggable: true,
      label: {show: true},
      edgeSymbol: ['none', 'none'],
      edgeSymbolSize: [4, 10],
      edgeLabel: {show: false},
      data: ppiData.value.nodes,
      links: ppiData.value.links,
      emphasis: {
        focus: 'adjacency',
        lineStyle: {
          width: 4,             // 优化点：高亮线宽适中
          color: '#00509d',     // 优化点：深蓝色主色调高亮
          opacity: 0.9
        },
        label: {
          color: '#333',        // 优化点：高亮字体颜色适应浅底色
          fontWeight: 'bold'
        }
      }
    }]
  };

  chartInstance.setOption(chartOption);

  chartInstance.on('click', (params) => {
    if (params.dataType === 'node') selectNode(params.name);
  });

  chartInstance.getZr().on('click', (e) => {
    if (selectedNode.value && !e.target) clearSelection();
  });

  window.addEventListener('resize', () => {
    if (chartInstance) chartInstance.resize();
  });
};

// 分页交互
const handleCurrentPageChange = (val) => {
  currentPage.value = val;
  loadLethalityTable(appliedSliderValue.value);
};

const handlePageSizeChange = (val) => {
  pageSize.value = val;
  currentPage.value = 1;
  loadLethalityTable(appliedSliderValue.value);
};

// 筛选交互
const handleGoClick = () => {
  if (hasData.value && currentSliderValue.value !== appliedSliderValue.value) {
    filteredGene.value = '';
    loadAllData(currentSliderValue.value);
  }
};

// 重置交互
const handleResetClick = () => {
  if (hasData.value && currentSliderValue.value !== -1.0) {
    filteredGene.value = '';
    loadAllData(-1.0);
  }
};

// 窗口resize防抖处理
const handleWindowResize = () => {
  clearTimeout(resizeDebounceTimer);
  resizeDebounceTimer = setTimeout(() => {
    if (chartInstance) chartInstance.resize();
  }, 100);
};

// 生命周期
onMounted(() => {
  loadAllData(-1.0);
  window.addEventListener('resize', handleWindowResize);
});

onUnmounted(() => {
  window.removeEventListener('resize', handleWindowResize);
  if (chartInstance) chartInstance.dispose();
  clearSelection();
  clearTimeout(resizeDebounceTimer);
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

/* 分页控件样式优化，已将弃用的 ::v-deep 升级为 Vue 3 标准的 :deep() */
:deep(.el-pagination) {
  margin-top: 15px;
}

:deep(.el-pagination__total) {
  color: #666;
}

:deep(.el-pager li) {
  margin: 0 5px;
}

:deep(.el-pager li.active) {
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
  margin-top: 25px;
  text-align: center;
  padding: 10px 0;
}

/* 按钮悬停效果 */
:deep(.el-button--primary) {
  transition: all 0.2s;
}

:deep(.el-button--primary:hover) {
  transform: translateY(-2px);
}

:deep(.el-button--success) {
  transition: all 0.2s;
}

:deep(.el-button--success:hover) {
  transform: translateY(-2px);
}

:deep(.el-tag) {
  transition: all 0.2s;
}

:deep(.el-tag:hover) {
  background-color: #e0efff;
  transform: translateY(-1px);
}
</style>