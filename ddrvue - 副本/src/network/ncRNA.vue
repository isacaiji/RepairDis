<template>
  <div class="bg-gray-50 min-h-screen flex flex-col" style="height: 100%; background-color: #f8f9fb;">
    <!-- 顶部控制区 -->
    <div class="mb-4 flex-shrink-0" style="width: 90%; margin: 0 auto; padding-top: 30px;">
      <div class="header" style="margin-bottom: 10px; text-align: left;">
        <h1 style="color: #003f88; margin-bottom: 5px; font-size: 24px; font-weight: 600;">ncRNA Target Regulatory Network</h1>
      </div>
      <el-divider></el-divider>
      <div class="sl-hunt">
        <p><strong>ncRNA Target Regulatory Network:</strong> ncRNA Target Regulatory Network
          is a systems-level framework embodying regulatory circuitry exerted by ncRNAs over cognate target molecules,
          with nodes as molecular entities and edges as empirically validated regulatory interactions.</p>
        <p><strong>Data Source:</strong> The data is derived from literature collected in the
          <a href="https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/php/index.php" target="_blank" style="color: #00509d;"> miRTarBase </a>
          database</p>
      </div>

      <!-- 搜索和控制区域卡片 -->
      <el-card class="search-card" style="margin-top: 20px; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px;">
        <div style="margin: 0 auto; display: flex; align-items: center; gap: 10px; width: 100%; padding: 10px; flex-wrap: wrap;">
          <span style="color: #333; font-weight: 500; white-space: nowrap;">Name:</span>
          <el-input
              v-model="targetNcRNA"
              placeholder="e.g., hsa-miR-15a-5p or BCL2"
              style="height: 40px; flex: 1; min-width: 250px;"
              clearable
              @clear="handleInputClear"
              @keyup.enter="handleNcRNALocate"
          />
          <el-button
              type="primary"
              @click="handleNcRNALocate"
              :disabled="!targetNcRNA.trim()"
              style="height: 40px; padding: 0 20px; background: linear-gradient(to right, #003f88, #00509d); border: none; font-size: 14px; white-space: nowrap;"
          >
            Locate
          </el-button>
          <el-button
              type="success"
              @click="handleExportTxt"
              :disabled="!hasData"
              :loading="isExporting"
              style="height: 40px; padding: 0 20px; background: linear-gradient(to right, #0077b6, #0096c7); border-color: #0077b6; font-size: 14px; white-space: nowrap;"
          >
            <ElIcon><Download /></ElIcon> Export Data
          </el-button>
        </div>

        <!-- 图例 -->
        <div style="margin-top: 15px; display: flex; gap: 25px; flex-wrap: wrap; padding: 10px 10px 0; background-color: #f9f9f9; border-radius: 8px; border: 1px solid #eee;">
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 16px; height: 16px; background-color: #dc3545; clip-path: polygon(50% 0%, 100% 50%, 50% 100%, 0% 50%);"></div>
            <span class="text-sm">ncRNA</span>
          </div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 16px; height: 16px; border-radius: 50%; background-color: #28a745;"></div>
            <span class="text-sm">Target Gene</span>
          </div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 20px; height: 2px; background-color: #6c757d;"></div>
            <span class="text-sm">ncRNA → Target Gene</span>
          </div>
        </div>
      </el-card>
    </div>

    <el-divider></el-divider>

    <!-- 图表 + 信息面板 -->
    <div class="flex-1 flex flex-col gap-4 p-0 md:p-4 overflow-hidden" style="padding-bottom: 10px">
      <div class="flex-1 flex flex-col md:flex-row gap-4" style="width: 98%; margin: 0 auto">
        <!-- 图表容器 -->
        <div
            id="ncRNA-chart-container"
            class="w-full md:w-3/4 border rounded overflow-hidden flex-1"
            :style="{ minHeight: '500px' }"
            v-if="hasData || isLoading"
        >
          <el-loading
              v-if="isLoading"
              target="#ncRNA-chart-container"
              text="Loading relationships..."
              background="rgba(255, 255, 255, 0.8)"
          />
          <el-alert
              v-else-if="loadError"
              type="error"
              :description="loadError"
              show-icon
              class="m-4"
          />
          <el-empty
              v-else-if="!hasData && !isLoading"
              description="Enter a name and click 'Locate' to view network"
              class="h-full flex items-center justify-center"
          />
        </div>

        <!-- 信息面板 -->
        <el-card class="w-full md:w-1/4 bg-white rounded-lg border p-4" v-if="hasData || isLoading">
          <div v-if="isLoading" class="flex justify-center items-center h-40">
            <el-loading-spinner size="large" />
          </div>

          <div v-else-if="queryInfo" class="space-y-3">
            <h3 style="color: #003f88; margin: 0; font-size: 16px; font-weight: 600;">{{ queryInfo.title }}</h3>
            <div class="text-sm">
              <p><span style="font-weight: bold; color: #333;">Total relationships:</span> {{ queryInfo.relationCount }}</p>
              <p><span style="font-weight: bold; color: #333;">Related targets:</span> {{ queryInfo.targetCount }}</p>
            </div>
            <el-divider style="margin: 10px 0;"></el-divider>
            <div v-if="queryInfo.targetGenes.length > 0">
              <p class="text-sm font-weight: bold;">Related Entities:</p>
              <div class="flex flex-wrap gap-2 mt-1">
                <el-tag
                    v-for="gene in queryInfo.targetGenes.slice(0, 8)"
                    :key="gene"
                    style="cursor: pointer; background-color: #f0f7ff; color: #1e88e5; border-color: #bfdbfe;"
                    @click="handleGeneRelocate(gene)"
                >
                  {{ gene }}
                </el-tag>
                <el-tag v-if="queryInfo.targetGenes.length > 8" type="info">+{{ queryInfo.targetGenes.length - 8 }} more</el-tag>
              </div>
            </div>
          </div>
        </el-card>
      </div>

      <el-divider></el-divider>

      <!-- 表格 -->
      <el-card class="flex-shrink-0 result-card" style="min-height: 400px; background-color: #fff; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px; padding: 20px;">
        <div v-if="hasData" class="mb-2 font-bold text-[#00509d]">
          Regulatory Relationships for {{ targetNcRNA }} (Total {{ total }} pairs)
        </div>

        <el-table
            :data="tableData"
            v-if="hasData"
            border
            max-height="500px"
            style="width: 100%;"
            :header-cell-style="{ 'background-color': '#f0f5ff', 'color': '#003a8c', 'font-size': '14px' }"
            :row-style="{ 'font-size': '14px' }"
            v-loading="isLoadingTable"
            element-loading-text="Loading table data..."
        >
          <el-table-column prop="mirTarBaseId" label="miRTarBase ID" min-width="120" show-overflow-tooltip />
          <el-table-column prop="ncRNA" label="ncRNA" min-width="120" show-overflow-tooltip />
          <el-table-column prop="ncRNASpecies" label="ncRNA Species" min-width="100" />
          <el-table-column prop="targetGene" label="Target Gene" min-width="120" show-overflow-tooltip />
          <el-table-column prop="targetGeneEntrezId" label="Entrez ID" min-width="100" />
          <el-table-column prop="targetGeneSpecies" label="Target Species" min-width="100" />
          <el-table-column prop="experiments" label="Experiments" min-width="150" show-overflow-tooltip />
          <el-table-column prop="supportType" label="Support Type" min-width="100" />
          <el-table-column prop="reference" label="PMID" min-width="100">
            <template #default="scope">
              <el-link
                  v-if="scope.row.reference"
                  type="primary"
                  :href="`https://pubmed.ncbi.nlm.nih.gov/${scope.row.reference}/`"
                  target="_blank"
                  style="color: #1150bd"
              >
                {{ scope.row.reference }}
              </el-link>
              <span v-else>—</span>
            </template>
          </el-table-column>
          <el-table-column label="Action" width="150" fixed="right">
            <template #default="scope">
              <el-button
                  type="text"
                  @click="highlightRelationship(scope.row)"
                  style="color: #1150bd"
              >
                Highlight in Network
              </el-button>
            </template>
          </el-table-column>
        </el-table>

        <!-- 分页 -->
        <div class="pagination-container flex justify-center mt-4" v-if="hasData && total > 0"
             style="padding-bottom: 20px;margin: 20px auto 0;">
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
            v-if="!hasData && !isLoading && !loadError"
            description="Enter a name and click 'Locate' to view relationships"
            class="py-4"
        />
      </el-card>
    </div>
  </div>
</template>

<script setup>
import {ref, onMounted, onUnmounted, nextTick} from 'vue';
import * as echarts from 'echarts';
import {ElMessage, ElIcon, ElLoading} from 'element-plus';
import {Download} from '@element-plus/icons-vue';

const apiBase = 'http://121.37.88.191:9016';

// 状态
const targetNcRNA = ref('');
const isLoading = ref(false);
const isLoadingTable = ref(false);
const loadError = ref('');
const isExporting = ref(false);
const hasData = ref(false);
const queryInfo = ref(null);
const currentPage = ref(1);
const pageSize = ref(10);
const total = ref(0);
const tableData = ref([]);
const regData = ref({nodes: [], links: []});

let chartInstance = null;

// 核心：定位
const handleNcRNALocate = async () => {
  if (!targetNcRNA.value.trim()) {
    ElMessage.warning('Please enter a name');
    return;
  }

  const inputName = targetNcRNA.value.trim();
  isLoading.value = true;
  loadError.value = '';
  hasData.value = false;

  if (chartInstance) {
    chartInstance.dispose();
    chartInstance = null;
  }

  try {
    const params = new URLSearchParams({name: inputName});
    const response = await fetch(`${apiBase}/api/network/ncrna?${params.toString()}`);
    if (!response.ok) {
      throw new Error(`Request failed (Status: ${response.status})`);
    }

    const data = await response.json();
    if (!data.links || data.links.length === 0) {
      throw new Error(`No relationships found for: ${inputName}`);
    }

    const searchMode = data.searchMode; // "ncRNA" 或 "gene"
    const queryTerm = data.query;       // 实际查询词（大小写可能不同）

    const allNodes = new Set();
    data.links.forEach(link => {
      allNodes.add(link.source);
      allNodes.add(link.target);
    });

    const nodes = Array.from(allNodes).map(name => {
      const isQueryNode = name.toLowerCase() === queryTerm.toLowerCase();

      if (isQueryNode) {
        if (searchMode === 'gene') {
          return {
            name,
            type: 'query_gene',
            symbol: 'circle',
            symbolSize: 22,
            itemStyle: {color: '#28a745'}
          };
        } else {
          return {
            name,
            type: 'query_ncRNA',
            symbol: 'diamond',
            symbolSize: 22,
            itemStyle: {color: '#dc3545'}
          };
        }
      }

      const isNcRNA = data.links.some(link => link.source === name);
      if (isNcRNA) {
        return {
          name,
          type: 'ncRNA',
          symbol: 'diamond',
          symbolSize: 18,
          itemStyle: {color: '#dc3545'}
        };
      } else {
        return {
          name,
          type: 'target_gene',
          symbol: 'circle',
          symbolSize: 18,
          itemStyle: {color: '#28a745'}
        };
      }
    });

    const links = data.links.map(link => ({
      source: link.source,
      target: link.target,
      lineStyle: {
        width: 2.5,
        curveness: 0.2,
        color: '#6c757d',
        opacity: 0.9
      }
    }));

    regData.value = {nodes, links};
    hasData.value = true;

    // 更新信息面板
    const nonQueryNodes = Array.from(allNodes).filter(n => n.toLowerCase() !== queryTerm.toLowerCase());
    queryInfo.value = {
      title: searchMode === 'gene' ? `Gene: ${queryTerm}` : `ncRNA: ${queryTerm}`,
      targetCount: nonQueryNodes.length,
      relationCount: links.length,
      targetGenes: nonQueryNodes
    };

    // 加载表格
    await loadTargetTable(inputName);

    await nextTick();
    initNcRNAEcharts();

  } catch (error) {
    loadError.value = error.message;
    ElMessage.error(loadError.value);
  } finally {
    isLoading.value = false;
  }
};

// 从信息面板重新定位
const handleGeneRelocate = (name) => {
  targetNcRNA.value = name;
  handleNcRNALocate();
};

// 加载表格
const loadTargetTable = async (name) => {
  isLoadingTable.value = true;
  try {
    const params = new URLSearchParams({
      name,
      pageNum: currentPage.value,
      pageSize: pageSize.value
    });
    const response = await fetch(`${apiBase}/api/network/ncrna/table?${params.toString()}`);
    if (!response.ok) throw new Error(`Table load failed (Status: ${response.status})`);

    const res = await response.json();
    tableData.value = res.list;
    total.value = res.total || 0;
  } catch (error) {
    ElMessage.error(`Table load failed: ${error.message}`);
  } finally {
    isLoadingTable.value = false;
  }
};

// 初始化图表
const initNcRNAEcharts = () => {
  const container = document.getElementById('ncRNA-chart-container');
  if (!container || !chartInstance) {
    chartInstance = echarts.init(container);
  }

  const option = {
    animationDuration: 1500,
    animationEasingUpdate: 'quinticInOut',
    tooltip: {
      trigger: 'item',
      formatter: (params) => {
        if (params.dataType === 'node') {
          let label = 'Unknown';
          if (params.data.type === 'query_ncRNA') label = 'Query ncRNA';
          else if (params.data.type === 'ncRNA') label = 'ncRNA';
          else if (params.data.type === 'query_gene') label = 'Query Gene';
          else if (params.data.type === 'target_gene') label = 'Target Gene';
          return `${params.name} (${label})`;
        } else {
          return `${params.data.source} → ${params.data.target}`;
        }
      }
    },
    series: [{
      type: 'graph',
      layout: 'force',
      force: {repulsion: 300, edgeLength: 150, gravity: 0.2},
      roam: true,
      draggable: true,
      label: {
        show: true,
        fontSize: 12,
        color: '#333',
        overflow: 'truncate',
        width: 60
      },
      edgeSymbol: ['none', 'arrow'],
      edgeSymbolSize: [4, 10],
      data: regData.value.nodes,
      links: regData.value.links,
      emphasis: {
        focus: 'adjacency',
        lineStyle: {width: 5},
        label: {fontSize: 14}
      }
    }]
  };

  chartInstance.setOption(option);

  chartInstance.on('click', (params) => {
    if (params.dataType === 'node' && ['target_gene', 'ncRNA'].includes(params.data.type)) {
      handleGeneRelocate(params.data.name);
    }
  });

  const handleResize = () => chartInstance?.resize();
  window.addEventListener('resize', handleResize);
  onUnmounted(() => window.removeEventListener('resize', handleResize));
};

// 高亮
const highlightRelationship = (row) => {
  if (!chartInstance) return;

  chartInstance.dispatchAction({type: 'downplay', seriesIndex: 0});

  const nodes = regData.value.nodes;
  const ncRNAIndex = nodes.findIndex(n => n.name === row.ncRNA);
  const geneIndex = nodes.findIndex(n => n.name === row.targetGene);
  const linkIndex = regData.value.links.findIndex(
      l => l.source === row.ncRNA && l.target === row.targetGene
  );

  if (ncRNAIndex !== -1) chartInstance.dispatchAction({type: 'highlight', seriesIndex: 0, dataIndex: ncRNAIndex});
  if (geneIndex !== -1) chartInstance.dispatchAction({type: 'highlight', seriesIndex: 0, dataIndex: geneIndex});
  if (linkIndex !== -1) chartInstance.dispatchAction({type: 'highlight', seriesIndex: 0, dataIndex: linkIndex});
};

// 其他方法
const handleGeneDetail = (name) => ElMessage.info(`Detail: ${name}`);
const handleInputClear = () => {
  // 重置所有状态
  hasData.value = false;
  queryInfo.value = null;
  tableData.value = [];
  total.value = 0;
  targetNcRNA.value = '';
  if (chartInstance) {
    chartInstance.dispose();
    chartInstance = null;
  }
};

const handleExportTxt = async () => {
  if (!hasData.value || !targetNcRNA.value.trim()) return;
  isExporting.value = true;
  try {
    const params = new URLSearchParams({name: targetNcRNA.value.trim()});
    const response = await fetch(`${apiBase}/api/network/ncrna/export-current?name=${targetNcRNA.value}`, {
      headers: {'Accept': 'text/tab-separated-values'}
    });
    if (!response.ok) throw new Error(`Export failed (Status: ${response.status})`);
    const blob = await response.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${targetNcRNA.value.trim()}_relationships.tsv`;
    a.click();
    URL.revokeObjectURL(url);
    ElMessage.success('Exported successfully');
  } catch (error) {
    ElMessage.error(`Export failed: ${error.message}`);
  } finally {
    isExporting.value = false;
  }
};

const handleCurrentPageChange = (val) => {
  currentPage.value = val;
  if (targetNcRNA.value.trim()) loadTargetTable(targetNcRNA.value.trim());
};
const handlePageSizeChange = (val) => {
  pageSize.value = val;
  currentPage.value = 1;
  if (targetNcRNA.value.trim()) loadTargetTable(targetNcRNA.value.trim());
};

onMounted(() => {
  chartInstance = null;
});
onUnmounted(() => {
  if (chartInstance) {
    chartInstance.dispose();
    chartInstance = null;
  }
});
</script>

<style scoped>
#ncRNA-chart-container {
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

::v-deep .el-button--success {
  transition: all 0.2s;
}

::v-deep .el-button--success:hover {
  transform: translateY(-2px);
}

::v-deep .el-tag {
  transition: all 0.2s;
}

::v-deep .el-tag:hover {
  background-color: #e0efff;
  transform: translateY(-1px);
}
</style>