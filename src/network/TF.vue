<template>
  <div class="bg-gray-50 min-h-screen flex flex-col" style="height: 100%; background-color: #f8f9fb;">
    <!-- 顶部控制区：聚焦基因搜索 -->
    <div class="mb-4 flex-shrink-0" style="width: 90%; margin: 0 auto; padding-top: 30px;">
      <div class="header" style="margin-bottom: 10px; text-align: left;">
        <h1 style="color: #003f88; margin-bottom: 5px; font-size: 24px; font-weight: 600;">TF Regulatory Network</h1>
      </div>
      <el-divider></el-divider>
      <div class="sl-hunt">
        <p><strong>TF Regulatory Network:</strong> Transcription Factor Regulatory Network
          is a systems-level architecture encapsulating transcriptional control exerted by TFs over target genes,
          with nodes as molecular players and edges as empirically substantiated regulatory events.</p>
        <p><strong>Data Source:</strong> The data is derived from literature collected in the
          <a href="https://www.grnpedia.org/trrust/" target="_blank" style="color: #00509d;"> TRRUST </a>
          database</p>
      </div>

      <!-- 搜索和控制区域卡片 -->
      <el-card class="search-card" style="margin-top: 20px; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px;">
        <!-- 核心：基因搜索框 + 按钮组 -->
        <div style="margin: 0 auto; display: flex; align-items: center; gap: 10px; width: 100%; padding: 10px; flex-wrap: wrap;">
          <span style="color: #333; font-weight: 500; white-space: nowrap;">Gene:</span>
          <el-input
              v-model="targetGene"
              placeholder="Enter gene symbol (e.g., BRCA1, CDKN1A)"
              style="height: 40px; flex: 1; min-width: 250px;"
              clearable
              @clear="handleInputClear"
              @keyup.enter="handleGeneLocate"
          />
          <el-button
              type="primary"
              @click="handleGeneLocate"
              :disabled="!targetGene.trim()"
              style="height: 40px; padding: 0 20px; background: linear-gradient(to right, #003f88, #00509d); border: none; font-size: 14px; white-space: nowrap;"
          >
            Locate Gene
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

        <!-- 图例：调控关系类型 -->
        <div style="margin-top: 15px; display: flex; gap: 25px; flex-wrap: wrap; padding: 10px 10px 0; background-color: #f9f9f9; border-radius: 8px; border: 1px solid #eee;">
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 20px; height: 2px; background-color: #28a745; border-radius: 1px;"></div>
            <span class="text-sm">Activation</span>
          </div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 20px; height: 2px; background-color: #dc3545; border-radius: 1px;"></div>
            <span class="text-sm">Repression</span>
          </div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 20px; height: 2px; background-color: #6c757d; border-radius: 1px;"></div>
            <span class="text-sm">Unknown</span>
          </div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 16px; height: 16px; border-radius: 50%; background-color: #1e88e5;"></div>
            <span class="text-sm">Query Gene</span>
          </div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 16px; height: 16px; border-radius: 50%; background-color: #ffb300;"></div>
            <span class="text-sm">Regulatory Gene</span>
          </div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <div style="width: 16px; height: 16px; border-radius: 50%; background-color: #ab47bc;"></div>
            <span class="text-sm">Target Gene</span>
          </div>
        </div>
      </el-card>
    </div>

    <el-divider></el-divider>

    <!-- 下方视图区：图表 + 表格 -->
    <div class="flex-1 flex flex-col gap-4 p-0 md:p-4 overflow-hidden" style="padding-bottom: 10px">
      <!-- 图表容器 + 信息面板 -->
      <div class="flex-1 flex flex-col md:flex-row gap-4" style="width: 98%; margin: 0 auto">
        <!-- 图表容器：仅在有数据时显示 -->
        <div
            id="reg-chart-container"
            class="w-full md:w-3/4 border rounded overflow-hidden flex-1"
            :style="{ minHeight: '500px' }"
            v-if="hasData || isLoading"
        >
          <el-loading
              v-if="isLoading"
              target="#reg-chart-container"
              text="Loading regulatory relationships..."
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
              v-else-if="!hasData && !isLoading"
              description="Please enter a gene and click 'Locate Gene' to view regulatory network"
              class="h-full flex items-center justify-center"
          />
        </div>

        <!-- 基因信息面板：显示查询基因的调控详情 -->
        <el-card class="w-full md:w-1/4 bg-white rounded-lg border p-4" v-if="hasData || isLoading">
          <div v-if="isLoading" class="flex justify-center items-center h-40">
            <el-loading-spinner size="large" />
          </div>

          <div v-else-if="queryGeneInfo" class="space-y-3">
            <h3 style="color: #003f88; margin: 0; font-size: 16px; font-weight: 600;">Gene: {{ targetGene }}</h3>
            <div class="text-sm">
              <p><span style="font-weight: bold; color: #333;">As TF regulating:</span> {{ queryGeneInfo.asTFCount }} genes</p>
              <p><span style="font-weight: bold; color: #333;">As target regulated by:</span> {{ queryGeneInfo.asTargetCount }} TFs</p>
              <p><span style="font-weight: bold; color: #333;">Activations:</span> {{ queryGeneInfo.activationCount }} relationships</p>
              <p><span style="font-weight: bold; color: #333;">Repressions:</span> {{ queryGeneInfo.repressionCount }} relationships</p>
              <p><span style="font-weight: bold; color: #333;">Unknown regulation:</span> {{ queryGeneInfo.unknownCount }} relationships</p>
            </div>
            <el-divider style="margin: 10px 0;"></el-divider>
            <div v-if="queryGeneInfo.relatedGenes.length > 0">
              <p class="text-sm font-weight: bold;">Related Genes:</p>
              <div class="flex flex-wrap gap-2 mt-1">
                <el-tag
                    v-for="gene in queryGeneInfo.relatedGenes.slice(0, 8)"
                    :key="gene"
                    style="cursor: pointer; background-color: #f0f7ff; color: #1e88e5; border-color: #bfdbfe;"
                    @click="handleGeneRelocate(gene)"
                >
                  {{ gene }}
                </el-tag>
                <el-tag v-if="queryGeneInfo.relatedGenes.length > 8" type="info">+{{ queryGeneInfo.relatedGenes.length - 8 }} more</el-tag>
              </div>
            </div>
          </div>
        </el-card>
      </div>

      <el-divider></el-divider>

      <!-- 调控关系表格 -->
      <el-card class="flex-shrink-0 result-card" style="min-height: 400px; background-color: #fff; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px; padding: 20px;">
        <div v-if="hasData" class="mb-2 font-bold text-[#00509d]">
          Regulatory Relationships for {{ targetGene }} (Total {{ total }} pairs)
        </div>

        <el-table
            ref="regTable"
            :data="tableData"
            v-if="hasData"
            :key="tableRefreshKey"
            border
            max-height="500px"
            style="width: 100%;"
            :header-cell-style="{ 'background-color': '#f0f5ff', 'color': '#003a8c', 'font-size': '14px' }"
            :row-style="{ 'font-size': '14px' }"
            v-loading="isLoadingTable"
            element-loading-text="Loading table data..."
        >
          <el-table-column prop="tf" label="TF" align="center" />
          <el-table-column prop="target" label="Target" align="center" />
          <el-table-column prop="regulationType" label="Regulation Type" align="center">
            <template #default="scope">
              <span :style="{
                color: scope.row.regulationType === 'Activation' ? '#28a745' :
                       scope.row.regulationType === 'Repression' ? '#dc3545' : '#6c757d',
                fontWeight: 'bold'
              }">
                {{ scope.row.regulationType }}
              </span>
            </template>
          </el-table-column>
          <el-table-column prop="evidence" label="Evidence Source (PMID)" align="center" />
          <el-table-column label="Action" align="center">
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

        <!-- 分页控件 -->
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
            v-if="!hasData && !isLoading && !loadError"
            description="Please enter a gene and click 'Locate Gene' to view regulatory relationships"
            class="py-4"
        />
      </el-card>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted, onUnmounted, nextTick } from 'vue';
import * as echarts from 'echarts';
import { ElMessage, ElIcon, ElLoading } from 'element-plus';
import { Download } from '@element-plus/icons-vue';

// 后端API基础地址
// const apiBase = 'http://localhost:9016';
const apiBase = 'http://121.37.88.191:9016';
// 状态管理
const targetGene = ref('');
const isLoading = ref(false);
const isLoadingTable = ref(false);
const loadError = ref('');
const isExporting = ref(false);
const hasData = ref(false);
const queryGeneInfo = ref(null);
const tableRefreshKey = ref(true);

// 分页状态
const currentPage = ref(1);
const pageSize = ref(10);
const total = ref(0);
const tableData = ref([]);

// 网络数据存储
const regData = ref({
  nodes: [],
  links: []
});

// 图表实例
let chartInstance = null;
let highlightedElement = null;

// 核心方法：定位基因并加载其调控网络
const handleGeneLocate = async () => {
  if (!targetGene.value.trim()) {
    ElMessage.warning('Please enter a gene symbol first');
    return;
  }

  const gene = targetGene.value.trim();
  isLoading.value = true;
  loadError.value = '';
  hasData.value = false;

  // 重置图表实例
  if (chartInstance) {
    chartInstance.dispose();
    chartInstance = null;
  }

  try {
    const response = await fetch(`${apiBase}/api/network/tf?gene=${gene}`);
    if (!response.ok) {
      if (response.status === 404) {
        throw new Error(`No regulatory relationships found for ${gene}`);
      }
      throw new Error(`Request failed (Status: ${response.status})`);
    }

    const data = await response.json();
    if (!data.links || data.links.length === 0) {
      throw new Error(`No regulatory relationships found for ${gene}`);
    }

    // 格式化节点和边 - 支持基因同时作为TF和target
    const allGenes = new Set();
    data.links.forEach(link => {
      allGenes.add(link.source);
      allGenes.add(link.target);
    });

    // 创建节点 - 区分查询基因、其他TF和靶基因
    const nodes = Array.from(allGenes).map(name => {
      if (name === gene) {
        return {
          name,
          type: 'query',
          symbolSize: 22,
          itemStyle: { color: '#1e88e5' }
        };
      } else if (data.links.some(link => link.source === name)) {
        return {
          name,
          type: 'tf',
          symbolSize: 18,
          itemStyle: { color: '#ffb300' }
        };
      } else {
        return {
          name,
          type: 'target',
          symbolSize: 18,
          itemStyle: { color: '#ab47bc' }
        };
      }
    });

    // 创建边 - 保持TF到target的方向
    const links = data.links.map(link => ({
      ...link,
      source: link.source,
      target: link.target,
      lineStyle: {
        width: 2.5,
        curveness: 0.2,
        color: getLinkColor(link.regulationType),
        opacity: 0.9
      }
    }));

    regData.value = { nodes, links };
    hasData.value = true;

    // 更新基因信息
    queryGeneInfo.value = {
      asTFCount: links.filter(l => l.source === gene).length,
      asTargetCount: links.filter(l => l.target === gene).length,
      activationCount: links.filter(l => l.regulationType === 'Activation').length,
      repressionCount: links.filter(l => l.regulationType === 'Repression').length,
      unknownCount: links.filter(l => l.regulationType === 'Unknown').length,
      relatedGenes: Array.from(allGenes).filter(name => name !== gene)
    };

    // 加载表格数据
    await loadRegulationTable(gene);

    // 确保DOM更新后再初始化图表
    await nextTick();
    initRegEcharts();

  } catch (error) {
    loadError.value = error.message;
    ElMessage.error(loadError.value);
  } finally {
    isLoading.value = false;
  }
};

// 获取边的颜色
const getLinkColor = (regulationType) => {
  switch (regulationType) {
    case 'Activation': return '#28a745';
    case 'Repression': return '#dc3545';
    default: return '#6c757d';
  }
};

// 加载表格数据
const loadRegulationTable = async (gene) => {
  isLoadingTable.value = true;
  try {
    const params = new URLSearchParams({
      gene: gene,
      pageNum: currentPage.value,
      pageSize: pageSize.value
    });
    const response = await fetch(`${apiBase}/api/network/tf/table?${params.toString()}`);
    if (!response.ok) throw new Error(`Table request failed (Status: ${response.status})`);

    const res = await response.json();
    tableData.value = res.list.map(item => ({
      tf: item.source,
      target: item.target,
      regulationType: item.regulationType || 'Unknown',
      evidence: item.evidence || 'No PMID'
    }));
    total.value = res.total || 0;
  } catch (error) {
    ElMessage.error(`Table load failed: ${error.message}`);
  } finally {
    isLoadingTable.value = false;
  }
};

// 初始化图表
const initRegEcharts = () => {
  const container = document.getElementById('reg-chart-container');
  if (!container) {
    ElMessage.error('Chart container not found');
    return;
  }

  // 销毁可能存在的旧实例
  if (chartInstance) {
    chartInstance.dispose();
  }

  // 初始化图表实例
  chartInstance = echarts.init(container);

  // 图表配置
  const option = {
    animationDuration: 1500,
    animationEasingUpdate: 'quinticInOut',
    tooltip: {
      trigger: 'item',
      formatter: params => {
        if (params.dataType === 'node') {
          let typeLabel = '';
          switch(params.data.type) {
            case 'query': typeLabel = 'Query Gene'; break;
            case 'tf': typeLabel = 'Regulatory TF'; break;
            case 'target': typeLabel = 'Target Gene'; break;
          }
          return `${params.name} (${typeLabel})`;
        } else {
          return `${params.data.source} → ${params.data.target}: ${params.data.regulationType}`;
        }
      }
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
        lineStyle: { width: 5 },
        label: { fontSize: 14 }
      }
    }]
  };

  // 设置配置项
  chartInstance.setOption(option);

  // 绑定点击事件
  chartInstance.on('click', (params) => {
    if (chartInstance && params.dataType === 'node' && params.data.name !== targetGene.value) {
      handleGeneRelocate(params.data.name);
    }
  });

  // 窗口resize适配
  const handleResize = () => {
    if (chartInstance) chartInstance.resize();
  };
  window.addEventListener('resize', handleResize);

  // 组件卸载时移除监听
  onUnmounted(() => {
    window.removeEventListener('resize', handleResize);
  });
};

// 高亮调控关系
const highlightRelationship = (relation) => {
  if (!chartInstance || !regData.value.links.length || !regData.value.nodes.length) {
    return;
  }

  // 清除之前的高亮
  chartInstance.dispatchAction({ type: 'downplay', seriesIndex: 0 });

  // 高亮TF节点
  const tfIndex = regData.value.nodes.findIndex(node => node.name === relation.tf);
  if (tfIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'highlight',
      seriesIndex: 0,
      dataIndex: tfIndex
    });
  }

  // 高亮靶基因节点
  const targetIndex = regData.value.nodes.findIndex(node => node.name === relation.target);
  if (targetIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'highlight',
      seriesIndex: 0,
      dataIndex: targetIndex
    });
  }

  // 高亮调控边
  const linkIndex = regData.value.links.findIndex(link =>
      link.source === relation.tf &&
      link.target === relation.target &&
      link.regulationType === relation.regulationType
  );

  if (linkIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'highlight',
      seriesIndex: 0,
      dataIndex: linkIndex
    });
  }

  highlightedElement = relation;
};

// 从信息面板重新定位基因
const handleGeneRelocate = (geneName) => {
  targetGene.value = geneName;
  handleGeneLocate();
};

// 输入框清空处理
const handleInputClear = () => {
  if (hasData.value) {
    hasData.value = false;
    regData.value = { nodes: [], links: [] };
    queryGeneInfo.value = null;
    tableData.value = [];
    total.value = 0;
    if (chartInstance) {
      chartInstance.dispose();
      chartInstance = null;
    }
  }
};

// 导出数据
const handleExportTxt = async () => {
  if (!hasData.value || !targetGene.value.trim()) return;

  isExporting.value = true;
  try {
    const response = await fetch(`${apiBase}/api/network/tf/export-current?name=${targetGene.value.trim()}`, {
      method: 'GET',
      headers: { 'Accept': 'text/tab-separated-values' }
    });

    if (!response.ok) throw new Error(`Export failed (Status: ${response.status})`);
    const blob = await response.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${targetGene.value.trim()}_regulatory_relationships.tsv`;
    a.click();
    URL.revokeObjectURL(url);

    ElMessage.success('Regulatory data exported successfully');
  } catch (error) {
    ElMessage.error(`Export failed: ${error.message}`);
  } finally {
    isExporting.value = false;
  }
};

// 分页交互
const handleCurrentPageChange = (val) => {
  currentPage.value = val;
  if (targetGene.value.trim()) {
    loadRegulationTable(targetGene.value.trim());
  }
};

const handlePageSizeChange = (val) => {
  pageSize.value = val;
  currentPage.value = 1;
  if (targetGene.value.trim()) {
    loadRegulationTable(targetGene.value.trim());
  }
};

// 生命周期
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
#reg-chart-container {
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