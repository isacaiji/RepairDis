<template>
  <div class="bg-gray-50 min-h-screen flex flex-col" style="height: 100%; background-color: #f8f9fb;">
    <!-- 顶部控制区：聚焦基因搜索 -->
    <div class="mb-4 flex-shrink-0" style="width: 90%; margin: 0 auto; padding-top: 30px;">
      <div class="header" style="margin-bottom: 10px; text-align: left;">
        <h1 style="color: #003f88; margin-bottom: 5px; font-size: 24px; font-weight: 600;">Cross-talk Network</h1>
      </div>
      <el-divider></el-divider>
      <div class="sl-hunt">
        <p><strong>Cross-talk Network:</strong> cross-talk network denotes an elaborate molecular or signaling
          circuitry wherein discrete pathways reciprocally interfere, modulate, or crosspond, thereby orchestrating
          physiological or pathological cascades in biological systems via synergistic or antagonistic crosstalk.</p>
        <p><strong>Data Source:</strong> The data is derived from literature collected in the
          <a href="https://www.genome.jp/kegg/" style="color: #00509d">KEGG</a> and
          <a href="http://signor.uniroma2.it/downloads.php" style="color: #00509d">SIGNOR </a> DataBase</p>
      </div>

      <!-- 搜索和控制区域卡片 -->
      <el-card class="search-card" style="margin-top: 20px; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px;">
        <!-- 基因搜索框 + 定位按钮 + 取消按钮 -->
        <div style="margin: 0 auto; display: flex; align-items: center; gap: 10px; width: 100%; padding: 10px; flex-wrap: wrap;">
          <span style="color: #333; font-weight: 500; white-space: nowrap;">Gene:</span>
          <el-input
              v-model="targetGene"
              placeholder="Enter gene symbol (e.g., BRCA1, TP53)"
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
          <!-- 取消按钮：仅在有选中基因时显示 -->
          <el-button
              type="warning"
              @click="handleCancel"
              v-if="targetGene.trim()"
              style="height: 40px; padding: 0 20px; font-size: 14px; white-space: nowrap;"
          >
            Cancel
          </el-button>
          <el-button
              type="success"
              @click="handleExport"
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
          <div v-for="(color, pathway) in pathwayColors" :key="pathway" style="display: flex; align-items: center; gap: 8px;">
            <div :style="{width: '16px', height: '16px', borderRadius: '50%', backgroundColor: color}"></div>
            <span class="text-sm">{{ pathway }}</span>
          </div>
        </div>
      </el-card>
    </div>

    <el-divider></el-divider>

    <!-- 下方视图区：图表 + 表格 -->
    <div class="flex-1 flex flex-col gap-4 p-0 md:p-4 overflow-hidden" style="padding-bottom: 10px">
      <!-- 图表容器 + 信息面板 -->
      <div class="flex-1 flex flex-col md:flex-row gap-4" style="width: 98%; margin: 0 auto">
        <!-- 图表容器 -->
        <div
            id="cross-chart-container"
            class="w-full md:w-3/4 border rounded overflow-hidden flex-1"
            :style="{ minHeight: '500px' }"
            v-if="hasData || isLoading"
        >
          <el-loading
              v-if="1"
              target="#cross-chart-container"
              text="Loading network..."
              background="rgba(255,255,255,0.8)"
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
              description="No cross-talk network data available"
              class="h-full flex items-center justify-center"
          />
        </div>
      </div>

      <el-divider></el-divider>

      <!-- 调控关系表格 -->
      <el-card class="flex-shrink-0 result-card" style="min-height: 400px; background-color: #fff; border: none; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08); border-radius: 10px; padding: 20px;">
        <div v-if="hasData" class="mb-2 font-bold text-[#00509d]">
          {{ targetGene ? `Cross-talk Relationships for ${targetGene}` : 'All Cross-talk Relationships' }}
          (Total {{ total }} pairs)
        </div>

        <el-table
            :data="paginatedData"
            v-if="hasData"
            border
            max-height="600px"
            style="width: 100%;"
            :header-cell-style="{ 'background-color': '#f0f5ff', 'color': '#003a8c', 'font-size': '14px' }"
            :row-style="{ 'font-size': '14px' }"
            v-loading="isLoadingTable"
            element-loading-text="Loading table data..."
        >
          <el-table-column prop="pathway" label="Source Pathway" align="center"/>
          <el-table-column prop="source" label="Source Gene" align="center"/>
          <el-table-column prop="targetPathway" label="Target Pathway" align="center"/>
          <el-table-column prop="target" label="Target Gene" align="center"/>
          <el-table-column prop="effect" label="Effect" align="center">
            <template #default="scope">
              <span :style="{
                color: scope.row.effect === 'Activation' ? '#28a745' :
                       scope.row.effect === 'Repression' ? '#dc3545' : '#6c757d',
                fontWeight: 'bold'
              }">
                {{ scope.row.effect }}
              </span>
            </template>
          </el-table-column>
          <el-table-column prop="pmid" label="PMID" align="center"/>
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
            description="No cross-talk relationships found"
            class="py-4"
        />
      </el-card>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted, onUnmounted, nextTick, watch } from 'vue';
import * as echarts from 'echarts';
import { ElMessage, ElIcon } from 'element-plus';
import { Download } from '@element-plus/icons-vue';
import crosstalkData from '@/assets/network/crosstalk_data.json'; // 本地数据

// 状态管理
const targetGene = ref('');
const isLoading = ref(false);
const isLoadingTable = ref(false);
const loadError = ref('');
const isExporting = ref(false);
const hasData = ref(false);
const networkStats = ref({
  totalGenes: 0,
  totalPathways: 0,
  totalInteractions: 0,
  activationCount: 0,
  repressionCount: 0,
  topGenes: []
});

// 分页状态
const currentPage = ref(1);
const pageSize = ref(10);
const total = ref(0);
const geneTableData = ref([]);
const paginatedData = ref([]);

// 网络数据存储
const regData = ref({
  nodes: [],
  links: []
});

// 图表实例
let chartInstance = null;

// 通路颜色
const pathwayColors = {
  "DNA repair": "#1f77b4",
  "p53 signaling": "#ff7f0e",
  "PI3K-AKT signaling": "#2ca02c",
  "TGF-beta signaling": "#d62728"
};

// 页面挂载时加载全量数据
onMounted(() => {
  loadFullNetwork();

});

// 加载全量网络数据
const loadFullNetwork = async () => {
  isLoading.value = true;
  isLoadingTable.value = true;
  loadError.value = '';
  hasData.value = false;

  try {
    // 处理所有调控关系边
    const regulationLinks = crosstalkData.edges
        .filter(l => l.type === 'regulation')
        .map(l => ({
          ...l,
          lineStyle: {
            width: 2,
            curveness: 0.2,
            color: l.effect === 'Activation' ? '#28a745' :
                l.effect === 'Repression' ? '#dc3545' : '#6c757d'
          }
        }));

    // 处理成员关系边
    const membershipLinks = crosstalkData.edges.filter(l => l.type === 'membership');

    // 提取所有节点
    const allNodes = crosstalkData.nodes.map(n => {
      let color = n.type === 'pathway'
          ? (pathwayColors[n.id] || '#28a745')
          : (n.pathway && pathwayColors[n.pathway[0]]) || '#ab47bc';

      return {
        name: n.id,
        type: n.type,
        pathway: n.pathway || [],
        symbolSize: n.type === 'pathway' ? 50 : 20,
        itemStyle: { color }
      };
    });

    // 保存网络数据
    regData.value = {
      nodes: allNodes,
      links: [...regulationLinks, ...membershipLinks]
    };
    hasData.value = true;

    // 计算网络统计信息
    const geneNodes = allNodes.filter(n => n.type === 'gene');
    const pathwayNodes = allNodes.filter(n => n.type === 'pathway');

    // 统计基因互作次数
    const geneInteractionCount = {};
    regulationLinks.forEach(link => {
      geneInteractionCount[link.source] = (geneInteractionCount[link.source] || 0) + 1;
      geneInteractionCount[link.target] = (geneInteractionCount[link.target] || 0) + 1;
    });

    // 获取互作最多的前5个基因
    const topGenes = Object.entries(geneInteractionCount)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 5)
        .map(item => item[0]);

    networkStats.value = {
      totalGenes: geneNodes.length,
      totalPathways: pathwayNodes.length,
      totalInteractions: regulationLinks.length,
      activationCount: regulationLinks.filter(l => l.effect === 'Activation').length,
      repressionCount: regulationLinks.filter(l => l.effect === 'Repression').length,
      topGenes
    };

    // 处理表格数据
    geneTableData.value = regulationLinks.map(l => {
      const sourceNode = crosstalkData.nodes.find(n => n.id === l.source);
      const targetNode = crosstalkData.nodes.find(n => n.id === l.target);
      return {
        pathway: sourceNode?.pathway?.[0] || '',
        targetPathway: targetNode?.pathway?.[0] || '',
        ...l
      };
    });
    total.value = geneTableData.value.length;
    updatePaginatedData();

    // 初始化图表
    await nextTick();
    initNetwork();

  } catch (error) {
    loadError.value = `Failed to load network: ${error.message}`;
    ElMessage.error(loadError.value);
  } finally {
    isLoading.value = false;
    isLoadingTable.value = false;
  }
};

// 定位特定基因
const handleGeneLocate = async () => {
  if (!targetGene.value.trim()) {
    ElMessage.warning('Please enter a gene symbol first');
    return;
  }

  const gene = targetGene.value.trim();
  isLoading.value = true;
  isLoadingTable.value = true;
  loadError.value = '';

  try {
    // 筛选包含目标基因的边
    const relevantLinks = crosstalkData.edges
        .filter(l => l.type === 'regulation' && (l.source === gene || l.target === gene));

    if (relevantLinks.length === 0) {
      throw new Error(`No cross-talk relationships found for ${gene}`);
    }

    // 提取相关节点
    const allNodes = new Set();
    relevantLinks.forEach(link => {
      allNodes.add(link.source);
      allNodes.add(link.target);
    });
    // 添加通路节点
    crosstalkData.nodes
        .filter(n => n.type === 'pathway')
        .forEach(n => allNodes.add(n.id));

    // 创建节点数据
    const nodes = Array.from(allNodes).map(name => {
      const originalNode = crosstalkData.nodes.find(n => n.id === name);
      if (!originalNode) return null;

      // 标记查询基因
      if (name === gene) {
        return {
          name,
          type: originalNode.type,
          pathway: originalNode.pathway || [],
          symbolSize: 22,
          itemStyle: { color: '#1e88e5' } // 查询基因特殊颜色
        };
      }

      // 普通节点
      let color = originalNode.type === 'pathway'
          ? (pathwayColors[originalNode.id] || '#28a745')
          : (originalNode.pathway && pathwayColors[originalNode.pathway[0]]) || '#ab47bc';

      return {
        name,
        type: originalNode.type,
        pathway: originalNode.pathway || [],
        symbolSize: originalNode.type === 'pathway' ? 50 : 20,
        itemStyle: { color }
      };
    }).filter(Boolean);

    // 创建边数据
    const regulationLinks = relevantLinks.map(l => ({
      ...l,
      lineStyle: {
        width: 2,
        curveness: 0.2,
        color: l.effect === 'Activation' ? '#28a745' :
            l.effect === 'Repression' ? '#dc3545' : '#6c757d'
      }
    }));

    // 添加成员关系边
    const membershipLinks = crosstalkData.edges
        .filter(l => l.type === 'membership' &&
            (allNodes.has(l.source) || allNodes.has(l.target)));

    regData.value = { nodes, links: [...regulationLinks, ...membershipLinks] };
    hasData.value = true;
    // 更新统计信息
    const relatedGenes = Array.from(allNodes)
        .filter(name => name !== gene && crosstalkData.nodes.find(n => n.id === name && n.type === 'gene'));

    networkStats.value = {
      ...networkStats.value,
      totalInteractions: regulationLinks.length,
      activationCount: regulationLinks.filter(l => l.effect === 'Activation').length,
      repressionCount: regulationLinks.filter(l => l.effect === 'Repression').length,
      topGenes: relatedGenes
    };

    // 处理表格数据
    geneTableData.value = regulationLinks.map(l => {
      const sourceNode = crosstalkData.nodes.find(n => n.id === l.source);
      const targetNode = crosstalkData.nodes.find(n => n.id === l.target);
      return {
        pathway: sourceNode?.pathway?.[0] || '',
        targetPathway: targetNode?.pathway?.[0] || '',
        ...l
      };
    });

    total.value = geneTableData.value.length;
    currentPage.value = 1; // 重置分页
    updatePaginatedData();
    // 初始化图表
    await nextTick();
    initNetwork();

  } catch (error) {
    loadError.value = error.message;
    ElMessage.error(loadError.value);
  } finally {
    isLoading.value = false;
    isLoadingTable.value = false;
  }
};

// 初始化网络图
const initNetwork = () => {
  const container = document.getElementById('cross-chart-container');
  if (!container) {
    ElMessage.error('Chart container not found');
    return;
  }

  // 销毁旧实例
  if (chartInstance) {
    chartInstance.dispose();
  }

  // 创建新实例
  chartInstance = echarts.init(container);

  // 图表配置
  const option = {
    animationDuration: 1500,
    animationEasingUpdate: 'quinticInOut',
    tooltip: {
      formatter: params => {
        if (params.dataType === 'node') {
          const pathways = params.data.pathway?.join(', ') || '';
          let typeLabel = params.data.type === 'pathway' ? 'Pathway' : 'Gene';
          if (params.data.name === targetGene.value) {
            typeLabel = 'Query Gene';
          }
          return `${params.data.name} (${typeLabel})${pathways ? `<br/>Pathways: ${pathways}` : ''}`;
        }
        return `${params.data.source} → ${params.data.target}: ${params.data.effect}`;
      }
    },
    series: [{
      type: 'graph',
      layout: 'force',
      force: { repulsion: 300, edgeLength: 150, gravity: 0.4 },
      roam: true,
      draggable: true,
      label: { show: true, fontSize: 12, color: '#333' },
      edgeSymbol: ['none', 'arrow'],
      edgeSymbolSize: [4, 10],
      data: regData.value.nodes,
      links: regData.value.links,
      emphasis: { focus: 'adjacency', lineStyle: { width: 5 } }
    }]
  };

  chartInstance.setOption(option);

  // 节点点击事件（切换基因）
  chartInstance.on('click', (params) => {
    if (params.dataType === 'node' && params.data.type === 'gene') {
      targetGene.value = params.data.name;
      handleGeneLocate();
    }
  });

  // 窗口大小适配
  const handleResize = () => chartInstance?.resize();
  window.addEventListener('resize', handleResize);
  onUnmounted(() => window.removeEventListener('resize', handleResize));
};

// 高亮关系
const highlightRelationship = (relation) => {
  if (!chartInstance || !regData.value.links.length || !regData.value.nodes.length) return;

  // 清除之前的高亮
  chartInstance.dispatchAction({ type: 'downplay', seriesIndex: 0 });

  // 高亮源节点
  const sourceIndex = regData.value.nodes.findIndex(node => node.name === relation.source);
  if (sourceIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'highlight',
      seriesIndex: 0,
      dataIndex: sourceIndex
    });
  }

  // 高亮目标节点
  const targetIndex = regData.value.nodes.findIndex(node => node.name === relation.target);
  if (targetIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'highlight',
      seriesIndex: 0,
      dataIndex: targetIndex
    });
  }

  // 高亮边
  const linkIndex = regData.value.links.findIndex(link =>
      link.source === relation.source &&
      link.target === relation.target &&
      link.effect === relation.effect
  );

  if (linkIndex !== -1) {
    chartInstance.dispatchAction({
      type: 'highlight',
      seriesIndex: 0,
      dataIndex: linkIndex
    });
  }
};

// 重新定位基因
const handleGeneRelocate = (geneName) => {
  targetGene.value = geneName;
  handleGeneLocate();
};

// 取消选中 - 重置为全量数据
const handleCancel = () => {
  targetGene.value = '';
  loadFullNetwork();
};

// 输入框清空处理 - 重置为全量数据
const handleInputClear = () => {
  if (targetGene.value.trim()) {
    handleCancel();
  }
};

// 导出数据
const handleExport = () => {
  if (!hasData.value) return;

  isExporting.value = true;
  try {
    const headers = ['pathway', 'source', 'target', 'effect', 'pmid'];
    const tsv = [headers.join('\t')].concat(
        geneTableData.value.map(r => [r.pathway, r.source, r.target, r.effect, r.pmid].join('\t'))
    ).join('\n');

    const fileName = targetGene.value
        ? `${targetGene.value.trim()}_cross_talk_data.tsv`
        : `full_cross_talk_network_data.tsv`;

    const blob = new Blob([tsv], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = fileName;
    a.click();
    URL.revokeObjectURL(url);

    ElMessage.success('Network data exported successfully');
  } catch (error) {
    ElMessage.error(`Export failed: ${error.message}`);
  } finally {
    isExporting.value = false;
  }
};

// 分页处理
const updatePaginatedData = () => {
  const start = (currentPage.value - 1) * pageSize.value;
  paginatedData.value = geneTableData.value.slice(start, start + pageSize.value);
};

const handleCurrentPageChange = (val) => {
  currentPage.value = val;
  updatePaginatedData();
};

const handlePageSizeChange = (val) => {
  pageSize.value = val;
  currentPage.value = 1;
  updatePaginatedData();
};

// 监听表格数据变化更新分页
watch(geneTableData, updatePaginatedData);

// 组件卸载时清理图表
onUnmounted(() => {
  if (chartInstance) {
    chartInstance.dispose();
    chartInstance = null;
  }
});
</script>

<style scoped>
#cross-chart-container {
  width: 100%;
  height: 100%;
}

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