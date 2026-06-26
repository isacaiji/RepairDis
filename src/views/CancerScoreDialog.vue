<template>
  <el-dialog
      :model-value="modelValue"
      @update:model-value="$emit('update:modelValue', $event)"
      :title="null"
      width="900px"
      destroy-on-close
      @opened="renderChart"
      class="custom-academic-dialog"
  >
    <template #header>
      <div class="custom-dialog-header">
        <div class="title-main">{{ geneName }} <span class="title-suffix">Distribution</span></div>
        <div class="title-sub">Cancer-specific gene-level MO-DDRweight across TCGA cancer types</div>
      </div>
    </template>

    <div class="dialog-body-layout">
      <div v-if="scoreData.length > 0" ref="chartRef" class="echarts-viewport"></div>
      <el-empty v-else description="No cancer score data available" :image-size="120" />

      <div v-if="scoreData.length > 0" class="chart-legend-box">
        <div class="legend-note">Values are displayed as MO_DDRweight x 100; the gene detail page shows the mean across available cancer types.</div>
        <div class="legend-item"><span class="dot danger"></span> High MO-DDRweight (&gt;=80)</div>
        <div class="legend-item"><span class="dot warning"></span> Moderate MO-DDRweight (50-80)</div>
        <div class="legend-item"><span class="dot success"></span> Low MO-DDRweight (&lt;50)</div>
      </div>
    </div>

    <template #footer>
      <div class="dialog-footer-minimal">
        <el-button @click="$emit('update:modelValue', false)" class="btn-close">Dismiss</el-button>
      </div>
    </template>
  </el-dialog>
</template>

<script setup>
import { ref, onUnmounted } from 'vue';
import * as echarts from 'echarts';

const props = defineProps({
  modelValue: Boolean,
  geneName: String,
  scoreData: {
    type: Array,
    default: () => []
  }
});

const emit = defineEmits(['update:modelValue']);
const chartRef = ref(null);
let myChart = null;

const renderChart = () => {
  if (!chartRef.value) return;
  myChart = echarts.init(chartRef.value);

  const displayData = props.scoreData
      .map(item => ({
        name: item.name,
        value: Number(item.value)
      }))
      .filter(item => Number.isFinite(item.value))
      .sort((a, b) => b.value - a.value);

  const option = {
    backgroundColor: 'transparent',
    tooltip: {
      trigger: 'axis',
      backgroundColor: 'rgba(255, 255, 255, 0.98)',
      borderColor: '#e2e8f0',
      borderWidth: 1,
      textStyle: { color: '#1e293b', fontSize: 12 },
      shadowColor: 'rgba(0, 0, 0, 0.05)',
      shadowBlur: 10,
      formatter: (params) => {
        const data = params[0];
        return `<div style="padding:4px 8px">
                  <div style="font-weight:bold;margin-bottom:4px">${data.name}</div>
                  <div style="color:${data.color.colorStops[0].color}">MO-DDRweight x 100: ${Number(data.value).toFixed(2)}</div>
                </div>`;
      }
    },
    grid: { left: '4%', right: '4%', bottom: '15%', top: '10%', containLabel: true },
    // 增加数据缩放滑块，方便查看 33 个癌种
    dataZoom: [{
      type: 'inside',
      start: 0,
      end: 100
    }, {
      type: 'slider',
      height: 18,
      bottom: 10,
      borderColor: 'transparent',
      backgroundColor: '#f1f5f9',
      fillerColor: 'rgba(56, 189, 248, 0.1)',
      handleStyle: { color: '#cbd5e1' },
      textStyle: { color: 'transparent' }
    }],
    xAxis: {
      type: 'category',
      data: displayData.map(d => d.name),
      axisLine: { lineStyle: { color: '#cbd5e1' } },
      axisLabel: {
        interval: 0,
        rotate: 45,
        fontSize: 11,
        color: '#64748b',
        fontWeight: 600
      },
      axisTick: { show: false }
    },
    yAxis: {
      type: 'value',
      min: 0,
      max: 100,
      axisLine: { show: false },
      axisLabel: { color: '#94a3b8', fontSize: 11 },
      splitLine: { lineStyle: { type: 'dashed', color: '#f1f5f9' } }
    },
    series: [{
      name: 'MO-DDRweight',
      type: 'bar',
      barWidth: '55%',
      itemStyle: {
        borderRadius: [4, 4, 0, 0],
        // 增加渐变色质感
        color: new echarts.graphic.LinearGradient(0, 0, 0, 1, [
          { offset: 0, color: '' }, // 动态填充
          { offset: 1, color: '' }
        ])
      },
      // 这里的根据数值动态上色逻辑进行了精修
      data: displayData.map(item => {
        let colorTop, colorBottom;
        if (item.value >= 80) { colorTop = '#ef4444'; colorBottom = '#fca5a5'; }
        else if (item.value >= 50) { colorTop = '#f59e0b'; colorBottom = '#fcd34d'; }
        else { colorTop = '#10b981'; colorBottom = '#6ee7b7'; }

        return {
          value: Number(item.value),
          itemStyle: {
            color: new echarts.graphic.LinearGradient(0, 0, 0, 1, [
              { offset: 0, color: colorTop },
              { offset: 1, color: colorBottom }
            ])
          }
        };
      }),
      emphasis: {
        itemStyle: { opacity: 0.8, shadowBlur: 10, shadowColor: 'rgba(0,0,0,0.1)' }
      }
    }]
  };

  myChart.setOption(option);
  window.addEventListener('resize', handleResize);
};

const handleResize = () => myChart && myChart.resize();

onUnmounted(() => {
  window.removeEventListener('resize', handleResize);
  if (myChart) myChart.dispose();
});
</script>

<style scoped>
/* 弹窗整体风格微调 */
:deep(.custom-academic-dialog) {
  border-radius: 16px;
  overflow: hidden;
  box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.15);
}

.custom-dialog-header {
  padding: 20px 0 10px;
}

.title-main {
  font-size: 24px;
  font-weight: 800;
  color: #0f172a;
  letter-spacing: -0.5px;
}

.title-suffix {
  color: #94a3b8;
  font-weight: 400;
}

.title-sub {
  font-size: 13px;
  color: #64748b;
  margin-top: 4px;
}

.dialog-body-layout {
  padding: 0 10px;
}

.echarts-viewport {
  width: 100%;
  height: 480px;
  margin-top: 10px;
}

/* 底部图例 */
.chart-legend-box {
  display: flex;
  flex-wrap: wrap;
  justify-content: center;
  gap: 30px;
  margin-top: 15px;
  padding: 12px;
  background: #f8fafc;
  border-radius: 8px;
}

.legend-note {
  flex-basis: 100%;
  text-align: center;
  font-size: 12px;
  color: #64748b;
}

.legend-item {
  font-size: 12px;
  font-weight: 700;
  color: #475569;
  display: flex;
  align-items: center;
  gap: 8px;
}

.dot { width: 8px; height: 8px; border-radius: 50%; }
.dot.danger { background: #ef4444; }
.dot.warning { background: #f59e0b; }
.dot.success { background: #10b981; }

.dialog-footer-minimal {
  padding-top: 10px;
}

.btn-close {
  border-radius: 8px;
  padding: 8px 20px;
  font-weight: 600;
}
</style>
