<template>
  <div>
    <!-- 搜索表单 -->
    <div style="margin: 20px; text-align: center;display: flex;justify-content: center">
      <div style="display: flex;justify-content: center;align-items: center">
        <span style="font-size: 18px">Gene: &nbsp;&nbsp;</span>
      </div>
      <el-select
          v-model="gene"
          filterable
          remote
          placeholder="please enter gene symbol"
          :remote-method="remoteMethod"
          :loading="loading"
          @change="getMethyByGene"
          style="width: 25%"
      >
        <el-option
            v-for="option in options"
            :key="option"
            :label="option"
            :value="option"
        />
      </el-select>
      <el-button type="primary" @click="getMethyByGene" style="margin-left: 20px">Click</el-button>
    </div>

    <hr />
    <div style="width: auto; height: 600px" id="methy" class="chart"></div>
  </div>
</template>

<script setup>
import { ref, onMounted, watch } from 'vue'
import axios from 'axios'
import * as echarts from 'echarts'

// API 配置
const baseURL = 'http://121.37.88.191:83'
const geneListURL = `$http://121.37.88.191:9016/api/genes/all`
const methyURL = `${baseURL}/dzx/methy`

// 响应式变量
const gene = ref('ATM')
const options = ref([])
const allGenes = ref([])
const loading = ref(false)
const methyList = ref([])

// 模糊搜索基因
const remoteMethod = (query) => {
  if (query !== '') {
    options.value = allGenes.value.filter((item) =>
        item.toLowerCase().includes(query.toLowerCase())
    )
  } else {
    options.value = []
  }
}

// 获取甲基化数据并渲染图表
const getMethyByGene = async () => {
  try {
    loading.value = true
    const res = await axios.get(`${methyURL}/${gene.value}`)
    methyList.value = res.data
    drawBoxplot()
  } catch (err) {
    console.error(err)
  } finally {
    loading.value = false
  }
}

// 初始化 ECharts 图表
const drawBoxplot = () => {
  const chartDom = document.getElementById('methy')
  const myChart = echarts.init(chartDom)
  const data = methyList.value.map(item => item.dataList)
  const cancers = methyList.value.map(item => item.cancer)

  myChart.setOption({
    dataset: [
      { source: data },
      {
        transform: {
          type: 'boxplot',
          config: {
            itemNameFormatter: param => cancers[param.value]
          }
        }
      },
      {
        fromDatasetIndex: 1,
        fromTransformResult: 1
      }
    ],
    title: {
      text: 'DNA Methylation',
      left: 'center'
    },
    tooltip: {
      trigger: 'item',
      axisPointer: {
        type: 'shadow'
      }
    },
    grid: {
      left: '10%',
      right: '10%',
      bottom: '10%'
    },
    xAxis: {
      type: 'category',
      name: 'Cancer Type',
      nameGap: 40,
      nameLocation: 'middle',
      nameTextStyle: {
        color: '#000',
        fontSize: 16
      },
      axisLabel: {
        interval: 0,
        rotate: 40
      }
    },
    yAxis: {
      type: 'value',
      name: 'Methylation Beta Value',
      nameRotate: 90,
      nameGap: 40,
      nameLocation: 'middle',
      nameTextStyle: {
        color: '#000',
        fontSize: 16
      },
      splitArea: {
        show: true
      }
    },
    series: [
      {
        name: 'boxplot',
        type: 'boxplot',
        datasetIndex: 1
      },
      {
        name: 'outlier',
        type: 'scatter',
        datasetIndex: 2
      }
    ]
  })
}

// 初次挂载时获取基因列表
onMounted(async () => {
  try {
    const res = await axios.get(geneListURL)
    allGenes.value = res.data
  } catch (err) {
    console.error(err)
  }
})
</script>

<style scoped>
.chart {
  width: 100%;
}
</style>
