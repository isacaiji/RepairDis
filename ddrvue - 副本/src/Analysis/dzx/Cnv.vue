<template>
  <div>
    <!-- 搜索表单 -->
    <div style="margin: 20px; text-align: center;display: flex;justify-content: center">
      <div style="display: flex;justify-content: center;align-items: center">
        <span style="font-size: 16px">Gene:&nbsp;&nbsp;</span>
      </div>
      <el-select
          v-model="gene"
          filterable
          remote
          placeholder="please enter gene symbol"
          :remote-method="remoteMethod"
          @change="getCnvData"
          style="width: 25%"
      >
        <el-option
            v-for="option in options"
            :key="option"
            :label="option"
            :value="option"
        />
      </el-select>
      <el-button type="primary" @click="getCnvData" style="margin-left: 10px">Click</el-button>
    </div>
    <hr />
    <div style="width: auto; height: 600px" id="cnv" class="chart" v-loading="loading"></div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue'
import * as echarts from 'echarts'
import { ElMessage } from 'element-plus'
import axios from 'axios'

// 状态
const allGenes = ref([])
const options = ref([])
const gene = ref('ATM')
const cnvList = ref([])
const loading = ref(false)
const baseURL = 'http://121.37.88.191:83/dzx'

// 远程搜索
const remoteMethod = (query) => {
  if (query !== '') {
    options.value = allGenes.value.filter((item) =>
        item.toLowerCase().includes(query.toLowerCase())
    )
  } else {
    options.value = []
  }
}

// 请求数据
const getCnvData = () => {
  if (!gene.value) return
  loading.value = true
  axios
      .get(`${baseURL}/cnv/${gene.value}`)
      .then((res) => {
        cnvList.value = res.data
        cnvInit()
      })
      .catch(() => {
      })
      .finally(() => {
        loading.value = false
      })
}

// 初始化 ECharts 图表
const cnvInit = () => {
  const data = []
  const cancers = []

  for (const cnv of cnvList.value) {
    data.push(cnv.dataList)
    cancers.push(cnv.cancer)
  }

  const chart = echarts.init(document.getElementById('cnv'))
  chart.setOption({
    dataset: [
      { source: data },
      {
        transform: {
          type: 'boxplot',
          config: {
            itemNameFormatter: function (param) {
              return cancers[param.value]
            }
          }
        }
      },
      {
        fromDatasetIndex: 1,
        fromTransformResult: 1
      }
    ],
    tooltip: {
      trigger: 'item',
      axisPointer: {
        type: 'shadow'
      }
    },
    grid: {
      left: '10%',
      right: '10%',
      bottom: '15%'
    },
    title: [
      {
        text: 'CNV',
        left: 'center'
      },
      {
        text: 'Cancer Type',
        left: 'center',
        top: '95%'
      }
    ],
    xAxis: {
      type: 'category',
      position: 'bottom',
      show: true,
      axisLabel: {
        interval: 0,
        rotate: 40
      },
      nameTextStyle: {
        color: '#000',
        fontSize: 16
      }
    },
    yAxis: {
      type: 'value',
      nameRotate: 90,
      nameGap: 30,
      nameLocation: 'middle',
      name: 'Copy Number Variance',
      nameTextStyle: {
        color: '#000',
        fontSize: 16
      }
    },
    series: [
      {
        name: 'boxplot',
        type: 'boxplot',
        datasetIndex: 1,
        width: '50%'
      },
      {
        name: 'outlier',
        type: 'scatter',
        datasetIndex: 2
      }
    ]
  })
}

// 初始获取基因数据
onMounted(() => {
  axios.get('http://121.37.88.191:9016/api/genes/all').then((res) => {
    allGenes.value = res.data
  })
})
</script>
