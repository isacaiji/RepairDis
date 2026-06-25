<template>
  <div>
    <!--搜索表单-->
    <div style="margin: 20px; text-align: center;display: flex;justify-content: center">
      <div style="display: flex;justify-content: center;align-items: center">
        <span style="font-size: 18px">Gene:&nbsp;&nbsp;</span>
      </div>

      <el-select
          v-model="gene"
          filterable
          remote
          placeholder="please enter gene symbol"
          :remote-method="remoteMethod"
          @change="logfcAndMut"
          style="width: 25%"
      >
        <el-option
            v-for="option in options"
            :key="option"
            :label="option"
            :value="option"
        />
      </el-select>
      <el-button type="primary" @click="logfcAndMut" style="margin-left: 10px">Click</el-button>
    </div>
    <hr />
    <div style="width: auto; height: 600px" id="mut" class="chart"></div>
  </div>
</template>

<script setup>
import { ref, onMounted, nextTick } from 'vue'
import * as echarts from 'echarts'
import axios from 'axios'

// 基因搜索相关
const gene = ref('ATM')
const allGenes = ref([])
const options = ref([])

// 突变数据
const logfcAndMutData = ref([])

// 模糊搜索函数
const remoteMethod = (query) => {
  if (query !== '') {
    options.value = allGenes.value.filter(item =>
        item.toLowerCase().includes(query.toLowerCase())
    )
  } else {
    options.value = []
  }
}

// 获取图表数据并初始化
const logfcAndMut = async () => {
  const { data } = await axios.get(`http://121.37.88.191:83/dzx/mut/${gene.value}`)
  logfcAndMutData.value = data
  await nextTick()
  mutInit()
}

// ECharts 初始化图表
const mutInit = () => {
  const cancers = ['Cancer']
  const muts = ['Mutation']
  const nomutes = ['No Mutation']

  logfcAndMutData.value.forEach(item => {
    cancers.push(item.cancer)
    muts.push(item.mut)
    nomutes.push(1 - item.mut)
  })

  const sets = []
  const title = []

  let x = 8, y = 16
  let a = 7, b = 25

  for (let i = 1; i < cancers.length; i++) {
    sets.push({
      type: 'pie',
      radius: '11%',
      center: [`${x}%`, `${y}%`],
      encode: {
        itemName: 'Cancer',
        value: cancers[i]
      }
    })
    title.push({
      subtext: cancers[i],
      left: `${a}%`,
      top: `${b}%`,
      textAlign: 'center'
    })

    x += 12
    a += 12
    if (i % 8 === 0) {
      x = 8
      y += 22
      a = 8
      b += 22
    }
  }

  echarts.init(document.getElementById('mut')).setOption({
    title,
    legend: {},
    tooltip: {},
    label: {
      position: 'outer',
      alignTo: 'none',
      bleedMargin: 5
    },
    dataset: {
      source: [cancers, muts, nomutes]
    },
    series: sets
  })
}

// 初始化获取基因列表
onMounted(async () => {
  const { data } = await axios.get('http://121.37.88.191:8989/rnamodifications/gene')
  allGenes.value = data
})
</script>

<style scoped>
.chart {
  margin: auto;
}
</style>
