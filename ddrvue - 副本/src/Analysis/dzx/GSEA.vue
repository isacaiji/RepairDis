<template>
  <div style="margin: auto; text-align: center">
    <!-- 搜索表单 -->
    <div style="margin: 20px;">
      <el-form :inline="true" class="demo-form-inline">
        <el-form-item label="All Gene">
          <el-switch
              v-model="allGene"
              active-color="#13ce66"
              inactive-color="#ff4949"
          />
        </el-form-item>

        <el-form-item label="Cancer">
          <el-select v-model="cancer" filterable placeholder="Cancer">
            <el-option
                v-for="item in cancers"
                :key="item"
                :label="item"
                :value="item"
            />
          </el-select>
        </el-form-item>

        <!-- all gene -->
        <el-form-item label="Gene" v-if="allGene">
          <el-input v-model="gene" placeholder="please enter gene symbol" />
        </el-form-item>

        <!-- longevity gene -->
        <el-form-item label="Gene" v-else>
          <el-input v-model="gene" placeholder="please enter gene symbol" />
        </el-form-item>

        <el-form-item>
          <el-button type="primary" @click="getPic">Click</el-button>
        </el-form-item>
      </el-form>
    </div>

    <hr />

    <!-- 结果展示 -->
    <el-empty description="Wait For Your Click" v-if="!show" style="height: 600px" />
    <div v-loading="loading" v-else>
      <el-image
          style="height: 600px"
          :src="'data:image/png;base64,' + img64"
          fit="contain"
      />
      <br />
      <el-button type="primary" round @click="downloadImg">Download</el-button>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue'
import axios from 'axios'

// API 配置
// const baseURL = 'http://121.37.88.191:9016'
const baseURL = 'http://121.37.88.191:8989'
const api = {
  rURL: baseURL + '/r',
  dzxURL: baseURL + 'dzx'
}

// 响应式数据
const show = ref(false)
const allGene = ref(false)
const gene = ref('ATM')
const cancer = ref('BLCA')
const img64 = ref('')
const loading = ref(false)
const allGenes = ref([])
const options = ref([])

const cancers = [
  'BLCA', 'BRCA', 'CHOL', 'COAD', 'ESCA',
  'HNSC', 'KICH', 'KIRC', 'KIRP', 'LIHC',
  'LUAD', 'LUSC', 'PRAD', 'READ', 'STAD',
  'THCA', 'UCEC'
]

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

// 获取图像
const getPic = async () => {
  show.value = true
  loading.value = true
  img64.value = ''
  try {
    const res = await axios.get(`${api.rURL}/gsea/${gene.value}/${cancer.value}`)
    img64.value = res.data
  } catch (err) {
    console.error(err)
  } finally {
    loading.value = false
  }
}

// 下载图像
const downloadImg = () => {
  const a = document.createElement('a')
  a.href = `${api.dzxURL}/down/gsea/${gene.value}/${cancer.value}`
  a.click()
}

// 初始化获取基因列表
onMounted(async () => {
  try {
    const res = await axios.get(`${api.rURL}/gene`)
    allGenes.value = res.data
  } catch (err) {
    console.error(err)
  }
})
</script>
