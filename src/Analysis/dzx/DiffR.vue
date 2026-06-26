<template>
  <div style="margin: auto; text-align: center">
    <!-- 搜索表单 -->
    <div style="margin: 20px;">
      <el-form :inline="true" class="demo-form-inline">
        <el-form-item label="All Gene">
          <el-switch
              v-model="allGene"
              active-color="#13ce66"
              inactive-color="#ff4949">
          </el-switch>
        </el-form-item>

        <!-- 输入框或远程选择框 -->
        <el-form-item label="Gene" v-if="allGene">
          <el-input v-model="gene" placeholder="please enter gene symbol" />
        </el-form-item>

        <el-form-item label="Gene" v-else>
          <el-select
              v-model="gene"
              filterable
              remote
              placeholder="please enter gene symbol"
              :remote-method="remoteMethod"
              :loading="loading">
            <el-option
                v-for="option in options"
                :key="option"
                :label="option"
                :value="option" />
          </el-select>
        </el-form-item>

        <el-form-item>
          <el-button type="primary" @click="getDiffPic" >Click</el-button>
        </el-form-item>

      </el-form>
    </div>

    <hr />

    <el-empty description="Wait For Your Click" v-if="!show" style="height: 600px" />
    <div v-loading="loading" v-else>
      <!-- 图片展示 -->
      <el-image
          style="height: 600px; width: 1100px"
          :src="'data:image/png;base64,' + img64"
          fit="contain"
      />
      <br />
      <el-button type="primary" round @click="downloadDiffImg">Download</el-button>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue'
import axios from 'axios'

// 数据状态
const show = ref(false)
const gene = ref('ATM')
const allGene = ref(false)
const loading = ref(false)
const allGenes = ref([])
const options = ref([])
const img64 = ref('')
const diffList = ref([])

const baseAPI = 'http://121.37.88.191:83'
const api = {
  rURL: baseAPI + '/r',
  dzxURL: baseAPI + '/dzx'
}

// 模糊搜索方法
const remoteMethod = (query) => {
  if (query !== '') {
    options.value = allGenes.value.filter((item) =>
        item.toLowerCase().includes(query.toLowerCase())
    )
  } else {
    options.value = []
  }
}

// 获取图片数据
const getDiffPic = async () => {
  show.value = true
  loading.value = true
  img64.value = ''
  try {
    const resImg = await axios.get(`${api.rURL}/diffpic/${gene.value}`)
    img64.value = resImg.data
    const resData = await axios.get(`${api.dzxURL}/diff/${gene.value}`)
    diffList.value = resData.data
  } catch (error) {
    console.error(error)
  } finally {
    loading.value = false
  }
}

// 下载图像
const downloadDiffImg = () => {
  const a = document.createElement('a')
  a.href = `${api.dzxURL}/diff/down/img/${gene.value}`
  a.click()
}

// 下载数据
const downloadDiffData = () => {
  const a = document.createElement('a')
  a.href = `${api.dzxURL}/diff/down/data/${gene.value}`
  a.click()
}

// 初始化基因数据
onMounted(async () => {
  try {
    const res = await axios.get(`${api.dzxURL}/api/genes/all`)
    allGenes.value = res.data
  } catch (err) {
    console.error(err)
  }
})
</script>
