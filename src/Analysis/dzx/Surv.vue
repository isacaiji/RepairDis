<template>
  <div style="margin: auto; text-align: center">
    <div style="margin: 20px;">
      <el-form :inline="true" class="demo-form-inline">
        <el-form-item label="All Gene">
          <el-switch v-model="allGene" active-color="#13ce66" inactive-color="#ff4949" />
        </el-form-item>

        <el-form-item label="Cancer">
          <el-select v-model="cancer" filterable placeholder="Cancer">
            <el-option v-for="item in cancers" :key="item" :label="item" :value="item" />
          </el-select>
        </el-form-item>

        <!-- All gene input -->
        <el-form-item label="Gene" v-if="allGene">
          <el-input v-model="gene" placeholder="Please enter gene symbol" />
        </el-form-item>

        <!-- Longevity gene select -->
        <el-form-item label="Gene" v-else>
          <el-select
              v-model="gene"
              filterable
              remote
              placeholder="Please enter gene symbol"
              :remote-method="remoteMethod"
              :loading="loading"
          >
            <el-option v-for="option in options" :key="option" :label="option" :value="option" />
          </el-select>
        </el-form-item>

        <el-form-item>
          <el-button type="primary" @click="getPic" style="margin-left: 10px">Click</el-button>
        </el-form-item>

      </el-form>
    </div>

    <hr />
    <el-empty description="Wait For Your Click" v-if="!show" style="height: 600px" />
    <div v-loading="loading" v-else>
      <el-image
          style="height: 64vh"
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

const gene = ref('ATM')
const cancer = ref('BLCA')
const img64 = ref('')
const show = ref(false)
const loading = ref(false)
const allGene = ref(false)
const options = ref([])
const allGenes = ref([])
const cancers = ref([])

const remoteMethod = (query) => {
  if (query !== '') {
    options.value = allGenes.value.filter((item) =>
        item.toLowerCase().includes(query.toLowerCase())
    )
  } else {
    options.value = []
  }
}

const getPic = async () => {
  show.value = true
  loading.value = true
  img64.value = ''
  try {
    const { data } = await axios.get(`http://121.37.88.191:83/r/surv/${gene.value}/${cancer.value}`)
    img64.value = data
  } finally {
    loading.value = false
  }
}

const downloadImg = () => {
  const a = document.createElement('a')
  a.href = `http://121.37.88.191:83/dzx/down/surv/${gene.value}/${cancer.value}`
  a.setAttribute('download', '')
  a.click()
}

onMounted(async () => {
  // 获取基因列表
  const { data } = await axios.get('http://121.37.88.191:9016/api/genes/all')
  allGenes.value = data

  cancers.value = [
    'BLCA', 'BRCA', 'CHOL', 'COAD', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'KIRP',
    'LIHC', 'LUAD', 'LUSC', 'PRAD', 'READ', 'STAD', 'THCA', 'UCEC'
  ]
})
</script>

<style scoped>
</style>
