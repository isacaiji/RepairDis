<template>
  <div>
    <div style="margin: 20px;">
      <el-form :inline="true" class="demo-form-inline">
        <span style="font-size: 16px">Gene:&nbsp;&nbsp;</span>
        <el-select
            v-model="gene"
            filterable
            remote
            placeholder="please enter gene symbol"
            :remote-method="remoteMethod"
            :loading="loading"
            style="width: 300px;margin-right: 10px"
        >
          <el-option
              v-for="option in options"
              :label="option"
              :value="option"
              :key="option"
          />
        </el-select>
        <el-button type="primary" @click="getReceptorPic">Click</el-button>
      </el-form>
    </div>

    <hr>

    <el-empty description="Wait For Your Click" v-show="!show" style="height: 600px"></el-empty>

    <div v-loading="loading" v-show="show">
      <div class="checkpoint-left" style="background-color: #ffffff; height: 600px;">
        <el-image
            style="height: 450px; margin-top: 50px"
            :src="'data:image/png/;base64,' + png64"
            fit="contain"
        />
        <br>
        <el-button type="primary" round @click="downloadDiffImg">Download</el-button>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue'
import axios from 'axios'
import jsPDF from "jspdf";

// 基础URL设置
const baseURL = 'http://121.37.88.191:8989'

// 响应式数据
const png64 = ref("")
const show = ref(false)
const allGenes = ref([])
const options = ref([])
const receptor = ref([])
const gene = ref('ATM')
const loading = ref(false)

// 模糊搜索
const remoteMethod = (query) => {
  if (query !== "") {
    options.value = allGenes.value.filter((item) =>
        item.toLowerCase().indexOf(query.toLowerCase()) > -1
    )
  } else {
    options.value = []
  }
}

// 获取受体图片
const getReceptorPic = () => {
  show.value = true
  loading.value = true
  png64.value = ""

  axios.get(`${baseURL}/r/immu/receptor/${gene.value}`)
      .then(imgRes => {
        png64.value = imgRes.data
      })
      .catch(err => {
        console.error('Error fetching data:', err)
      })
      .finally(() => {
        loading.value = false
      })
}

// 下载图片
const downloadDiffImg = () => {
  let doc = new jsPDF({
    orientation: 'landscape',
    unit: 'px',
    format: [1886, 1028]
  })

  doc.addImage('data:image/png;base64,' + png64.value, 'PNG', 0, 0)
  doc.save(`${gene.value}_receptor_image.pdf`)
}

// 初始化获取基因列表
onMounted(() => {
  axios.get(`http://121.37.88.191:9016/api/genes/all`)
      .then(res => {
        allGenes.value = res.data
      })
      .catch(err => {
        console.error('Error fetching genes:', err)
      })
})
</script>

<style scoped>
</style>