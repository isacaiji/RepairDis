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
        <el-button type="primary" @click="getCheckpointPic">Click</el-button>
      </el-form>
    </div>

    <hr>

    <!-- 1. 增大空状态容器高度，和图片容器匹配 -->
    <el-empty description="Wait For Your Click" v-show="!show" style="height: 800px"></el-empty>

    <div v-loading="loading" v-show="show">
      <!-- 2. 增大图片容器高度，宽度铺满 -->
      <div class="checkpoint-left" style="background-color: #ffffff; height: 800px; width: 100%;">
        <!-- 3. 调整图片尺寸，增加预览功能 -->
        <el-image
            style="height: 700px; width: 100%; margin-top: 20px;"
            :src="'data:image/png/;base64,' + png64"
            fit="contain"
            :preview-src-list="['data:image/png/;base64,' + png64]"
        />
        <br>
        <el-button type="primary" round @click="downloadDiffImg">Download</el-button>
      </div>
    </div>
  </div>
</template>

<script setup>
import {ref, onMounted} from 'vue'
import axios from 'axios'
import jsPDF from "jspdf";

// 基础URL设置
const baseURL = 'http://121.37.88.191:8989'
// const baseURL = 'http://121.37.88.191:83/dzx/cnv/'

// 响应式数据
const png64 = ref("")
const show = ref(false)
const allGenes = ref([])
const options = ref([])
const checkpoint = ref([])
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

// 获取检查点图片
const getCheckpointPic = () => {
  show.value = true
  loading.value = true
  png64.value = ""

  axios.get(`${baseURL}/r/immu/checkpoint/${gene.value}`)
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

// 下载图片（优化：根据图片真实尺寸生成PDF，避免拉伸/压缩）
const downloadDiffImg = () => {
  // 创建临时图片对象获取真实尺寸
  const img = new Image()
  img.src = 'data:image/png;base64,' + png64.value

  img.onload = () => {
    // 根据图片尺寸自动设置PDF方向和大小
    const isLandscape = img.width > img.height
    const doc = new jsPDF({
      orientation: isLandscape ? 'landscape' : 'portrait',
      unit: 'px',
      format: [img.width, img.height]
    })

    // 按真实尺寸添加图片，避免变形
    doc.addImage(img, 'PNG', 0, 0, img.width, img.height)
    doc.save(`${gene.value}_checkpoint_image.pdf`)
  }
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
/* 优化图片容器样式，防止溢出 */
.checkpoint-left {
  overflow: auto;
}

/* 强制保持图片比例，避免拉伸变形 */
.el-image {
  object-fit: contain !important;
}
</style>