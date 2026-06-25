<template>
  <div id="main-container" style="min-height: 750px;width: 100%;background-color: #f9f9f9">
    <!-- 标题区域 -->
    <div style="min-height: 50px;background-color: #f9f9f9;display: flex;align-items: center;">
      <div class="font-container" style="display: flex;flex: 1;">
        <span style="font-size: 25px;font-weight: 600;color: #003f88;margin-left: 35px">{{selectedGene}}</span>
        <span style="font-size: 25px;font-weight: 600;margin-left: 10px">Evolution Analysis Result</span>
      </div>
    </div>
    <el-divider style="margin-top: -10px;background: linear-gradient(90deg, rgba(0,63,136,0) 0%, rgba(0,63,136,0.3) 50%, rgba(0,63,136,0) 100%);height: 1px;"></el-divider>

    <!-- 主要内容区域 -->
    <div class="selected-content" style="max-width: 1200px;margin: 0 auto;padding: 20px;">
      <!-- 控制区域 -->
      <div class="content-wrapper" style="padding: 25px;border-radius: 8px;">
        <div class="input-group">
          <el-switch
              v-model="cr"
              class="switch-cr"
              active-value="C"
              inactive-value="R"
              style="--el-switch-on-color: #003f88; --el-switch-off-color: #003f88;"
              size="large"
              inline-prompt
              active-text="Circular"
              inactive-text="Rectangular"
              @change = "getResult"
          ></el-switch>

          <el-select
              v-model="selectedGene"
              placeholder="Select a gene"
              class="select-style"
              style="width: 220px;"
          >
            <el-option
                v-for="gene in geneList"
                :key="gene"
                :label="gene"
                :value="gene"
            ></el-option>
          </el-select>

          <el-button @click="getResult" class="el-button--primary button-style">Get Result</el-button>
          <el-button @click="download" class="el-button--success button-style">Download</el-button>
        </div>
      </div>

      <!-- 图片展示区域 - 限制最大显示比例 -->
      <div v-if="imageUrl" class="image-container content-wrapper">
        <div class="image-scroll-wrapper">
          <img
              :src="imageUrl"
              alt="Tree Image"
              class="responsive-image"
              :style="{
              maxWidth: '80%',  // 缩小到原来的80%
              maxHeight: '600px',  // 限制最大高度
              objectFit: 'contain',
              margin: '0 auto',  // 居中显示
              display: 'block'
            }"
          />
        </div>
      </div>

      <!-- 错误信息展示区域 -->
      <div v-if="errorMessage" class="error-area">
        <p>{{ errorMessage }}</p>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import {ref, onMounted} from 'vue';
import axios from 'axios';

// 切换状态 Circle或者Rectangle   active是C,inactive是R
const cr = ref("C");

// 存储选中的基因
const selectedGene = ref('ATM');
// 存储图片的 URL
const imageUrl = ref('');
// 存储错误信息
const errorMessage = ref('');

// 定义基因列表
const geneList = ref();

const treeListUrl = "http://121.37.88.191:9016/api/genes/all";
const r = axios.get(treeListUrl).then((res)=>{
  geneList.value = res.data;
})

// 搜索函数
const getResult = async () => {
  if (selectedGene.value) {
    try {
      // 构建图片的 URL
      const url = `http://121.37.88.191:9016/analysis/tree/${selectedGene.value}-${cr.value}`;
      // 发送请求检查图片是否存在
      const response = await axios.get(url, { responseType: 'blob' });
      if (response.status === 200) {
        // 创建临时 URL 用于显示图片
        imageUrl.value = URL.createObjectURL(response.data);
        errorMessage.value = '';
      } else {
        imageUrl.value = '';
        errorMessage.value = 'Image not found. Please try another gene.';
      }
    } catch (error) {
      imageUrl.value = '';
      errorMessage.value = 'Image not found. Please try another gene.';
    }
  } else {
    errorMessage.value = 'Please select a gene.';
  }
};

// 下载函数
const download = async () => {
  if (selectedGene.value) {
    try {
      // 构建图片的 URL
      const url = `http://121.37.88.191:9016/analysis/tree/${selectedGene.value}-${cr.value}`;
      // 发送请求获取图片数据
      const response = await axios.get(url, { responseType: 'blob' });
      if (response.status === 200) {
        // 创建 Blob 对象
        const blob = new Blob([response.data], { type: response.headers['content-type'] });
        // 创建下载链接
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = `${selectedGene.value}.png`;
        // 触发下载
        link.click();
        // 释放 URL 对象
        URL.revokeObjectURL(link.href);
        errorMessage.value = '';
      } else {
        errorMessage.value = 'Image not found. Cannot download.';
      }
    } catch (error) {
      errorMessage.value = 'Error downloading image. Please try again.';
    }
  } else {
    errorMessage.value = 'Please select a gene first.';
  }
};

onMounted(()=>{
  getResult();
})
</script>

<style scoped>
/* 原有样式保持不变 */
.input-group {
  display: flex;
  justify-content: center;
  align-items: center;
  gap: 15px;
  margin-bottom: 0;
  margin-top: 0;
  width: 100%;
  margin-left: 0;
}

.select-style {
  border-radius: 4px;
  font-size: 16px;
  outline: none;
  transition: all 0.3s ease;
  border-color: #e2e8f0;
}

.select-style:focus {
  border-color: #003f88;
  box-shadow: 0 0 0 3px rgba(0, 63, 136, 0.1);
}

.button-style {
  padding: 8px 20px;
  font-size: 16px;
  border-radius: 4px;
  cursor: pointer;
  transition: all 0.2s ease;
}

::v-deep .el-button--primary {
  background-color: #003f88;
  border-color: #003f88;
}

::v-deep .el-button--primary:hover {
  background-color: #002855;
  border-color: #002855;
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
}

::v-deep .el-button--success {
  background-color: #1150bd;
  border-color: #1150bd;
}

::v-deep .el-button--success:hover {
  background-color: #003a8c;
  border-color: #003a8c;
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
}

/* 图片容器优化 */
.image-container {
  margin: 25px auto;
  padding: 30px;
  text-align: center;
}

/* 滚动容器 - 当图片过大时出现滚动条 */
.image-scroll-wrapper {
  overflow: auto;
  max-height: 700px; /* 固定最大高度 */
  padding: 10px;
  margin: 0 auto;
  border: 1px solid #e2e8f0;
  border-radius: 4px;
  background-color: #f9f9f9;
}

.responsive-image {
  transition: all 0.3s ease;
}

.error-area {
  color: #dc2626;
  margin: 150px auto;
  text-align: center;
  font-size: 16px;
}

.content-wrapper {
  background-color: #ffffff;
  border-radius: 8px;
  box-shadow: 0 3px 15px rgba(0, 0, 0, 0.07);
  margin-bottom: 25px;
}
</style>