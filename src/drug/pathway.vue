<template>
  <div id="pathway">
    <div class="pathway-title">
      <!-- 药物靶点通路模块标题简介 -->
    </div>
    <div class="pathway-select">
      <!-- 搜索框 -->
      <div style="text-align: left;height: 200px; ">
        <el-form :inline="true" class="demo-form-inline" >
          <span style="margin-left: 65px;">Gene List:&nbsp&nbsp</span>
          <div style="margin-left: 60px;width: 800px;height: 150px">
            <el-input
                type="textarea"
                :rows="5"
                v-model="genes"
                placeholder="Please enter gene symbols, separated by commas"
                style="width: 800px;height: 150px"
            >
            </el-input>
          </div>
          <div style="margin-left: 60px">
            <span @click="clicksample" style="color: #3a8ee6; cursor: pointer;">Sample</span>
            <el-button type="primary" @click="getPathwayPic" style="margin-left: 100px;margin-top: 10px" >Click</el-button>
          </div>
        </el-form>
      </div>
    </div>
    <el-divider></el-divider>
    <div class="pathway-container">
      <div style="min-height: 640px">
        <el-empty description="Wait For Your Click" v-show="!show" style="height: 600px"></el-empty>
        <div v-loading="loading" v-show="show">
          <el-image
              style="width: 800px;height: 600px"
              :src="'data:image/png;base64,' + png64"
              fit="contain">
          </el-image>
          <br style="margin-top: 40px">
          <el-button type="primary" round @click="downloadDiffImg()">Download</el-button>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref } from 'vue';
import jsPDF from "jspdf";
import { ElMessage } from 'element-plus';
import axios from 'axios';

// 基础URL
const baseURL = 'http://121.37.88.191:8989';

// 响应式数据
const png64 = ref("");
const show = ref(false);
const options = ref([]);
const genes = ref('');
const loading = ref(false);
const allGene = ref(false);
const diffList = ref([]);

// 获取药物通路图
const getPathwayPic = () => {
  show.value = true;
  loading.value = true;
  png64.value = "";

  // 将输入的多个基因分隔开来
  const geneList = genes.value.split(/[, \n]+/).map(gene => gene.trim()).filter(gene => gene);

  if (geneList.length > 200) {
    ElMessage({
      message: 'The number of genes cannot exceed 200, so please re-enter.',
      type: 'warning'
    });
    return;
  }

  // 输出输入的多基因信息
  console.log(geneList)
  axios.get(`${baseURL}/r/drugpathway/${geneList}`).then(res => {
    png64.value = res.data;
    loading.value = false;
  });
};

// 下载数据
const downloadDiffData = () => {
  let headerRow = ['Gene', 'Drug', 'p.value', 'cor'];
  let csvContent = headerRow.join(",") + "\n";
  let csvData = diffList.value.map(row => {
    let { id, ...rowData } = row;
    return Object.values(rowData).join(",");
  });
  csvContent += csvData.join("\n");

  let blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
  let link = document.createElement("a");
  let filename = genes.value.replace(/,/g, '_') + '_drug.csv';
  if (link.download !== undefined) {
    let url = URL.createObjectURL(blob);
    link.setAttribute("href", url);
    link.setAttribute("download", filename);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  }
};

// 下载图片
const downloadDiffImg = () => {
  if (!png64.value){
    ElMessage({
      message:'There are no analysis results on the current page',
      type:'warning'
    });
    return;
  }

  let doc = new jsPDF({
    orientation: 'landscape',
    unit: 'px',
    format: [857, 2572]
  });

  doc.addImage('data:image/png;base64,' + png64.value, 'PNG', 0, 0);
  let filename = genes.value.replace(/,/g, '_') + '_drug.pdf';
  doc.save(filename);
};

// 点击示例
const clicksample = () => {
  genes.value = 'EZH2 TP53 BRCA1';
};
</script>

<style scoped>
#pathway{
  width: 1100px;
  margin: 0 auto;
  border-style: solid;
  position: relative;
  border-color: rgba(168,168,168,0);
  background-color: #ffffff;
}
</style>