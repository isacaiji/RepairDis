<template>
  <div class="drug">
    <div class="drug-body">
      <div>
        <el-form :inline="true" class="demo-form-inline">
          <span style="font-size: 16px">Gene:&nbsp&nbsp</span>
          <el-select
              v-model="gene"
              filterable
              remote
              placeholder="please enter gene symbol"
              :remote-method="remoteMethod"
              :loading="loading"
              style="width: 300px;margin-right: 10px">
            <el-option
                v-for="option in options"
                :label="option"
                :value="option"
                :key="option">
            </el-option>
          </el-select>
          <el-button type="primary" @click="getcheckpointPic">Click</el-button>
        </el-form>
      </div>

      <hr style="margin-bottom: 40px;margin-top: 30px">

      <el-empty description="Wait For Your Click" v-show="!show" style="height: 600px"></el-empty>

      <div v-loading="loading" v-show="show">
        <el-image
            style="width: 1100px"
            :src="'data:image/png/;base64,' + png64"
            fit="contain">
        </el-image>
        <br>
        <el-button type="primary" round @click="downloadDiffImg" style="margin-top: 30px">Download</el-button>
        <el-divider></el-divider>
        <el-table
            :stripe="true"
            :data="diffList"
            height="400"
            border
            style="width: 100%;">
          <el-table-column
              prop="gene"
              label="Gene">
          </el-table-column>
          <el-table-column
              prop="drugname"
              label="Drug">
          </el-table-column>
          <el-table-column
              prop="value"
              label="p.value">
          </el-table-column>
          <el-table-column
              prop="cor"
              label="cor">
          </el-table-column>
        </el-table>
        <br style="margin-top: 40px">
        <el-button type="primary" round @click="downloadDiffData">Download</el-button>
        <br>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue';
import jsPDF from "jspdf";
import axios from 'axios';

// 基础URL
const baseURL = 'http://121.37.88.191:8989';

// 响应式数据
const png64 = ref("");
const show = ref(false);
const allGenes = ref([]);
const options = ref([]);
const gene = ref('ATM');
const loading = ref(false);
const allGene = ref(false);
const diffList = ref([]);

// 模糊搜索
const remoteMethod = (query) => {
  if (query !== "") {
    options.value = allGenes.value.filter((item) => {
      return item.toLowerCase().indexOf(query.toLowerCase()) > -1
    });
  } else {
    options.value = [];
  }
};

// 获取免疫检查点图
const getcheckpointPic = () => {
  show.value = true;
  loading.value = true;
  png64.value = "";
  axios.get(`${baseURL}/r/drug/${gene.value}`).then(res => {
    png64.value = res.data;
    axios.get(`${baseURL}/data/drug/${gene.value}`).then(res => {
      diffList.value = res.data;
    });
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
  let filename = gene.value + '_drug.csv';
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
  let doc = new jsPDF({
    orientation: 'landscape',
    unit: 'px',
    format: [857, 2572]
  });

  doc.addImage('data:image/png;base64,' + png64.value, 'PNG', 0, 0);
  let filename = gene.value + '_drug.pdf';
  doc.save(filename);
};

// 挂载时获取数据
onMounted(() => {
  axios.get(`http://121.37.88.191:9016/api/genes/all`).then(res => {
    allGenes.value = res.data;
  });
});
</script>

<style scoped>
.drug{
  width: 1100px;
  margin: 0 auto;
  border-style: solid;
  border-color: rgba(168,168,168,0);
  background-color: #ffffff;
}
.drug-body{
  position: relative;
  width: 100%;
  height: 90%;
  background: #ffffff;
}
</style>