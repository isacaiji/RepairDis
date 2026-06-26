<template>
  <div class="gene-search-container">
    <!-- 搜索表单 -->
    <div class="search-form-wrapper">
      <div class="search-form">
        <el-select
            v-model="searchQuery"
            filterable
            remote
            placeholder="Enter gene name"
            :remote-method="remoteMethod"
            :loading="loading"
            @change="searchGenes"
        >
          <el-option
              v-for="option in options"
              :label="option"
              :value="option"
              :key="option">
            {{option}}
          </el-option>
        </el-select>
        <el-button type="primary" @click="searchGenes">
          Search
        </el-button>
      </div>
      <!-- 添加搜索例子 -->
      <div class="search-examples">
        <span>eg: </span>
        <span @click="fillSearchQuery('TP53')">TP53</span>
        <span>、</span>
        <span @click="fillSearchQuery('BRCA1')">BRCA1</span>
      </div>
    </div>
    <!-- 显示搜索结果表格 -->
    <div class="empty-res" v-if="searchResults.length === 0">
      <el-empty description="No DATA" style="background-color: white;width: 1200px;height: 750px;margin: 0 auto;" />
    </div>
    <div class="result-wrapper" v-if="searchResults.length > 0">
      <div class="result-header">
        <div class="search-content">Search for: {{ searchQuery }}</div>
        <div class="result-count">Total: {{ searchResults.length }} items</div>
      </div>
      <div class="gene-table-container">
        <div class="table-wrapper">
          <el-table :data="searchResults" stripe border highlight-current-row>
            <el-table-column prop="geneName" label="Gene" align="center"></el-table-column>
            <el-table-column prop="ensembl" label="Ensembl ID" align="center"></el-table-column>
            <el-table-column prop="pathway" label="pathway" align="center"></el-table-column>
            <el-table-column prop="pmid" label="PMID" align="center"></el-table-column>
            <el-table-column label="Operation" >
              <template #default="{ row }">
                <el-button type="custom-button" @click="goToGenedetail(row.id)">View Details</el-button>
              </template>
            </el-table-column>
          </el-table>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import {ref, onMounted} from 'vue';
import {useRouter, useRoute} from 'vue-router';
import axios from 'axios';

// 搜索关键词
const searchQuery = ref('');
// 存储搜索结果
const searchResults = ref([]);
const router = useRouter();
const route = useRoute();
// 模糊搜索
const allGenes = ref([]);
const options = ref([]);
const loading = ref(false);

// 执行搜索
const searchGenes = async () => {
  try {
    const response = await axios.get(`http://121.37.88.191:9016/api/genes/search?query=${searchQuery.value}`);
    searchResults.value = response.data;
  } catch (error) {
    console.error('Error searching genes:', error);
  }
};

function goToGenedetail(id) {
  router.push({path: '/detail', query: {id: id}});
}

// 填充搜索框
const fillSearchQuery = (example) => {
  searchQuery.value = example;
  searchGenes();
};

// 在组件挂载时获取路由传递的查询参数并执行搜索
onMounted(() => {
  axios.get('http://121.37.88.191:9016/api/genes/all').then(res => {
    allGenes.value = res.data;
  });
});

const remoteMethod = (query) => {
  // 如果用户输入内容了，就发请求拿数据，远程搜索模糊查询
  if (query !== "") {
    loading.value = true; // 开始拿数据喽
    options.value = allGenes.value.filter((item) => {
      // 大于-1说明只要有就行，不论是不是开头也好，中间也好，或者结尾也好
      return item.toLowerCase().indexOf(query.toLowerCase()) > -1
    })
    loading.value = false // 拿到数据
  } else {
    options.value = [];
  }
}
</script>

<style scoped>
.gene-search-container {
  padding: 40px;
  max-width: 1400px;
  margin: 0 auto;
  background-color: #f9f9f9; /* 更柔和的背景色 */
  border-radius: 8px;
}

.search-form-wrapper {
  display: flex;
  flex-direction: column;
  align-items: center;
  margin-bottom: 20px;
}

.search-form {
  display: flex;
  width: 60%;
  gap: 10px;
  position: relative;
}

.search-form .el-input {
  flex: 1;
}

.search-form .el-input .el-input__inner {
  background-color: #fff;
  border: 2px solid #002855; /* 使用主题色 */
  color: #333;
  padding: 10px;
  border-radius: 5px;
  outline: none;
  transition: border-color 0.3s ease, box-shadow 0.3s ease;
}

.search-form .el-input .el-input__inner::placeholder {
  color: #aaa;
}

.search-form .el-input .el-input__inner:hover {
  border-color: #0056b3;
  box-shadow: 0 0 5px rgba(0, 123, 255, 0.3);
}

.search-form .el-input .el-input__inner:focus {
  border-color: #0056b3;
  box-shadow: 0 0 5px rgba(0, 123, 255, 0.5);
}

.search-form .el-button {
  background-color: #002855; /* 使用主题色 */
  border-color: #002855;
  color: #fff;
  padding: 0 20px;
  border-radius: 5px;
  transition: background-color 0.3s ease;
}

.search-form .el-button:hover {
  background-color: #0056b3;
  border-color: #0056b3;
}

.search-examples {
  margin-top: 10px;
  text-align: center;
}

.search-examples span {
  cursor: pointer;
  color: #1a80ed;
}

.search-examples span:first-child {
  cursor: default;
  color: #333;
}

.result-wrapper {
  box-shadow: 0 0 15px rgba(0, 0, 0, 0.15);
  border-radius: 8px;
  overflow: hidden;
  background-color: #fff;
  margin-top: 20px;
}

.result-header {
  display: flex;
  justify-content: space-between;
  padding: 10px;
  background-color: #f0f0f0;
  border-bottom: 1px solid #ccc;
}

.search-content {
  font-weight: bold;
}

.result-count {
  font-weight: bold;
}

.el-table {
  font-size: 14px;
}

.detail-card {
  background-color: #fff;
  box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
  border-radius: 4px;
  padding: 20px;
}

.gene-table-container {
  padding: 20px;
  max-width: 1200px;
  margin: 0 auto;
}

.table-wrapper {
  box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
  border-radius: 4px;
  overflow: hidden;
  min-height: 500px;
}

.el-table {
  font-size: 14px;
}

.el-button--custom-button {
  background-color: #007BFF;
  border: none;
  color: white;
  padding: 8px 16px;
  text-align: center;
  text-decoration: none;
  display: inline-block;
  font-size: 14px;
  border-radius: 4px;
  cursor: pointer;
  transition: background-color 0.3s ease, box-shadow 0.3s ease;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
}

.el-button--custom-button:hover {
  background-color: #66b0ff;
  box-shadow: 0 4px 8px rgba(0, 0, 0, 0.3);
}

.el-button--custom-button:active {
  background-color: #0056b3;
  box-shadow: 0 1px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}

.search-suggestions {
  position: absolute;
  top: 100%;
  left: 0;
  width: 100%;
  background-color: white;
  border: 1px solid #ccc;
  border-top: none;
  border-radius: 0 0 4px 4px;
  z-index: 10;
}

.search-suggestions .el-option {
  padding: 8px 12px;
  cursor: pointer;
}

.search-suggestions .el-option:hover {
  background-color: #f0f0f0;
}
</style>