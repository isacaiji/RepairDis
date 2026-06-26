<template>
  <div class="main-container">
    <!-- 搜索框 -->
    <div class="search-container">
      <el-input
          type="text"
          placeholder="Please enter the name"
          v-model="geneQuery"
          class="gene-input"
          @keyup.enter="fetchStructure"
          style="width: 40%;margin-right: 10px;min-height: 40px;"
          @focus="focus"
          @blur="blur"
      ></el-input>
      <el-button
          @click="fetchStructure"
          style="background-color: #1179f4;border-radius: 8px;min-width: 100px;min-height: 40px;"
      >
        <span style="color: white;font-size: 15px">Search</span>
      </el-button>
      <Transition name="slide-down">
        <ul v-if="showSuggestions" class="suggestions-list">
          <li v-for="suggestion in filteredGenes" :key="suggestion" @click="selectSuggestion(suggestion)">{{ suggestion }}</li>
        </ul>
      </Transition>
    </div>
    <!-- 3D 图像 -->
    <div style="height: 500px;display: flex;justify-content: center;align-items: center">
      <div class="mol-container" ref="container"></div>
      <!-- 错误信息 -->
      <div class="err" v-if="errmsg.length > 0 ">
        <p  style="text-align: center;font-size: 20px;margin-top: 200px;color: red">{{ errmsg }}</p>
      </div>
    </div>

  </div>
</template>

<script setup lang="ts">
// 接收蛋白质名称，从后端获取数据
import * as $3Dmol from '3dmol/build/3Dmol';
import {computed, onMounted, ref, watch} from 'vue';
import axios from 'axios';
import geneList from "@/components/geneList"

// 搜索栏
// 定义基因名称数据
const genes = geneList;
// 输入框内容
const geneQuery = ref('');
// 过滤后的基因列表
const filteredGenes = ref<string[]>([]);
// 是否显示建议列表
const showSuggestions = ref(false);

// 输入框点击时需要显示的数据
function getSelectedData() {
  if (geneQuery.value) {
    filteredGenes.value = genes.filter((gene) =>
        gene.toLowerCase().includes(geneQuery.value.toLowerCase())
    );
    showSuggestions.value = filteredGenes.value.length > 1;
  } else {
    filteredGenes.value = genes;
    showSuggestions.value = true;
  }
}

// 输入框点击事件
function focus() {
  getSelectedData();
}

// 输入框失焦
function blur() {
  setTimeout(() => {
    showSuggestions.value = false;
  }, 200);
}

// 选择建议项
function selectSuggestion(suggestion: string) {
  geneQuery.value = suggestion;
  showSuggestions.value = false;
}

watch(geneQuery, (newValue) => {
  if (newValue) {
    getSelectedData();
  } else {
    showSuggestions.value = false;
  }
});

const container = ref(null);
// 参数
const config = { backgroundColor: 'white' };

// 新增：是否显示结构信息
const showInfo = ref(false);
// 新增：结构信息文本
const infoText = ref('');

function createViewers(data: string) {
  let viewer = $3Dmol.createViewer(container.value, config);
  viewer.addModel(data, 'pdb');
  viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
  viewer.zoomTo();
  viewer.render();
  viewer.zoom(1.2, 2000);
}

// 请求地址
const url = computed(() => `http://121.37.88.191:9016/evolution/proteins/${geneQuery.value}`);
const errmsg = ref('');

//
function fetchStructure() {
  axios({
    method: 'get',
    url: url.value,
    responseType: 'text',
  })
      .then((res) => {
        createViewers(res.data);
      })
      .catch((error) => {
        errmsg.value = 'The Protein Structure is NOT FOUND!';
        console.log(error);
      });
}
</script>

<style scoped>
.mol-container {
  width: 100%;
  height: 500px;
  position: relative;
  margin-top: 20px;
}

.suggestions-list {
  list-style-type: none;
  padding: 0;
  margin: 0;
  border: 1px solid #ccc;
  width: 35.5%;
  max-height: 200px;
  overflow-y: auto;
  position: absolute;
  background-color: white;
  z-index: 1;
  left: 435px;
}

.suggestions-list li {
  padding: 8px;
  cursor: pointer;
}

.suggestions-list li:hover {
  background-color: #f0f0f0;
}

/* 定义 slide-down 过渡动画 */
.slide-down-enter-active,
.slide-down-leave-active {
  transition: all 0.3s ease;
}

.slide-down-enter-from,
.slide-down-leave-to {
  transform: translateY(-10px);
  opacity: 0;
}

</style>