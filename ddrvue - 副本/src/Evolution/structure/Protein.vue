<template>
  <div class="main-container">
    <!-- 3D 图像 -->
    <div style="height: 250px;width: 500px;display: flex;justify-content: center;align-items: center">
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

const props = defineProps({
  name:String,
})

const container = ref(null);
// 参数
const config = { backgroundColor: 'white' };

function createViewers(data: string) {
  let viewer = $3Dmol.createViewer(container.value, config);
  viewer.addModel(data, 'pdb');
  viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
  viewer.zoomTo();
  viewer.render();
  viewer.zoom(1.2, 2000);
}
// 121.37.88.191
// const apiBase = ref('http://localhost:9016');
const apiBase = ref('http://121.37.88.191:9016');

//请求地址
const url = computed(() => `${apiBase.value}/proteins/${props.name}`);
const errmsg = ref('');

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
onMounted(()=>{
  fetchStructure();
})

</script>

<style scoped>
.mol-container {
  width: 100%;
  height: 100%;
  position: relative;
  margin-top: 20px;
}
</style>