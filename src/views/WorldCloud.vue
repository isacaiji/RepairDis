<template>
  <div>
    <svg :width="width" :height="height">
      <a class="fontA" v-for="(tag, index) in tags" :key="`tag-${index}`">
        <text
            :id="tag.id"
            :x="tag.x"
            :y="tag.y"
            :data-ensemble-id="tag.ID"
            :font-size="12 * (600 / (600 - tag.z))"
            :fill-opacity="(400 + tag.z) / 600"
            @mousemove="listenerMove"
            @mouseout="listenerOut"
            @click="clickToPage"
        >
          {{ tag.text }}
        </text>
      </a>
    </svg>
    <el-dialog
        v-model="dialogVisible"
        width="80%"
    >
      <div style="border-style: solid; padding: 30px; border-radius: 20px; margin: 100px">
        <ListCard ref="result" :id="id" :basicShow="dialogVisible"/>
        <div style="text-align: center; margin-top: 20px">
          <span slot="footer" class="dialog-footer" style="position:relative;margin-top: 20px">
            <el-button @click="dialogVisible = false">close</el-button>
          </span>
        </div>
      </div>
    </el-dialog>
  </div>
</template>

<script setup >
import {ref, computed, onMounted,watch} from 'vue';
import ListCard from "@/summary/Listcardl.vue";
import axios from "axios";
import {useRouter} from "vue-router";

// 定义 props
const props = defineProps({
  genes: {
    type: Array,
    default: () => []
  }
});
// 定义响应式数据

const dialogVisible = ref(false);
const ensembleID = ref('');
const id = ref('');
const width = ref(500);
const height = ref(400);
const tagsNum = ref(0);
const RADIUS = ref(200);
const speedX = ref(Math.PI / 360 / 1.5 / 2);
const speedY = ref(Math.PI / 360 / 1.5 / 2);
const tags = ref([]);
const timer = ref(null);
const result = ref(null);
const dialog = ref(null);

const router = useRouter();
// 计算属性
const CX = computed(() => width.value / 2);
const CY = computed(() => height.value / 2);

// 初始化数据
const initData = () => {
  let newTags = [];
  tagsNum.value = props.genes.length;
  for (let i = 0; i < props.genes.length; i++) {
    let tag = {};
    let k = -1 + (2 * (i + 1) - 1) / tagsNum.value;
    let a = Math.acos(k);
    let b = a * Math.sqrt(tagsNum.value * Math.PI);
    tag.text = props.genes[i].geneName;
    tag.ID = props.genes[i].id.toString();
    tag.x = CX.value + RADIUS.value * Math.sin(a) * Math.cos(b);
    tag.y = CY.value + RADIUS.value * Math.sin(a) * Math.sin(b);
    tag.z = RADIUS.value * Math.cos(a);
    tag.id = i;
    newTags.push(tag);
  }
  tags.value = newTags;
};

// X 轴旋转函数
const rotateX = (angleX) => {
  const cos = Math.cos(angleX);
  const sin = Math.sin(angleX);
  for (let tag of tags.value) {
    const y1 = (tag.y - CY.value) * cos - tag.z * sin + CY.value;
    const z1 = tag.z * cos + (tag.y - CY.value) * sin;
    tag.y = y1;
    tag.z = z1;
  }
};

// Y 轴旋转函数
const rotateY = (angleY) => {
  const cos = Math.cos(angleY);
  const sin = Math.sin(angleY);
  for (let tag of tags.value) {
    const x1 = (tag.x - CX.value) * cos - tag.z * sin + CX.value;
    const z1 = tag.z * cos + (tag.x - CX.value) * sin;
    tag.x = x1;
    tag.z = z1;
  }
};

// 运行标签动画
const runTags = () => {
  if (timer.value) {
    clearInterval(timer.value);
    timer.value = null;
  }
  timer.value = setInterval(() => {
    rotateX(speedX.value);
    rotateY(speedY.value);
  }, 17);
};

// 鼠标移动事件处理
const listenerMove = (e) => {
  if (e.target.id) {
    clearInterval(timer.value);
  }
};

// 鼠标移出事件处理
const listenerOut = (e) => {
  if (e.target.id) {
    runTags();
  }
};

// 点击事件处理
// const clickToPage = (e) => {
//   id.value = e.target.getAttribute("data-ensemble-id");
//   console.log(id.value);
//   dialogVisible.value = true;
//   result.value.getGeneData(id.value);
//   // axios.get(`http://localhost:8082/api/genes/${id.value}`);
// };

function clickToPage(e){
  id.value = e.target.getAttribute("data-ensemble-id");
  router.push({ path: '/detail', query: { id: id.value } });
}

// 监听 genes 的变化
watch(() => props.genes, (newVal) => {
  if (Array.isArray(newVal) && newVal.length > 0) {
    initData();
    runTags();
  }
}, { immediate: true });

// 组件挂载时运行标签动画
onMounted(() => {
  if (Array.isArray(props.genes) && props.genes.length > 0) {
    initData();
    runTags();
  } else {
    console.error("genes is not a valid array or is empty");
  }
});
</script>

<style scoped>
.fontA {
  fill: #607b8f;
  /*font-weight: bold;*/
}
.fontA:hover {
  fill: rgba(227, 66, 66, 0.92);
  cursor: pointer;
}
</style>