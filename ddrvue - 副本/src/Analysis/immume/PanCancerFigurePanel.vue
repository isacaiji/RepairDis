<template>
  <div>
    <div style="margin: 20px;">
      <el-form :inline="true" class="demo-form-inline">
        <span style="font-size: 16px">Cancer:&nbsp;&nbsp;</span>
        <el-select
            v-model="selectedCancer"
            filterable
            placeholder="please select cancer type"
            style="width: 300px;margin-right: 10px"
        >
          <el-option
              v-for="item in cancers"
              :key="item"
              :label="item"
              :value="item"
          />
        </el-select>
        <el-button type="primary" @click="showFigure">Click</el-button>
        <span class="analysis-title">{{ title }}</span>
      </el-form>
    </div>

    <hr>

    <el-empty description="Wait For Your Click" v-show="!show" style="height: 800px"></el-empty>

    <div v-show="show" class="figure-panel">
      <el-image
          class="figure-image"
          :src="imageSrc"
          fit="contain"
          :preview-src-list="[imageSrc]"
      />
      <br>
      <el-button type="primary" round @click="openPdf">Download</el-button>
    </div>
  </div>
</template>

<script setup>
import {computed, ref} from 'vue'

const props = defineProps({
  title: {
    type: String,
    required: true
  },
  moduleDir: {
    type: String,
    required: true
  },
  filePattern: {
    type: String,
    required: true
  },
})

const base = '/repairdis/cancer-immune'
const selectedCancer = ref('LUAD')
const show = ref(false)

const cancers = [
  'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC',
  'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV',
  'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA',
  'THYM', 'UCEC', 'UCS', 'UVM'
]

const fileBase = computed(() => props.filePattern.replace('{cancer}', selectedCancer.value))

const imageSrc = computed(() =>
    `${base}/${props.moduleDir}/${selectedCancer.value}/${fileBase.value}.png`
)

const pdfSrc = computed(() =>
    `${base}/${props.moduleDir}/${selectedCancer.value}/${fileBase.value}.pdf`
)

const showFigure = () => {
  show.value = true
}

const openPdf = () => {
  window.open(pdfSrc.value, '_blank', 'noopener')
}
</script>

<style scoped>
.analysis-title {
  display: inline-flex;
  align-items: center;
  min-height: 32px;
  margin-left: 16px;
  color: #073763;
  font-weight: 600;
}

.figure-panel {
  height: 800px;
  width: 100%;
  padding: 20px;
  box-sizing: border-box;
  overflow: auto;
  text-align: center;
  background-color: #ffffff;
}

.figure-image {
  height: 700px;
  width: 100%;
}
</style>
