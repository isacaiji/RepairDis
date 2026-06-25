<template>
  <div class="analysis-container">
    <!-- 菜单栏 -->
    <div class="main-menu-container">
      <div class="main-menu">
        <div
            :class="{ 'menu-item': true, 'active-menu': activeIndex === 'base' }"
            @click="handleMenuClick('base')"
        >
          Base Analysis
        </div>
        <div
            :class="{ 'menu-item': true, 'active-menu': activeIndex === 'immune' }"
            @click="handleMenuClick('immune')"
        >
          Immune Analysis
        </div>
      </div>
    </div>

    <div class="custom-tabs">
      <div
          v-for="(tab, index) in currentTabs"
          :key="index"
          :class="{ 'active-tab': activeTabIndex === index, 'tab-item': true }"
          @click="activeTabIndex = index"
      >
        {{ tab.label }}
      </div>
    </div>

    <!-- 内容展示区域 -->
    <div class="tab-content">
      <component :is="currentTabs[activeTabIndex]?.component" />
    </div>
  </div>
</template>

<script setup>
import {computed, ref} from 'vue'

// Import your components
import Cnv from '@/Analysis/dzx/Cnv.vue'
import DiffR from '@/Analysis/dzx/DiffR.vue'
import Methy from '@/Analysis/dzx/Methy.vue'
import Mut from '@/Analysis/dzx/Mut.vue'
import Surv from '@/Analysis/dzx/Surv.vue'
import Expression from '@/Analysis/dzx/Expression.vue'
import CLUSTER from "@/Analysis/dzx/CLUSTER.vue";

//immune
import Checkpoint from '@/Analysis/immume/checkpoint.vue'
import Chemokine from '@/Analysis/immume/chemokine.vue'
import Immucell from '@/Analysis/immume/immucell.vue'
import Immunescore from '@/Analysis/immume/immunescore.vue'
import Immuinhibitor from '@/Analysis/immume/immuinhibitor.vue'
import Immustimulator from '@/Analysis/immume/immustimulator.vue'
import Receptor from '@/Analysis/immume/receptor.vue'
import PanCancerEstimate from '@/Analysis/immume/PanCancerEstimate.vue'
import PanCancerCellType from '@/Analysis/immume/PanCancerCellType.vue'
import PanCancerSuppression from '@/Analysis/immume/PanCancerSuppression.vue'
import PanCancerExclusion from '@/Analysis/immume/PanCancerExclusion.vue'
import PanCancerCheckpoint from '@/Analysis/immume/PanCancerCheckpoint.vue'
import PanCancerTide from '@/Analysis/immume/PanCancerTide.vue'
import GSEA from '@/Analysis/dzx/GSEA.vue'

import Pathway from "@/drug/pathway.vue";
import Drug from "@/drug/drug.vue";

const menuTabs = {
  immune: [
    {label: 'DDR-state TME', component: PanCancerEstimate},
    {label: 'DDR-state cell', component: PanCancerCellType},
    {label: 'DDR-state suppression', component: PanCancerSuppression},
    {label: 'DDR-state exclusion', component: PanCancerExclusion},
    {label: 'DDR-state checkpoint', component: PanCancerCheckpoint},
    {label: 'DDR-state TIDE', component: PanCancerTide},
    {label: 'Checkpoint', component: Checkpoint},
    {label: 'Chemokine', component: Chemokine},
    // {label: 'Immu cell', component: Immucell},
    // {label: 'Immune score', component: Immunescore},
    {label: 'Immune inhibitor', component: Immuinhibitor},
    {label: 'Immune stimulator', component: Immustimulator},
    {label: 'Receptor', component: Receptor},
    {label: 'ssGSEA', component: GSEA}
  ],
  base: [
    {label: 'Copy Number Variation', component: Cnv},
    // {label: 'mRNA Expression', component: DiffR},
    {label: 'DNA Methylation', component: Methy},
    {label: 'Mutation', component: Mut},
    {label: 'Survival', component: Surv},
    {label: 'Expression', component: Expression},
    {label: 'CLUSTER', component: CLUSTER}
  ]
};

const activeIndex = ref('base'); // 默认选中base
const activeTabIndex = ref(0);// 默认选中第一个标签页

// 根据当前激活的菜单，动态计算显示的标签页
const currentTabs = computed(() => {
  return menuTabs[activeIndex.value] || []
});

// 菜单点击事件
const handleMenuClick = (key) => {
  activeIndex.value = key
  activeTabIndex.value = 0
};
</script>

<style scoped>
.analysis-container {
  margin: 0 auto;
  background-color: #f8f9fb;
}

/* 主菜单容器 */
.main-menu-container {
  background-color: white;
  border-radius: 12px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  margin: 20px;
  padding: 16px;
}

/* 主菜单样式 */
.main-menu {
  display: flex;
  gap: 20px;
  border-bottom: 2px solid #e6f0ff;
  padding-bottom: 16px;
}

/* 菜单项样式 */
.menu-item {
  padding: 12px 30px;
  border-radius: 8px;
  cursor: pointer;
  font-size: 18px;
  font-weight: 600;
  color: #666;
  transition: all 0.3s ease;
  position: relative;
}

/* 菜单项悬停效果 */
.menu-item:hover {
  color: #00509d;
  background-color: #f0f5ff;
}

/* 激活菜单项样式 */
.active-menu {
  color: #003f88;
  background-color: #e6f0ff;
}

/* 激活菜单项底部指示器 */
.active-menu::after {
  content: '';
  position: absolute;
  bottom: -18px;
  left: 50%;
  transform: translateX(-50%);
  width: 80%;
  height: 4px;
  background: linear-gradient(to right, #003f88, #00509d);
  border-radius: 2px;
}

/* 子标签栏样式 */
.custom-tabs {
  display: flex;
  flex-wrap: wrap;
  justify-content: flex-start;
  margin: 0 20px 20px;
  gap: 10px;
}

.tab-item {
  padding: 10px 20px;
  border-radius: 8px;
  background: #f0f5ff;
  color: #00509d;
  cursor: pointer;
  font-weight: 500;
  font-size: 14px;
  transition: all 0.3s;
}

.tab-item:hover:not(.active-tab) {
  background: #e6f0ff;
}

.active-tab {
  background: linear-gradient(to right, #003f88, #00509d);
  color: white;
}

/* 内容区域 */
.tab-content {
  border-top: 1px solid #e6f0ff;
  padding: 20px;
  background-color: white;
  margin: 0 20px;
  border-radius: 0 0 12px 12px;
  min-height: 600px;
}
</style>
