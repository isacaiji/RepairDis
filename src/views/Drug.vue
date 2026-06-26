<template>
  <div class="drug-main-container">
    <!-- 模仿 Analysis 组件风格的 Tab 栏 -->
    <div class="main-menu-container">
      <div class="main-menu">
        <div
            :class="{ 'menu-item': true, 'active-menu': activeTab === 'overview' }"
            @click="handleMenuClick('overview')"
        >
          Drug Overview
        </div>
        <div
            :class="{ 'menu-item': true, 'active-menu': activeTab === 'relevance' }"
            @click="handleMenuClick('relevance')"
        >
          Drug Relevance
        </div>
        <div
            :class="{ 'menu-item': true, 'active-menu': activeTab === 'pathway' }"
            @click="handleMenuClick('pathway')"
        >
          Drug Pathway
        </div>
      </div>
    </div>

    <!-- 内容展示区域 -->
    <div class="tab-content">
      <component
          :is="activeTab === 'overview' ? DrugOverview :
              activeTab === 'relevance' ? MainDrug : Pathway"
      />
    </div>
  </div>
</template>

<script setup>
import { ref } from 'vue';
// 导入三个子组件
import DrugOverview from '@/drug/drug.vue';
import MainDrug from '@/drug/MainDrug.vue';
import Pathway from '@/drug/Pathway.vue';

// 激活的Tab（默认选中第一个）
const activeTab = ref('overview');

// Tab切换事件
const handleMenuClick = (key) => {
  activeTab.value = key;
  window.scrollTo({ top: 0, behavior: 'smooth' });
};
</script>

<style scoped>
.drug-main-container {
  margin: 0 auto;
  background-color: #f8f9fb;
  min-height: 100vh;
  padding: 0 50px;
}

/* 主菜单容器（完全模仿 Analysis 组件） */
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

/* 激活菜单项底部渐变指示器（核心风格） */
.active-menu::after {
  content: '';
  position: absolute;
  bottom: -18px; /* 刚好对齐底部边框 */
  left: 50%;
  transform: translateX(-50%);
  width: 80%;
  height: 4px;
  background: linear-gradient(to right, #003f88, #00509d);
  border-radius: 2px;
}

/* 内容区域（模仿 Analysis 组件） */
.tab-content {
  border-top: 1px solid #e6f0ff;
  padding: 20px;
  background-color: white;
  margin: 0 20px;
  border-radius: 0 0 12px 12px;
  min-height: 800px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
}
</style>