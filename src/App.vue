<template>
  <div id="app">
    <div class="app-container">
      <el-container class="flex-container">
        <el-header class="da-head" style="position: relative; height: 55px">
          <!-- LOGO区域 -->
          <div class="LOGO-container">
            <router-link to="/home" class="brand-link" aria-label="RepairDis home">
              <img src="@/assets/logo/repairdis-brand.png" alt="RepairDis" class="brand-icon">
              <span class="brand-word">
                <span class="brand-word-main">Repair</span><span class="brand-word-accent">Dis</span>
              </span>
            </router-link>
          </div>

          <!-- 导航菜单容器 -->
          <div class="da-wrap">
            <div
                v-for="(menu, index) in fixedOrderMenuList"
                :key="`menu-${index}`"
                :class="['item1', { 'has-submenu': menu.type === 'dropdown' }]"
            >
              <!-- 1. 普通导航项（无下拉） -->
              <template v-if="menu.type === 'normal'">
                <router-link :to="menu.path" class="flex-content">
                  <!-- 新增active判断：当前路由包含菜单项path时激活 -->
                  <div class="menu" :class="{ 'active': route.path.includes(menu.path) && menu.path !== '/' }">
                    <span> {{ menu.title }}</span>
                  </div>
                </router-link>
              </template>

              <!-- 2. 带下拉的导航项 -->
              <template v-else-if="menu.type === 'dropdown'">
                <div
                    class="menu"
                    @mouseenter="showDropdownIndex = index"
                    @mouseleave="showDropdownIndex = -1"

                :class="{ 'active': menu.subList.some(item => route.path.includes(item.path)) }"
                >
                <span> {{ menu.title }}</span>
                <div class="submenu" v-if="showDropdownIndex === index">
                  <router-link
                      v-for="(subItem, subIndex) in menu.subList"
                      :key="`sub-${index}-${subIndex}`"
                      :to="subItem.path"
                      class="submenu-item"
                      @click="showDropdownIndex = -1"
                  :class="{ 'sub-active': route.path === subItem.path }"
                  >
                  {{ subItem.label }}
                  </router-link>
                </div>
            </div>
</template>

<!-- 3. 特殊项（Contact us，无路由跳转） -->
<template v-else-if="menu.type === 'special'">
  <div class="menu" @click="sendEmail">
    <span style="font-size: 18px"><i class="el-icon-message"></i> {{ menu.title }}</span>
  </div>
</template>
</div>
</div>
</el-header>

<!-- 主内容区域 -->
<el-main style="width: 1400px;min-height: 800px; margin:0 auto;">
<router-view></router-view>
</el-main>

<!-- 页脚 -->
<el-footer style="height: auto; margin: -25px 0 0; padding: 0">
  <Footer/>
</el-footer>
</el-container>
</div>
</div>
</template>

<script setup>
import { ref, watch } from 'vue';
import { useRoute, useRouter } from 'vue-router';
import Footer from "@/views/footer.vue"

// 核心：固定顺序的导航菜单配置
const fixedOrderMenuList = ref([
  { type: 'normal', title: 'Home', path: '/home' },
  { type: 'normal', title: 'Browser', path: '/browser' },
  { type: 'normal', title: 'Molecular landscape', path: '/analysis' },
  { type: 'normal', title: 'Drug', path: '/drug'},
  {
    type: 'dropdown',
    title: 'Network',
    subList: [
      { path: '/network/ppi', label: 'PPI Network' },
      { path: '/network/cross', label: 'Cross-talk Network' },
      { path: '/network/TF', label: 'TF Network' },
      { path: '/network/ncRNA',label: 'ncRNA Network' },
    ]
  },
  {
    type: 'dropdown',
    title: 'Synthetic-lethality',
    subList: [
      { path: '/sl/network', label: 'Interaction' }
    ]
  },
  { type: 'normal', title: 'Evolution', path: '/evolution' },
  { type: 'normal', title: 'Help' , path: '/help' },
  { type: 'special', title: 'Contact us' }
]);

// 下拉菜单控制状态
const showDropdownIndex = ref(-1);

// 路由相关
const router = useRouter();
const route = useRoute();
const color = ref('');
const textcolor = ref('');

// 发送邮件方法
const sendEmail = () => {
  const email = "lguo@njupt.edu.cn";
  const subject = "RepairDis Inquiry";
  window.location.href = `mailto:${email}?subject=${encodeURIComponent(subject)}`;
};

// 监听路由变化：切换路由时关闭所有下拉
watch(
    () => route.path,
    (to) => {
      showDropdownIndex.value = -1;

      if (to.startsWith('/evolution')) {
        color.value = '#ffffff';
        textcolor.value = '#4f5154';
      } else {
        color.value = '';
        textcolor.value = '#ffffff';
      }
      if (to === '/network') {
        router.push('/network/ppi');
      }
    }
);
</script>

<style>
#app {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  position: relative;
  min-height: 750px;
  text-align: center;
  color: #2c3e50;
  background-color: #e3e8ec;
}
::-webkit-scrollbar {
  display: none;
}

body,
div,
dl,
dt,
dd,
ul,
ol,
li,
h1,
h2,
h3,
h4,
h5,
h6,
pre,
form,
fieldset,
legend,
input,
textarea,
button,
p,
blockquote,
th,
td {
  margin: 0;
}

body {
  text-align: center;
  font-family: Helvetica Neue, Helvetica, Arial, Microsoft Yahei, Hiragino Sans GB, Heiti SC, WenQuanYi Micro Hei, sans-serif;
}

li {
  list-style: none;
}

a {
  text-decoration: none;
  color: white;
}

img {
  border: none;
}

router-link {
  color: #e1e1e1;
}

.da-head {
  position: relative;
  width: 1360px;
  margin: 0 auto -20px auto;
  background: #002855;
}

.LOGO-container {
  height: 55px;
  width: 205px;
  flex-shrink: 0;
  display: flex;
  align-items: center;
}

.brand-link {
  height: 55px;
  display: inline-flex;
  align-items: center;
  gap: 8px;
  padding: 0 14px 0 12px;
  color: #ffffff;
  text-decoration: none;
  border-radius: 0 14px 14px 0;
  background: linear-gradient(90deg, rgba(255,255,255,0.12), rgba(255,255,255,0.03));
  transition: background 0.2s ease, transform 0.2s ease;
}

.brand-link:hover {
  background: linear-gradient(90deg, rgba(0,168,166,0.28), rgba(255,255,255,0.08));
  transform: translateX(2px);
}

.brand-icon {
  height: 32px;
  width: 32px;
  object-fit: contain;
  filter: none;
}

.brand-word {
  font-family: Georgia, 'Times New Roman', serif;
  font-size: 24px;
  font-weight: 800;
  letter-spacing: -0.3px;
  line-height: 1;
  white-space: nowrap;
}

.brand-word-main {
  color: #ffffff;
}

.brand-word-accent {
  color: #17b6ae;
}

.da-wrap {
  position: relative;
  display: flex;
  justify-content: flex-end;
  margin-top: -55px;
  margin-left: 205px;
  width: calc(100% - 205px);
}

.item1 {
  padding: 0 3px;
  color: #e1e1e1;
  position: relative;
}

.menu {
  width: 100%;
  height: 55px;
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 0 12px;
  text-align: center;
  position: relative;
  font-weight: bold;
  font-size: 18px;
  line-height: 55px;
  cursor: pointer;
  flex-direction: column;
  transition: all 0.3s ease; /* 平滑过渡效果 */
}

/* 新增激活状态样式 */
.menu.active {
  background-color: white; /* 背景变白 */
  color: #002855; /* 字体颜色设为主题色 */
}

/* 二级菜单容器样式 */
.submenu {
  position: absolute;
  top: 51px;
  width: auto;
  min-width: 150px;
  background: linear-gradient(180deg, #002855 0%, #003f88 100%);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
  z-index: 1000;
  border-radius: 0 0 8px 8px;
  background: #fff;
  line-height: 45px;
}

/* 二级菜单项 */
.submenu-item {
  display: block;
  color: #333;
  font-size: 14px;
  text-align: center;
  transition: background 0.2s;
  padding: 2px 0;
}

/* 子菜单激活状态 */
.submenu-item.sub-active {
  background: #183f88;
  color: #fff;
}

.submenu-item:hover {
  background: #183f88;
  color: #fff;
}

.menu:hover {
  font-size: 18px;
  box-shadow: 1px 1px 1px rgb(158, 163, 168);
  border: 1px solid white;
}

.app-container {
  position: relative;
}

.flex-container {
  margin-bottom: 0;
}
</style>
