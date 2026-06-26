<template>
  <div class="homepage">
    <div class="hero-section">
      <div class="hero-overlay"></div>

      <div class="hero-content">
        <header class="header">
          <h1 class="main-title">RepairDis</h1>
          <div class="tagline">
            A Curated Multi-omics Database for DNA Damage Repair Analysis in Cancer Research
          </div>
        </header>

        <!-- 搜索区 + 左侧 Logo -->
        <div class="search-row">
          <div class="search-container">
            <div class="search-form">
              <el-input
                  v-model="searchQuery"
                  placeholder="Enter gene name (e.g., TP53, BRCA1)"
                  @keyup.enter="searchGenes"
                  clearable
                  class="search-input"
              ></el-input>

              <el-button type="primary" @click="searchGenes" class="search-btn">
                <el-icon><Search /></el-icon>
                <span>Search</span>
              </el-button>
            </div>

            <div class="search-examples">
              <span>Examples: </span>
              <span @click="fillSearchQuery('TP53')" class="example-tag">TP53</span>
              <span @click="fillSearchQuery('BRCA1')" class="example-tag">BRCA1</span>
              <span @click="fillSearchQuery('ATM')" class="example-tag">ATM</span>
            </div>
          </div>
        </div>
      </div>
    </div>

    <div class="main-content">
      <div class="intro-section">
        <el-card class="intro-card">
          <h2 class="section-title">About RepairDis</h2>
          <p class="intro-text">
            RepairDis is a comprehensive online platform for DNA damage repair
            (DDR)-related cancer research. It integrates core analytical functions
            including DDR gene cluster analysis, survival prognosis assessment,
            drug-target association prediction, and pathway visualization. These
            features enable systematic exploration of DDR molecular characteristics,
            regulatory networks, and their roles in cancer progression across diverse
            cancer types. As a user-friendly one-stop solution, it streamlines research
            workflows, supporting in-depth investigations into DDR mechanisms and precision therapy development.
          </p>
        </el-card>

        <el-card class="data-card">
          <h2 class="section-title">Database Statistics</h2>
          <div class="data-stats">
            <div class="stat-item">
              <div class="stat-value">236</div>
              <div class="stat-label">Repair Genes</div>
            </div>
            <div class="stat-item">
              <div class="stat-value">44</div>
              <div class="stat-label">Species</div>
            </div>
          </div>
        </el-card>
      </div>

      <div class="features-section">
        <h2 class="section-title">Analytical Modules</h2>
        <div class="features-grid">
          <el-card class="feature-card" @click="navigateTo('/analysis')">
            <div class="card-inner">
              <div class="card-image">
                <img src="@/assets/Analytical%20Modules/molecular.png" alt="Molecular" class="module-image">
              </div>
              <div class="card-content">
                <h3 class="card-title">Molecular Landscape</h3>
                <p class="card-desc">Integrate multi-omics data to analyze DDR gene expression, mutation patterns</p>
              </div>
              <div class="card-action">
                <span>Explore</span><el-icon class="action-icon"><ArrowRight /></el-icon>
              </div>
            </div>
          </el-card>

          <el-card class="feature-card" @click="navigateTo('/network/ppi')">
            <div class="card-inner">
              <div class="card-image">
                <img src="@/assets/Analytical%20Modules/network.png" alt="Network" class="module-image">
              </div>
              <div class="card-content">
                <h3 class="card-title">Interaction Network</h3>
                <p class="card-desc">Visualize protein-protein, gene-gene, and pathway interactions</p>
              </div>
              <div class="card-action">
                <span>Explore</span><el-icon class="action-icon"><ArrowRight /></el-icon>
              </div>
            </div>
          </el-card>

          <el-card class="feature-card" @click="navigateTo('/drug')">
            <div class="card-inner">
              <div class="card-image">
                <img src="@/assets/Analytical%20Modules/drug.png" alt="Drug" class="module-image">
              </div>
              <div class="card-content">
                <h3 class="card-title">Drug Analysis</h3>
                <p class="card-desc">Predict potential drugs targeting DDR pathways and analyze relationships</p>
              </div>
              <div class="card-action">
                <span>Explore</span><el-icon class="action-icon"><ArrowRight /></el-icon>
              </div>
            </div>
          </el-card>

          <el-card class="feature-card" @click="navigateTo('/evolution')">
            <div class="card-inner">
              <div class="card-image">
                <img src="@/assets/Analytical%20Modules/evolution.png" alt="Evolution analysis" class="module-image">
              </div>
              <div class="card-content">
                <h3 class="card-title">Evolution Analysis</h3>
                <p class="card-desc">Explore cross-species conservation of DDR genes using phylogenetic trees</p>
              </div>
              <div class="card-action">
                <span>Explore</span><el-icon class="action-icon"><ArrowRight /></el-icon>
              </div>
            </div>
          </el-card>
        </div>
      </div>

      <div class="highlight-section">
        <h2 class="section-title">DNA Repair Pathways with Key Therapeutic Targets</h2>
        <el-card class="visual-card no-title highlight-card">
          <div class="highlight-image-container">
            <img :src="ddrCancer" alt="DNA Repair Pathways" class="highlight-img">
          </div>
        </el-card>
      </div>

      <div class="visualization-section">
        <el-card class="visual-card no-title">
          <div class="carousel-container">
            <div class="carousel-wrapper">
              <button
                  v-if="carouselImages.length > 1"
                  class="carousel-btn carousel-btn-left"
                  @click="prevSlide"
                  :disabled="currentIndex === 0"
              >
                <el-icon><ArrowLeft /></el-icon>
              </button>

              <div class="carousel-main">
                <transition name="carousel-fade" mode="out-in">
                  <div class="carousel-item" :key="currentIndex">
                    <div class="carousel-image-wrapper-large">
                      <img :src="carouselImages[currentIndex].url" class="carousel-img">
                    </div>
                    <div class="carousel-caption">
                      {{ carouselImages[currentIndex].desc }}
                    </div>
                  </div>
                </transition>
              </div>

              <button
                  v-if="carouselImages.length > 1"
                  class="carousel-btn carousel-btn-right"
                  @click="nextSlide"
                  :disabled="currentIndex === carouselImages.length - 1"
              >
                <el-icon><ArrowRight /></el-icon>
              </button>
            </div>

            <div class="carousel-indicators" v-if="carouselImages.length > 1">
              <button
                  v-for="(item, index) in carouselImages"
                  :key="index"
                  class="indicator-dot"
                  :class="{ active: index === currentIndex }"
                  @click="goToSlide(index)"
              ></button>
            </div>
          </div>
        </el-card>

        <el-card class="visual-card no-title">
          <div class="wordcloud-container">
            <WorldCloud v-if="allGeneData" :genes="allGeneData" class="wordcloud"></WorldCloud>
            <div v-else class="loading-text">Loading gene distribution...</div>
          </div>
        </el-card>
      </div>
    </div>
  </div>
</template>

<script setup>
import { onMounted, ref } from 'vue';
import { useRouter } from 'vue-router';
import axios from 'axios';
import WorldCloud from "@/views/WorldCloud.vue";
import { ArrowRight, Search, ArrowLeft } from "@element-plus/icons-vue";

// 导入图片
import ddrCancer from '@/assets/ddr1.png';
import ddrPathway from '@/assets/ddr2.png';

const carouselImages = ref([
  {
    url: ddrPathway,
    desc: 'DNA Repair Failure in BRCA-Mutant Breast Cancer'
  }
]);

const currentIndex = ref(0);

const prevSlide = () => {
  if (currentIndex.value > 0) currentIndex.value--;
};

const nextSlide = () => {
  if (currentIndex.value < carouselImages.value.length - 1) currentIndex.value++;
};

const goToSlide = (index) => {
  currentIndex.value = index;
};

const searchQuery = ref('');
const router = useRouter();

const searchGenes = async () => {
  if (!searchQuery.value.trim()) return;
  try {
    const response = await axios.get(`http://121.37.88.191:9016/api/genes/search?query=${searchQuery.value}`);
    if (response.data && response.data.length > 0) {
      await router.push({ path: '/detail', query: { id: response.data[0].id } });
    }
  } catch (error) {
    console.error('Error:', error);
  }
};

const fillSearchQuery = (example) => {
  searchQuery.value = example;
};

const navigateTo = (path) => {
  router.push(path);
  window.scrollTo(0, 0);
};

const allGeneData = ref(null);

onMounted(() => {
  axios.get('http://121.37.88.191:9016/api/genes').then(res => {
    allGeneData.value = res.data;
  });
});
</script>

<style scoped>
.homepage {
  font-family: 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  color: #333;
  background-color: #f9fafb;
  min-height: 100vh;
  display: flex;
  flex-direction: column;
}

/* Hero Section */
.hero-section {
  position: relative;
  height: 550px;
  background-image: url('../assets/back.jpg');
  background-size: cover;
  background-position: center;
  overflow: hidden;
}

.hero-overlay {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  background: linear-gradient(135deg, rgba(0, 40, 85, 0.85), rgba(0, 63, 136, 0.8));
}

.hero-content {
  position: relative;
  height: 100%;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  padding: 0 20px;
  max-width: 1200px;
  margin: 0 auto;
}

.header {
  text-align: center;
  margin-bottom: 26px;
}

.main-title {
  font-size: 64px;
  font-weight: 700;
  color: #ffffff;
  margin-bottom: 15px;
  letter-spacing: 1px;
}

.tagline {
  font-size: 18px;
  color: #b3d1ff;
  font-style: italic;
}

/* Search */
.search-row {
  width: 100%;
  max-width: 920px;
  display: flex;
  align-items: flex-start;
  justify-content: center;
  position: relative;
}

.search-container {
  flex: 0 1 800px;
  max-width: 800px;
  margin: 0 auto;
  background-color: rgba(255, 255, 255, 0.95);
  padding: 20px;
  border-radius: 10px;
  box-shadow: 0 10px 30px rgba(0, 0, 0, 0.15);
}

.search-form {
  display: flex;
  width: 100%;
  gap: 10px;
}

.search-input {
  flex: 1;
  height: 50px !important;
  font-size: 16px;
}

.search-btn {
  height: 50px !important;
  padding: 0 25px;
  font-size: 16px;
  background: linear-gradient(to right, #003f88, #00509d);
  border: none;
  transition: all 0.3s ease;
  color: white;
}

.search-btn:hover {
  background: linear-gradient(to right, #002855, #003f88);
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(0, 63, 136, 0.3);
}

.search-examples {
  margin-top: 15px;
  color: #495057;
  font-size: 14px;
}

.example-tag {
  color: #00509d;
  font-weight: 500;
  cursor: pointer;
  padding: 0 8px;
  margin: 0 4px;
  border-radius: 4px;
  transition: all 0.2s;
}

.example-tag:hover {
  background-color: #e6f0ff;
  text-decoration: underline;
}

/* Main Content */
.main-content {
  max-width: 1400px;
  margin: 0 auto;
  padding: 70px 20px;
  flex: 1;
}

.section-title {
  color: #003f88;
  font-size: 30px;
  margin-bottom: 25px;
  padding-bottom: 12px;
  border-bottom: 2px solid #e0e7ff;
  font-weight: 600;
}

/* Intro Section */
.intro-section {
  display: flex;
  gap: 30px;
  margin-bottom: 70px;
  align-items: stretch;
}

.intro-card {
  flex: 2;
  border-radius: 12px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  border: none;
  padding: 35px;
  background-color: #ffffff;
  display: flex;
  flex-direction: column;
}

.data-card {
  flex: 1;
  border-radius: 12px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  border: none;
  padding: 35px;
  background-color: #ffffff;
  display: flex;
  flex-direction: column;
}

.intro-text {
  font-size: 17px;
  line-height: 1.8;
  color: #555;
  text-align: justify;
  flex: 1;
}

.data-stats {
  display: flex;
  gap: 20px;
  margin-top: 10px;
  flex: 1;
  align-items: center;
  justify-content: center;
}

.stat-item {
  text-align: center;
  padding: 40px 20px;
  background-color: #f0f5ff;
  border-radius: 8px;
  transition: transform 0.3s ease, box-shadow 0.3s ease;
  flex: 1;
  max-width: 150px;
}

.stat-item:hover {
  transform: translateY(-5px);
  box-shadow: 0 6px 15px rgba(0, 63, 136, 0.1);
}

.stat-value {
  font-size: 36px;
  font-weight: 700;
  color: #003f88;
  margin-bottom: 5px;
}

.stat-label {
  font-size: 15px;
  color: #666;
}

/* Features Section */
.features-section {
  margin-bottom: 70px;
}

.features-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
  gap: 30px;
  margin-top: 30px;
}

.feature-card {
  border-radius: 12px;
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.05);
  border: none;
  overflow: hidden;
  transition: all 0.3s ease;
  cursor: pointer;
  background-color: #ffffff;
}

.feature-card:hover {
  transform: translateY(-8px);
  box-shadow: 0 12px 20px rgba(0, 63, 136, 0.1);
}

.card-inner {
  padding: 30px;
  display: flex;
  flex-direction: column;
  height: 100%;
  text-align: center;
}

.card-image {
  margin: 15px 0;
  height: 140px;
  overflow: hidden;
  border-radius: 8px;
  background-color: #f8f9fa;
  display: flex;
  align-items: center;
  justify-content: center;
}

.module-image {
  max-width: 100%;
  max-height: 100%;
  object-fit: cover;
  transition: transform 0.3s ease;
}

.feature-card:hover .module-image {
  transform: scale(1.05);
}

.card-content {
  flex: 1;
  margin-bottom: 15px;
}

.card-title {
  color: #003f88;
  font-size: 20px;
  margin-bottom: 10px;
  font-weight: 600;
}

.card-desc {
  color: #666;
  font-size: 14px;
  line-height: 1.6;
}

.card-action {
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 8px;
  padding: 10px 0 0;
  color: #00509d;
  font-weight: 500;
  font-size: 14px;
}

.action-icon {
  transition: transform 0.3s ease;
}

.feature-card:hover .action-icon {
  transform: translateX(5px);
}

/* Highlight 区域 */
.highlight-section {
  margin-bottom: 60px;
}

.highlight-card {
  border-radius: 12px;
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.05);
  border: none;
  background-color: #ffffff;
  padding: 10px;
}

.highlight-image-container {
  width: 100%;
  background-color: #ffffff;
  border-radius: 8px;
  overflow: hidden;
  text-align: center;
  padding: 10px 0;
}

.highlight-img {
  width: 100%;
  height: auto;
  object-fit: contain;
  display: block;
  margin: 0 auto;
}

/* Visualization Section */
.visualization-section {
  display: flex;
  gap: 30px;
  margin-bottom: 70px;
  align-items: stretch;
}

.visual-card {
  flex: 1;
  border-radius: 12px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  border: none;
  padding: 35px;
  background-color: #ffffff;
  min-height: 580px;
  display: flex;
  flex-direction: column;
}

.visual-card.no-title {
  padding-top: 20px;
}

/* 轮播区域 */
.carousel-container {
  flex: 1;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  padding: 10px;
}

.carousel-wrapper {
  width: 100%;
  display: flex;
  align-items: center;
  position: relative;
}

.carousel-main {
  flex: 1;
  overflow: hidden;
  position: relative;
  min-height: 480px;
}

.carousel-item {
  width: 100%;
  height: 100%;
}

.carousel-image-wrapper-large {
  width: 100%;
  height: 480px;
  display: flex;
  justify-content: center;
  align-items: center;
  background-color: #ffffff;
  border-radius: 6px;
  overflow: hidden;
}

.carousel-img {
  max-width: 100%;
  max-height: 100%;
  object-fit: contain;
}

.carousel-caption {
  text-align: center;
  padding: 15px 0 0;
  color: #555;
  font-size: 15px;
  font-weight: 500;
}

/* 轮播按钮 */
.carousel-btn {
  width: 36px;
  height: 36px;
  border-radius: 4px;
  background-color: #f0f5ff;
  color: #00509d;
  border: 1px solid #dbeafe;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 14px;
  cursor: pointer;
  z-index: 5;
  transition: all 0.2s ease;
  flex-shrink: 0;
}

.carousel-btn:hover:not(:disabled) {
  background-color: #e0e7ff;
  border-color: #bfdbfe;
}

.carousel-btn:disabled {
  opacity: 0.3;
  cursor: not-allowed;
}

.carousel-btn-left {
  margin-right: 15px;
}

.carousel-btn-right {
  margin-left: 15px;
}

/* 词云 */
.wordcloud-container {
  flex: 1;
  display: flex;
  justify-content: center;
  align-items: center;
  padding: 20px;
}

.wordcloud {
  width: 100%;
  height: 100%;
  min-height: 480px;
}

.loading-text {
  color: #666;
  font-size: 16px;
  text-align: center;
}

/* 指示器 */
.carousel-indicators {
  display: flex;
  justify-content: center;
  gap: 8px;
  margin-top: 15px;
}

.indicator-dot {
  width: 8px;
  height: 8px;
  border-radius: 50%;
  background-color: #e0e7ff;
  border: none;
  cursor: pointer;
  transition: all 0.2s ease;
}

.indicator-dot.active {
  background-color: #00509d;
  width: 24px;
  border-radius: 4px;
}

.carousel-fade-enter-active,
.carousel-fade-leave-active {
  transition: opacity 0.3s ease;
}

.carousel-fade-enter-from,
.carousel-fade-leave-to {
  opacity: 0;
}
</style>
