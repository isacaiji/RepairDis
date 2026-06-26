<template>
  <div class="analysis-container">
    <div class="module-header">
      <div>
        <h2>Evolution Analysis</h2>
        <p>Explore cross-species conservation of DDR genes using phylogenetic trees.</p>
      </div>
    </div>

    <div class="tree-content">
      <Evolution />
    </div>
  </div>
</template>

<script setup>
import Evolution from '@/views/Tree.vue';
</script>

<style scoped>
.analysis-container {
  max-width: 1200px;
  margin: 40px auto;
  padding: 28px;
  border: 1px solid rgba(0, 40, 85, 0.18);
  border-radius: 18px;
  background: linear-gradient(180deg, #ffffff 0%, #f6fafc 100%);
  box-shadow: 0 16px 40px rgba(0, 40, 85, 0.08);
}

.module-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 20px;
  text-align: left;
}

.module-header h2 {
  margin: 0 0 6px;
  color: #082b52;
  font-size: 26px;
  font-weight: 800;
}

.module-header p {
  margin: 0;
  color: #5f7287;
  font-size: 14px;
}

.tree-content {
  border-top: 1px solid #dce6ee;
  padding-top: 22px;
}
</style>
