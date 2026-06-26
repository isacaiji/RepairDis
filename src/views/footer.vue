<template>
  <div class="site-footer">
    <div class="footer-inner">
      <p>Key Laboratory of Flexible Electronics (LoFE) & Institute of Advanced Materials (IAM)</p>
      <p>Nanjing University of Posts and Telecommunications, Nanjing, 210023, China</p>
      <p>Jiangsu Key Laboratory for Molecular and Medical Biotechnology, School of Life Science, Nanjing Normal University, Nanjing, 210024, China</p>
      <p>© All right Reversed. ICP9036104</p>
    </div>
  </div>
</template>

<script setup>
const currentYear = new Date().getFullYear()
</script>

<style scoped>
.site-footer {
  position: relative;
  width: 100%;
  min-height: 90px;
  background: #002855;
  color: rgba(255, 255, 255, 0.86);
  border-top: 3px solid #17b6ae;
}

.footer-inner {
  max-width: 1360px;
  margin: 0 auto;
  padding: 14px 24px 12px;
  text-align: center;
}

p {
  margin: 4px 0;
  font-size: 12px;
  line-height: 1.4;
  color: rgba(255, 255, 255, 0.76);
}

a {
  color: #b8fff7;
  font-weight: 600;
  text-decoration: none;
}

a:hover {
  color: #ffffff;
  text-decoration: underline;
}

@media (max-width: 900px) {
  .footer-inner {
    padding: 14px 16px;
  }
}
</style>
