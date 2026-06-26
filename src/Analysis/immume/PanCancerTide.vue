<template>
  <PanCancerFigurePanel
      title="TIDE-predicted immunotherapy response"
      module-dir="K_TIDE_response"
      file-pattern="{cancer}_TIDE_score_response_combined"
  />
</template>

<script setup>
import PanCancerFigurePanel from './PanCancerFigurePanel.vue'
</script>
