<template>
  <div>
    <el-card class="box-card">
      <div style="text-align: center">
        <span style="font-size: 50px">Contact Information</span>
      </div>
      <hr />
      <div style="text-align: left">
        <p>
          Name: Zhang yan
          <hr />
          Address: Nanjing University of Posts and Telecommunications
          <hr />
          Email: 2939439845@qq.com
        </p>
      </div>
    </el-card>
  </div>
</template>

<script setup>
import { ElCard } from 'element-plus';
</script>

<style scoped>
.box-card {
  margin-top: 100px;
  margin-bottom: 100px;
  border-width: 2px;
  border-color: #333745;
  border-radius: 30px;
}
</style>