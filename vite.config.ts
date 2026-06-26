import { defineConfig } from 'vite';
import vue from '@vitejs/plugin-vue';
import path from 'path';
import vueJsx from '@vitejs/plugin-vue-jsx';
import api from "./src/api";

const aliasPath = path.resolve(__dirname, 'src');
console.log('Resolved alias path:', aliasPath); // 打印解析后的路径

export default defineConfig({
  base: './',
  plugins: [
      vue(),
      vueJsx()
  ],
  server: {
    proxy: {
      '/api': {
        target: api.myURL, // 后端服务器地址
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, '')
      }
    },
  },
  resolve: {
    //配置@为./src
    alias: {
      '@': path.resolve(__dirname, 'src')
    }
  }
});