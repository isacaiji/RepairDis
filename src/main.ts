//引入createApp用于创建应用
import {createApp} from 'vue'
//引入APP根组件
import APP from './App.vue'
//引入路由
import router from './router'
import ElementPlus from 'element-plus';
import 'element-plus/dist/index.css';

const app =  createApp(APP);
app.use(router);
app.use(ElementPlus);
app.mount('#app');