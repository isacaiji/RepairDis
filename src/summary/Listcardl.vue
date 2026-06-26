<template>
  <div class="list-detail-body">
    <h1 style="text-align: center; margin: 10px;">{{geneData.geneName}} Gene Information</h1>
    <div class="list--detail-body-table" v-show="basicShow">
      <table class="list--detail-body-tablebody" >
        <tbody>
        <tr >
          <th>Gene:</th>
          <td style="color: #3a8ee6;"><div style="margin-left: 20px">{{geneData.geneName}}</div></td>
        </tr>
        <tr style="background-color: #f6f6f6">
          <th>Accession:</th>
          <td><div style="margin-left: 20px">{{geneData.accession}}</div></td>
        </tr>
        <tr>
          <th>Description:</th>
          <td><div style="margin-left: 20px">{{geneData.description}}</div></td>
        </tr>
<!--        <tr >-->
<!--          <th>EnsembleID:</th>-->
<!--          <td style="color: #3a8ee6;">-->
<!--            <button style="color:#3a8ee6;margin-left: 16px; border: 0;font-size:16px;background-color: #ffffff"  :underline="false" @click="genecard" >-->
<!--              {{this.rnadata.ensembl}}-->
<!--            </button>-->
<!--          </td>-->
<!--        </tr>-->
<!--        <tr style="background-color: #f6f6f6">-->
<!--          <th>Function:</th>-->
<!--          <td><div style="margin-left: 20px">{{this.rnadata.function}}</div></td>-->
<!--        </tr>-->
<!--        <tr>-->
<!--          <th>Target:</th>-->
<!--          <td><div style="margin-left: 20px">{{this.rnadata.target}}</div></td>-->
<!--        </tr>-->
<!--        <tr style="background-color: #f6f6f6">-->
<!--          <th>Location:</th>-->
<!--          <td><div style="margin-left: 20px">{{this.rnadata.location}}</div></td>-->
<!--        </tr>-->
<!--        <tr>-->
<!--          <th>Title:</th>-->
<!--          <td>-->
<!--            <div style="margin-left: 20px">-->
<!--              {{this.rnadata.title}}-->
<!--            </div>-->

<!--          </td>-->
<!--        </tr>-->
<!--        <tr style="background-color: #f6f6f6">-->
<!--          <th>Method:</th>-->
<!--          <td><div style="margin-left: 20px">{{this.rnadata.method}}</div></td>-->
<!--        </tr>-->
<!--        <tr>-->
<!--          <th>Cells:</th>-->
<!--          <td><div style="margin-left: 20px">{{this.rnadata.cells}}</div></td>-->
<!--        </tr>-->
        <tr style="background-color: #f6f6f6">
          <th>proteinSequences:</th>
          <td><div style="margin-left: 20px">{{geneData.proteinSequences}}</div></td>
        </tr>
<!--        <tr>-->
<!--          <th>PMID:</th>-->
<!--          <td>-->
<!--            <button style="color:#3a8ee6;margin-left: 16px; border: 0;font-size:16px;background-color: #ffffff"  :underline="false" @click="clickview" >-->
<!--              {{this.rnadata.pmid}}-->
<!--            </button>-->
<!--          </td>-->
<!--        </tr>-->
        </tbody>
      </table>


    </div>
  </div>

</template>

<script setup>
import {onMounted, ref} from 'vue'
import axios from "axios";
//基因信息，为json
const geneData = ref([]);
//定义props
const props = defineProps({
  id: {
    type: [String, Number],
    required: true
  },
  basicShow: {
    type: Boolean,
    default: false
  }
});

onMounted(()=> {
  getGeneData(props.id);
});
//获取单个基因信息
function getGeneData(id){
  axios.get("http://121.37.88.191:9016/api/genes/" + id).then((res)=>{
    geneData.value = res.data;
    console.log("===================");
  })
}

//打开PMID
function clickView(){
  window.open('https://pubmed.ncbi.nlm.nih.gov/' + geneData.value.pmid + '/')
}
//打开geneCard
function geneCard(){
  const genecardsUrl = `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${geneData.value.geneName}&keywords=${geneData.value.geneName}`;
  window.open(genecardsUrl, '_blank');
}

defineExpose({
  getGeneData,
  geneCard,
});
</script>

<style scoped>
body, div, dl, dt, dd, ul, ol, li, h1, h2, h3, h4, h5, h6, pre, form, fieldset, legend, input, textarea, button, p, blockquote, th, td
{margin: 0;}
body{text-align: center;font-family: Helvetica Neue,Helvetica,Arial,Microsoft Yahei,Hiragino Sans GB,Heiti SC,WenQuanYi Micro Hei,sans-serif;}
li{ list-style: none;}
a{text-decoration: none;
  color: #5a5a5a
}
button:hover{
  text-decoration: underline; /* 添加下划线 */
}
img{border: none;}
.list-detail-body{
  position: relative;
  /*display: flex;*/
  /*flex-direction: column;*/
  /*justify-content: space-between;*/
  left: 0;

  width: 100%;
  height: 45vw;
}
.list--detail-body-table{
  position: relative;
  width: 100%;
  height: 90%;
  border-radius: 5px;
  background-color: #ffffff;
  margin-bottom: 100px;
}
/*表格在div块中的位置*/
.list--detail-body-tablebody{
  position: relative;
  top: 6%;
  left: 3.4%;
  width: 90%;
  height: 98%;
  border: 1px solid #c5c5c5;
  border-collapse: collapse;/*表格合成单一线框*/
}
th{
  text-align:left;
  text-indent: 10px;
  font-size: 16px;
  width: 20%;
  font-weight: initial;
  border: 1px solid #c5c5c5;
  color: #6c6c6c;
}
td{
  width: 70%;
  border: 1px solid #c5c5c5;
  font-size: 16px;
  text-align: left;
  /*text-indent: 7%;*/
  color: #2f2f2f;
}
/*.el-button-viewpaper{*/
/*  position: relative;*/
/*  top:63px;*/
/*  color: #4a91da;*/
/*  border-radius: 7px;*/
/*}*/
/*.el-button-viewpaper:focus{*/
/*  background: #ffffff;*/
/*}*/
</style>
