// 基础URL配置
const baseURL = 'http://121.37.88.191:83';
const myURL = 'http://121.37.88.191:9016';

// API端点
const xdlURL = baseURL + '/xdl';
const dzxURL = baseURL + '/dzx';
const rURL = baseURL + '/r';

const cancers = [
    'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM',
    'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC',
    'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD',
    'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM'
] as const;

// 导出类型
export type CancerType = typeof cancers[number];

// 导出配置
export default {
    dzxURL,
    xdlURL,
    rURL,
    cancers,
    myURL,
} as const;
