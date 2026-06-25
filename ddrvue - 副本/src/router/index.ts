import {createRouter, createWebHistory, RouteRecordRaw} from 'vue-router';

const routes :Array<RouteRecordRaw> = [
    {
        path: '/',
        redirect: '/home'
    },
    {
        path: '/home',
        name: 'Home',
        component: () => import('@/views/Home.vue')
    },
    {
        path: '/browser',
        name: 'Browser',
        component: () => import('@/views/Browser.vue')
    },
    {
        path: '/analysis',
        name: 'Analysis',
        component: () => import('@/views/Analysis.vue'),
        // children: [
        //     // 免疫相关组件（7个）
        //     {
        //         path: 'immu/checkpoint',
        //         name: 'ImmuCheckpoint',
        //         component: () => import('@/Analysis/immume/checkpoint.vue')
        //     },
        //     {
        //         path: 'immu/chemokine',
        //         name: 'ImmuChemokine',
        //         component: () => import('@/Analysis/immume/chemokine.vue')
        //     },
        //     {
        //         path: 'immu/immucell',
        //         name: 'ImmuCell',
        //         component: () => import('@/Analysis/immume/immucell.vue')
        //     },
        //     {
        //         path: 'immu/immunescore',
        //         name: 'ImmuneScore',
        //         component: () => import('@/Analysis/immume/immunescore.vue')
        //     },
        //     {
        //         path: 'immu/immuinhibitor',
        //         name: 'ImmuInhibitor',
        //         component: () => import('@/Analysis/immume/immuinhibitor.vue')
        //     },
        //     {
        //         path: 'immu/immustimulator',
        //         name: 'ImmuStimulator',
        //         component: () => import('@/Analysis/immume/immustimulator.vue')
        //     },
        //     {
        //         path: 'immu/receptor',
        //         name: 'ImmuReceptor',
        //         component: () => import('@/Analysis/immume/receptor.vue')
        //     },
        //     // dzx
        //     {
        //         path: 'dzx/cnv',
        //         name: 'DzxCnv',
        //         component: () => import('@/Analysis/dzx/Cnv.vue')
        //     },
        //     {
        //         path: 'dzx/diff',
        //         name: 'DzxDiff',
        //         component: () => import('@/Analysis/dzx/DiffR.vue')
        //     },
        //     {
        //         path: 'dzx/gsea',
        //         name: 'DzxGsea',
        //         component: () => import('@/Analysis/dzx/GSEA.vue')
        //     },
        //     {
        //         path: 'dzx/methy',
        //         name: 'DzxMethy',
        //         component: () => import('@/Analysis/dzx/Methy.vue')
        //     },
        //     {
        //         path: 'dzx/mut',
        //         name: 'DzxMut',
        //         component: () => import('@/Analysis/dzx/Mut.vue')
        //     },
        //     {
        //         path: 'dzx/surv',
        //         name: 'DzxSurv',
        //         component: () => import('@/Analysis/dzx/Surv.vue')
        //     }
        // ]
    },
    {
        path: '/drug',
        name: 'Drug',
        component: () => import('@/views/Drug.vue'),
    },
    {
        path: '/sl',
        name: 'SL',
        component: () => import('@/views/SL.vue'),
        children: [
            {
                path: 'network',
                name: 'SLNetwork',
                component: () => import('@/sl/SLNetwork.vue')
            },
        ]
    },
    {
        path: '/network',
        name: 'Network',
        redirect: '/network/ppi',
        children: [
            {
                path: 'ppi',
                name: 'PpiNetwork',
                component: () => import('@/network/ppi.vue')
            },
            {
                path: 'tf',
                name: 'TF',
                component: () => import('@/network/TF.vue')
            },
            {
                path: 'ncrna',
                name: 'ncRNA',
                component: () => import('@/network/ncRNA.vue')
            },
            {
                path: 'cross',
                name: 'Cross',
                component: () => import('@/network/Cross.vue')
            }
        ]
    },
    {
        path: '/evolution',
        name: 'Evolution',
        component : () => import('@/Evolution/EvolutionAll.vue'),
    },
    {
        path: '/help',
        name: 'Help',
        component: () => import('@/views/Help.vue')
    },
    {
        path: '/contract_us',
        name: 'ContractUs',
        component: () => import('@/views/ContractUs.vue')
    },
    {
        path:'/detail',
        name:'Detail',
        component: () => import('@/views/GeneDetail.vue')
    }
];

const router = createRouter({
    history: createWebHistory('/DDRAD'),
    routes
});

export default router;
