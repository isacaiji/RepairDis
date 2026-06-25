package com.bishe.ddr_springboot.entity;

import lombok.Data;
import java.util.List;
import java.util.Map;

/**
 * ncRNA靶向调控网络合并实体类
 * 兼具数据库映射和接口返回功能
 */
@Data
public class NcRNATarget {
    // ------------------------------
    // 1. 数据库映射字段
    // ------------------------------
    private String mirTarBaseId;         // miRTarBase ID
    private String ncRNA;                // ncRNA名称
    private String ncRNASpecies;         // ncRNA物种
    private String targetGene;           // 靶基因名称
    private String targetGeneEntrezId;   // 靶基因Entrez ID
    private String targetGeneSpecies;    // 靶基因物种
    private String experiments;          // 实验方法
    private String supportType;          // 支持类型
    private String reference;           // 参考文献PMID
    private String type;                 // ncRNA类型（miRNA/lncRNA等）

    // ------------------------------
    // 2. 接口返回字段（网络数据）
    // ------------------------------
    private List<Node> nodes;            // 网络节点
    private List<Link> links;            // 网络连接
    private List<Map<String, String>> tableData; // 表格数据
    private Integer total;               // 总数（分页用）

    /**
     * 节点内部类（前端可视化用）
     */
    @Data
    public static class Node {
        private String name;             // 节点名称（ncRNA或基因名）
        private String type;             // 类型：query_ncRNA/ncRNA/target_gene
        private String symbol;           // 图形符号
        private Integer symbolSize;      // 图形大小
        private Map<String, String> itemStyle; // 样式（颜色等）
    }

    /**
     * 边内部类（前端可视化用）
     */
    @Data
    public static class Link {
        private String source;           // 源节点（ncRNA）
        private String target;           // 目标节点（靶基因）
        private Map<String, Object> lineStyle; // 边样式
    }
}