package com.bishe.ddr_springboot.entity;

import lombok.Data;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * 合并后的SL实体类（兼具数据库映射和接口返回功能）
 */
@Data
public class SLResponse {
    // ------------------------------
    // 1. 数据库映射字段（对应SLNetwork）
    // ------------------------------
    private String geneA;             // 数据库字段：基因A
    private String geneB;             // 数据库字段：基因B
    private String geminiSensitive;   // 数据库字段：敏感性评分（字符串）
    private String cellLine;          // 数据库字段：细胞系

    // ------------------------------
    // 2. 接口返回字段（对应SLResponse）
    // ------------------------------
    private List<Node> nodes;         // 网络节点
    private List<Link> links;         // 网络连接
    private Map<String, Double> sensitivityStats; // 统计信息
    private Double maxSensitivity;    // 最大敏感性
    private Double minSensitivity;    // 最小敏感性
    private List<Map<String, String>> originalRecords; // 原始记录
    private Set<String> qualifiedGeneNames; // 符合条件的基因

    // 内部类：节点和连接（前端展示用）
    @Data
    public static class Node {
        private String name;
    }

    @Data
    public static class Link {
        private String source;
        private String target;
        private Double sensitivity; // 转换后的Double类型
        private String cellLine;
    }
}