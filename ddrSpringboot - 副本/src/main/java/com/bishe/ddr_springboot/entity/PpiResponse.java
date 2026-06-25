package com.bishe.ddr_springboot.entity;

import lombok.Data;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * PPI网络数据返回实体
 * 包含前端所需的所有字段：nodes、links、node1Counts、maxCount、minCount、originalRecords
 */
@Data
public class PpiResponse {

    private List<Node> nodes;

    private List<Link> links;

    private Map<String, Integer> node1Counts;

    private Integer maxCount;

    private Integer minCount = 0;

    // 存储筛选后的原始数据记录（供前端表格展示）
    private List<Map<String, String>> originalRecords;

    // 新增：存储符合阈值条件的节点名称集合（用于导出时筛选TSV行）
    private Set<String> qualifiedNodeNames;

    @Data
    public static class Node {
        private String name; // 蛋白质ID（与前端node.name对应）
    }

    @Data
    public static class Link {
        private String source; // 源节点名称（对应Node.name）
        private String target; // 目标节点名称（对应Node.name）
    }
}
