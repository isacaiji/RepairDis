package com.bishe.ddr_springboot.entity;

import lombok.Data;

/**
 * 药物-靶点关系实体类
 */
@Data
public class DrugTarget {
    private String drugbankId;         // DrugBank ID
    private String drugName;           // 药物名称
    private String drugType;           // 药物类型
    private String approved;           // 批准状态
    private String targetName;         // 靶点名称
    private String organism;           // 物种
    private String action;             // 作用方式
    private String geneName;           // 基因名称
    private String uniprotId;          // UniProt ID
    private String polypeptideName;    // 多肽名称
    private String polypeptideId;      // 多肽ID
}