package com.bishe.ddr_springboot.entity;  // 包路径改为entity

import lombok.Data;

/**
 * TF调控关系实体类（映射Excel中的数据）
 */
@Data
public class TFRegulation {
    // 转录因子（TF）名称
    private String source;
    // 靶基因名称
    private String target;
    // 调控类型（Activation/Repression/Unknown）
    private String regulationType;
    // 证据来源（如PMID）
    private String evidence;
}