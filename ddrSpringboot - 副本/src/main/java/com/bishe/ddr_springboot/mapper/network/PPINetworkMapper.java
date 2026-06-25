package com.bishe.ddr_springboot.mapper.network;

import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;
import java.util.List;
import java.util.Map;

/**
 * PPI网络数据访问层（仅静态SQL，无动态标签）
 */
@Mapper
public interface PPINetworkMapper {

    /**
     * 统计所有node1的出现次数
     */
    @Select("SELECT node1 AS name, COUNT(*) AS totalCount FROM ppi GROUP BY node1")
    List<Map<String, Object>> selectNode1Counts();

    /**
     * 查询所有互作关系（node1-node2）
     */
    @Select("SELECT node1 AS source, node2 AS target FROM ppi")
    List<Map<String, String>> selectAllInteractions();

    /**
     * 查询所有原始记录（用于表格展示）
     */
    @Select("SELECT " +
            "node1, node2, " +
            "node1_string_id, node2_string_id, " +
            "combined_score " +
            "FROM ppi")
    List<Map<String, String>> selectAllOriginalRecords();

    /**
     * 查询所有记录（用于全量导出）
     */
    @Select("SELECT " +
            "node1, node2, " +
            "node1_string_id, node2_string_id, " +
            "neighborhood_on_chromosome, gene_fusion, phylogenetic_cooccurrence, " +
            "homology, coexpression, experimentally_determined_interaction, " +
            "database_annotated, automated_textmining, combined_score " +
            "FROM ppi")
    List<Map<String, Object>> selectAllForExport();

    /**
     * 根据节点名称导出相关记录（用于当前查询结果导出）
     * 包含该节点作为node1或node2的所有记录
     */
    @Select("SELECT " +
            "node1, node2, " +
            "node1_string_id, node2_string_id, " +
            "neighborhood_on_chromosome, gene_fusion, phylogenetic_cooccurrence, " +
            "homology, coexpression, experimentally_determined_interaction, " +
            "database_annotated, automated_textmining, combined_score " +
            "FROM ppi " +
            "WHERE node1 = #{name} OR node2 = #{name}")
    List<Map<String, Object>> selectForExport(@Param("name") String name);
}