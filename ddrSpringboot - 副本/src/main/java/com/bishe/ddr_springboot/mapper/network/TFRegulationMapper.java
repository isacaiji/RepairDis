package com.bishe.ddr_springboot.mapper.network;

import com.bishe.ddr_springboot.entity.TFRegulation;
import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;

import java.util.List;
import java.util.Map;

/**
 * TF调控关系数据访问层（对应表格：TF,Target,Mode of Regulation,References (PMID)）
 */
@Mapper
public interface TFRegulationMapper {

    /**
     * 查询所有TF调控关系
     */
    @Select("SELECT " +
            "TF AS source, " +
            "Target AS target, " +
            "`Mode of Regulation` AS regulationType, " +
            "`References (PMID)` AS evidence " +
            "FROM TF")
    List<TFRegulation> selectAll();

    /**
     * 根据基因查询相关的调控关系
     * 包含该基因作为TF（source）或靶基因（target）的情况
     */
    @Select("SELECT " +
            "TF AS source, " +
            "Target AS target, " +
            "`Mode of Regulation` AS regulationType, " +
            "`References (PMID)` AS evidence " +
            "FROM TF " +
            "WHERE TF = #{gene} OR Target = #{gene}")
    List<TFRegulation> selectByGene(@Param("gene") String gene);

    /**
     * 根据TF名称查询其调控的靶基因
     */
    @Select("SELECT " +
            "TF AS source, " +
            "Target AS target, " +
            "`Mode of Regulation` AS regulationType, " +
            "`References (PMID)` AS evidence " +
            "FROM TF " +
            "WHERE TF = #{tfName}")
    List<TFRegulation> selectByTF(@Param("tfName") String tfName);

    /**
     * 根据靶基因查询调控它的TF
     */
    @Select("SELECT " +
            "TF AS source, " +
            "Target AS target, " +
            "`Mode of Regulation` AS regulationType, " +
            "`References (PMID)` AS evidence " +
            "FROM TF " +
            "WHERE Target = #{targetGene}")
    List<TFRegulation> selectByTargetGene(@Param("targetGene") String targetGene);

    /**
     * 查询所有记录用于导出（保留原始字段名）
     */
    @Select("SELECT " +
            "TF, Target, `Mode of Regulation`, `References (PMID)` " +
            "FROM TF")
    List<Map<String, Object>> selectAllForExport();
}