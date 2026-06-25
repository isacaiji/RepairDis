package com.bishe.ddr_springboot.mapper.network;

import com.bishe.ddr_springboot.entity.NcRNATarget;
import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;

import java.util.List;
import java.util.Map;

@Mapper
public interface NcRNATargetMapper {

    String SELECT_FIELDS =
            "`miRTarBase ID` AS mirTarBaseId, " +
                    "`miRNA` AS ncRNA, " +
                    "`Species (miRNA)` AS ncRNASpecies, " +
                    "`Target Gene` AS targetGene, " +
                    "`Target Gene (Entrez ID)` AS targetGeneEntrezId, " +
                    "`Species (Target Gene)` AS targetGeneSpecies, " +
                    "Experiments AS experiments, " +
                    "`Support Type` AS supportType, " +
                    "`References (PMID)` AS reference";

    // 全量（用于缓存等，可选）
    @Select("SELECT " + SELECT_FIELDS + " FROM miRNA")
    List<NcRNATarget> selectAll();

    // 按 ncRNA 查询
    @Select("SELECT " + SELECT_FIELDS + " FROM miRNA WHERE `miRNA` = #{name}")
    List<NcRNATarget> selectByNcRNAName(@Param("name") String name);

    @Select("SELECT " + SELECT_FIELDS + " FROM miRNA WHERE `miRNA` = #{name} LIMIT #{start}, #{pageSize}")
    List<NcRNATarget> selectByNcRNAnamePage(
            @Param("name") String name,
            @Param("start") Integer start,
            @Param("pageSize") Integer pageSize
    );

    @Select("SELECT COUNT(*) FROM miRNA WHERE `miRNA` = #{name}")
    Integer selectCountByNcRNAname(@Param("name") String name);

    // 按靶基因查询（新增）
    @Select("SELECT " + SELECT_FIELDS + " FROM miRNA WHERE `Target Gene` = #{targetGene}")
    List<NcRNATarget> selectByTargetGene(@Param("targetGene") String targetGene);

    @Select("SELECT " + SELECT_FIELDS + " FROM miRNA WHERE `Target Gene` = #{targetGene} LIMIT #{start}, #{pageSize}")
    List<NcRNATarget> selectByTargetGenePage(
            @Param("targetGene") String targetGene,
            @Param("start") Integer start,
            @Param("pageSize") Integer pageSize
    );

    @Select("SELECT COUNT(*) FROM miRNA WHERE `Target Gene` = #{targetGene}")
    Integer selectCountByTargetGene(@Param("targetGene") String targetGene);

    // 导出用（返回原始列，Map）
    @Select("SELECT * FROM miRNA WHERE `miRNA` = #{name} OR `Target Gene` = #{name}")
    List<Map<String, Object>> selectForExport(@Param("name") String name);

    // 新增：全量导出用（返回所有记录的原始字段）
    @Select("SELECT * FROM miRNA")
    List<Map<String, Object>> selectAllForExport();
}