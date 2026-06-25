package com.bishe.ddr_springboot.mapper;

import com.bishe.ddr_springboot.entity.DrugTarget;
import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;

import java.util.List;
import java.util.Map;

@Mapper
public interface DrugTargetMapper {

    String SELECT_FIELDS =
            "drugbank_id AS drugbankId, " +
                    "drug_name AS drugName, " +
                    "drug_type AS drugType, " +
                    "approved, " +
                    "target_name AS targetName, " +
                    "organism, " +
                    "action, " +
                    "gene_name AS geneName, " +
                    "uniprot_id AS uniprotId, " +
                    "polypeptide_name AS polypeptideName, " +
                    "polypeptide_id AS polypeptideId";

    /**
     * 查询所有药物靶点关系
     */
    @Select("SELECT " + SELECT_FIELDS + " FROM drug_targets")
    List<DrugTarget> selectAll();

    /**
     * 分页查询所有药物靶点
     */
    @Select("SELECT " + SELECT_FIELDS + " FROM drug_targets LIMIT #{startIndex}, #{pageSize}")
    List<DrugTarget> selectByPage(
            @Param("startIndex") Integer startIndex,
            @Param("pageSize") Integer pageSize
    );

    /**
     * 查询总记录数
     */
    @Select("SELECT COUNT(*) FROM drug_targets")
    Integer selectTotalCount();

    /**
     * 根据药物名称查询
     */
    @Select("SELECT " + SELECT_FIELDS + " FROM drug_targets WHERE drug_name = #{drugName}")
    List<DrugTarget> selectByDrugName(@Param("drugName") String drugName);

    /**
     * 根据基因名称查询
     */
    @Select("SELECT " + SELECT_FIELDS + " FROM drug_targets WHERE gene_name = #{geneName}")
    List<DrugTarget> selectByGeneName(@Param("geneName") String geneName);

    /**
     * 分页查询（按药物名称）
     */
    @Select("SELECT " + SELECT_FIELDS + " FROM drug_targets WHERE drug_name = #{drugName} LIMIT #{start}, #{pageSize}")
    List<DrugTarget> selectByDrugNamePage(
            @Param("drugName") String drugName,
            @Param("start") Integer start,
            @Param("pageSize") Integer pageSize
    );

    /**
     * 分页查询（按基因名称）
     */
    @Select("SELECT " + SELECT_FIELDS + " FROM drug_targets WHERE gene_name = #{geneName} LIMIT #{start}, #{pageSize}")
    List<DrugTarget> selectByGeneNamePage(
            @Param("geneName") String geneName,
            @Param("start") Integer start,
            @Param("pageSize") Integer pageSize
    );

    /**
     * 统计药物相关记录数
     */
    @Select("SELECT COUNT(*) FROM drug_targets WHERE drug_name = #{drugName}")
    Integer selectCountByDrugName(@Param("drugName") String drugName);

    /**
     * 统计基因相关记录数
     */
    @Select("SELECT COUNT(*) FROM drug_targets WHERE gene_name = #{geneName}")
    Integer selectCountByGeneName(@Param("geneName") String geneName);

    /**
     * 多条件查询（分页）
     */
    @Select("SELECT " + SELECT_FIELDS + " FROM drug_targets " +
            "WHERE (drug_name LIKE CONCAT('%', #{drugName}, '%') OR #{drugName} = '') " +
            "AND (gene_name LIKE CONCAT('%', #{geneName}, '%') OR #{geneName} = '') " +
            "AND (approved = #{approved} OR #{approved} = '') " +
            "LIMIT #{startIndex}, #{pageSize}")
    List<DrugTarget> selectByConditions(
            @Param("drugName") String drugName,
            @Param("geneName") String geneName,
            @Param("approved") String approved,
            @Param("startIndex") Integer startIndex,
            @Param("pageSize") Integer pageSize
    );

    /**
     * 多条件查询总记录数
     */
    @Select("SELECT COUNT(*) FROM drug_targets " +
            "WHERE (drug_name LIKE CONCAT('%', #{drugName}, '%') OR #{drugName} = '') " +
            "AND (gene_name LIKE CONCAT('%', #{geneName}, '%') OR #{geneName} = '') " +
            "AND (approved = #{approved} OR #{approved} = '')")
    Integer selectCountByConditions(
            @Param("drugName") String drugName,
            @Param("geneName") String geneName,
            @Param("approved") String approved
    );

    /**
     * 导出用（原始字段，按名称查询）
     */
    @Select("SELECT * FROM drug_targets WHERE drug_name = #{name} OR gene_name = #{name}")
    List<Map<String, Object>> selectForExport(@Param("name") String name);

    /**
     * 全量导出用（原始字段）
     */
    @Select("SELECT * FROM drug_targets")
    List<Map<String, Object>> selectAllForExport();

    /**
     * 按条件导出用（原始字段）
     */
    @Select("SELECT * FROM drug_targets " +
            "WHERE (drug_name LIKE CONCAT('%', #{drugName}, '%') OR #{drugName} = '') " +
            "AND (gene_name LIKE CONCAT('%', #{geneName}, '%') OR #{geneName} = '') " +
            "AND (approved = #{approved} OR #{approved} = '')")
    List<Map<String, Object>> selectByConditionsForExport(
            @Param("drugName") String drugName,
            @Param("geneName") String geneName,
            @Param("approved") String approved
    );

    /**
     * 药物名称模糊联想
     **/
    @Select("SELECT DISTINCT TRIM(drug_name) AS drugName " +
            "FROM drug_targets " +
            "WHERE drug_name LIKE CONCAT('%', #{keyword}, '%') " +
            "ORDER BY TRIM(drug_name) ASC")
    List<String> selectDrugNameSuggestions(@Param("keyword") String keyword);

    /**
     * 基因名称模糊联想
     */
    @Select("SELECT DISTINCT TRIM(gene_name) AS geneName " +
            "FROM drug_targets " +
            "WHERE gene_name LIKE CONCAT('%', #{keyword}, '%') " +
            "ORDER BY TRIM(gene_name) ASC")
    List<String> selectGeneNameSuggestions(@Param("keyword") String keyword);
}