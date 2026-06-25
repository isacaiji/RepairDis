package com.bishe.ddr_springboot.mapper;

import com.bishe.ddr_springboot.entity.Gene;
import com.bishe.ddr_springboot.entity.GeneSummary;
import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Param;
import org.apache.ibatis.annotations.Select;
import java.util.List;

@Mapper
public interface GeneMapper {
    /**
     * 分页查询基因摘要
     */
    @Select("SELECT id, geneName, ensembl AS ensemblId, pathway, pmid FROM ddrgenes LIMIT #{startIndex}, #{pageSize}")
    List<GeneSummary> findGeneSummariesByPage(
            @Param("startIndex") Integer startIndex,
            @Param("pageSize") Integer pageSize
    );

    /**
     * 查询所有基因摘要（全量）
     */
    @Select("SELECT id, geneName, ensembl AS ensemblId, pathway, pmid FROM ddrgenes")
    List<GeneSummary> findAllGeneSummaries();

    /**
     * 查询总条数
     */
    @Select("SELECT COUNT(*) FROM ddrgenes")
    Integer findGeneTotalCount();

    /**
     * 查询所有基因完整信息
     */
    @Select("SELECT * FROM ddrgenes")
    List<Gene> findAllGenes();

    /**
     * 查询所有基因名称
     */
    @Select("SELECT geneName FROM ddrgenes")
    List<String> findAllGeneNames();

    /**
     * 根据ID查询单个基因完整信息
     */
    @Select("SELECT * FROM ddrgenes WHERE id = #{id}")
    Gene findGeneById(@Param("id") Integer id);

    /**
     * 根据名称查询单个基因完整信息
     */
    @Select("SELECT * FROM ddrgenes WHERE geneName = #{name}")
    Gene findGeneByName(@Param("name") String name);

    /**
     * 根据关键词模糊查询基因
     */
    @Select("SELECT * FROM ddrgenes WHERE geneName LIKE CONCAT('%', #{query}, '%')")
    List<Gene> findGenesByQuery(@Param("query") String query);

}
