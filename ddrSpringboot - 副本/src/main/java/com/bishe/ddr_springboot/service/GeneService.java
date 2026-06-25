package com.bishe.ddr_springboot.service;

import com.bishe.ddr_springboot.entity.CancerScore;
import com.bishe.ddr_springboot.entity.Gene;
import com.bishe.ddr_springboot.entity.GeneSummary;
import com.bishe.ddr_springboot.entity.PageResult;
import java.util.List;

public interface GeneService {
    /**
     * 分页查询基因摘要
     */
    PageResult<GeneSummary> getGeneSummariesByPage(Integer pageNum, Integer pageSize);

    /**
     * 查询所有基因摘要（全量）
     */
    List<GeneSummary> getGeneSummaries();

    /**
     * 查询所有基因完整信息
     */
    List<Gene> getGenes();

    /**
     * 根据ID查询单个基因完整信息
     */
    Gene getGeneById(Integer id);

    /**
     * 根据名称查询单个基因完整信息
     */
    Gene getGeneByName(String name);

    /**
     * 根据关键词搜索基因
     */
    List<Gene> getGenesByQuery(String query);

    /**
     * 查询所有基因名称
     */
    List<String> getAllGeneNames();

    /**
     * 根据基因名称查询分数
     */
    CancerScore getCancerScoresByGeneName(String geneName);
}
