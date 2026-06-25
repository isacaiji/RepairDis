package com.bishe.ddr_springboot.service;

import com.bishe.ddr_springboot.entity.CancerScore;
import com.bishe.ddr_springboot.entity.Gene;
import com.bishe.ddr_springboot.entity.GeneSummary;
import com.bishe.ddr_springboot.entity.PageResult;
import com.bishe.ddr_springboot.mapper.GeneMapper;
import org.springframework.stereotype.Service;
import java.util.List;

@Service
public class GeneServiceImpl implements GeneService {

    private final GeneMapper geneMapper;
    private final MoDdrWeightService moDdrWeightService;

    public GeneServiceImpl(GeneMapper geneMapper, MoDdrWeightService moDdrWeightService) {
        this.geneMapper = geneMapper;
        this.moDdrWeightService = moDdrWeightService;
    }

    @Override
    public PageResult<GeneSummary> getGeneSummariesByPage(Integer pageNum, Integer pageSize) {
        // 参数校验
        if (pageNum == null || pageNum < 1) pageNum = 1;
        if (pageSize == null || pageSize < 1 || pageSize > 100) pageSize = 15;

        // 计算起始索引
        Integer startIndex = (pageNum - 1) * pageSize;

        // 查询数据
        List<GeneSummary> summaryList = geneMapper.findGeneSummariesByPage(startIndex, pageSize);
        Integer total = geneMapper.findGeneTotalCount();

        // 封装结果
        return new PageResult<>(summaryList, total, pageNum, pageSize);
    }

    @Override
    public List<GeneSummary> getGeneSummaries() {
        return geneMapper.findAllGeneSummaries();
    }

    @Override
    public List<Gene> getGenes() {
        return geneMapper.findAllGenes();
    }

    @Override
    public Gene getGeneById(Integer id) {
        Gene gene = geneMapper.findGeneById(id);
        if (gene != null && gene.getGeneName() != null) {
            Double score = moDdrWeightService.getTotalScoreByGene(gene.getGeneName());
            gene.setMeanMoDdrWeight(score);

        }
        return gene;
    }

    @Override
    public Gene getGeneByName(String name) {
        Gene gene = geneMapper.findGeneByName(name);
        if (gene != null && gene.getGeneName() != null) {
            Double score = moDdrWeightService.getTotalScoreByGene(gene.getGeneName());
            gene.setMeanMoDdrWeight(score);
        }
        return gene;
    }

    @Override
    public List<Gene> getGenesByQuery(String query) {
        return geneMapper.findGenesByQuery(query);
    }

    @Override
    public List<String> getAllGeneNames() {
        return geneMapper.findAllGeneNames();
    }

    @Override
    public CancerScore getCancerScoresByGeneName(String geneName) {
        return moDdrWeightService.getCancerScoreByGene(geneName);
    }

}
