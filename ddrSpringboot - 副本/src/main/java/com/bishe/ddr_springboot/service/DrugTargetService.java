package com.bishe.ddr_springboot.service;

import com.bishe.ddr_springboot.entity.DrugTarget;
import com.bishe.ddr_springboot.entity.PageResult;

import java.util.List;
import java.util.Map;

/**
 * 药物靶点关系服务接口（独立服务，不依赖NetworkService）
 */
public interface DrugTargetService {

    /**
     * 查询所有药物靶点关系
     */
    List<DrugTarget> getAllDrugTargets();

    /**
     * 分页查询药物靶点关系
     * @param pageNum 页码
     * @param pageSize 每页条数
     * @return 分页结果
     */
    PageResult<DrugTarget> getDrugTargetsByPage(Integer pageNum, Integer pageSize);

    /**
     * 根据药物名称查询
     */
    List<DrugTarget> getByDrugName(String drugName);

    /**
     * 根据基因名称查询
     */
    List<DrugTarget> getByGeneName(String geneName);

    /**
     * 条件筛选查询（支持多参数组合）
     */
    PageResult<DrugTarget> queryByConditions(Map<String, Object> params);

    /**
     * 导出所有数据为TSV格式
     */
    byte[] exportAll();

    /**
     * 按条件导出数据为TSV格式
     */
    byte[] exportByConditions(Map<String, Object> params);

    /**
     * @param params 包含查询参数（如name：药物名或基因名）
     * @return 网络数据（节点、连接等）
     */
    Map<String, Object> getNetworkData(Map<String, Object> params);

    List<String> getDrugNameSuggestions(String keyword);
    List<String> getGeneNameSuggestions(String keyword);
}