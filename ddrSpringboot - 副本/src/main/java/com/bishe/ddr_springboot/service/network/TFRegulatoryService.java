package com.bishe.ddr_springboot.service.network;

import com.bishe.ddr_springboot.entity.TFRegulation;
import java.util.List;
import java.util.Map;

/**
 * TF调控关系服务接口（继承通用网络服务接口）
 */
public interface TFRegulatoryService extends NetworkService {

    /**
     * 实现通用接口：获取网络拓扑数据（节点和连接）
     * @param params 筛选参数（如基因名称、调控类型等）
     * @return 网络数据（包含节点、连接等）
     */
    @Override
    Map<String, Object> getNetworkData(Map<String, Object> params);

    /**
     * 实现通用接口：获取表格分页数据
     * @param params 分页参数（pageNum, pageSize）和筛选参数
     * @return 分页结果（total, list, pageNum, pageSize）
     */
    @Override
    Map<String, Object> getTableData(Map<String, Object> params);

    /**
     * 实现通用接口：导出全部数据
     * @return TSV格式字节数组
     */
    @Override
    byte[] exportAllData();

    /**
     * 实现通用接口：导出选中节点数据
     * @return TSV格式字节数组
     */
    @Override
    byte[] exportCurrentData(String name, Double s);

    // 以下为TF调控关系特有的扩展方法
    /**
     * 根据基因查询相关调控关系（作为TF或靶基因）
     * @param gene 基因名称
     * @return 调控关系列表
     */
    List<TFRegulation> getRegulationsByGene(String gene);

    /**
     * 根据TF名称查询其调控的靶基因
     * @param tfName TF名称
     * @return 调控关系列表
     */
    List<TFRegulation> getRegulationsByTF(String tfName);

    /**
     * 根据靶基因查询调控它的TF
     * @param targetGene 靶基因名称
     * @return 调控关系列表
     */
    List<TFRegulation> getRegulatorsOfTarget(String targetGene);
}