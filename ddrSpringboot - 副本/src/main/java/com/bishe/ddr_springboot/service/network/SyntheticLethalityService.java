package com.bishe.ddr_springboot.service.network;

import java.util.Map;

/**
 * 合成致死网络服务接口（继承通用接口，保持一致性）
 */
public interface SyntheticLethalityService extends NetworkService {

    /**
     * 获取网络拓扑数据（节点和连接等）
     * @param params 筛选参数（如sensitivityThreshold）
     * @return 网络响应数据，包含节点、连接、统计信息等
     */
    Map<String, Object> getNetworkData(Map<String, Object> params);

    /**
     * 获取表格分页数据
     * @param params 分页参数（pageNum, pageSize）和筛选参数（sensitivityThreshold, gene）
     * @return 分页结果，包含总条数、当前页数据、页码和页大小
     */
    Map<String, Object> getTableData(Map<String, Object> params);

    byte[] exportAllData();

    /**
     * 导出选中节点相关的数据
     * @return TSV格式字节数组
     */
    byte[] exportCurrentData(String name, Double s);
}