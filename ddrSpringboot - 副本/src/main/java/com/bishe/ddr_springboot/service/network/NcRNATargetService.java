package com.bishe.ddr_springboot.service.network;

import java.util.List;
import java.util.Map;

/**
 * ncRNA靶向关系服务接口（遵循通用网络服务规范）
 */
public interface NcRNATargetService extends NetworkService {

    /**
     * 获取网络拓扑数据（节点和连接）
     * @param params 筛选参数（包含ncRNA名称等）
     * @return 网络数据（nodes、links等）
     */
    @Override
    Map<String, Object> getNetworkData(Map<String, Object> params);

    /**
     * 获取表格分页数据
     * @param params 分页参数（pageNum、pageSize）和筛选参数
     * @return 分页结果（total、list、pageNum、pageSize）
     */
    @Override
    Map<String, Object> getTableData(Map<String, Object> params);

    /**
     * 导出全部符合条件的数据
     *
     * @return TSV格式字节数组
     */
    @Override
    byte[] exportAllData();

    /**
     * 导出选中节点相关的数据
     *
     *
     * @return TSV格式字节数组
     */
    @Override
    byte[] exportCurrentData(String name, Double s);
}