package com.bishe.ddr_springboot.service.network;

import java.util.Map;

public interface NetworkService {
    // 获取网络拓扑数据（节点和连接）
    Map<String, Object> getNetworkData(Map<String, Object> params);

    // 获取表格分页数据
    Map<String, Object> getTableData(Map<String, Object> params);

    // 导出全部数据
    byte[] exportAllData();

    // 导出选中节点数据
    byte[] exportCurrentData(String name, Double s);
}
