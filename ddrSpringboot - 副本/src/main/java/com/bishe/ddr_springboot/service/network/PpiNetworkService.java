package com.bishe.ddr_springboot.service.network;

import com.bishe.ddr_springboot.entity.PpiResponse;
import java.util.List;
import java.util.Map;

/**
 * PPI网络服务接口
 */
public interface PpiNetworkService extends NetworkService {

    /**
     * 获取网络数据（含节点和连接）
     * @param params 筛选参数（如degreeThreshold）
     * @return 网络响应实体
     */
    Map<String, Object> getNetworkData(Map<String, Object> params);

    /**
     * 获取分页表格数据
     * @param params 分页参数（pageNum, pageSize）和筛选参数
     * @return 分页结果（total, list, pageNum, pageSize）
     */
    Map<String, Object> getTableData(Map<String, Object> params);

    /**
     * 全量导出数据
     * @return TSV格式字节数组
     */
    byte[] exportAllData();

    /**
     * 导出选中节点数据
     * @return TSV格式字节数组
     */
    byte[] exportCurrentData(String name, Double s);
}