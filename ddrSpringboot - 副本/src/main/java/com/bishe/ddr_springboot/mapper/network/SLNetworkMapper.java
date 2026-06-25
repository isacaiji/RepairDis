package com.bishe.ddr_springboot.mapper.network;

import com.bishe.ddr_springboot.entity.SLResponse;
import org.apache.ibatis.annotations.Mapper;
import org.apache.ibatis.annotations.Select;
import java.util.List;

/**
 * 数据访问层：操作合并实体SLResponse（数据库映射部分）
 */
@Mapper
public interface SLNetworkMapper {

    /**
     * 查询全量数据（仅映射数据库字段：geneA、geneB、geminiSensitive、cellLine）
     */
    // 修改 SLNetworkMapper 的 @Select 语句：
    @Select("SELECT " +
            "`Gene A` AS geneA, " +
            "`Gene B` AS geneB, " +
            "`GEMINI sensitive` AS geminiSensitive, " +
            "`Cell line` AS cellLine " +
            "FROM SLnetwork")
    List<SLResponse> selectAll();
}