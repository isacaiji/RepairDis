package com.bishe.ddr_springboot.entity;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import java.util.List;

/**
 * 分页结果封装类，统一返回分页数据格式
 */
@Data
@NoArgsConstructor
@AllArgsConstructor
public class PageResult<T> {
    // 当前页数据列表
    private List<T> list;
    // 总数据条数
    private Integer total;
    // 当前页码
    private Integer pageNum;
    // 每页条数
    private Integer pageSize;
}
