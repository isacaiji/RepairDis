package com.bishe.ddr_springboot.service.network;

import com.bishe.ddr_springboot.entity.PpiResponse;
import com.bishe.ddr_springboot.mapper.network.PPINetworkMapper;
import org.springframework.stereotype.Service;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Service
public class PpiNetworkServiceImpl implements PpiNetworkService {

    private final PPINetworkMapper ppiNetworkMapper;

    public PpiNetworkServiceImpl(PPINetworkMapper ppiNetworkMapper) {
        this.ppiNetworkMapper = ppiNetworkMapper;
    }

    // 安全解析数字参数（兼容 Long/Integer）
    private long parseNumber(Object value, long defaultValue) {
        if (value == null) return defaultValue;
        if (value instanceof Number) return ((Number) value).longValue();
        try {
            return Long.parseLong(value.toString());
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }

    @Override
    public Map<String, Object> getNetworkData(Map<String, Object> params) {
        // 解析阈值参数
        int degreeThreshold = (int) parseNumber(params.get("degreeThreshold"), 0);

        // 查询节点计数和所有连接
        List<Map<String, Object>> node1Counts = ppiNetworkMapper.selectNode1Counts();
        List<Map<String, String>> allInteractions = ppiNetworkMapper.selectAllInteractions();

        // 筛选符合阈值的节点（修复Long转Integer问题）
        Set<String> qualifiedNodes = node1Counts.stream()
                .filter(countMap -> parseNumber(countMap.get("totalCount"), 0) >= degreeThreshold)
                .map(countMap -> countMap.get("name").toString())
                .collect(Collectors.toSet());

        // 筛选关联连接
        List<PpiResponse.Link> filteredLinks = allInteractions.stream()
                .filter(linkMap -> qualifiedNodes.contains(linkMap.get("source")))
                .map(linkMap -> {
                    PpiResponse.Link link = new PpiResponse.Link();
                    link.setSource(linkMap.get("source"));
                    link.setTarget(linkMap.get("target"));
                    return link;
                })
                .collect(Collectors.toList());

        // 构建节点列表
        Set<String> allNodeNames = filteredLinks.stream()
                .flatMap(link -> Stream.of(link.getSource(), link.getTarget()))
                .collect(Collectors.toSet());

        List<PpiResponse.Node> nodes = allNodeNames.stream()
                .map(name -> {
                    PpiResponse.Node node = new PpiResponse.Node();
                    node.setName(name);
                    return node;
                })
                .collect(Collectors.toList());

        // 计算最大计数（使用Long避免类型问题）
        long maxCount = node1Counts.stream()
                .mapToLong(countMap -> parseNumber(countMap.get("totalCount"), 0))
                .max()
                .orElse(0);

        // 封装结果
        Map<String, Object> result = new HashMap<>();
        result.put("nodes", nodes);
        result.put("links", filteredLinks);
        result.put("node1Counts", node1Counts);
        result.put("maxCount", maxCount);
        result.put("minCount", 0L);

        return result;
    }

    @Override
    public Map<String, Object> getTableData(Map<String, Object> params) {
        // 解析参数
        int pageNum = (int) parseNumber(params.get("pageNum"), 1);
        int pageSize = (int) parseNumber(params.get("pageSize"), 10);
        int degreeThreshold = (int) parseNumber(params.get("degreeThreshold"), 0);
        String selectedNode = params.getOrDefault("selectedNode", "").toString().trim();

        // 获取符合阈值的节点
        Set<String> qualifiedNodes = ppiNetworkMapper.selectNode1Counts().stream()
                .filter(countMap -> parseNumber(countMap.get("totalCount"), 0) >= degreeThreshold)
                .map(countMap -> countMap.get("name").toString())
                .collect(Collectors.toSet());

        // 查询所有原始记录并筛选
        List<Map<String, String>> allRecords = ppiNetworkMapper.selectAllOriginalRecords();
        List<Map<String, String>> filteredRecords = allRecords.stream()
                .filter(record -> qualifiedNodes.contains(record.get("node1")))
                .filter(record -> selectedNode.isEmpty()
                        || selectedNode.equals(record.get("node1"))
                        || selectedNode.equals(record.get("node2")))
                .collect(Collectors.toList());

        // 分页处理
        int total = filteredRecords.size();
        int startIndex = Math.max(0, (pageNum - 1) * pageSize);
        int endIndex = Math.min(startIndex + pageSize, total);
        List<Map<String, String>> pageRecords = startIndex < total
                ? filteredRecords.subList(startIndex, endIndex)
                : Collections.emptyList();

        // 封装分页结果
        Map<String, Object> pageResult = new HashMap<>();
        pageResult.put("total", total);
        pageResult.put("list", pageRecords);
        pageResult.put("pageNum", pageNum);
        pageResult.put("pageSize", pageSize);

        return pageResult;
    }

    @Override
    public byte[] exportAllData() {
        // 导出全表
        List<Map<String, Object>> all = ppiNetworkMapper.selectAllForExport();
        return convertToTsv(all);
    }

    @Override
    public byte[] exportCurrentData(String name, Double s) {
        // 导出当前查询结果
        List<Map<String, Object>> current = ppiNetworkMapper.selectForExport(name);
        return convertToTsv(current);
    }

    /**
     * 转换为TSV格式
     */
    private byte[] convertToTsv(List<Map<String, Object>> data) {
        if (data.isEmpty()) {
            return "No data to export".getBytes(StandardCharsets.UTF_8);
        }

        Set<String> headers = data.get(0).keySet();
        StringBuilder tsv = new StringBuilder(String.join("\t", headers) + "\n");

        for (Map<String, Object> row : data) {
            List<String> values = headers.stream()
                    .map(key -> row.get(key) != null ? row.get(key).toString() : "")
                    .collect(Collectors.toList());
            tsv.append(String.join("\t", values)).append("\n");
        }

        return tsv.toString().getBytes(StandardCharsets.UTF_8);
    }
}