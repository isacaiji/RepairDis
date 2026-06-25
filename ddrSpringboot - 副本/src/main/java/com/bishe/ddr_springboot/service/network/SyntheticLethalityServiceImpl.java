package com.bishe.ddr_springboot.service.network;

import com.bishe.ddr_springboot.entity.SLResponse;
import com.bishe.ddr_springboot.mapper.network.SLNetworkMapper;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Service
@Slf4j
public class SyntheticLethalityServiceImpl implements SyntheticLethalityService {

    private final SLNetworkMapper slNetworkMapper;
    // 缓存全量数据（使用合并实体的数据库字段）
    private List<SLResponse> allSLData;

    public SyntheticLethalityServiceImpl(SLNetworkMapper slNetworkMapper) {
        this.slNetworkMapper = slNetworkMapper;
    }

    @PostConstruct
    public void init() {
        allSLData = slNetworkMapper.selectAll();
        log.info("SL数据初始化完成，共{}条记录", allSLData.size());
    }

    /**
     * 实现通用接口：获取网络拓扑数据（返回Map，包含节点、连接等）
     */
    @Override
    public Map<String, Object> getNetworkData(Map<String, Object> params) {
        // 解析阈值参数（默认-1）
        double threshold = params.get("sensitivityThreshold") != null
                ? Double.parseDouble(params.get("sensitivityThreshold").toString())
                : -1;

        // 1. 筛选符合条件的数据
        List<SLResponse> filteredData = allSLData.stream()
                .filter(record -> {
                    try {
                        double sensitivity = Double.parseDouble(record.getGeminiSensitive());
                        return sensitivity < threshold;
                    } catch (NumberFormatException e) {
                        return false;
                    }
                })
                .collect(Collectors.toList());

        // 2. 提取节点
        Set<String> geneSet = filteredData.stream()
                .flatMap(record -> Stream.of(record.getGeneA(), record.getGeneB()))
                .collect(Collectors.toSet());
        List<SLResponse.Node> nodes = geneSet.stream()
                .map(gene -> {
                    SLResponse.Node node = new SLResponse.Node();
                    node.setName(gene);
                    return node;
                })
                .collect(Collectors.toList());

        // 3. 构建连接
        List<SLResponse.Link> links = filteredData.stream()
                .map(record -> {
                    SLResponse.Link link = new SLResponse.Link();
                    link.setSource(record.getGeneA());
                    link.setTarget(record.getGeneB());
                    link.setSensitivity(Double.parseDouble(record.getGeminiSensitive()));
                    link.setCellLine(record.getCellLine());
                    return link;
                })
                .collect(Collectors.toList());

        // 4. 计算统计值
        List<Double> allSensitivities = allSLData.stream()
                .map(record -> Double.parseDouble(record.getGeminiSensitive()))
                .collect(Collectors.toList());
        double max = allSensitivities.stream().mapToDouble(Double::doubleValue).max().orElse(0.0);
        double min = allSensitivities.stream().mapToDouble(Double::doubleValue).min().orElse(0.0);
        double avg = allSensitivities.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);

        // 5. 组装原始记录
        List<Map<String, String>> originalRecords = filteredData.stream()
                .map(record -> {
                    Map<String, String> map = new HashMap<>();
                    map.put("geneA", record.getGeneA());
                    map.put("geneB", record.getGeneB());
                    map.put("geminiSensitive", record.getGeminiSensitive());
                    map.put("cellLine", record.getCellLine());
                    return map;
                })
                .collect(Collectors.toList());

        // 6. 封装为统一Map格式返回（兼容Java 8，用HashMap初始化sensitivityStats）
        Map<String, Object> result = new HashMap<>();
        result.put("nodes", nodes);
        result.put("links", links);
        result.put("maxSensitivity", max);
        result.put("minSensitivity", min);

        // Java 8兼容：手动创建HashMap并添加统计数据
        Map<String, Double> sensitivityStats = new HashMap<>();
        sensitivityStats.put("max", max);
        sensitivityStats.put("min", min);
        sensitivityStats.put("average", avg);
        result.put("sensitivityStats", sensitivityStats);

        result.put("originalRecords", originalRecords);
        result.put("qualifiedGeneNames", geneSet);

        return result;
    }

    /**
     * 实现通用接口：获取表格分页数据
     */
    @Override
    public Map<String, Object> getTableData(Map<String, Object> params) {
        // 解析参数
        double threshold = Double.parseDouble(params.getOrDefault("threshold", "-1.0").toString());
        int pageNum = Integer.parseInt(params.getOrDefault("pageNum", "1").toString());
        int pageSize = Integer.parseInt(params.getOrDefault("pageSize", "10").toString());
        String gene = params.get("gene") != null ? params.get("gene").toString().trim() : "";

        // 参数校验
        pageNum = Math.max(pageNum, 1);
        pageSize = Math.max(Math.min(pageSize, 100), 10);
        int offset = (pageNum - 1) * pageSize;

        // 1. 阈值筛选
        List<SLResponse> thresholdData = allSLData.stream()
                .filter(record -> {
                    try {
                        return Double.parseDouble(record.getGeminiSensitive()) < threshold;
                    } catch (NumberFormatException e) {
                        return false;
                    }
                })
                .collect(Collectors.toList());

        // 2. 基因筛选
        List<SLResponse> filteredData = thresholdData.stream()
                .filter(record -> {
                    if (gene.isEmpty()) return true;
                    return gene.equals(record.getGeneA()) || gene.equals(record.getGeneB());
                })
                .collect(Collectors.toList());

        // 3. 分页处理
        int total = filteredData.size();
        int start = Math.min(offset, total);
        int end = Math.min(start + pageSize, total);
        List<SLResponse> pageData = start < end ? filteredData.subList(start, end) : new ArrayList<>();

        // 4. 转换表格格式
        List<Map<String, String>> tableList = pageData.stream()
                .map(record -> {
                    Map<String, String> row = new HashMap<>();
                    row.put("geneA", record.getGeneA());
                    row.put("geneB", record.getGeneB());
                    row.put("geminiSensitive", record.getGeminiSensitive());
                    row.put("cellLine", record.getCellLine());
                    return row;
                })
                .collect(Collectors.toList());

        // 5. 封装分页结果（统一格式）
        Map<String, Object> result = new HashMap<>();
        result.put("total", total);
        result.put("list", tableList);
        result.put("pageNum", pageNum);
        result.put("pageSize", pageSize);

        return result;
    }

    @Override
    public byte[] exportAllData() {
        List<SLResponse> slResponses = slNetworkMapper.selectAll();
        return buildTsv(slResponses);
    }

    @Override
    public byte[] exportCurrentData(String name, Double s) {
        // 1. 校验基因名称必填
        if (name == null || name.trim().isEmpty()) {
            log.warn("导出当前数据失败：基因名称为空");
            return "基因名称不能为空".getBytes(StandardCharsets.UTF_8);
        }
        String targetGene = name.trim();

        // 2. 处理阈值参数：优先使用传入的s，无值则默认0.5（与getTableData保持一致）
        double threshold = (s != null) ? s : 0.5;

        // 3. 筛选逻辑：包含目标基因且符合阈值条件的记录
        List<SLResponse> currentData = allSLData.stream()
                .filter(record -> {
                    try {
                        // 解析敏感性值并与阈值比较
                        double sensitivity = Double.parseDouble(record.getGeminiSensitive());
                        boolean thresholdPass = sensitivity < threshold;

                        // 检查是否包含目标基因（geneA或geneB）
                        boolean geneMatch = targetGene.equals(record.getGeneA())
                                || targetGene.equals(record.getGeneB());

                        return thresholdPass && geneMatch;
                    } catch (NumberFormatException e) {
                        log.error("解析敏感性值失败（记录：{}）", record, e);
                        return false; // 跳过解析失败的记录
                    }
                })
                .collect(Collectors.toList());

        // 4. 处理无数据场景
        if (currentData.isEmpty()) {
            String message = String.format("未找到基因【%s】在阈值【%s】下的相关数据", targetGene, threshold);
            log.info(message);
            return message.getBytes(StandardCharsets.UTF_8);
        }

        // 5. 构建TSV并返回
        return buildTsv(currentData);
    }

    /**
     * 构建TSV格式数据
     */
    private byte[] buildTsv(List<SLResponse> data) {
        StringBuilder tsv = new StringBuilder();
        tsv.append("GeneA\tGeneB\tGEMINI sensitive\tCell line\n");
        data.forEach(record -> {
            tsv.append(record.getGeneA()).append("\t")
                    .append(record.getGeneB()).append("\t")
                    .append(record.getGeminiSensitive()).append("\t")
                    .append(record.getCellLine()).append("\n");
        });
        return tsv.toString().getBytes(StandardCharsets.UTF_8);
    }
}