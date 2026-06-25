package com.bishe.ddr_springboot.service.network;

import com.bishe.ddr_springboot.entity.NcRNATarget;
import com.bishe.ddr_springboot.mapper.network.NcRNATargetMapper;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.util.StringUtils;

import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

@Service
@Slf4j
public class NcRNATargetServiceImpl implements NcRNATargetService {

    private final NcRNATargetMapper ncRNATargetMapper;

    public NcRNATargetServiceImpl(NcRNATargetMapper ncRNATargetMapper) {
        this.ncRNATargetMapper = ncRNATargetMapper;
    }

    @Override
    public Map<String, Object> getNetworkData(Map<String, Object> params) {
        String query = params.getOrDefault("name", "").toString().trim();
        if (!StringUtils.hasText(query)) {
            log.warn("查询名称为空，返回空网络数据");
            return Collections.emptyMap();
        }

        // 分别查询 ncRNA 和 靶基因
        List<NcRNATarget> byNcRNA = ncRNATargetMapper.selectByNcRNAName(query);
        List<NcRNATarget> byGene = ncRNATargetMapper.selectByTargetGene(query);

        List<NcRNATarget> finalList;
        String searchMode;

        if (!byNcRNA.isEmpty() && byGene.isEmpty()) {
            finalList = byNcRNA;
            searchMode = "ncRNA";
        } else if (byNcRNA.isEmpty() && !byGene.isEmpty()) {
            finalList = byGene;
            searchMode = "gene";
        } else if (!byNcRNA.isEmpty() && !byGene.isEmpty()) {
            // 合并（去重可选，此处简单合并）
            finalList = new ArrayList<>(byNcRNA);
            finalList.addAll(byGene);
            searchMode = "both";
        } else {
            log.info("未查询到与 {} 相关的调控关系", query);
            return Collections.emptyMap();
        }

        // 提取所有节点名称
        Set<String> ncRNANames = new HashSet<>();
        Set<String> targetGeneNames = new HashSet<>();
        for (NcRNATarget t : finalList) {
            if (StringUtils.hasText(t.getNcRNA())) {
                ncRNANames.add(t.getNcRNA().trim());
            }
            if (StringUtils.hasText(t.getTargetGene())) {
                targetGeneNames.add(t.getTargetGene().trim());
            }
        }

        Set<String> nodeNames = new HashSet<>();
        nodeNames.add(query); // 查询词本身
        nodeNames.addAll(ncRNANames);
        nodeNames.addAll(targetGeneNames);

        // 构建节点
        List<Map<String, Object>> nodes = new ArrayList<>();
        for (String name : nodeNames) {
            Map<String, Object> node = new HashMap<>();
            node.put("name", name);

            if (name.equalsIgnoreCase(query)) {
                if ("ncRNA".equals(searchMode) || ("both".equals(searchMode) && ncRNANames.contains(name))) {
                    node.put("type", "query_ncRNA");
                    node.put("symbol", "triangle");
                    node.put("symbolSize", 22);
                    node.put("itemStyle", Collections.singletonMap("color", "#dc3545"));
                } else {
                    node.put("type", "query_gene");
                    node.put("symbol", "circle");
                    node.put("symbolSize", 22);
                    node.put("itemStyle", Collections.singletonMap("color", "#28a745"));
                }
            } else if (ncRNANames.contains(name)) {
                node.put("type", "ncRNA");
                node.put("symbol", "triangle");
                node.put("symbolSize", 18);
                node.put("itemStyle", Collections.singletonMap("color", "#dc3545"));
            } else {
                node.put("type", "target_gene");
                node.put("symbol", "circle");
                node.put("symbolSize", 18);
                node.put("itemStyle", Collections.singletonMap("color", "#28a745"));
            }
            nodes.add(node);
        }

        // 构建边
        List<Map<String, Object>> links = new ArrayList<>();
        for (NcRNATarget t : finalList) {
            if (StringUtils.hasText(t.getNcRNA()) && StringUtils.hasText(t.getTargetGene())) {
                Map<String, Object> link = new HashMap<>();
                link.put("source", t.getNcRNA().trim());
                link.put("target", t.getTargetGene().trim());

                Map<String, Object> lineStyle = new HashMap<>();
                lineStyle.put("width", 2.5);
                lineStyle.put("curveness", 0.2);
                lineStyle.put("color", "#6c757d");
                link.put("lineStyle", lineStyle);

                links.add(link);
            }
        }

        Map<String, Object> result = new HashMap<>();
        result.put("nodes", nodes);
        result.put("links", links);
        result.put("searchMode", searchMode);
        result.put("query", query);
        return result;
    }

    @Override
    public Map<String, Object> getTableData(Map<String, Object> params) {
        String query = params.getOrDefault("name", "").toString().trim();
        int pageNum = parseToInt(params.get("pageNum"), 1);
        int pageSize = parseToInt(params.get("pageSize"), 10);
        pageNum = Math.max(pageNum, 1);
        pageSize = Math.max(Math.min(pageSize, 100), 10);
        int startIndex = (pageNum - 1) * pageSize;

        List<NcRNATarget> list;
        int total;

        int countNcRNA = ncRNATargetMapper.selectCountByNcRNAname(query);
        int countGene = ncRNATargetMapper.selectCountByTargetGene(query);

        if (countNcRNA > 0 && countGene == 0) {
            list = ncRNATargetMapper.selectByNcRNAnamePage(query, startIndex, pageSize);
            total = countNcRNA;
        } else if (countNcRNA == 0 && countGene > 0) {
            list = ncRNATargetMapper.selectByTargetGenePage(query, startIndex, pageSize);
            total = countGene;
        } else if (countNcRNA > 0 && countGene > 0) {
            // 简单处理：优先返回 ncRNA 的结果（或可合并，但分页复杂）
            list = ncRNATargetMapper.selectByNcRNAnamePage(query, startIndex, pageSize);
            total = countNcRNA;
        } else {
            int finalPageNum = pageNum;
            int finalPageSize = pageSize;
            return new HashMap<String, Object>() {{
                put("total", 0);
                put("list", new ArrayList<Map<String, String>>());
                put("pageNum", finalPageNum);
                put("pageSize", finalPageSize);
            }};
        }

        List<Map<String, String>> tableList = new ArrayList<>();
        for (NcRNATarget t : list) {
            Map<String, String> row = new HashMap<>();
            row.put("mirTarBaseId", trimValue(t.getMirTarBaseId()));
            row.put("ncRNA", trimValue(t.getNcRNA()));
            row.put("ncRNASpecies", trimValue(t.getNcRNASpecies()));
            row.put("targetGene", trimValue(t.getTargetGene()));
            row.put("targetGeneEntrezId", trimValue(t.getTargetGeneEntrezId()));
            row.put("targetGeneSpecies", trimValue(t.getTargetGeneSpecies()));
            row.put("experiments", trimValue(t.getExperiments()));
            row.put("supportType", trimValue(t.getSupportType()));
            row.put("reference", trimValue(t.getReference()));
            tableList.add(row);
        }

        Map<String, Object> pageResult = new HashMap<>();
        pageResult.put("total", total);
        pageResult.put("list", tableList);
        pageResult.put("pageNum", pageNum);
        pageResult.put("pageSize", pageSize);
        return pageResult;
    }

    @Override
    public byte[] exportAllData() {
        // 导出全表
        List<Map<String, Object>> all = ncRNATargetMapper.selectAllForExport();
        return buildTsv(all);
    }

    @Override
    public byte[] exportCurrentData(String name,Double s) {
        // 导出当前查询结果
        List<Map<String, Object>> current = ncRNATargetMapper.selectForExport(name);
        return buildTsv(current);
    }

    private byte[] buildTsv(List<Map<String, Object>> data) {
        if (data.isEmpty()) {
            return "No data to export".getBytes(StandardCharsets.UTF_8);
        }

        String[] headers = {
                "miRTarBase ID", "ncRNA", "Species (ncRNA)",
                "Target Gene", "Target Gene (Entrez ID)", "Species (Target Gene)",
                "Experiments", "Support Type", "References (PMID)"
        };

        StringBuilder tsv = new StringBuilder();
        tsv.append(String.join("\t", headers)).append("\n");

        for (Map<String, Object> row : data) {
            List<String> values = new ArrayList<>();
            for (String header : headers) {
                String field = mapHeaderToField(header);
                Object val = row.get(field);
                values.add(trimValue(val));
            }
            tsv.append(String.join("\t", values)).append("\n");
        }

        return tsv.toString().getBytes(StandardCharsets.UTF_8);
    }

    private String mapHeaderToField(String header) {
        switch (header) {
            case "miRTarBase ID": return "miRTarBase ID"; // 对应 Map 中的 "miRTarBase ID"
            case "ncRNA": return "miRNA"; // 对应 Map 中的 "miRNA"（ncRNA 的实际存储键）
            case "Species (ncRNA)": return "Species (miRNA)"; // 对应 Map 中的 "Species (miRNA)"
            case "Target Gene": return "Target Gene"; // 对应 Map 中的 "Target Gene"
            case "Target Gene (Entrez ID)": return "Target Gene (Entrez ID)"; // 对应实际键
            case "Species (Target Gene)": return "Species (Target Gene)"; // 对应实际键
            case "Experiments": return "Experiments"; // 对应实际键
            case "Support Type": return "Support Type"; // 对应实际键
            case "References (PMID)": return "References (PMID)"; // 注意：原数据中可能是 "References (PMID)"，需确认
            default: return "";
        }
    }

    private String trimValue(Object value) {
        return value != null ? value.toString().trim() : "";
    }

    private int parseToInt(Object param, int defaultValue) {
        if (param == null) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(param.toString().trim());
        } catch (NumberFormatException e) {
            log.warn("参数转换为int失败，使用默认值: {}", defaultValue, e);
            return defaultValue;
        }
    }
}