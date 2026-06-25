package com.bishe.ddr_springboot.service.network;

import com.bishe.ddr_springboot.entity.TFRegulation;
import com.bishe.ddr_springboot.mapper.network.TFRegulationMapper;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.util.StringUtils;

import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

@Service
@Slf4j
public class TFRegulatoryServiceImpl implements TFRegulatoryService {

    private final TFRegulationMapper tfRegulationMapper;

    // 构造方法注入 Mapper（替代 ExcelReaderUtil）
    public TFRegulatoryServiceImpl(TFRegulationMapper tfRegulationMapper) {
        this.tfRegulationMapper = tfRegulationMapper;
        log.info("TFRegulatoryServiceImpl initialized with database mapper");
    }

    /**
     * 根据基因查询相关调控关系（作为TF或靶基因）— 从数据库查询
     */
    @Override
    public List<TFRegulation> getRegulationsByGene(String gene) {
        if (!StringUtils.hasText(gene)) {
            return new ArrayList<>();
        }
        return tfRegulationMapper.selectByGene(gene.trim());
    }

    /**
     * 分页查询（内部方法，复用）
     */
    private Map<String, Object> getRegulationsByPage(String gene, int pageNum, int pageSize) {
        pageNum = Math.max(pageNum, 1);
        pageSize = Math.max(Math.min(pageSize, 100), 10);

        List<TFRegulation> allRegs = getRegulationsByGene(gene);
        int total = allRegs.size();

        if (total == 0) {
            int finalPageNum1 = pageNum;
            int finalPageSize1 = pageSize;
            return new HashMap<String, Object>() {{
                put("total", 0);
                put("list", new ArrayList<>());
                put("pageNum", finalPageNum1);
                put("pageSize", finalPageSize1);
            }};
        }

        int fromIndex = (pageNum - 1) * pageSize;
        int toIndex = Math.min(fromIndex + pageSize, total);
        List<TFRegulation> pageData = allRegs.subList(fromIndex, toIndex);

        int finalPageNum = pageNum;
        int finalPageSize = pageSize;
        return new HashMap<String, Object>() {{
            put("total", total);
            put("list", pageData);
            put("pageNum", finalPageNum);
            put("pageSize", finalPageSize);
        }};
    }

    // === 添加以下两个方法到 TFRegulatoryServiceImpl 类中 ===

    @Override
    public List<TFRegulation> getRegulationsByTF(String tfName) {
        if (!StringUtils.hasText(tfName)) {
            return new ArrayList<>();
        }
        return tfRegulationMapper.selectByTF(tfName.trim());
    }

    @Override
    public List<TFRegulation> getRegulatorsOfTarget(String targetGene) {
        if (!StringUtils.hasText(targetGene)) {
            return new ArrayList<>();
        }
        return tfRegulationMapper.selectByTargetGene(targetGene.trim());
    }

    // === 新导出方法（调用 Mapper）===

    @Override
    public byte[] exportAllData() {
        // 调用 Mapper 全表导出
        List<Map<String, Object>> allData = tfRegulationMapper.selectAllForExport();
        return buildTsvFromMap(allData);
    }

    @Override
    public byte[] exportCurrentData(String name, Double s) {
        if (!StringUtils.hasText(name)) {
            return "Gene name is required".getBytes(StandardCharsets.UTF_8);
        }
        // 先查 TFRegulation 列表，再转为 Map（复用 buildTsv 逻辑）
        List<TFRegulation> regulations = tfRegulationMapper.selectByGene(name.trim());
        return buildTsv(regulations);
    }

    // === 网络和表格方法（保持不变，但数据来自数据库）===

    @Override
    public Map<String, Object> getNetworkData(Map<String, Object> params) {
        String gene = params.getOrDefault("gene", "").toString().trim();
        List<TFRegulation> regulations = getRegulationsByGene(gene);

        Set<String> nodeNames = new HashSet<>();
        for (TFRegulation reg : regulations) {
            if (reg != null) {
                nodeNames.add(reg.getSource());
                nodeNames.add(reg.getTarget());
            }
        }

        List<Map<String, String>> nodes = nodeNames.stream()
                .map(name -> {
                    Map<String, String> node = new HashMap<>();
                    node.put("name", name);
                    boolean isTF = regulations.stream()
                            .anyMatch(reg -> reg != null && name.equals(reg.getSource()));
                    node.put("type", isTF ? "TF" : "Target");
                    return node;
                })
                .collect(Collectors.toList());

        List<Map<String, String>> links = regulations.stream()
                .filter(Objects::nonNull)
                .map(reg -> {
                    Map<String, String> link = new HashMap<>();
                    link.put("source", reg.getSource());
                    link.put("target", reg.getTarget());
                    link.put("regulationType", reg.getRegulationType());
                    return link;
                })
                .collect(Collectors.toList());

        Map<String, Object> result = new HashMap<>();
        result.put("nodes", nodes);
        result.put("links", links);
        result.put("total", regulations.size());
        return result;
    }

    @Override
    public Map<String, Object> getTableData(Map<String, Object> params) {
        String gene = params.getOrDefault("gene", "").toString().trim();
        int pageNum = parseToInt(params.get("pageNum"), 1);
        int pageSize = parseToInt(params.get("pageSize"), 10);
        return getRegulationsByPage(gene, pageNum, pageSize);
    }

    // === TSV 构建方法 ===

    /**
     * 从 TFRegulation 列表构建 TSV
     */
    private byte[] buildTsv(List<TFRegulation> data) {
        StringBuilder tsv = new StringBuilder();
        tsv.append("TF\tTarget\tMode of Regulation\tReferences (PMID)\n");
        for (TFRegulation reg : data) {
            if (reg == null) continue;
            tsv.append(safeTrim(reg.getSource())).append("\t")
                    .append(safeTrim(reg.getTarget())).append("\t")
                    .append(safeTrim(reg.getRegulationType())).append("\t")
                    .append(safeTrim(reg.getEvidence())).append("\n");
        }
        return tsv.toString().getBytes(StandardCharsets.UTF_8);
    }

    /**
     * 从 Map 列表构建 TSV（用于全表导出）
     */
    private byte[] buildTsvFromMap(List<Map<String, Object>> data) {
        StringBuilder tsv = new StringBuilder();
        tsv.append("TF\tTarget\tMode of Regulation\tReferences (PMID)\n");
        for (Map<String, Object> row : data) {
            tsv.append(safeTrim(row.get("TF"))).append("\t")
                    .append(safeTrim(row.get("Target"))).append("\t")
                    .append(safeTrim(row.get("Mode of Regulation"))).append("\t")
                    .append(safeTrim(row.get("References (PMID)"))).append("\n");
        }
        return tsv.toString().getBytes(StandardCharsets.UTF_8);
    }

    private String safeTrim(Object value) {
        return value != null ? value.toString().trim() : "";
    }

    private int parseToInt(Object param, int defaultValue) {
        if (param == null) return defaultValue;
        try {
            return Integer.parseInt(param.toString().trim());
        } catch (NumberFormatException e) {
            log.warn("参数转换失败，使用默认值: {}", defaultValue);
            return defaultValue;
        }
    }
}