package com.bishe.ddr_springboot.service;

import com.bishe.ddr_springboot.entity.DrugTarget;
import com.bishe.ddr_springboot.entity.PageResult;
import com.bishe.ddr_springboot.mapper.DrugTargetMapper;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.util.StringUtils;

import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

@Service
@Slf4j
public class DrugTargetServiceImpl implements DrugTargetService {

    private final DrugTargetMapper drugTargetMapper;

    public DrugTargetServiceImpl(DrugTargetMapper drugTargetMapper) {
        this.drugTargetMapper = drugTargetMapper;
    }

    @Override
    public List<DrugTarget> getAllDrugTargets() {
        log.info("查询所有药物靶点关系");
        return drugTargetMapper.selectAll();
    }

    @Override
    public PageResult<DrugTarget> getDrugTargetsByPage(Integer pageNum, Integer pageSize) {
        pageNum = Math.max(pageNum == null ? 1 : pageNum, 1);
        pageSize = Math.max(Math.min(pageSize == null ? 10 : pageSize, 100), 10);
        int startIndex = (pageNum - 1) * pageSize;

        log.info("分页查询药物靶点关系：页码={}, 每页条数={}", pageNum, pageSize);
        List<DrugTarget> list = drugTargetMapper.selectByPage(startIndex, pageSize);
        int total = drugTargetMapper.selectTotalCount();

        return new PageResult<>(list, total, pageNum, pageSize);
    }

    @Override
    public List<DrugTarget> getByDrugName(String drugName) {
        if (!StringUtils.hasText(drugName)) {
            log.warn("查询药物靶点时，药物名称为空");
            return Collections.emptyList();
        }
        String trimName = drugName.trim();
        log.info("根据药物名称查询靶点：{}", trimName);
        List<DrugTarget> result = drugTargetMapper.selectByDrugName(trimName);
        log.info("药物[{}]查询到{}条靶点记录", trimName, result.size());
        return result;
    }

    @Override
    public List<DrugTarget> getByGeneName(String geneName) {
        if (!StringUtils.hasText(geneName)) {
            log.warn("查询药物靶点时，基因名称为空");
            return Collections.emptyList();
        }
        String trimName = geneName.trim();
        log.info("根据基因名称查询靶点：{}", trimName);
        List<DrugTarget> result = drugTargetMapper.selectByGeneName(trimName);
        log.info("基因[{}]查询到{}条靶点记录", trimName, result.size());
        return result;
    }

    @Override
    public PageResult<DrugTarget> queryByConditions(Map<String, Object> params) {
        // 解析分页参数
        int pageNum = params.get("pageNum") != null ? (Integer) params.get("pageNum") : 1;
        int pageSize = params.get("pageSize") != null ? (Integer) params.get("pageSize") : 10;
        pageNum = Math.max(pageNum, 1);
        pageSize = Math.max(Math.min(pageSize, 100), 10);
        int startIndex = (pageNum - 1) * pageSize;

        // 构建查询条件
        String drugName = params.get("drugName") != null ? params.get("drugName").toString().trim() : "";
        String geneName = params.get("geneName") != null ? params.get("geneName").toString().trim() : "";
        String approved = params.get("approved") != null ? params.get("approved").toString().trim() : "";

        log.info("多条件查询药物靶点：药物={}, 基因={}, 审批状态={}, 页码={}, 每页条数={}",
                drugName, geneName, approved, pageNum, pageSize);

        // 执行查询
        List<DrugTarget> list = drugTargetMapper.selectByConditions(
                drugName, geneName, approved, startIndex, pageSize);
        int total = drugTargetMapper.selectCountByConditions(drugName, geneName, approved);

        return new PageResult<>(list, total, pageNum, pageSize);
    }

    /**
     * 新增：实现药物-靶点网络数据构建
     */
    @Override
    public Map<String, Object> getNetworkData(Map<String, Object> params) {
        // 1. 处理查询参数（支持无参数场景）
        String query = params != null ? params.getOrDefault("name", "").toString().trim() : "";
        boolean isFullData = !StringUtils.hasText(query); // 无参数或参数为空时查询全量数据

        List<DrugTarget> targetList;
        if (isFullData) {
            // 无参数：查询全部药物-靶点关系
            targetList = drugTargetMapper.selectAll(); // 需要在DrugTargetMapper中新增全量查询方法
            log.info("查询全量药物靶点网络数据，共{}条记录", targetList.size());
        } else {
            // 有参数：按原有逻辑查询（药物或基因匹配）
            List<DrugTarget> byDrug = drugTargetMapper.selectByDrugName(query);
            List<DrugTarget> byGene = drugTargetMapper.selectByGeneName(query);
            Set<DrugTarget> combined = new HashSet<>(byDrug);
            combined.addAll(byGene);
            targetList = new ArrayList<>(combined);
            log.info("查询[{}]相关药物靶点关系，共{}条记录", query, targetList.size());
        }

        // 2. 处理无结果场景
        if (targetList.isEmpty()) {
            log.info("未查询到{}药物靶点关系", isFullData ? "全量" : "与[" + query + "]相关的");
            return Collections.emptyMap();
        }

        // 3. 提取节点（药物和基因）
        Set<String> drugNodes = targetList.stream()
                .map(DrugTarget::getDrugName)
                .filter(StringUtils::hasText)
                .map(String::trim)
                .collect(Collectors.toSet());
        Set<String> geneNodes = targetList.stream()
                .map(DrugTarget::getGeneName)
                .filter(StringUtils::hasText)
                .map(String::trim)
                .collect(Collectors.toSet());

        // 4. 构建节点列表（区分类型和样式）
        List<Map<String, Object>> nodes = new ArrayList<>();
        // 添加药物节点
        nodes.addAll(drugNodes.stream()
                .map(drug -> {
                    Map<String, Object> node = new HashMap<>();
                    node.put("name", drug);
                    node.put("type", "drug");
                    node.put("symbol", "diamond");
                    node.put("symbolSize", 20);
                    node.put("itemStyle", Collections.singletonMap("color", "#1E90FF"));
                    return node;
                })
                .collect(Collectors.toList()));
        // 添加基因节点
        nodes.addAll(geneNodes.stream()
                .map(gene -> {
                    Map<String, Object> node = new HashMap<>();
                    node.put("name", gene);
                    node.put("type", "gene");
                    node.put("symbol", "circle");
                    node.put("symbolSize", 18);
                    node.put("itemStyle", Collections.singletonMap("color", "#32CD32"));
                    return node;
                })
                .collect(Collectors.toList()));

        // 5. 构建连接（药物-基因）
        List<Map<String, Object>> links = targetList.stream()
                .filter(t -> StringUtils.hasText(t.getDrugName()) && StringUtils.hasText(t.getGeneName()))
                .map(t -> {
                    Map<String, Object> link = new HashMap<>();
                    link.put("source", t.getDrugName().trim());
                    link.put("target", t.getGeneName().trim());
                    link.put("action", t.getAction() != null ? t.getAction().trim() : "unknown");

                    Map<String, Object> lineStyle = new HashMap<>();
                    lineStyle.put("width", 2);
                    lineStyle.put("curveness", 0.1);
                    link.put("lineStyle", lineStyle);
                    return link;
                })
                .collect(Collectors.toList());

        // 6. 封装结果（增加全量数据标记）
        Map<String, Object> result = new HashMap<>();
        result.put("nodes", nodes);
        result.put("links", links);
        result.put("total", links.size());
        result.put("query", isFullData ? "all" : query); // 全量查询时标记为"all"
        result.put("isFullData", isFullData); // 明确标记是否为全量数据
        log.info("药物靶点网络构建完成：节点数={}, 连接数={}", nodes.size(), links.size());
        return result;
    }

    public Map<String, Object> getTableData(Map<String, Object> params) {
        String query = params.getOrDefault("name", "").toString().trim();
        int pageNum = parseToInt(params.get("pageNum"), 1);
        int pageSize = parseToInt(params.get("pageSize"), 10);
        pageNum = Math.max(pageNum, 1);
        pageSize = Math.max(Math.min(pageSize, 100), 10);
        int startIndex = (pageNum - 1) * pageSize;

        List<DrugTarget> list;
        int total;

        // 优先按药物查询，再按基因查询
        int countByDrug = drugTargetMapper.selectCountByDrugName(query);
        int countByGene = drugTargetMapper.selectCountByGeneName(query);

        if (countByDrug > 0) {
            list = drugTargetMapper.selectByDrugNamePage(query, startIndex, pageSize);
            total = countByDrug;
        } else if (countByGene > 0) {
            list = drugTargetMapper.selectByGeneNamePage(query, startIndex, pageSize);
            total = countByGene;
        } else {
            log.info("无符合条件的药物靶点记录：{}", query);
            return buildEmptyPageResult(pageNum, pageSize);
        }

        // 转换为表格数据
        List<Map<String, String>> tableList = list.stream()
                .map(target -> {
                    Map<String, String> row = new HashMap<>();
                    row.put("drugbankId", trimValue(target.getDrugbankId()));
                    row.put("drugName", trimValue(target.getDrugName()));
                    row.put("drugType", trimValue(target.getDrugType()));
                    row.put("approved", trimValue(target.getApproved()));
                    row.put("targetName", trimValue(target.getTargetName()));
                    row.put("organism", trimValue(target.getOrganism()));
                    row.put("action", trimValue(target.getAction()));
                    row.put("geneName", trimValue(target.getGeneName()));
                    row.put("uniprotId", trimValue(target.getUniprotId()));
                    row.put("polypeptideName", trimValue(target.getPolypeptideName()));
                    row.put("polypeptideId", trimValue(target.getPolypeptideId()));
                    return row;
                })
                .collect(Collectors.toList());

        Map<String, Object> pageResult = new HashMap<>();
        pageResult.put("total", total);
        pageResult.put("list", tableList);
        pageResult.put("pageNum", pageNum);
        pageResult.put("pageSize", pageSize);
        return pageResult;
    }

    @Override
    public byte[] exportAll() {
        log.info("导出所有药物靶点数据");
        List<Map<String, Object>> allData = drugTargetMapper.selectAllForExport();
        return buildTsv(allData);
    }

    @Override
    public byte[] exportByConditions(Map<String, Object> params) {
        String drugName = params.get("drugName") != null ? params.get("drugName").toString().trim() : "";
        String geneName = params.get("geneName") != null ? params.get("geneName").toString().trim() : "";
        String approved = params.get("approved") != null ? params.get("approved").toString().trim() : "";

        log.info("按条件导出药物靶点数据：药物={}, 基因={}, 审批状态={}", drugName, geneName, approved);
        List<Map<String, Object>> data = drugTargetMapper.selectByConditionsForExport(drugName, geneName, approved);
        return buildTsv(data);
    }

    public byte[] exportCurrentData(String name, Double s) {
        if (!StringUtils.hasText(name)) {
            log.warn("导出当前药物靶点数据失败：名称为空");
            return "名称不能为空".getBytes(StandardCharsets.UTF_8);
        }
        String trimName = name.trim();
        log.info("导出与[{}]相关的药物靶点数据", trimName);
        List<Map<String, Object>> data = drugTargetMapper.selectForExport(trimName);
        return buildTsv(data);
    }

    /**
     * 构建TSV格式数据
     */
    private byte[] buildTsv(List<Map<String, Object>> data) {
        if (data.isEmpty()) {
            return "No data to export".getBytes(StandardCharsets.UTF_8);
        }

        String[] headers = {
                "drugbank_id", "drug_name", "drug_type", "approved",
                "target_name", "organism", "action", "gene_name",
                "uniprot_id", "polypeptide_name", "polypeptide_id"
        };

        StringBuilder tsv = new StringBuilder(String.join("\t", headers) + "\n");

        for (Map<String, Object> row : data) {
            List<String> values = Arrays.stream(headers)
                    .map(header -> {
                        Object value = row.get(header);
                        return value != null ? value.toString().trim() : "";
                    })
                    .collect(Collectors.toList());
            tsv.append(String.join("\t", values)).append("\n");
        }

        return tsv.toString().getBytes(StandardCharsets.UTF_8);
    }

    private String trimValue(String value) {
        return StringUtils.hasText(value) ? value.trim() : "";
    }

    private int parseToInt(Object param, int defaultValue) {
        if (param == null) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(param.toString());
        } catch (NumberFormatException e) {
            log.warn("参数解析失败，使用默认值：{}", defaultValue, e);
            return defaultValue;
        }
    }

    private Map<String, Object> buildEmptyPageResult(int pageNum, int pageSize) {
        Map<String, Object> result = new HashMap<>();
        result.put("total", 0);
        result.put("list", Collections.emptyList());
        result.put("pageNum", pageNum);
        result.put("pageSize", pageSize);
        return result;
    }
    @Override
    public List<String> getDrugNameSuggestions(String keyword) {
        if (!StringUtils.hasText(keyword)) {
            log.warn("药物联想搜索关键词为空");
            return Collections.emptyList();
        }
        String trimKeyword = keyword.trim();
        log.info("获取药物名称联想建议：关键词={}", trimKeyword);
        List<String> suggestions = drugTargetMapper.selectDrugNameSuggestions(trimKeyword);
        log.info("药物联想建议查询结果：{}条", suggestions.size());
        return suggestions;
    }

    @Override
    public List<String> getGeneNameSuggestions(String keyword) {
        if (!StringUtils.hasText(keyword)) {
            log.warn("基因联想搜索关键词为空");
            return Collections.emptyList();
        }
        String trimKeyword = keyword.trim();
        log.info("获取基因名称联想建议：关键词={}", trimKeyword);
        List<String> suggestions = drugTargetMapper.selectGeneNameSuggestions(trimKeyword);
        log.info("基因联想建议查询结果：{}条", suggestions.size());
        return suggestions;
    }
}