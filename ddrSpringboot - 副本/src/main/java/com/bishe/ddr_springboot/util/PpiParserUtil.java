package com.bishe.ddr_springboot.util;

import com.bishe.ddr_springboot.entity.PpiResponse;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.util.StringUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.stream.Collectors;

/**
 * PPI文件解析工具类（读取 resources/static/network/ppi/string_interactions.tsv）
 */
public class PpiParserUtil {
    private static final Logger log = LoggerFactory.getLogger(PpiParserUtil.class);

    // 定义 resources 下的 PPI 文件路径（固定为 string_interactions.tsv）
    private static final String PPI_RESOURCE_PATH = "static/network/ppi/string_interactions.tsv";
    // 支持的分隔符（TSV文件默认制表符，兼容逗号以防格式错误）
    private static final String[] SEPARATORS = {"\t", ","};


    // 对外提供文件路径的方法（供导出接口使用）
    public static String getPpiResourcePath() {
        return PPI_RESOURCE_PATH;
    }

    /**
     * 解析并筛选PPI资源文件（根据阈值筛选Node1出现次数）
     * @param degreeThreshold 筛选阈值（Node1出现次数≥此值才保留）
     * @return 筛选后的PPI响应数据
     */
    public static PpiResponse parseAndFilterPpiResourceFile(Integer degreeThreshold) {
        // 1. 初始化容器（存储原始数据和统计信息）
        Set<String> allNodeNameSet = new HashSet<>();
        List<OriginalRecord> originalRecords = new ArrayList<>(); // 存储原始记录
        Map<String, Integer> node1CountMap = new HashMap<>();

        // 2. 读取并解析文件
        Resource ppiResource = new ClassPathResource(PPI_RESOURCE_PATH);
        if (!ppiResource.exists()) {
            log.error("PPI文件不存在于 resources 路径：{}", PPI_RESOURCE_PATH);
            throw new RuntimeException("PPI文件缺失，请检查 resources/static/network/ppi 目录");
        }

        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(ppiResource.getInputStream(), "UTF-8")
        )) {
            String line;
            int lineNum = 0;
            String separator = null;

            while ((line = br.readLine()) != null) {
                lineNum++;
                line = line.trim();

                // 跳过空行和注释行
                if (StringUtils.isEmpty(line) || line.startsWith("#")) {
                    continue;
                }

                // 确定分隔符
                if (separator == null) {
                    for (String sep : SEPARATORS) {
                        if (line.contains(sep)) {
                            separator = sep;
                            log.info("识别PPI文件分隔符：{}（第{}行）", sep, lineNum);
                            break;
                        }
                    }
                    if (separator == null) {
                        log.error("PPI文件格式错误（第{}行）：非TSV/CSV格式", lineNum);
                        throw new RuntimeException("文件格式错误，需为制表符分隔的TSV文件");
                    }
                }

                // 分割行数据
                String[] parts = line.split(separator);
                if (parts.length < 2) {
                    log.warn("跳过无效行（第{}行）：字段数量不足", lineNum);
                    continue;
                }

                String node1 = parts[0].trim();
                String node2 = parts[1].trim();

                // 跳过空节点
                if (StringUtils.isEmpty(node1) || StringUtils.isEmpty(node2)) {
                    log.warn("跳过无效行（第{}行）：蛋白质ID为空", lineNum);
                    continue;
                }

                // 保存原始记录（包含所有字段，供前端表格展示）
                originalRecords.add(new OriginalRecord(parts));

                // 统计所有节点和Node1出现次数
                allNodeNameSet.add(node1);
                allNodeNameSet.add(node2);
                node1CountMap.put(node1, node1CountMap.getOrDefault(node1, 0) + 1);
            }

            // 3. 根据阈值筛选Node1
            Set<String> qualifiedNode1Set = node1CountMap.entrySet().stream()
                    .filter(entry -> entry.getValue() >= degreeThreshold) // 核心筛选条件
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toSet());

            // 4. 筛选关联的节点（保留与合格Node1相关的所有节点）
            Set<String> filteredNodeSet = new HashSet<>();
            List<PpiResponse.Link> filteredLinks = new ArrayList<>();
            List<Map<String, String>> filteredRecords = new ArrayList<>();

            for (OriginalRecord record : originalRecords) {
                String node1 = record.parts[0].trim();
                String node2 = record.parts[1].trim();

                // 只保留合格Node1的记录
                if (qualifiedNode1Set.contains(node1)) {
                    // 保存符合条件的原始记录（转换为Map供前端表格使用）
                    Map<String, String> recordMap = new HashMap<>();
                    recordMap.put("node1", node1);
                    recordMap.put("node2", node2);
                    // 添加其他字段（如得分、来源等，根据实际TSV字段调整）
                    for (int i = 2; i < record.parts.length; i++) {
                        recordMap.put("field" + i, record.parts[i].trim());
                    }
                    filteredRecords.add(recordMap);

                    // 保存关联的节点和边
                    filteredNodeSet.add(node1);
                    filteredNodeSet.add(node2);

                    PpiResponse.Link link = new PpiResponse.Link();
                    link.setSource(node1);
                    link.setTarget(node2);
                    filteredLinks.add(link);
                }
            }

            // 5. 处理筛选结果
            List<PpiResponse.Node> filteredNodes = filteredNodeSet.stream()
                    .map(name -> {
                        PpiResponse.Node node = new PpiResponse.Node();
                        node.setName(name);
                        return node;
                    })
                    .collect(Collectors.toList());

            // 计算筛选后的统计值
            Integer maxCount = qualifiedNode1Set.isEmpty() ? 0 :
                    node1CountMap.entrySet().stream()
                            .filter(entry -> qualifiedNode1Set.contains(entry.getKey()))
                            .map(Map.Entry::getValue)
                            .max(Integer::compare)
                            .orElse(0);

            // 6. 封装返回数据
            PpiResponse response = new PpiResponse();
            response.setNodes(filteredNodes);
            response.setLinks(filteredLinks);
            response.setNode1Counts(node1CountMap); // 保留完整计数，供前端Slider范围使用
            response.setMaxCount(node1CountMap.values().stream().max(Integer::compare).orElse(0)); // 原始最大值
            response.setMinCount(0);
            response.setOriginalRecords(filteredRecords); // 筛选后的原始记录（供前端表格）

            log.info("PPI筛选完成：阈值={}，保留节点数={}，连接数={}，记录数={}",
                    degreeThreshold, filteredNodes.size(), filteredLinks.size(), filteredRecords.size());
            return response;

        } catch (IOException e) {
            log.error("读取PPI文件流失败", e);
            throw new RuntimeException("文件读取错误：" + e.getMessage());
        }
    }

    /**
     * 内部辅助类：存储原始行记录
     */
    private static class OriginalRecord {
        String[] parts;
        OriginalRecord(String[] parts) {
            this.parts = parts;
        }
    }

    /**
     * 兼容原有方法：不筛选，返回全部数据（调用带默认阈值0的筛选方法）
     */
    public static PpiResponse parsePpiResourceFile() {
        return parseAndFilterPpiResourceFile(0); // 阈值0表示不筛选
    }
}
