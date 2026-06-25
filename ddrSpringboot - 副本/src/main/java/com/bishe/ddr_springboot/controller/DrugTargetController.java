package com.bishe.ddr_springboot.controller;

import com.bishe.ddr_springboot.entity.DrugTarget;
import com.bishe.ddr_springboot.entity.PageResult;
import com.bishe.ddr_springboot.service.DrugTargetService;
import lombok.extern.slf4j.Slf4j;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import javax.annotation.Resource;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * 药物靶点关系控制器
 */
@Slf4j
@RestController
@RequestMapping("/api/drug-target")
public class DrugTargetController {

    @Resource
    private DrugTargetService drugTargetService;

    /**
     * 获取所有药物靶点关系
     */
    @GetMapping("/all")
    public List<DrugTarget> getAll() {
        return drugTargetService.getAllDrugTargets();
    }

    /**
     * 分页查询
     */
    @GetMapping("/page")
    public PageResult<DrugTarget> getByPage(
            @RequestParam(required = false) Integer pageNum,
            @RequestParam(required = false) Integer pageSize) {
        return drugTargetService.getDrugTargetsByPage(pageNum, pageSize);
    }

    /**
     * 按药物名称查询
     */
    @GetMapping("/by-drug")
    public List<DrugTarget> getByDrugName(@RequestParam String drugName) {
        return drugTargetService.getByDrugName(drugName);
    }

    /**
     * 按基因名称查询
     */
    @GetMapping("/by-gene")
    public List<DrugTarget> getByGeneName(@RequestParam String geneName) {
        return drugTargetService.getByGeneName(geneName);
    }

    /**
     * 多条件查询
     */
    @PostMapping("/query")
    public PageResult<DrugTarget> queryByConditions(@RequestBody Map<String, Object> params) {
        return drugTargetService.queryByConditions(params);
    }

    /**
     * 导出全部数据
     */
    @GetMapping("/export/all")
    public ResponseEntity<byte[]> exportAll() {
        byte[] data = drugTargetService.exportAll();
        return buildExportResponse(data, "all_drug_targets.tsv");
    }

    /**
     * 按条件导出数据
     */
    @PostMapping("/export/conditions")
    public ResponseEntity<byte[]> exportByConditions(@RequestBody Map<String, Object> params) {
        byte[] data = drugTargetService.exportByConditions(params);
        return buildExportResponse(data, "filtered_drug_targets.tsv");
    }

    /**
     * 获取药物-靶点网络数据
     */
    @GetMapping("/network")
    public ResponseEntity<Map<String, Object>> getNetwork(
            @RequestParam Map<String, Object> params) {
        try {
            Map<String, Object> networkData = drugTargetService.getNetworkData(params);
            return ResponseEntity.ok(networkData);
        } catch (Exception e) {
            log.error("获取药物-靶点网络失败", e);
            return ResponseEntity.status(500).body(Collections.singletonMap("error", e.getMessage()));
        }
    }

    /**
     * 药物名称联想建议（去掉limit参数）
     */
    @GetMapping("/suggest/drug")
    public List<String> getDrugSuggestions(@RequestParam String keyword) {
        try {
            return drugTargetService.getDrugNameSuggestions(keyword);
        } catch (Exception e) {
            log.error("获取药物联想建议失败", e);
            return Collections.emptyList();
        }
    }

    /**
     * 基因名称联想建议（去掉limit参数）
     */
    @GetMapping("/suggest/gene")
    public List<String> getGeneSuggestions(@RequestParam String keyword) {
        try {
            return drugTargetService.getGeneNameSuggestions(keyword);
        } catch (Exception e) {
            log.error("获取基因联想建议失败", e);
            return Collections.emptyList();
        }
    }

    /**
     * 构建导出响应
     */
    private ResponseEntity<byte[]> buildExportResponse(byte[] data, String fileName) {
        HttpHeaders headers = new HttpHeaders();
        headers.add("Content-Disposition", "attachment; filename=\"" + fileName + "\"");
        headers.add("Content-Type", "text/tab-separated-values; charset=UTF-8");
        return new ResponseEntity<>(data, headers, HttpStatus.OK);
    }
}