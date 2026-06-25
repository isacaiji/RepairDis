package com.bishe.ddr_springboot.controller.network;

import com.bishe.ddr_springboot.service.network.NetworkService;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpStatus;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.util.StringUtils;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestParam;

import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.Map;

public abstract class AbstractNetworkController {
    protected abstract NetworkService getNetworkService();

    /**
     * 获取网络数据
     */
    @GetMapping
    public ResponseEntity<?> getNetworkData(@RequestParam Map<String, Object> params) {
        try {
            Map<String, Object> result = getNetworkService().getNetworkData(params);
            return ResponseEntity.ok(result);
        } catch (IllegalArgumentException e) {
            return ResponseEntity.badRequest().body(e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body("获取网络数据失败：" + e.getMessage());
        }
    }

    /**
     * 获取表格分页数据
     */
    @GetMapping("/table")
    public ResponseEntity<?> getTableData(@RequestParam Map<String, Object> params) {
        try {
            // 解析分页参数
            int pageNum = params.containsKey("pageNum") ? Integer.parseInt(params.get("pageNum").toString()) : 1;
            int pageSize = params.containsKey("pageSize") ? Integer.parseInt(params.get("pageSize").toString()) : 10;
            params.put("pageNum", pageNum);
            params.put("pageSize", pageSize);

            Map<String, Object> result = getNetworkService().getTableData(params);
            return ResponseEntity.ok(result);
        } catch (NumberFormatException e) {
            return ResponseEntity.badRequest().body("分页参数格式错误");
        } catch (Exception e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body("获取表格数据失败：" + e.getMessage());
        }
    }

    /**
     * 导出全库数据（所有记录）
     */
    @GetMapping("/export")
    public ResponseEntity<byte[]> exportAllData() { // 无需参数
        try {
            // 注意：这里需要强转，或扩展 NetworkService 接口
            byte[] data = getNetworkService().exportAllData();
            String fileName = getExportFileName() + "_full";
            return createExportResponse(data, fileName);
        } catch (Exception e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body(("Export all failed: " + e.getMessage()).getBytes(StandardCharsets.UTF_8));
        }
    }
    /**
     * 导出当前查询条件下的数据（即“当前结果”）
     */
    @GetMapping("/export-current") // 建议改名，避免 confusion
    public ResponseEntity<byte[]> exportCurrentData(
            @RequestParam String name,
            @RequestParam(required = false) Double degreeThreshold) {
        if (!StringUtils.hasText(name)) {
            return ResponseEntity.badRequest()
                    .body("参数 'name' 不能为空".getBytes(StandardCharsets.UTF_8));
        }
        try {
            byte[] data = getNetworkService().exportCurrentData(name.trim(),degreeThreshold);
            String fileName = getExportFileName() + "_for_" + name.trim();
            return createExportResponse(data, fileName);
        } catch (Exception e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body(("导出当前数据失败: " + e.getMessage()).getBytes(StandardCharsets.UTF_8));
        }
    }

    // 创建导出响应
    protected ResponseEntity<byte[]> createExportResponse(byte[] data, String baseFileName) {
        HttpHeaders headers = new HttpHeaders();
        headers.setContentType(MediaType.parseMediaType("text/tab-separated-values"));
        String fileName = baseFileName + ".tsv";

        String encodedFileName = new String(
                fileName.getBytes(StandardCharsets.UTF_8),
                StandardCharsets.ISO_8859_1
        );
        headers.setContentDispositionFormData("attachment", encodedFileName);
        return new ResponseEntity<>(data, headers, HttpStatus.OK);
    }

    // 获取导出文件名（由子类实现）
    protected abstract String getExportFileName();
}
