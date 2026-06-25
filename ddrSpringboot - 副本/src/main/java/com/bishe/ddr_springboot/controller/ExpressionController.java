package com.bishe.ddr_springboot.controller;

import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpStatus;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

@RestController
@RequestMapping("/expression")
public class ExpressionController {

    /**
     * 1) 返回 PNG 图片 —— 前端展示使用
     * GET /expression/img?gene=TP53
     */
    @GetMapping("/img")
    public ResponseEntity<byte[]> getGenePng(@RequestParam String gene) {

        String path = "static/expression/" + gene + ".png";
        Resource resource = new ClassPathResource(path);

        try {
            if (!resource.exists()) {
                return ResponseEntity.status(HttpStatus.NOT_FOUND)
                        .body(("PNG not found: " + gene).getBytes("UTF-8"));
            }

            byte[] imageBytes = readFileToBytes(resource);

            HttpHeaders headers = new HttpHeaders();
            headers.setContentType(MediaType.IMAGE_PNG);

            return new ResponseEntity<>(imageBytes, headers, HttpStatus.OK);

        } catch (Exception e) {
            e.printStackTrace();
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(null);
        }
    }


    /**
     * 2) 返回 PDF 文件 —— 前端下载使用
     * GET /expression/pdf?gene=TP53
     */
    @GetMapping("/pdf")
    public ResponseEntity<byte[]> getGenePdf(@RequestParam String gene) {

        String path = "static/expression/" + gene + ".pdf";
        Resource resource = new ClassPathResource(path);

        try {
            if (!resource.exists()) {
                return ResponseEntity.status(HttpStatus.NOT_FOUND)
                        .body(("PDF not found: " + gene).getBytes("UTF-8"));
            }

            byte[] pdfBytes = readFileToBytes(resource);

            HttpHeaders headers = new HttpHeaders();
            headers.setContentType(MediaType.APPLICATION_PDF);
            headers.setContentDispositionFormData(gene + ".pdf", gene + ".pdf");

            return new ResponseEntity<>(pdfBytes, headers, HttpStatus.OK);

        } catch (Exception e) {
            e.printStackTrace();
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(null);
        }
    }


    // 工具方法：Java8 兼容的文件读取方式
    private byte[] readFileToBytes(Resource resource) throws IOException {

        InputStream inputStream = resource.getInputStream();
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();

        byte[] data = new byte[1024];
        int nRead;

        while ((nRead = inputStream.read(data, 0, data.length)) != -1) {
            buffer.write(data, 0, nRead);
        }

        buffer.flush();
        return buffer.toByteArray();
    }
}
