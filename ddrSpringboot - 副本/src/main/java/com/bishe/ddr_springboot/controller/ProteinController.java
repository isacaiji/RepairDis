package com.bishe.ddr_springboot.controller;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.io.Resource;
import org.springframework.core.io.ResourceLoader;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import java.io.IOException;
import java.util.Objects;

@RestController
@RequestMapping("/proteins")
public class ProteinController {

    private final ResourceLoader resourceLoader; // 注入 ResourceLoader

    @Autowired
    public ProteinController(ResourceLoader resourceLoader) {
        this.resourceLoader = resourceLoader;
    }

    @GetMapping("/{filename}")
    public ResponseEntity<Resource> proteinStructure(
            @PathVariable String filename) throws IOException {

        // 使用 classpath: 前缀加载类路径下的资源，路径直接用 / 分隔（Spring 自动适配系统）
        String resourcePath = "classpath:static/protein/" + filename + ".pdb";
        Resource resource = resourceLoader.getResource(resourcePath);

        if (!resource.exists()) {
            return ResponseEntity.notFound().build();
        }

        HttpHeaders headers = new HttpHeaders();
        headers.add(HttpHeaders.CONTENT_DISPOSITION, "attachment; filename=" + resource.getFilename());
        headers.setContentType(MediaType.parseMediaType("chemical/x-pdb"));

        return ResponseEntity.ok()
                .headers(headers)
                .contentLength(resource.contentLength())
                .contentType(Objects.requireNonNull(headers.getContentType()))
                .body(resource);
    }
}