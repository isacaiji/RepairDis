package com.bishe.ddr_springboot.controller;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.rendering.PDFRenderer;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import javax.annotation.PostConstruct;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.Base64;

@RestController
@RequestMapping("/cluster")
public class ClusterController {

    // ================================
    //       缓存结构（带 TTL）
    // ================================
    private static class CacheEntry {
        List<String> images;
        long expireAt;

        CacheEntry(List<String> images, long expireAt) {
            this.images = images;
            this.expireAt = expireAt;
        }
    }

    // 30 分钟缓存时间（毫秒）
    private static final long CACHE_TTL = 30 * 60 * 1000;

    // 缓存 Map（线程安全）
    private static final ConcurrentHashMap<String, CacheEntry> CACHE = new ConcurrentHashMap<>();


    // 启动时运行一个定时任务：每 10 分钟清理过期缓存
    @PostConstruct
    public void initCleaner() {
        ScheduledExecutorService scheduler = Executors.newSingleThreadScheduledExecutor();
        scheduler.scheduleAtFixedRate(() -> {
            long now = System.currentTimeMillis();
            CACHE.entrySet().removeIf(e -> e.getValue().expireAt < now);
        }, 10, 10, TimeUnit.MINUTES);
    }


    @GetMapping("/consensus/images")
    public ResponseEntity<?> getConsensusImages(@RequestParam String cancer) {
        int[] pages = {1, 6};  // 第 2 页、第 7 页
        return loadPagesWithCache(cancer, "consensus.pdf", pages);
    }


    @GetMapping("/icl/images")
    public ResponseEntity<?> getICLImages(@RequestParam String cancer) {
        int[] pages = {2}; // 第 3 页
        return loadPagesWithCache(cancer, "icl.pdf", pages);
    }

    // ================================
    //       新增：survival.pdf 接口
    // ================================
    @GetMapping("/survival/images")
    public ResponseEntity<?> getSurvivalImages(@RequestParam String cancer) {
        int[] pages = {1}; // PDF 索引从 0 开始，1 对应第二页
        return loadPagesWithCache(cancer, "survival.pdf", pages);
    }


    /**
     * 带自动过期缓存的 PDF 请求
     */
    private ResponseEntity<?> loadPagesWithCache(String cancer, String filename, int[] pageIndexes) {
        try {
            String cacheKey = cancer + "_" + filename;
            long now = System.currentTimeMillis();

            // 1. 缓存检查（不变）
            CacheEntry entry = CACHE.get(cacheKey);
            if (entry != null && entry.expireAt > now) {
                return ResponseEntity.ok(entry.images);
            }

            // 2. 加载JAR包内资源（核心修改：使用ClassPathResource）
            String folderName = cancer + "-0909";
            // 资源路径：直接写classpath下的相对路径（无需拼接user.dir），Spring自动识别
            String resourcePath = "static/CLUSTER/" + folderName + "/" + filename;
            ClassPathResource resource = new ClassPathResource(resourcePath);

            // 检查资源是否存在（兼容本地文件夹和JAR包环境）
            if (!resource.exists()) {
                return ResponseEntity.status(404).body("File not found: " + filename);
            }

            // 关键：通过getResourceAsStream()获取资源输入流（而非File对象）
            // PDDocument支持通过InputStream加载文件，无需File对象
            InputStream inputStream = resource.getInputStream();
            PDDocument document = PDDocument.load(inputStream); // 改用InputStream加载
            PDFRenderer renderer = new PDFRenderer(document);

            List<String> images = new ArrayList<>();
            // PDF渲染逻辑（完全不变）
            for (int pageIndex : pageIndexes) {
                if (pageIndex < 0 || pageIndex >= document.getNumberOfPages()) {
                    images.add("");
                    continue;
                }
                BufferedImage image = renderer.renderImageWithDPI(pageIndex, 120);
                ByteArrayOutputStream baos = new ByteArrayOutputStream();
                ImageIO.write(image, "png", baos);
                String base64 = "data:image/png;base64," +
                        Base64.getEncoder().encodeToString(baos.toByteArray());
                images.add(base64);
            }

            document.close();
            inputStream.close(); // 关闭输入流

            // 3. 更新缓存（不变）
            CACHE.put(cacheKey, new CacheEntry(images, now + CACHE_TTL));
            return ResponseEntity.ok(images);

        } catch (Exception e) {
            return ResponseEntity.status(500).body("Error: " + e.getMessage());
        }
    }

}