package com.bishe.ddr_springboot.controller;

import com.bishe.ddr_springboot.entity.CancerScore;
import com.bishe.ddr_springboot.entity.Gene;
import com.bishe.ddr_springboot.entity.GeneSummary;
import com.bishe.ddr_springboot.entity.PageResult;
import com.bishe.ddr_springboot.service.GeneService;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.rendering.PDFRenderer;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.http.*;
import org.springframework.web.bind.annotation.*;

import javax.imageio.ImageIO;
import javax.servlet.http.HttpServletResponse;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.List;

@RestController
@RequestMapping("/api/genes")
public class GeneController {

    private final GeneService geneService;
    // 热图文件存储路径（相对于resources目录）
    private static final String HEATMUTATION_DIR = "static/heatmutations/";

    public GeneController(GeneService geneService) {
        this.geneService = geneService;
    }

    /**
     * 分页查询基因摘要
     */
    @GetMapping("/summary")
    public ResponseEntity<PageResult<GeneSummary>> getGeneSummariesByPage(
            @RequestParam(required = false, defaultValue = "1") Integer pageNum,
            @RequestParam(required = false, defaultValue = "15") Integer pageSize) {
        PageResult<GeneSummary> pageResult = geneService.getGeneSummariesByPage(pageNum, pageSize);
        return new ResponseEntity<>(pageResult, HttpStatus.OK);
    }

    /**
     * 查询全量基因摘要
     */
    @GetMapping("/summary/all")
    public ResponseEntity<List<GeneSummary>> getAllGeneSummaries() {
        List<GeneSummary> summaries = geneService.getGeneSummaries();
        return new ResponseEntity<>(summaries, HttpStatus.OK);
    }

    /**
     * 获取所有基因完整信息
     */
    @GetMapping
    public ResponseEntity<List<Gene>> getGenes() {
        List<Gene> genes = geneService.getGenes();
        return new ResponseEntity<>(genes, HttpStatus.OK);
    }

    /**
     * 获取所有基因名称
     */
    @GetMapping("/all")
    public List<String> getGeneNames() {
        return geneService.getAllGeneNames();
    }

    /**
     * 根据ID获取单个基因详情
     */
    @GetMapping("/{id}")
    public ResponseEntity<Gene> getGeneById(@PathVariable Integer id) {
        Gene gene = geneService.getGeneById(id);
        if (gene != null) {
            return new ResponseEntity<>(gene, HttpStatus.OK);
        } else {
            return new ResponseEntity<>(HttpStatus.NOT_FOUND);
        }
    }

    /**
     * 根据名称获取单个基因详情
     */
    @GetMapping("/name/{name}")
    public ResponseEntity<Gene> getGeneByName(@PathVariable String name) {
        Gene gene = geneService.getGeneByName(name);
        if (gene != null) {
            return new ResponseEntity<>(gene, HttpStatus.OK);
        } else {
            return new ResponseEntity<>(HttpStatus.NOT_FOUND);
        }
    }

    /**
     * 根据关键词搜索基因
     */
    @GetMapping("/search")
    public ResponseEntity<List<Gene>> getGenesByQuery(@RequestParam String query) {
        List<Gene> genes = geneService.getGenesByQuery(query);
        return new ResponseEntity<>(genes, HttpStatus.OK);
    }


    // 2. 新增图片转换接口
    @GetMapping("/{geneName}/heatmutation-image")
    public ResponseEntity<byte[]> getHeatmapAsImage(@PathVariable String geneName, HttpServletResponse response) {
        response.setHeader("Access-Control-Allow-Origin", "*");

        try {
            // 1. 读取PDF文件
            String pdfFileName = geneName + "_lollipop (1).pdf";
            String pdfFilePath = "static/heatmutations/" + pdfFileName;
            Resource pdfResource = new ClassPathResource(pdfFilePath);

            if (!pdfResource.exists()) {
                return ResponseEntity.status(HttpStatus.NOT_FOUND)
                        .body(("No PDF file found for " + geneName).getBytes());
            }

            // 2. 使用PDFBox将PDF转为PNG图片
            PDDocument document = PDDocument.load(pdfResource.getInputStream());
            PDFRenderer renderer = new PDFRenderer(document);

            // 只转换第一页（如果是多页PDF）
            BufferedImage image = renderer.renderImage(0, 1.5f); // 1.5f为缩放比例

            // 3. 将图片转为字节数组
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            ImageIO.write(image, "png", baos);
            byte[] imageBytes = baos.toByteArray();

            // 4. 关闭资源
            document.close();
            baos.close();

            // 5. 设置图片响应头
            HttpHeaders headers = new HttpHeaders();
            headers.setContentType(MediaType.IMAGE_PNG);
            headers.setContentLength(imageBytes.length);

            return new ResponseEntity<>(imageBytes, headers, HttpStatus.OK);

        } catch (IOException e) {
            e.printStackTrace();
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body(("Failed to convert PDF to image: " + e.getMessage()).getBytes());
        }
    }
    /**
     * 新增：获取基因的突变热图PDF文件
     * 接口：GET /api/genes/{geneName}/heatmutation
     *
     * @param geneName 基因名称
     * @return PDF文件的字节流响应
     */
    @GetMapping("/{geneName}/heatmutation")
    public ResponseEntity<byte[]> getHeatMutationFile(@PathVariable String geneName, HttpServletResponse response) {
        // 允许跨域（根据实际情况调整Origin）
        response.setHeader("Access-Control-Allow-Origin", "*");
        response.setHeader("Access-Control-Allow-Methods", "GET");

        try {
            String fileName = geneName + "_lollipop.pdf";
            String filePath = HEATMUTATION_DIR + fileName;
            Resource resource = new ClassPathResource(filePath);

            if (!resource.exists()) {
                String errorMsg = "Heat mutation file not found for gene: " + geneName;
                return ResponseEntity.status(HttpStatus.NOT_FOUND)
                        .body(errorMsg.getBytes(StandardCharsets.UTF_8));
            }

            // 读取文件（Java 8兼容方式）
            try (InputStream inputStream = resource.getInputStream();
                 ByteArrayOutputStream buffer = new ByteArrayOutputStream()) {

                int nRead;
                byte[] data = new byte[1024];
                while ((nRead = inputStream.read(data, 0, data.length)) != -1) {
                    buffer.write(data, 0, nRead);
                }
                byte[] fileContent = buffer.toByteArray();

                // 关键修改：设置为inline，让浏览器优先预览
                HttpHeaders headers = new HttpHeaders();
                headers.setContentType(MediaType.APPLICATION_PDF);
                headers.setContentLength(fileContent.length);

                // 重点：使用inline而非form-data，避免触发下载工具
                String encodedFileName = new String(
                        fileName.getBytes(StandardCharsets.UTF_8),
                        StandardCharsets.ISO_8859_1
                );
                headers.setContentDisposition(
                        ContentDisposition.inline()  //  inline表示在浏览器内打开
                                .filename(encodedFileName, StandardCharsets.ISO_8859_1)
                                .build()
                );

                return new ResponseEntity<>(fileContent, headers, HttpStatus.OK);
            }

        } catch (IOException e) {
            String errorMsg = "Failed to read heat mutation file: " + e.getMessage();
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
                    .body(errorMsg.getBytes(StandardCharsets.UTF_8));
        }
    }

    @GetMapping("/{geneName}/score")
    public ResponseEntity<CancerScore> getCancerScores(@PathVariable String geneName) {
        CancerScore scores = geneService.getCancerScoresByGeneName(geneName);
        if (scores != null) {
            return new ResponseEntity<>(scores, HttpStatus.OK);
        } else {
            return new ResponseEntity<>(HttpStatus.NOT_FOUND);
        }
    }
}