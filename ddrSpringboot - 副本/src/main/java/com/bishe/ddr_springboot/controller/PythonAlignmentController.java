package com.ysy.data.controller;

import org.springframework.core.io.InputStreamResource;
import org.springframework.http.*;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

@RestController
@CrossOrigin
@RequestMapping("/msa")
public class PythonAlignmentController {

  private static final String UPLOAD_DIR = "/shujupan/adi/fasta/uploads/";
  private static final String RESULT_DIR = "/shujupan/adi/fasta/results/";
  private static final String PYTHON_SCRIPT_PATH = "/shujupan/adi/fasta/geneAligned.py"; // 修改为你的实际路径（相对或绝对）

  /**
   * 接收上传文件，保存并调用 Python 脚本处理
   */
  @PostMapping("/upload")
  public ResponseEntity<?> handleUpload(@RequestParam("file") MultipartFile file) {
    if (file.isEmpty()) {
      return ResponseEntity.badRequest().body("上传文件为空");
    }

    try {
      // 创建上传目录
      File uploadDir = new File(UPLOAD_DIR);
      if (!uploadDir.exists()) uploadDir.mkdirs();

      // 保存上传的文件
      String originalFilename = file.getOriginalFilename();
      File savedFile = new File(UPLOAD_DIR + originalFilename);
      file.transferTo(savedFile);

      // 调用 Python 脚本进行处理
      String resultLog = runPython(savedFile.getAbsolutePath(), RESULT_DIR);

      // 输出文件名
      String baseName = originalFilename.substring(0, originalFilename.lastIndexOf('.'));
      String resultFileName = baseName + "Aligned.fasta";

      Map<String, String> result = new HashMap<>();
      result.put("filename", resultFileName);
      return ResponseEntity.ok(result);

    } catch (Exception e) {
      return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
        .body("处理失败：" + e.getMessage());
    }
  }

  /**
   * 下载处理后生成的结果文件
   */
  @GetMapping("/download/{filename}")
  public ResponseEntity<?> downloadFile(@PathVariable("filename") String filename) {
    File file = new File(RESULT_DIR + filename);
    if (!file.exists()) {
      return ResponseEntity.status(HttpStatus.NOT_FOUND).body("文件不存在: " + filename);
    }

    try {
      InputStreamResource resource = new InputStreamResource(new FileInputStream(file));
      HttpHeaders headers = new HttpHeaders();
      headers.add(HttpHeaders.CONTENT_DISPOSITION, "attachment; filename=" + file.getName());

      return ResponseEntity.ok()
        .headers(headers)
        .contentLength(file.length())
        .contentType(MediaType.APPLICATION_OCTET_STREAM)
        .body(resource);

    } catch (IOException e) {
      return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR)
        .body("下载出错：" + e.getMessage());
    }
  }

  /**
   * Java 调用 Python 脚本
   */
  private String runPython(String inputPath, String outputDir) throws IOException, InterruptedException {
    ProcessBuilder builder = new ProcessBuilder(
      "python3", PYTHON_SCRIPT_PATH, inputPath, outputDir
    );
    builder.directory(new File(".")); // 项目根目录
    builder.redirectErrorStream(true);

    Process process = builder.start();

    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
    StringBuilder output = new StringBuilder();
    String line;
    while ((line = reader.readLine()) != null) {
      output.append(line).append("\n");
    }

    int exitCode = process.waitFor();
    if (exitCode != 0) {
      throw new RuntimeException("Python脚本执行失败：\n" + output);
    }

    return output.toString();
  }
}
