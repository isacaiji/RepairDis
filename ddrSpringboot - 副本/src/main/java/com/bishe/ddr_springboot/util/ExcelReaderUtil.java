package com.bishe.ddr_springboot.util;

import com.alibaba.excel.EasyExcel;
import com.bishe.ddr_springboot.entity.TFRegulation;  // 导入entity包下的类
import org.springframework.core.io.ClassPathResource;
import org.springframework.stereotype.Component;

import java.io.InputStream;
import java.util.List;

@Component
public class ExcelReaderUtil {

    /**
     * 读取resources/static/tf/TF.xlsx中的调控关系数据
     */
    public List<TFRegulation> readTFRegulations() {
        try {
            ClassPathResource resource = new ClassPathResource("static/tf/TF.xlsx");
            InputStream inputStream = resource.getInputStream();

            // 解析Excel并映射到entity类
            List<TFRegulation> regulations = EasyExcel.read(inputStream)
                    .head(TFRegulation.class)
                    .sheet()
                    .doReadSync();

            inputStream.close();
            return regulations;
        } catch (Exception e) {
            throw new RuntimeException("读取TF.xlsx失败：" + e.getMessage());
        }
    }
}