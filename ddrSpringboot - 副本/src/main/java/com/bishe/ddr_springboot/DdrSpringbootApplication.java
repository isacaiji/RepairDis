package com.bishe.ddr_springboot;

import org.mybatis.spring.annotation.MapperScan;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

@SpringBootApplication
@MapperScan("com.bishe.ddr_springboot.mapper") // 扫描包
public class DdrSpringbootApplication {

    public static void main(String[] args) {
        SpringApplication.run(DdrSpringbootApplication.class, args);
    }

}
