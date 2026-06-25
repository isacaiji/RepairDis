package com.bishe.ddr_springboot.controller.network;

import com.bishe.ddr_springboot.service.network.NetworkService;
import com.bishe.ddr_springboot.service.network.SyntheticLethalityService;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import javax.annotation.Resource;

/**
 * 合成致死网络
 * 路径：/api/synthetic-lethality
 */
@RestController
@RequestMapping("/api/network/sl")
public class SyntheticLethalityController extends AbstractNetworkController {

    @Resource
    private SyntheticLethalityService syntheticLethalityService; // 需要实现NetworkService接口

    @Override
    protected NetworkService getNetworkService() {
        return syntheticLethalityService;
    }

    @Override
    protected String getExportFileName() {
        return "synthetic_lethality_data";
    }
}