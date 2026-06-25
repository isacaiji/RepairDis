package com.bishe.ddr_springboot.controller.network;

import com.bishe.ddr_springboot.service.network.NetworkService;
import com.bishe.ddr_springboot.service.network.TFRegulatoryService;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import javax.annotation.Resource;

/**
 * TF调控网络控制器（继承抽象网络控制器）
 * 路径：/api/tf-regulatory
 */
@RestController
@RequestMapping("/api/network/tf")
public class TFRegulatoryController extends AbstractNetworkController {

    @Resource
    private TFRegulatoryService tfRegulatoryService; // 需要实现NetworkService接口

    @Override
    protected NetworkService getNetworkService() {
        return tfRegulatoryService;
    }

    @Override
    protected String getExportFileName() {
        return "tf_regulatory_data";
    }
}