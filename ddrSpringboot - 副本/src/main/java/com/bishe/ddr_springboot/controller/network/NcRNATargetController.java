package com.bishe.ddr_springboot.controller.network;

import com.bishe.ddr_springboot.service.network.NetworkService;
import com.bishe.ddr_springboot.service.network.NcRNATargetService;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import javax.annotation.Resource;

/**
 * ncRNA靶向关系控制器（继承抽象网络控制器，复用通用接口）
 */
@RestController
@RequestMapping("/api/network/ncrna")
public class NcRNATargetController extends AbstractNetworkController {

    @Resource
    private NcRNATargetService ncRNATargetService;

    // 提供具体的服务实现
    @Override
    protected NetworkService getNetworkService() {
        return ncRNATargetService;
    }

    // 导出文件默认名称
    @Override
    protected String getExportFileName() {
        return "ncrna_target_relationships";
    }
}