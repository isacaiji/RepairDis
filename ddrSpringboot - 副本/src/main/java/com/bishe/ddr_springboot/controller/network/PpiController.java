package com.bishe.ddr_springboot.controller.network;

import com.bishe.ddr_springboot.service.network.NetworkService;
import com.bishe.ddr_springboot.service.network.PpiNetworkService;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import javax.annotation.Resource;

@RestController
@RequestMapping("/api/network/ppi")
public class PpiController extends AbstractNetworkController {

    @Resource
    private PpiNetworkService ppiNetworkService;

    @Override
    protected NetworkService getNetworkService() {
        return ppiNetworkService;
    }

    @Override
    protected String getExportFileName() {
        return "ppi_network_data";
    }
}