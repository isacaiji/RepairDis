package com.bishe.ddr_springboot.entity;

import lombok.Data;

@Data
public class CancerScore {
    private String gene;
    private Double totalScore;
    private String scoreMethod;
    private Double acc;
    private Double blca;
    private Double brca;
    private Double cesc;
    private Double chol;
    private Double coad;
    private Double dlbc;
    private Double esca;
    private Double gbm;
    private Double hnsc;
    private Double kich;
    private Double kirc;
    private Double kirp;
    private Double laml;
    private Double lgg;
    private Double lihc;
    private Double luad;
    private Double lusc;
    private Double meso;
    private Double ov;
    private Double paad;
    private Double pcpg;
    private Double prad;
    private Double read;
    private Double sarc;
    private Double skcm;
    private Double stad;
    private Double tgct;
    private Double thca;
    private Double thym;
    private Double ucec;
    private Double ucs;
    private Double uvm;
}
