package com.bishe.ddr_springboot.entity;

import javax.persistence.*;
import com.fasterxml.jackson.annotation.JsonProperty;

@Entity
@Table(name = "ddrgenes")
public class Gene {
    @Id
    @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Integer id;
    private String geneName;
    private String pathway;
    private String Ensembl;
    private String Function;
    private String Title;
    private String Abstract;
    private String Year;
    private String PMID;
    @Transient
    private Double drfsScore;
    @Transient
    private Double meanMoDdrWeight;


    public Integer getId() {
        return id;
    }

    public void setId(Integer id) {
        this.id = id;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public String getPathway() {
        return pathway;
    }

    public void setPathway(String pathway) {
        this.pathway = pathway;
    }

    public String getEnsembl() {
        return Ensembl;
    }

    public void setEnsembl(String ensembl) {
        Ensembl = ensembl;
    }

    public String getFunction() {
        return Function;
    }

    public void setFunction(String function) {
        Function = function;
    }

    public String getTitle() {
        return Title;
    }

    public void setTitle(String title) {
        Title = title;
    }

    public String getAbstract() {
        return Abstract;
    }

    public void setAbstract(String anAbstract) {
        Abstract = anAbstract;
    }

    public String getYear() {
        return Year;
    }

    public void setYear(String year) {
        Year = year;
    }

    public String getPMID() {
        return PMID;
    }

    public void setPMID(String PMID) {
        this.PMID = PMID;
    }

    public Double getDrfsScore() {return drfsScore;}

    public void setDrfsScore(Double drfsScore) {
        this.drfsScore = drfsScore;
        this.meanMoDdrWeight = drfsScore;
    }

    public Double getMeanMoDdrWeight() {return meanMoDdrWeight;}

    public void setMeanMoDdrWeight(Double meanMoDdrWeight) {
        this.meanMoDdrWeight = meanMoDdrWeight;
        this.drfsScore = meanMoDdrWeight;
    }
}
