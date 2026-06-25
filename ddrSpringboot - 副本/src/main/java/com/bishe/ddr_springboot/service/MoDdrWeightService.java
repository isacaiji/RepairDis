package com.bishe.ddr_springboot.service;

import com.bishe.ddr_springboot.entity.CancerScore;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

@Service
public class MoDdrWeightService {

    private static final Logger log = LoggerFactory.getLogger(MoDdrWeightService.class);
    private static final String SCORE_RESOURCE = "data/mo_ddrweight_by_cancer.csv";
    private static final String SCORE_METHOD =
            "Gene-level MO-DDRweight from pan-cancer evidence fusion, displayed as MO_DDRweight x 100; mean score is the average across available cancer types.";

    private final Map<String, CancerScore> scoresByGene = new ConcurrentHashMap<>();

    @PostConstruct
    public void loadScores() {
        Map<String, Map<String, Double>> geneCancerScores = new LinkedHashMap<>();

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(
                new ClassPathResource(SCORE_RESOURCE).getInputStream(), StandardCharsets.UTF_8))) {

            String headerLine = reader.readLine();
            if (headerLine == null) {
                log.warn("MO-DDRweight resource is empty: {}", SCORE_RESOURCE);
                return;
            }

            Map<String, Integer> header = indexHeader(parseCsvLine(headerLine));
            Integer cancerIndex = header.get("Cancer");
            Integer geneIndex = header.get("Gene");
            Integer scoreIndex = header.get("MO_DDRweight");

            if (cancerIndex == null || geneIndex == null || scoreIndex == null) {
                log.warn("MO-DDRweight resource missing required columns: {}", SCORE_RESOURCE);
                return;
            }

            String line;
            while ((line = reader.readLine()) != null) {
                List<String> fields = parseCsvLine(line);
                if (fields.size() <= Math.max(scoreIndex, Math.max(cancerIndex, geneIndex))) {
                    continue;
                }

                String cancer = normalize(fields.get(cancerIndex));
                String gene = normalize(fields.get(geneIndex));
                Double weight = parseDouble(fields.get(scoreIndex));
                if (cancer.isEmpty() || gene.isEmpty() || weight == null) {
                    continue;
                }

                geneCancerScores
                        .computeIfAbsent(gene, key -> new LinkedHashMap<>())
                        .put(cancer, round2(weight * 100.0));
            }

            for (Map.Entry<String, Map<String, Double>> entry : geneCancerScores.entrySet()) {
                CancerScore score = buildCancerScore(entry.getKey(), entry.getValue());
                scoresByGene.put(entry.getKey(), score);
            }

            log.info("Loaded MO-DDRweight scores for {} genes from {}", scoresByGene.size(), SCORE_RESOURCE);
        } catch (IOException e) {
            log.warn("Failed to load MO-DDRweight resource: {}", SCORE_RESOURCE, e);
        }
    }

    public CancerScore getCancerScoreByGene(String geneName) {
        if (geneName == null) {
            return null;
        }
        return scoresByGene.get(normalize(geneName));
    }

    public Double getTotalScoreByGene(String geneName) {
        CancerScore score = getCancerScoreByGene(geneName);
        return score == null ? null : score.getTotalScore();
    }

    private CancerScore buildCancerScore(String gene, Map<String, Double> cancerScores) {
        CancerScore score = new CancerScore();
        score.setGene(gene);
        score.setScoreMethod(SCORE_METHOD);

        double sum = 0.0;
        int n = 0;
        for (Map.Entry<String, Double> entry : cancerScores.entrySet()) {
            Double value = entry.getValue();
            if (value == null) {
                continue;
            }
            setCancerScore(score, entry.getKey(), value);
            sum += value;
            n++;
        }

        score.setTotalScore(n == 0 ? null : round2(sum / n));
        return score;
    }

    private void setCancerScore(CancerScore score, String cancer, Double value) {
        switch (cancer) {
            case "ACC": score.setAcc(value); break;
            case "BLCA": score.setBlca(value); break;
            case "BRCA": score.setBrca(value); break;
            case "CESC": score.setCesc(value); break;
            case "CHOL": score.setChol(value); break;
            case "COAD": score.setCoad(value); break;
            case "DLBC": score.setDlbc(value); break;
            case "ESCA": score.setEsca(value); break;
            case "GBM": score.setGbm(value); break;
            case "HNSC": score.setHnsc(value); break;
            case "KICH": score.setKich(value); break;
            case "KIRC": score.setKirc(value); break;
            case "KIRP": score.setKirp(value); break;
            case "LAML": score.setLaml(value); break;
            case "LGG": score.setLgg(value); break;
            case "LIHC": score.setLihc(value); break;
            case "LUAD": score.setLuad(value); break;
            case "LUSC": score.setLusc(value); break;
            case "MESO": score.setMeso(value); break;
            case "OV": score.setOv(value); break;
            case "PAAD": score.setPaad(value); break;
            case "PCPG": score.setPcpg(value); break;
            case "PRAD": score.setPrad(value); break;
            case "READ": score.setRead(value); break;
            case "SARC": score.setSarc(value); break;
            case "SKCM": score.setSkcm(value); break;
            case "STAD": score.setStad(value); break;
            case "TGCT": score.setTgct(value); break;
            case "THCA": score.setThca(value); break;
            case "THYM": score.setThym(value); break;
            case "UCEC": score.setUcec(value); break;
            case "UCS": score.setUcs(value); break;
            case "UVM": score.setUvm(value); break;
            default: break;
        }
    }

    private Map<String, Integer> indexHeader(List<String> columns) {
        Map<String, Integer> index = new HashMap<>();
        for (int i = 0; i < columns.size(); i++) {
            index.put(columns.get(i), i);
        }
        return index;
    }

    private String normalize(String value) {
        return value == null ? "" : value.trim().toUpperCase(Locale.ROOT);
    }

    private Double parseDouble(String value) {
        if (value == null || value.trim().isEmpty()) {
            return null;
        }
        try {
            return Double.parseDouble(value.trim());
        } catch (NumberFormatException e) {
            return null;
        }
    }

    private double round2(double value) {
        return Math.round(value * 100.0) / 100.0;
    }

    private List<String> parseCsvLine(String line) {
        List<String> fields = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean inQuotes = false;

        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (ch == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    current.append('"');
                    i++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else if (ch == ',' && !inQuotes) {
                fields.add(current.toString());
                current.setLength(0);
            } else {
                current.append(ch);
            }
        }
        fields.add(current.toString());
        return fields;
    }
}
