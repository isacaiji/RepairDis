package com.bishe.ddr_springboot.repository;

import com.bishe.ddr_springboot.entity.Gene;
import org.springframework.data.jpa.repository.JpaRepository;

public interface GeneRepository extends JpaRepository<Gene, Integer> {
}