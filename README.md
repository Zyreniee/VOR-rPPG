# VOR: Variance-Optimized rPPG

[![DOI](https://zenodo.org/badge/DOI/FILL_AFTER_ZENODO.svg)](https://doi.org/10.5281/zenodo.FILL_AFTER_ZENODO)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

VOR (Variance-Optimized rPPG) is a novel multi-modal framework for robust remote photoplethysmography. It fuses multiple extraction algorithms (CHROM, POS, Green) across multiple regions of interest (ROIs) through a risk-aware spectral clustering engine.

## Core Features
- **Sigmoid Confidence Weighting**: Mathematically penalizes spectral artifacts and "harmonic traps".
- **Spectral Clustering**: Greedy fusion of multi-ROI candidates to establish consensus.
- **Temporal Trust Engine**: Stability analysis utilizing historical memory and method agreement logic.

## Academic Paper
The formal mathematical and algorithmic specification of the VOR engine is available in the `paper/` directory of this repository. 
- **Title:** VOR: A Variance-Optimized Multi-Modal Framework for Robust Remote Photoplethysmography
- **Author:** Yusuf Gunes (*Independent Researcher*)
- **ORCID:** [0009-0003-0173-9536](https://orcid.org/0009-0003-0173-9536)

## Implementation Details
This repository contains the C++20 reference implementation of the VOR component, specifically targeting deterministic, high-performance extraction for medical and telehealth applications. It is dependency-free by design to ensure broad portability.

*Note: This repository specifically contains the VOR extraction engine. The complete medical-grade triage infrastructure, which integrates VOR within a Vulkan-accelerated Heart-on-a-chip pipeline, remains part of the primary proprietary system.*

## Citation
If you use this algorithm in your research or product, please cite the Zenodo DOI provided above.
