# VOR: Variance-Optimized rPPG Specification

This document provides a formal mathematical and algorithmic specification for the VOR (Variance-Optimized rPPG) engine, designed to facilitate open-source implementations and clinical validation.

## 1. Core Logic
VOR treats remote photoplethysmography as a multi-modal consensus problem. It extracts candidates from multiple sources (Green, CHROM, POS) and regions of interest (ROIs), then filters them through a risk-aware trust engine.

## 2. Mathematical Formalism

### 2.1 Candidate Confidence ($\phi$)
For each spectral peak $f_i$, the confidence $\phi(f_i)$ is computed as:
$$ \phi(f_i) = \sigma\left(\frac{SNR_i - 3}{2}\right) \times (1 - \rho_b) \times (1 - \rho_h) \times \frac{1}{\sqrt{rank}} $$

Where:
- **$\sigma(z)$**: Sigmoid function $(1 + e^{-z})^{-1}$.
- **$SNR_i$**: Signal-to-noise ratio in dB.
- **$\rho_b$ (Boundary Risk)**: $ \max(0, 1 - \frac{|f_i - f_{center}|}{\Delta_{band}}) $, where $f_{center}$ is the center of the physiological band.
- **$\rho_h$ (Harmonic Risk)**: Probability that $f_i$ is a harmonic ($2k \times$) or sub-harmonic of another dominant peak.

### 2.2 Temporal Trust Engine
The final trust score $C_{cal}$ is recalibrated every cycle:
$$ C_{cal} = C_{raw} \times \sqrt{T} \times \sqrt{M} \times R $$

- **$T$ (Temporal Consistency)**: $\exp(-\sigma_{history}^2 / \sigma_{limit}^2)$, measuring stability over the last $N$ windows.
- **$M$ (Method Agreement)**: Variance between top-1 candidates of different extractors (POS vs. CHROM vs. Green).
- **$R$ (Resolution Quality)**: Ratio of the segment length to the minimum required for target frequency resolution.

## 3. Implementation Guidelines

### 3.1 Pre-processing
- **Resampling**: All ROI signals must be resampled to a uniform sampling frequency ($f_s$) before spectral analysis.
- **Detrending**: Use a linear or polynomial detrend (e.g., Tarvainen) to remove baseline illumination drift.

### 3.2 Fusion (Greedy Clustering)
Join all candidates $(f, \phi)$ in a pool and cluster them with a radius of 6 BPM.
Select the cluster with the highest aggregate score:
$$ Score(K) = \left( \sum_{j \in K} \phi_j \right) \cdot \ln(1 + |K|) \cdot \text{Tightness} $$

## 4. Open Source & License
This specification and the reference implementation in `VOREngine.cpp` are released under the **MIT License**.
