# VOR-rPPG: Variance-Optimized Multi-Modal Signal Fusion for Robust Remote Photoplethysmography

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19312766.svg)](https://doi.org/10.5281/zenodo.19312766)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**VOR-rPPG** is an ultra-high-performance C++20 research framework designed for medical-grade, non-contact pulse estimation from standard RGB video sequences.

## Core Innovations

- **Multi-Modal Sigma-Fusion**: Synthesizes pulse candidates from POS (Plane-Orthogonal-to-Skin), CHROM (Chrominance), and Green-channel extractors using a greedy spectral clustering strategy.
- **Risk-Aware Sigmoid Weighting ($\phi$)**: A novel confidence-scoring system that dynamically penalizes signals near physiological boundaries and those exhibiting harmonic interference.
- **Temporal Stability Engine**: A state-machine-driven trust mechanism ($T$) that suppresses jitter and ensures longitudinal consistency.
- **Zero-Copy Architecture**: Leveraging `std::span`, `std::atomic`, and modern C++20 primitives for <0.2ms processing latency per ROI.

## Build Requirements
- C++20 compatible compiler (GCC 11+, Clang 13+, MSVC 19.30+)
- CMake 3.15+

## Quick Start
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./vor_benchmark
```

## Citation
If you use this framework in your research,please cite the official preprint:

```bibtex
@software{vor_rppg_2026,
  author       = {Gön, Yusuf},
  title        = {VOR-rPPG: Variance-Optimized Remote Photoplethysmography Core},
  institution  = {Eskisehir Sabiha Gokcen Mesleki ve Teknik Anadolu Lisesi},
  month        = apr,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.19312766},
  url          = {https://doi.org/10.5281/zenodo.19312766}
}
```

## License
- Software: MIT License
- Publication/Data: CC-BY 4.0
Copyright (c) 2026 Yusuf Gön
