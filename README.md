# VOR-rPPG Core

Variance-Optimized rPPG (VOR) is a high-performance C++20 engine for multi-modal signal fusion in remote photoplethysmography. It implements risk-aware spectral clustering to mitigate motion artifacts and harmonic interference.

## Features
- **Multi-Modal Fusion**: Greedy clustering of candidates from POS, CHROM, and Green extractors.
- **Risk Annotation**: Sigmoid-based confidence weighting $(\phi)$ penalizing boundary and harmonic risks.
- **Temporal Trust**: Stability state machine for robust heart rate estimation in dynamic environments.

## Build Requirements
- C++20 compatible compiler (GCC 11+, Clang 13+, MSVC 19.30+)
- CMake 3.15+

## Building the Library
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Example Usage
A minimal benchmark demonstrating how to process multi-ROI RGB signals:
```bash
./vor_benchmark
```
Refer to `examples/benchmark.cpp` for implementation details.

## License
MIT License - Copyright (c) 2026 Yusuf Gön