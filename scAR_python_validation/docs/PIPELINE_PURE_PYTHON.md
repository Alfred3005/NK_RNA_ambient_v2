# 🛡️ Pipeline Memory: scAR Python Validation (V20-Native)

## Metadata
- **Objective**: Establish a high-fidelity, pure Python pipeline for NK transcriptomics.
- **Starting Point**: 1,582 donor-level H5AD files from `data/processed/scar_denoised/`.
- **Status**: Initialization Phase.

---

## 📜 Decision Log

### [2026-04-27] - Pipeline Isolation
- **Decision**: Created `scAR_python_validation/` to separate from `V20_CLEAN_ANALYSIS` (legacy/hybrid).
- **Rationale**: Prevent cross-contamination of metadata and gene signatures between the R-based and Python-native flows.

---

## 🧪 Data Lineage
1. **Raw Source**: `Datasets/131224_full_dataset.h5ad` (80GB).
2. **Denoising**: scAR (massive orchestrator run) -> `data/processed/scar_denoised/*.h5ad`.
3. **Current State**: Ready for consolidation and QC.

---

## 🚧 Progress Tracker
- [x] Directory structure setup.
- [/] Initializing memory log.
- [ ] Consolidation of 1,582 donor files.
- [ ] Adaptive QC (ddqc).
