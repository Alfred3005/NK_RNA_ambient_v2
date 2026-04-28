# 🛡️ Pipeline Memory: scAR Python Validation (V20-Native)

## Metadata
- **Objective**: Establish a high-fidelity, pure Python pipeline for NK transcriptomics.
- **Starting Point**: 1,582 donor-level H5AD files from `data/processed/scar_denoised/`.
- **Status**: Analysis Phase (Gold Standard Created).

---

## 📜 Decision Log

- **Decision**: Created `scAR_python_validation/` to separate from `V20_CLEAN_ANALYSIS` (legacy/hybrid).
- **Rationale**: Prevent cross-contamination of metadata and gene signatures between the R-based and Python-native flows.

### [2026-04-28] - Gold Standard Establishment
- **Decision**: Applied strict lineage purification (`NK_score > T_CELL_score` & `B_CELL_score < 0.1`) and mass critical filtering (`n_cells >= 200` per donor).
- **Rationale**: Mitigate noise from small donor fragments and heterotypic contamination, aligning with the reference workflow's donor count (502) while preserving the "rescue" benefits of the Python flow (547 donors final).

---

## 🧪 Data Lineage
1. **Raw Source**: `Datasets/131224_full_dataset.h5ad` (80GB).
2. **Denoising**: scAR (massive orchestrator run) -> `data/processed/scar_denoised/*.h5ad`.
3. **Current State**: Gold Standard created with 204,478 cells across 547 donors.

---

## 🚧 Progress Tracker
- [x] Directory structure setup.
- [x] Initializing memory log.
- [x] Consolidation of 1,582 donor files.
- [x] Adaptive QC (ddqc).
- [x] Gold Standard Purification (200 cell threshold).
- [ ] Differential Expression (Adult vs Old).
