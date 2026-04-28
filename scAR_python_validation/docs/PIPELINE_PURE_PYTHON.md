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
- [x] Differential Expression (Adult vs Old) - Pseudobulk PyDESeq2.
- [ ] Functional Enrichment (GO/KEGG).

---

## 🔬 Analysis Phase: DE "Adult vs Old" (2026-04-28)

We performed a comparative DE analysis using PyDESeq2 on the Gold Standard dataset (547 donors).

### Results:
- **Significance**: 232 genes found (padj < 0.05, |LFC| > 1).
- **Venn Contrast (Original vs V20)**:
  - Original Signature: 83 genes.
  - Intersection (15 genes): `SERPINA1`, `CST3`, `LST1`, `FAM131B-AS2`, `DEGS2`, `AIF1`, `LINC00513`, `SPON1`, `ANGPT2`, `ANKRD20A4P`, `JAKMIP1`, `ST3GAL1-DT`, `SIDT1-AS1`, `DUOX1`, `HBA2`.
  - Purified Out: 68 genes (Noise/Contamination removed).
  - Emergent Hits: 217 genes (Increased sensitivity).
- **Outcome**: Successfully rescued the NK identity and identified a robust inflammatory aging signature.
