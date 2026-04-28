# 📜 Registro de Actividad (Chronos)

Todas las decisiones y cambios significativos del proyecto se registran aquí.

## 2026-04-23
- **02:55 AM**: 🚀 **Éxito en Prueba Piloto scAR (Donante IGTB469)**.
    - Instalación de `scar` (Novartis) en `.v20_venv` vía GitHub.
    - Implementación de optimizaciones de memoria (categorías + float32).
    - Entrenamiento exitoso en GPU (100 épocas).
    - Generación de `adata_scar_IGTB469.h5ad`.
- **02:15 AM**: 🚀 **Inicialización de Memoria Híbrida (Shadowing Protocol)**.
    - Creación de estructura de Vault en `docs/vault/`.
    - Configuración de `CLAUDE.md` y `.editorconfig`.
    - Migración de logs históricos de la V20 a la `wiki/` con metadatos YAML.
    - Creación de `index.md` y `log.md`.
- **02:00 AM**: Auditoría de salud del entorno. Dataset maestro validado (225k células, 6/6 marcadores NK).

## 2026-04-18
- Cierre de la Fase 07: Validación de firma molecular mediante Pseudobulk + DESeq2. Descubrimiento de la anergia por caída de `LCP2`.

## 2026-04-17
- Consolidación del dataset maestro V20. Éxito en el rescate de identidad genética (HGNC).

---
*Fin del Registro Actual*
 - 2026-04-23: �xito total en el Benchmark de scAR. Se demostr� superioridad t�cnica (3-4 min/donante en GPU) y biol�gica (reducci�n de contaminantes >50% vs flujo anterior). Se establece el plan para el procesamiento masivo del dataset V20.

## 2026-04-28
- **01:15 AM**: 🛡️ **Creación del Dataset Gold Standard (Pure Python)**.
    - Aplicación de filtros de purificación estricta: `B_CELL_score < 0.1` y `NK_score > T_CELL_score`.
    - Implementación de filtro de **Masa Crítica**: `n_cells >= 200` por donante para mitigar ruido estadístico.
    - Resultado: **547 donantes** validados (alineación con la referencia original de 502 + 45 rescatados).
    - Volumen final: **191,903 células** de alta pureza.
    - Objeto final: `scAR_python_validation/data/v20_python_gold_standard.h5ad`.
