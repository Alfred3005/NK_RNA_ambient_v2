# 🧬 NK Aging Pipeline: The V20 Watershed (Era Fénix)

Este repositorio contiene la arquitectura bioinformática avanzada para el rescate y análisis de transcriptómica de células NK humanas, enfocada en el estudio del envejecimiento celular.

## 🚀 El Hito V20 (Abril 2026)
Tras enfrentar una crisis de integridad de datos por archivos masivos (~80GB) y contaminación de RNA ambiental, implementamos el **Protocolo Fénix**. Este proceso rescató ~130,000 células NK con un estándar de pureza genética y funcional sin precedentes.

### 🏁 Estado de la Misión
- **Dataset Maestro**: `V20_CLEAN_ANALYSIS/data/nk_v20_singlets.h5ad`
- **Volumen**: 196,091 células NK purificadas.
*   **Resultados Clave**: Hemos demostrado que el 84% de la firma de envejecimiento previa era ruido técnico, redefiniendo la narrativa hacia un modelo de **Anergia Intínseca** y **Estrés Oxidativo**.

---

## 🏛️ Santo y Seña (Estructura del Proyecto)

- `V20_CLEAN_ANALYSIS/`: **Directorio Activo**. Contiene todo el flujo de trabajo purificado.
    - `scripts/`: Pipeline numerado de 01 a 10 (scvi, SOLO, ddqc, PyDESeq2).
    - `docs/memory_logs/`: **Bitácora Maestra**. Reportes de validación y memos de descubrimiento.
- `legacy_scripts/`: Archivo histórico de la "Era Monster" (herramientas de rescate de fragmentos).
- `results/`: Figuras filtradas y tablas de DE finales.

---

## 🤖 AI Onboarding (Para Claude Code / Antigravity)
Si eres una IA encargada de auditar o continuar este proyecto, lee los siguientes documentos en este orden:
1. `V20_CLEAN_ANALYSIS/docs/memory_logs/WALKTHROUGH_MASTER.md`: Resumen técnico del pipeline.
2. `V20_CLEAN_ANALYSIS/docs/memory_logs/MEMO_DESCUBRIMIENTO.md`: La nueva base científica.
3. `V20_CLEAN_ANALYSIS/scripts/MASTER_SEQUENCE.md`: Lógica de ejecución de los scripts.

---

## 🛠️ Requerimientos resumidos
- `scanpy`, `scvi-py`, `pydeseq2`, `anndata`.
- Ambiente recomendado: `.v20_venv` (Python 3.10+).

**Investigador Principal:** Alfred3005  
**Asistente de IA:** Antigravity (Advanced Agentic Coding Agent)
