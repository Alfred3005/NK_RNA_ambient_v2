# Plan de Implementación: Árbol de Decisiones Técnicas (Tubería Bioinformática NK)

Este documento detalla la estructura y el contenido que tendrá el informe final sobre la evolución metodológica del proyecto `NK_pipeline_RNA_ambient`. El objetivo es documentar las rutas, justificaciones y variaciones en la "n" (células/donantes y genes resultantes) a lo largo del proceso.

## User Review Required

> [!IMPORTANT]
> Revisa si los números y las justificaciones plasmadas en la sección de "Cambios Propuestos" reflejan fielmente la historia clínica de los datos según tu memoria del proyecto. 
> Especialmente la justificación del paso V3 -> V4 (decisión administrativa para el enriquecimiento de vías).

## Open Questions

> [!WARNING]
> En la documentación de V3-Perfect se mencionan **182 genes** significativos, pero en el resumen de V4 se cita que V3 tenía **190 genes**. ¿Deseas que unifique esta cifra a 182, 190, o prefieres que mencione un rango?

## Proposed Changes

Crearemos un artefacto Markdown detallado (`TECHNICAL_DECISION_TREE.md`) que contendrá:

### 1. Contexto y Justificación del Proyecto
Una sección introductoria sobre la necesidad de rescatar un dataset masivo (~80GB, ~427k células) de células NK y purificarlo de RNA ambiental y dobletes.

### 2. Diagrama del Árbol de Decisiones (Mermaid)
Un flujograma visual `mermaid` que mostrará las siguientes bifurcaciones clave:
- **Raw Data** -> Decisión: Herramientas en R vs Nativo en Python.
- **Ambient RNA Removal** -> Decisión: A nivel de dataset vs A nivel de Donante (RAM constraints).
- **V2-Clean** -> Pseudobulk estándar vs Filtrado Declarativo.
- **V3-Perfect** -> MLE estándar vs LFC Shrinkage (Controversia de Falsos Positivos).
- **V4-Clean** -> Genes Totales vs Filtrado Ribosomal (RPS/RPL) para Enriquecimiento.

### 3. Cronología y Justificación Táctica (Niveles Táctico y Refinamiento)

Se redactará una narrativa cronológica que explicará los nodos del árbol:

*   **Punto de Partida (Era Monster/V20 Watershed)**
    *   **N**: ~427k células brutas.
    *   **Decisión**: Migrar a una tubería 100% nativa en Python.
    *   **Justificación**: Evitar la fricción de conversión entre objetos `AnnData` y `Seurat`. El tamaño monstruoso de la matriz hacía inviable el procesamiento en R.
*   **Purificación de RNA Ambiental (scAR)**
    *   **Decisión**: Aplicar la corrección a nivel de donante en lugar de hacerlo a nivel global del dataset.
    *   **Justificación**: Limitaciones de hardware (64GB RAM) y mejor resolución biológica.
    *   **N Resultante**: "Gold Standard" de 196,091 células NK purificadas.
*   **V2-Clean (El Primer Filtrado Declarativo)**
    *   **N**: 191,903 células de 547 donantes en Pseudobulk (PyDESeq2).
    *   **Decisión**: Excluir genes de Inmunoglobulinas (IGH/IGK/IGL) y TCRs.
    *   **Justificación**: Eliminar el ruido técnico preventivamente.
    *   **N Resultante (Firma)**: 229 genes significativos ($padj < 0.05, |LFC| > 1$).
*   **V3-Perfect (La Corrección del "Efecto Espejismo")**
    *   **Decisión**: Implementar LFC Shrinkage (estimación Bayesiana).
    *   **Justificación**: La alta "n" de donantes (547) inflaba el poder estadístico y generaba falsos positivos en genes con conteos bajísimos (low counts). El Shrinkage colapsó esos ruidos a cero.
    *   **N Resultante (Firma)**: ~182 genes significativos.
*   **V4-Clean (Alineamiento con Tesis y Enriquecimiento)**
    *   **Decisión**: Filtrar exhaustivamente 45 genes ribosomales (RPS/RPL).
    *   **Justificación**: Decisión administrativa y de dirección de Tesis. Los genes ribosomales mostraban una subexpresión masiva y robusta, pero sesgarían los análisis posteriores de enriquecimiento (GSEA/ORA) ocultando otros mecanismos del envejecimiento.
    *   **N Resultante (Firma)**: 140 genes de alta confianza.

## Verification Plan

### Manual Verification
- Te presentaré el documento renderizado en los artefactos de la conversación para que lo puedas revisar visualmente y exportar a tus notas.
- Podrás validar que el diagrama de decisiones (Mermaid) sea claro y estéticamente profesional, para que puedas integrarlo fácilmente en tu redacción de tesis o repositorios.
