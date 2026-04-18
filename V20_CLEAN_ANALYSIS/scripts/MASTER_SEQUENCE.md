# 📜 Master Sequence: V20 Pipeline Flow

Esta es la secuencia de ejecución lógica para reproducir el análisis de purificación NK.

| Orden | Script | Función Principal |
| :--- | :--- | :--- |
| **01** | `01-diagnostic-report.py` | Auditoría inicial del estado de los archivos h5ad. |
| **04** | `04-adaptive-qc.py` | Implementación de `ddqc`. Determina umbrales de cortes por cluster. |
| **05** | `05-preprocessing-metadata.py` | Armonización de grupos de edad y metadatos del donante. |
| **06** | `06-doublet-removal-solo.py` | Entrenamiento de scvi y ejecución de `SOLO` para limpiar agregados. |
| **07** | `07-evaluate-ambient-rna.py` | Generación de visualizaciones de contaminación ambiental (Violines). |
| **08** | `08-explore-legacy.py` | Exploración del dataset original de 34GB (Modo backed 'r'). |
| **09** | `09-comparative-de-analysis.py` | Comparativa Wilcoxon entre Legacy vs V20. |
| **10** | `10-pseudobulk-pydeseq2.py` | **Motor Estadístico Final**. Análisis PyDESeq2 por donante. |

## 🛠️ Scripts de Utilidad
- `append_genes_flat.py`: Herramienta para generar tablas de genes en formato Markdown para reportes.
- `nb_dump.txt`: Código extraído de los notebooks originales para referencia rápida.

---
*Nota: Los scripts 02 y 03 fueron etapas intermedias de diagnóstico durante el proceso de rescate.*
