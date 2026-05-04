# Walkthrough: Implementación y Validación de V3-Perfect

En esta sesión, elevamos el rigor estadístico del pipeline NK para abordar las controversias de falsos positivos en grandes cohortes.

## 1. Configuración del Entorno de Ejecución
Se optimizó el entorno WSL para permitir la ejecución de `PyDESeq2` sin depender de Conda, evitando conflictos con drivers NVIDIA:
- Instalación de dependencias en el entorno virtual del proyecto (`.venv_wsl`).
- Resolución de shebangs y rutas absolutas heredadas de migraciones de directorios anteriores.

## 2. Desarrollo del Pipeline "Perfect"
Se creó el script [10-pydeseq2-pseudobulk-perfect.py](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v3_perfect/scripts/10-pydeseq2-pseudobulk-perfect.py) con las siguientes mejoras:
- **LFC Shrinkage**: Implementación de `stat_res.lfc_shrink(coeff="age_group[T.old]")` para corregir la inflación de Fold Changes en genes de bajo conteo.
- **Purity Enforcement**: Mantenimiento de los filtros declarativos para genes IG y TCR.

## 3. Depuración y Correcciones Técnicas
Durante la ejecución, se resolvieron los siguientes errores críticos de la API de PyDESeq2:
- Corrección del argumento `coeff` (anteriormente `coef`).
- Identificación del string de coeficiente correcto generado por la fórmula interna (`age_group[T.old]`).

## 4. Auditoría de Resultados (V2 vs V3)
Se ejecutó un script de auditoría [compare_v2_v3.py](file:///c:/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v3_perfect/scripts/compare_v2_v3.py) para cuantificar el impacto:
- **Reducción de Ruido**: Se descartaron 47 genes que eran "artefactos" de baja expresión base.
- **Robustez Biológica**: Se validó que los **45 genes ribosomales** clave siguen siendo significativos tras la contracción bayesiana.

## 5. Documentación Final
Se generó un informe técnico detallado [REPORT_V3_PERFECT_DE_ANALYSIS.md](file:///C:/Users/PREDATOR/.gemini/antigravity/brain/1994367d-5869-4db8-b857-9d7f561b8110/REPORT_V3_PERFECT_DE_ANALYSIS.md) que explica matemáticamente por qué estos 182 genes constituyen la firma de oro definitiva para la tesis.

---

### Siguiente Paso
Proceder al **Análisis de Enriquecimiento Funcional (GSEA)** con la firma de 182 genes validada.
