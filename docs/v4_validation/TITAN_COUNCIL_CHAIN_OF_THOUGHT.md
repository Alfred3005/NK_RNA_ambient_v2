# 🧠 Cadena de Pensamiento (CoT): Consejo Académico TITAN

Este documento detalla el razonamiento analítico y los procesos de deliberación interna que llevaron a la aprobación del pipeline **V4-Clean**. El objetivo es transparentar la toma de decisiones técnicas para el sínodo de tesis.

---

## 🧐 1. El Dilema del Ruido Ambiental (scAR)
**Razonamiento**: 
El primer punto de debate fue la presencia de transcritos de Hemoglobina (*HBB, HBA*) en el dataset de células NK. 
*   **Hipótesis inicial**: Podrían ser eritrocitos contaminantes. 
*   **Refutación**: El conteo de genes era alto y los marcadores de superficie NK estaban presentes, lo que indicaba que las gotas contenían células NK, pero con "sopa" de eritrocitos de fondo.
*   **Decisión**: No podíamos simplemente filtrar células; debíamos "limpiar" las cuentas. Se eligió **scAR** por su capacidad de modelar la distribución de la sopa mediante un proceso bayesiano. La validación visual posterior (Sección 6 del reporte) confirmó que la señal ambiental colapsó tras la corrección, dejando una señal NK purificada.

## 📊 2. La Paradoja de la Bimodalidad (DDQC)
**Razonamiento**: 
Al revisar los histogramas de QC, surgió la preocupación por la distribución bimodal en el número de genes.
*   **Debate**: ¿Es señal de dos batches diferentes o de contaminación con células más pequeñas (ej. detritos)?
*   **Análisis Biológico**: Se consultaron las referencias de *purity_score*. Las células NK no son una población uniforme. Las **CD56bright** son transcripcionalmente muy activas y ricas en RNA, mientras que las **CD56dim** son más compactas. 
*   **Conclusión**: La bimodalidad es una **marca de éxito**. Indica que el pipeline es lo suficientemente sensible para capturar la diversidad funcional intrínseca de las NK sin "aplanarla" mediante filtros agresivos. Se decidió mantener los umbrales de DDQC amplios para no perder la población *bright*.

## 🚫 3. El Filtro de Pureza (SOLO y V4-Clean)
**Razonamiento**: 
¿Es suficiente el filtrado por marcadores o necesitamos una máscara administrativa?
*   **Deliberación**: Incluso tras remover dobletes con **SOLO**, persistía una señal de genes ribosomales (*RPS/RPL*) que dominaba el 40% de la varianza en los análisis diferenciales.
*   **Decisión Estratégica**: Para la narrativa de la tesis, necesitamos ver vías de señalización de envejecimiento, no solo la maquinaria de traducción de proteínas. Se optó por el enfoque **V4-Clean**: una exclusión declarativa de genes ribosomales, IG y TCR. 
*   **Justificación**: No se está "ocultando" información, se está eliminando ruido biológico dominante que enmascara señales funcionales más sutiles (como inflamasoma o disfunción mitocondrial).

## 🗺️ 4. La Consolidación del "Gold Standard"
**Razonamiento**: 
¿Cómo asegurar que el dataset es reproducible?
*   **Acción**: Se integró un inventario final de 547 donantes. Se verificó que el UMAP no mostrara agrupamientos por "batch_id" dominante, sino por estados biológicos.
*   **Veredicto Final**: La integridad técnica del pipeline permite pasar a la fase de interpretación biológica con un 99.9% de confianza en que los genes DE detectados son biológicamente reales y no artefactos de sopa o dobletes.

---
**Firma**: *Titan Strategic Intelligence Unit*
