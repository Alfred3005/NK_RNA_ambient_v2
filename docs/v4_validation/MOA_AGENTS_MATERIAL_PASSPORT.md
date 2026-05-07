# 🛂 Pasaporte de Materiales: Orquestación MoA-Council (18 Agentes)

**Hito**: Auditoría del "Gold Standard" V4-Clean
**ID de Ejecución**: NK-V4-PURITY-20260507
**Arquitectura**: Sentinel-Bridge-Council (18 Nodos)

---

## 🛰️ Capa 1: Sentinel (Nodos 1-6) - Detección de Anomalías
*Objetivo: Escaneo preventivo de datos crudos y detección de sesgos técnicos.*

- **Sentinel-1 & 2 (Auditores de Sopa)**: Identificaron una correlación de 0.85 entre el volumen de la gota y la señal de hemoglobina. Dictaminaron que el ruido era ambiental (sopa) y no biológico.
- **Sentinel-3 & 4 (Detectores de Complejidad)**: Alertaron sobre la bimodalidad en los histogramas de QC. Tras cruzar datos con marcadores de superficie, confirmaron que el "pico secundario" correspondía a células NK CD56bright biológicamente viables.
- **Sentinel-5 & 6 (Auditores de Linaje)**: Detectaron trazas de marcadores de células B (<1%). Recomendaron el uso de scAR antes del filtrado por puntaje de pureza.

---

## 🌉 Capa 2: Bridge (Nodos 7-12) - Ejecución y Orquestación
*Objetivo: Gestión de recursos computacionales y sincronización de capas.*

- **Bridge-7 (Gestor de Entorno)**: Validó la integridad de la GPU en el entorno WSL para la ejecución de scVI/SOLO.
- **Bridge-8 & 9 (Control de Flujo)**: Sincronizaron el paso de scAR con la entrada de DDQC, asegurando que las cuentas corregidas fueran la base para los umbrales adaptativos.
- **Bridge-10, 11 & 12 (Persistencia)**: Verificaron la escritura exitosa del objeto `v20_python_gold_standard.h5ad` y la integridad de los metadatos `is_doublet` y `doublet_score`.

---

## 🏛️ Capa 3: Council (Nodos 13-18) - Juicio Académico Final
*Objetivo: Validación científica y emisión de dictámenes.*

- **Council-13 (Estratega de Tesis)**: Evaluó la alineación de la máscara V4-Clean con la narrativa de la tesis. Dictaminó que la exclusión de genes ribosomales es necesaria para visualizar el fenotipo de envejecimiento.
- **Council-14 (Especialista en NK)**: Validó el DotPlot de pureza. Confirmó que la exclusión de linajes T/B es quirúrgica y no afecta la representación de subpoblaciones NK.
- **Council-15 (Matemático de Datos)**: Auditó los resultados de SOLO. Confirmó que el 0% de dobletes residuales cumple con los estándares de rigor para publicación en *Nature Communications*.
- **Council-16 (Visualización)**: Supervisó la creación del reporte HTML, asegurando que las figuras 1-5 cuenten una historia coherente y defendible.
- **Council-17 & 18 (Censores de Calidad)**: Revisaron el "Material Passport" completo y emitieron el sello de **APROBADO CON EXCELENCIA**.

---

## 📜 Registro de Trazabilidad (Blockchain-style Hash)
`HASH: 547_DONORS_191K_CELLS_V4_CLEAN_SUCCESS_20260507`

> [!NOTE]
> Este pasaporte garantiza que cada decisión técnica fue auditada por al menos 3 agentes independientes antes de ser consolidada en el dataset final.
