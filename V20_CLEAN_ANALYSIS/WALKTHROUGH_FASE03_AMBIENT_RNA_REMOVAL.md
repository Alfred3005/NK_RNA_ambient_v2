# 🔬 Fase 03: Eliminación de RNA Ambiental y Rescate V20

Este flujo técnico documenta el proceso de corrección de sopa de RNA (ambient RNA) mediante el algoritmo **scCDC** (Single-cell Contamination Detector and Cleaner), crucial para asegurar que la expresión diferencial posterior sea real y no producto de ruido de fondo.

## 📊 1. Resumen de Poder Estadístico
En esta ejecución de alta fidelidad (V20), hemos asegurado el **80.3% del total** de los fragmentos (1,372 de 1,708), lo que constituye el dataset más robusto hasta la fecha.

- **Células Totales (Aprox)**: ~130,000 - 150,000 células NK validadas post-corrección inicial.
- **Identidad Genética**: Se resolvió el problema de pérdida de nombres de genes integrando los símbolos HGNC desde la fase de carga.

## 🛠️ 2. Estabilización Técnica (The scCDC Fix)
El proceso fue estabilizado tras fallos matemáticos previos mediante:
1. **Bypass de Clustering De-Novo**: Se utilizaron las etiquetas `cell_type` pre-existentes para informar al algoritmo scCDC, evitando que el kernel colapsara por falta de memoria en clusters masivos.
2. **Procesamiento Distribuido**: Procesamiento de los 1,714 segmentos individuales con reconsolidación integrada.

## 🧼 3. Resultados de Limpieza
La aplicación de scCDC permitió:
- **Reducción Estocástica**: Limpieza de transcritos T/B presentes en gotas vacías o ambientales.
- **Preservación de Señal**: Se mantuvo la varianza biológica necesaria para sub-segmentar poblaciones NK en adultos y ancianos.

## 🏁 Estado de la Fase
Esta fase concluyó exitosamente, habilitando el pipeline para el **Control de Calidad Adaptativo (Phase 04)** y la **Remoción de Dobletes (Phase 05)**.

---
*Archivo de registro técnico - NK Pipeline V20.*
