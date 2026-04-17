# 🕰️ CHRONOS MEMORY: The Clean V20 Transition

Este documento sirve como ancla de memoria para el asistente Antigravity y el Usuario en el nuevo espacio de trabajo.

## 🏺 Legacy (The "Monster" Era)
- **Periodo**: Marzo - Abril 2026.
- **Error Crítico**: Pérdida de identidad génica (índices numéricos en lugar de HGNC) durante la integración de archivos h5ad de ~80GB.
- **Resultado**: El dataset maestro original de 161,460 células quedó inutilizable para análisis biológico.

## 🚀 The Rescue (Fénix V20)
- **Hito**: 17-Abr-2026.
- **Solución**: Implementación de mapeo forzado `feature_name` -> `index` durante el rescate de fragmentos.
- **Checkpoint Maestro**: `data/nk_v20_master.h5ad`. Una consolidación de ~130,000 células con 100% de genes HGNC válidos.

## 🎯 Vision (The Clean Path)
- **Enfoque**: Abandono total de datasets "monstruos" no procesados.
- **Prioridad**: Calidad y estadística robusta sobre cantidad absoluta.
- **Siguiente Paso**: Fase 04 - Control de Calidad Adaptativo (ddqc).

---
*Este documento protege el contexto de la Tesis NK contra la entropía de los cambios de carpeta.*
