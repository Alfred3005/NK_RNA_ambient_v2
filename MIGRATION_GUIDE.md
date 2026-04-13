# GUÍA DE MIGRACIÓN: NK Thesis Phoenix 🧬🔥

Sigue estos pasos para mover los archivos de este "Checkpoint" a tu nuevo entorno de trabajo.

## 1. Mapeo de Carpetas
Copia los archivos de `volatile-photon/migration_checkpoint/` a `NK_pipeline_RNA_ambient/` siguiendo este mapa:

| Archivo de Origen | Destino en Proyecto | Función |
| :--- | :--- | :--- |
| `docs/core_docs.md` | `docs/` | Diseño y Protocolo |
| `scripts/01-segmentation-rapids.py` | `scripts/` | Segmentación 80GB (Plan A/B) |
| `scripts/02-ambient-correction.R` | `scripts/` | Corrección scCDC |
| `scripts/03-custom-labeling.py` | `scripts/` | Etiquetado Edad y HVG |
| `config_files.md` | `./` (Raíz) | .gitignore y README |

## 2. Checklist de Despegue
- [ ] Mover el dataset de 80GB a `data/raw/131224_full_dataset.h5ad`.
- [ ] Ejecutar `pip install cudf scanpy psutil` en tu entorno virtual.
- [ ] Iniciar la nueva sesión de Antigravity apuntando a la nueva carpeta.

---
> [!NOTE]
> **Plan B Activado**: El script `01` ahora incluye la opción `--by-donor` para procesar por donador en caso de que un dataset sea demasiado grande para tu RAM.
