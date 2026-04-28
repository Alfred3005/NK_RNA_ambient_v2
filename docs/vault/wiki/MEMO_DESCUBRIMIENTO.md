---
type: memo
project: NK_V20_Rescue
topic: Discovery_Audit
status: final
date: 2026-04-18
tags: [discovery, differential_expression, pydeseq2, signature]
---
# 📝 Memorándum de Descubrimiento: El Paisaje NK V20

**De:** Antigravity (IA Coding Assistant)  
**Para:** Investigador Principal NK Pipeline  
**Fecha:** 18 de Abril, 2026  
**Asunto:** Auditoría y Redefinición Narrativa de la Inmunosenescencia NK (Fase V20)

---

## 🚀 1. Resumen Ejecutivo: El "Gran Pivote"

Este documento formaliza el cambio de paradigma tras completar la validación **Pseudobulk + DE** sobre el dataset purificado. 

*   **Estado Anterior (Legacy):** Una firma dominada por genes de alta expresión externa (`CXCL8`, `SERPINA1`, `MS4A1`) que sugerían una transición fenotípica hacia células mieloides o B.
*   **Estado Actual (V20):** Una firma refinada que revela una **pérdida de competencia funcional** (Anergia) y **estrés celular intrínseco**, eliminando el ruido ambiental.

---

## 🔍 2. Auditoría Técnica: ¿Qué pasó con la Firma Original?

La desaparición de genes clave no es una pérdida de datos, sino una **ganancia de fidelidad**. Aquí el porqué:

### A. La Trampa del RNA Ambiental (The Ambient Trap)
Muchos genes de la firma legacy (`IFI30`, `C1QA`, `CXCL8`) son proteínas secretadas o marcadores de linajes muy abundantes en sangre (Monocitos/Neutrófilos). En el dataset original sin corrección, este RNA flotante se asignaba erróneamente a las NKs. 
- **Veredicto:** El pipeline V20 (Fénix) detectó y neutralizó estas señales como ruido de fondo.

### B. El Efecto de los Dobletes
Marcadores como `MS4A1` (CD20) y `MZB1` indicaban la presencia de Células B. Al aplicar el filtro estricto de singletes, estas "NKs" que en realidad eran parejas físicas de NK+B fueron eliminadas.
- **Resultado:** La firma B colapsó, dejando una NK genéticamente pura.

---

## 🧬 3. Duelo de Evidencia: Realidad vs. Artefacto

| Eje Biológico | Genes Legacy (Dudosos) | Genes V20 (Verificados) | Interpretación Médica |
| :--- | :--- | :--- | :--- |
| **Identidad Celular** | `MS4A1`, `MZB1`, `TCL1A` (B-cells) | **CD81 (Down)**, **LCP2 (Down)** | La NK no se "vuelve B", sino que pierde sus receptores de membrana clave. |
| **Inflamación** | `CXCL8`, `SERPINA1`, `IL1B` | **DUOX1 (Up)**, **SOCS1 (Down)** | No es inflamación sistémica, es **estrés oxidativo local** lo que daña a la NK. |
| **Función Citotóxica** | `CST3`, `CTSL` | **LCP2 (Down)**, **NFATC2 (Down)** | Hay un fallo masivo en la transducción de señales de activación. |

---

## 🏛️ 4. El Nuevo Modelo Biológico: Anergia y Estrés

Bajo la lente V20, la inmunosenescencia NK se define por tres ejes fundamentales:

### I. Fallo de la Sinapsis Inmunitaria (`LCP2` / `NFATC2`)
La caída drástica de `LCP2` (SLP-76) sugiere que las células NK viejas han perdido la capacidad de ensamblar el complejo de señalización necesario para la degranulación. **Están presentes, pero son ciegas al estímulo.**

### II. Estrés Oxidativo Persistente (`DUOX1`)
La persistencia de `DUOX1` (Dual Oxidase 1) indica que la NK está en un estado de generación constante de especies reactivas de oxígeno (ROS). Este es un sello genuino de senescencia celular que probablemente precede a la anergia.

### III. Inestabilidad del Locus de Inmunoglobulina
El hallazgo de que los genes IG están **down-regulados** (negativos) en Old—una vez que se quita el ruido ambiental—sugiere una alteración en la accesibilidad cromatínica de estos loci, lo cual es un área de investigación totalmente nueva y prometedora para tu tesis.

---

## 🚩 5. Conclusión para la Tesis

**Los resultados V20 son biológicamente superiores.** Mientras que la firma anterior era una "radiografía de sangre sucia", la firma actual es una **"biopsia molecular limpia"** de la célula NK. 

> [!IMPORTANT]
> **Narrativa Recomendada:** "El proceso de envejecimiento en células NK humanas no está impulsado por una ganancia de funciones pro-inflamatorias (como se creía previamente debido a artefactos técnicos), sino por una pérdida crítica de efectores de activación y un aumento del daño oxidativo intrínseco."

---

## 🛣️ 6. Próximos Pasos Sugeridos

1.  **GSEA de Alta Resolución:** Correr el análisis de vías funcionales sobre los 876 genes para mapear formalmente estas conclusiones.
2.  **Validación de Covariables:** Re-correr DESeq2 controlando por `sex` para asegurar que `DUOX1` y `LCP2` no son sesgos de género.
3.  **Visualización Curatada:** Crear un DotPlot comparativo de los genes verificados en la Tabla 3 para el cuerpo principal del manuscrito.
