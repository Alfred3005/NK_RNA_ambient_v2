# Reporte de Pureza de Linaje: Enfoque en la Exclusión de Células B

Este documento describe detalladamente la estrategia biológica y computacional utilizada en nuestra tubería de análisis de secuenciación de ARN de célula única (scRNA-seq) para identificar, evaluar y excluir células pertenecientes al linaje de células B. Esta depuración es fundamental para garantizar que los análisis transcriptómicos posteriores (como los de expresión diferencial y abundancia de subtipos de células NK) se realicen sobre una población celular pura y libre de contaminación heterotípica.

---

## 1. El Linaje de Células B (Biología y Ontogenia)

Las **células B (linocitos B)** constituyen uno de los pilares de la inmunidad adaptativa humoral, cuya función principal es la producción de anticuerpos específicos de antígeno (inmunoglobulinas) y la presentación de antígenos a las células T colaboradoras.

### Ontogenia y Maduración
* **Origen:** Se originan a partir de células madre hematopoyéticas (HSC) en la médula ósea.
* **Compromiso de Linaje:** El factor de transcripción maestro **PAX5** es el regulador crítico que induce el programa transcripcional de células B (activando genes como `CD19`, `CD79A` y `CD79B`) mientras reprime activamente otros linajes (como Notch1 para células T o receptores mieloides).
* **Fases de Desarrollo:** Pasan por fases ordenadas de reordenamiento de los genes de inmunoglobulinas:
  1. **Pro-B:** Estadio inicial, inicio de recombinación de cadena pesada.
  2. **Pre-B:** Expresión del receptor de células B pre-inmaduro.
  3. **Célula B Inmadura:** Expresión de IgM en superficie. Migran de la médula ósea al bazo y nódulos linfáticos.
  4. **Célula B Madura (Naive):** Expresan IgM e IgD en superficie. Listas para reconocer antígenos.
* **Fase Efectora (Célula Plasmática):** Tras la activación antigénica y la cooperación con células T, las células B proliferan en los centros germinales (donde sufren hipermutación somática y cambio de isotipo) y se diferencian en:
  * **Células B de Memoria:** Células de vida larga listas para una respuesta secundaria rápida.
  * **Células Plasmáticas:** Fábricas de anticuerpos terminalmente diferenciadas. Durante este proceso, las células plasmáticas pierden la expresión de marcadores clásicos de superficie (como CD19 y CD20) y regulan fuertemente al alza la maquinaria de secreción proteica (genes como `MZB1`, `JCHAIN` y `SDC1` / CD138).

---

## 2. Estrategia de Purificación y Exclusión de Células B en el Pipeline NK

En proyectos de transcriptómica de células NK, la contaminación por células B es un desafío técnico relevante. Las células B abundan en sangre periférica (PBMC) y tejidos linfoides, y los métodos de aislamiento físico (como la selección magnética o FACS) pueden arrastrar células contaminantes o fragmentos celulares ("ambient RNA").

Para solucionar esto, nuestro pipeline implementa una **estrategia de purificación en dos pasos** a nivel bioinformático:

### Paso A: Puntuación de Linaje (Lineage Scoring)
Utilizando la biblioteca `Scanpy` (`sc.tl.score_genes`), calculamos puntuaciones específicas para cada célula basadas en un panel de genes de alta confianza para cada linaje principal (NK, Células T y Células B):

* **Marcadores de Células B utilizados para la puntuación (`B_CELL_score`):**
  * `CD19`: Marcador pan-B por excelencia, expresado desde Pro-B hasta antes de la célula plasmática.
  * `MS4A1` (CD20): Expresado en células B maduras, regulado a la baja en células plasmáticas.
  * `MZB1`: Proteína específica de células B de zona marginal y células plasmáticas. **Este marcador es clave** porque permite detectar células plasmáticas secretoras de anticuerpos que ya no expresan CD19 ni CD20.

### Paso B: Umbral de Exclusión Estricto
Para depurar el dataset maestro de células NK y obtener el "Gold Standard", aplicamos una máscara booleana estricta:

```python
# 1. B-Cell Purge (Eliminación de cualquier residuo del linaje B)
mask_b = adata.obs['B_CELL_score'] < 0.1

# 2. Dominancia de Identidad (Asegurar que sea NK y no T)
mask_nk = adata.obs['NK_score'] > adata.obs['T_CELL_score']

# Filtrado final
adata_clean = adata[mask_b & mask_nk].copy()
```

Cualquier célula con un `B_CELL_score >= 0.1` es descartada inmediatamente. Este umbral es extremadamente sensible y elimina no solo células B intactas, sino también dobletes celulares (p. ej., NK-B) y células NK con altos niveles de ruido transcriptómico de fondo (ambient RNA) procedente de células B lisadas durante el ensayo.

---

## 3. Diccionario de Marcadores de Linaje B

A continuación, se detalla la función fisiológica y el rol analítico de los marcadores clave del linaje B evaluados en la fase de control de calidad y validación de pureza:

| Símbolo del Gen | Nombre Común | Función Fisiológica | Rol en el Pipeline de Exclusión |
| :--- | :--- | :--- | :--- |
| **`CD19`** | CD19 | Glicoproteína transmembrana correceptora del receptor de células B (BCR). Facilita la transducción de señales de activación. | **Marcador Core.** Utilizado en el cálculo del `B_CELL_score` para identificar células B maduras y de memoria. |
| **`MS4A1`** | CD20 | Canal de calcio de membrana activado por ligando, crucial para la diferenciación y activación celular. | **Marcador Core.** Utilizado en el cálculo del `B_CELL_score` para identificar linfocitos B en fases maduras. |
| **`MZB1`** | MZB1 | Chaperona del retículo endoplásmico que ayuda al ensamblaje y secreción de anticuerpos (IgM). | **Marcador Core.** Indispensable para detectar células plasmáticas que carecen de marcadores de superficie tradicionales (`CD19`/`CD20`). |
| **`CD79A`** | Ig-alpha (CD79a) | Componente heterodimérico de señalización unido de forma no covalente a la inmunoglobulina de superficie. | **Marcador de Validación.** Evaluado en el control post-limpieza (`validation_B_2.ipynb`) para confirmar expresión nula en el dataset purificado. |
| **`CD79B`** | Ig-beta (CD79b) | Segundo componente del heterodímero de señalización del BCR (junto con CD79A). | **Marcador de Validación.** Evaluado para asegurar que no hay remanentes del complejo BCR activo. |
| **`IGHM`** | Cadena Mu de Ig | Región constante de la cadena pesada de la inmunoglobulina M (IgM), primera clase de anticuerpo expresada. | **Marcador de Validación.** Detecta contaminación por células B naive e inmaduras y anticuerpos libres solubles. |
| **`CD74`** | Cadena gamma de MHC-II | Chaperona implicada en el plegamiento y transporte de moléculas del Complejo Mayor de Histocompatibilidad Clase II. | **Marcador de Validación / APC.** Las células B, como células presentadoras de antígenos, lo expresan densamente. |
| **`HLA-DRA`** / **`HLA-DRB1`** | MHC Clase II HLA-DR | Moléculas heterodiméricas de superficie encargadas de presentar antígenos peptídicos exógenos a linfocitos T CD4+. | **Marcadores de Validación / APC.** Controlados para evitar contaminación por células inmunológicas presentadoras profesionales. |
| **`PAX5`** | PAX5 / BSAP | Factor de transcripción master comprometedor del linaje B. Regula epigenética y transcripcionalmente la identidad B. | **Marcador de Validación.** Su ausencia confirma la nula diferenciación o desdiferenciación hacia el linaje B en nuestras células NK. |
| **`JCHAIN`** | Cadena de Unión (J Chain) | Glicoproteína ácida indispensable para la polimerización de inmunoglobulinas secretadas (IgA dimérica e IgM pentamérica). | **Marcador de Validación.** Altamente específico de células plasmáticas secretoras maduras. |
| **`CD22`** | CD22 / Siglec-2 | Receptor de adhesión celular que se une a ácidos siálicos y actúa como regulador negativo de la señalización del BCR. | **Marcador de Validación.** Específico de células B maduras. |
| **`TCL1A`** | TCL1A | Coactivador de la ruta AKT implicado en la supervivencia y proliferación linfocitaria temprana. | **Marcador de Validación.** Específico de precursores y células B naive. |

---

## 4. Validación de Pureza y Resultados Obtenidos

La eficacia de esta metodología de exclusión se verificó visualmente mediante dotplots y estadísticos de distribución detallados en el script `scAR_python_validation/scripts/purify_and_audit.py` y el notebook de validación `validation_B_2.ipynb`:

1. **Expresión Residual Prácticamente Nula:** Tras aplicar el filtro `B_CELL_score < 0.1`, los marcadores canónicos como `CD19`, `MS4A1` (CD20) y `CD79A` alcanzaron valores medios de expresión cercanos a **0.00** en el dataset "Gold Standard".
2. **Depuración de Dobletes y Ambient RNA:** Se redujo drásticamente el ruido de fondo provocado por la lisis de células B durante el procesamiento de muestras, el cual suele adherirse a las membranas de las células NK (falsa expresión debida a ARN libre o "ambient RNA" capturado en las gotitas de microfluídica).
3. **Consolidación del Dataset Gold Standard:** El dataset consolidado resultante consta de **191,903 células de alta pureza** provenientes de **547 donantes** inmunológicamente bien definidos. La exclusión del linaje B nos permitió aislar la señal transcriptómica pura de envejecimiento inmunológico (adulto vs viejo) en células NK de forma limpia, libre del sesgo que causa la mezcla celular de PBMC.
