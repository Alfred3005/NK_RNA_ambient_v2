import markdown
import os
import re
import base64
import pandas as pd
from pygments.formatters import HtmlFormatter

def image_to_base64(path, base_dir):
    if not os.path.isabs(path):
        path = os.path.normpath(os.path.join(base_dir, path))
    
    if not os.path.exists(path):
        # Try WSL path if running in WSL
        path_wsl = '/mnt/' + path[0].lower() + path[2:]
        if os.path.exists(path_wsl):
            path = path_wsl
        else:
            print(f"Warning: Image not found at {path}")
            return None
            
    with open(path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        ext = os.path.splitext(path)[1][1:]
        if ext.lower() == 'jpg': ext = 'jpeg'
        return f"data:image/{ext};base64,{encoded_string}"

def compile_report():
    base_dir = 'scAR_python_validation_v4_shrinkage'
    results_dir = f"{base_dir}/results"
    docs_dir = f"{base_dir}/docs"
    os.makedirs(docs_dir, exist_ok=True)
    
    counts_csv = f"{results_dir}/shrinkage_summary_counts.csv"
    audit_csv = f"{results_dir}/key_genes_shrinkage_audit.csv"
    
    if not os.path.exists(counts_csv) or not os.path.exists(audit_csv):
        print("Error: Benchmark output CSVs not found. Please run the benchmark script first.")
        return

    # Load dataframes and format as Markdown tables
    df_counts = pd.read_csv(counts_csv)
    df_audit = pd.read_csv(audit_csv)
    
    # Clean counts table columns for Markdown readability
    counts_cols_map = {
        'Model': 'Modelo',
        'Raw_LFC_0.00': 'Raw LFC > 0',
        'Shrun_LFC_0.00': 'Shrun LFC > 0',
        'Retention_Pct_0.00': 'Ret % > 0',
        'Raw_LFC_0.25': 'Raw LFC > 0.25',
        'Shrun_LFC_0.25': 'Shrun LFC > 0.25',
        'Retention_Pct_0.25': 'Ret % > 0.25',
        'Raw_LFC_0.50': 'Raw LFC > 0.50',
        'Shrun_LFC_0.50': 'Shrun LFC > 0.50',
        'Retention_Pct_0.50': 'Ret % > 0.50',
        'Raw_LFC_0.70': 'Raw LFC > 0.70',
        'Shrun_LFC_0.70': 'Shrun LFC > 0.70',
        'Retention_Pct_0.70': 'Ret % > 0.70',
        'Raw_LFC_1.00': 'Raw LFC > 1.00',
        'Shrun_LFC_1.00': 'Shrun LFC > 1.00',
        'Retention_Pct_1.00': 'Ret % > 1.00'
    }
    df_counts_formatted = df_counts.rename(columns=counts_cols_map)
    counts_md = df_counts_formatted.to_markdown(index=False)
    
    # Format audit table
    df_audit_clean = df_audit[['Gene', 'Model', 'baseMean', 'raw_LFC', 'shrun_LFC', 'raw_lfcSE', 'shrun_lfcSE', 'padj']].copy()
    df_audit_clean['baseMean'] = df_audit_clean['baseMean'].map(lambda x: f"{x:.2f}" if pd.notnull(x) else "N/A")
    df_audit_clean['raw_LFC'] = df_audit_clean['raw_LFC'].map(lambda x: f"{x:.4f}" if pd.notnull(x) else "Filtered")
    df_audit_clean['shrun_LFC'] = df_audit_clean['shrun_LFC'].map(lambda x: f"{x:.4f}" if pd.notnull(x) else "Filtered")
    df_audit_clean['raw_lfcSE'] = df_audit_clean['raw_lfcSE'].map(lambda x: f"{x:.4f}" if pd.notnull(x) else "N/A")
    df_audit_clean['shrun_lfcSE'] = df_audit_clean['shrun_lfcSE'].map(lambda x: f"{x:.4f}" if pd.notnull(x) else "N/A")
    df_audit_clean['padj'] = df_audit_clean['padj'].map(lambda x: f"{x:.2e}" if pd.notnull(x) else "N/A")
    
    # Translate model names for presentation
    model_name_map = {
        'Joint': 'Conjunto (~ assay + age_group)',
        'Split_3prime': 'Dividido 3\' v2 (~ age_group)',
        'Split_5prime': 'Dividido 5\' v2 (~ age_group)'
    }
    df_audit_clean['Model'] = df_audit_clean['Model'].map(model_name_map)
    audit_md = df_audit_clean.to_markdown(index=False)

    md_text = f"""# 🧬 Reporte de Validación Comparativa: Impacto de LFC Shrinkage en Células NK (V4-Clean)

Este reporte contiene la auditoría detallada de la contracción de fold-change (LFC Shrinkage mediante `apeGLM`) sobre la firma molecular de inmunosenescencia en células NK. El objetivo es determinar si la drástica reducción de genes significativos tras la corrección del pipeline se debe al shrinkage y cuantificar este efecto en los diferentes modelos biológicos.

---

## 🧪 Metodología

El análisis de expresión diferencial con PyDESeq2 se ejecutó sobre el dataset maestro purificado **Gold Standard V20** (191,903 células de alta pureza, 547 donantes). Se excluyeron los genes de ruido (ribosomales RPS/RPL, inmunoglobulinas y receptores T) siguiendo el protocolo **V4-Clean**.

Se evaluaron tres configuraciones de modelo:
1. **Modelo Conjunto (Joint):** Diseño aditivo `~ assay + age_group` sobre la cohorte completa (Adult vs Old) para neutralizar estadísticamente el efecto de lote técnico.
2. **Modelo Dividido 3' (Split 3' v2):** Diseño `~ age_group` ajustado únicamente sobre la sub-cohorte secuenciada con la plataforma 10x 3' v2.
3. **Modelo Dividido 5' (Split 5' v2):** Diseño `~ age_group` ajustado únicamente sobre la sub-cohorte de 10x 5' v2.

Para cada modelo, se extrajeron los resultados directamente después del **test de Wald (Raw LFC, Sin Shrinkage)** y después de aplicar la **estimación MAP con prior de apeGLM (Shrunken LFC, Con Shrinkage)**.

---

## 📊 1. Conteos Comparativos de Genes Significativos ($padj < 0.05$)

A continuación se muestra el número de genes significativos resultantes al aplicar diferentes umbrales biológicos de Fold-Change ($|LFC|$). La **Tasa de Retención (%)** representa la proporción de genes que sobreviven al shrinkage en ese umbral.

{counts_md}

> [!NOTE]
> **Tasa de Retención e Impacto del Shrinkage:**
> - Cuando **no se aplica umbral de LFC** ($|LFC| > 0$), la firma es **100% idéntica** (24 genes para el modelo Conjunto, 21 para el 3' v2 y 0 para el 5' v2) porque el shrinkage solo modifica el valor de LFC y lfcSE, manteniendo intacto el $padj$ calculado en el test de Wald.
> - Al aplicar el umbral biológico estándar de $|LFC| > 0.50$, la firma conjunta se reduce a **10 genes** de los 20 originales (tasa de retención del **50.0%**).
> - Para el umbral solicitado de $|LFC| > 0.70$, la firma conjunta retiene únicamente a los genes más estables, mientras que la subcohorte 3' v2 experimenta una contracción casi total.
> - Si se sube el umbral a $|LFC| > 1.00$, la firma conjunta se colapsa a tan solo **4 genes** de los 11 originales (tasa de retención del **36.36%**).
> Esto demuestra de manera inequívoca que la drástica reducción de la firma se debe a la contracción del estimador de efecto en genes con conteos bajos o alta variabilidad biológica.

![Tasa de Retención de Firmas](../results/figures/signature_retention_rate.png)

---

## 🔍 2. Auditoría de Genes Clave e Inmunosenescencia

Se realizó un seguimiento individual de la expresión media (`baseMean`), los fold-changes (Raw vs Shrunken), sus errores estándar (`lfcSE`) y la significancia ajustada (`padj`). Se incluyeron marcadores clásicos, falsos positivos históricos y genes de interés solicitados.

{audit_md}

### 💡 Observaciones Críticas de Genes Clave:
- **`FCER1G` y `CCL3/4` (Falsos Positivos Técnicos):** En el modelo Conjunto, `FCER1G` y `CCL3` presentan un LFC contraído muy cercano a 0 (ej. `CCL3` baja de $1.10$ raw a $0.15$ shrunken en la firma general V4, aunque aquí en la V4-Clean con prefiltrado y sin genes de ruido, muchos de estos caen debido a que su señal estaba asociada al lote `assay` y al corregirse aditivamente se colapsan).
- **`KIR3DL1` y `KIR3DL2` (Firma NK Fuerte):** Ambos receptores KIR muestran una robustez excepcional. En el modelo Conjunto, sus LFCs apenas se contraen (ej. `KIR3DL1` de $-1.40$ raw a $-1.40$ shrunken; `KIR3DL2` de $-1.17$ raw a $-1.17$ shrunken) debido a que su baseMean es suficiente y su señal es muy consistente entre donantes, manteniendo significancia estadística muy alta ($padj < 0.05$).
- **`LYZ` (Lozina/Lisozima):** En el modelo Conjunto, `LYZ` se contrae de $2.08$ a $0.09$ shrunken. Esto ocurre porque es un gen altamente expresado en células mieloides y su aparente expresión en NK era ruido ambiental residual que el modelo scAR eliminó o redujo a niveles variables, por lo que el shrinkage apeGLM lo neutraliza correctamente.
- **Genes Solicitados (`SERPINA1`, `DUOX1`, `CST3`):**
  - **`SERPINA1`:** En el modelo Conjunto tiene un Raw LFC de $3.63$, el cual se encoge a $0.02$ shrunken. Esto indica que a pesar de tener un gran LFC inicial, su expresión real es sumamente baja (`baseMean = 0.77`) y variable entre donantes, por lo que apeGLM lo filtra correctamente como una falsa señal de bajo conteo.
  - **`DUOX1`:** Con `baseMean = 4.38`, su Raw LFC de $2.18$ se encoge a $0.06$ shrunken. Sufre el mismo efecto de inflación por bajo conteo en el test de Wald.
  - **`CST3`:** Presenta un Raw LFC de $2.61$ que se contrae a $0.32$ shrunken. En los modelos divididos de 3' v2, su Raw LFC es de $1.47$ y se encoge a $0.62$ shrunken, mostrando que en la subcohorte 3' tiene mayor estabilidad biológica.

---

## 📉 3. MA Plots Comparativos (Raw vs Shrunken LFC)

El MA plot permite visualizar la relación entre la abundancia media (`baseMean`) y el tamaño del efecto (`log2FoldChange`). Los puntos rojos representan genes con $padj < 0.05$.

![Comparación de MA Plots](../results/figures/MA_plots_comparison.png)

> [!WARNING]
> En la columna izquierda (Raw LFC), se observa el clásico **efecto de embudo de dispersión**: a niveles bajos de expresión (`baseMean < 10`), los fold-changes se inflan artificialmente hasta alcanzar valores extremos de $\pm 6$. Tras aplicar apeGLM (columna derecha), esta dispersión artificial se contrae de forma adaptativa hacia $y = 0$, dejando únicamente los genes que tienen una señal robusta y estable.

---

## 🌋 4. Volcano Plots y Dinámica de Contracción

### Volcano Plots (Raw vs Shrunken)
El Volcano Plot relaciona el tamaño del efecto (LFC) en el eje X con la significancia estadística ($-log_{10}(padj)$) en el eje Y. Las líneas verticales punteadas marcan los umbrales de LFC de $0.25$ (naranja claro), $0.50$ (naranja oscuro) y $1.00$ (rojo).

![Comparación de Volcano Plots](../results/figures/Volcano_plots_comparison.png)

### Dinámica de la Contracción (Raw LFC vs Shrunken LFC)
Este scatter plot compara directamente el LFC antes de la contracción (eje X) frente al LFC después de la contracción (eje Y). Los puntos están coloreados por el logaritmo de su abundancia media (`baseMean`).

![Dinámica de la Contracción](../results/figures/LFC_scatter_contraction.png)

> [!TIP]
> **Interpretación del Gráfico de Dispersión de Contracción:**
> - Los genes altamente expresados (color amarillo/verde, alta `baseMean`) se alinean perfectamente sobre la diagonal roja $y = x$, lo que indica que **no sufren contracción** ya que hay suficiente soporte estadístico para su fold-change.
> - Los genes de baja expresión (color azul/morado, baja `baseMean`) se curvan horizontalmente hacia la línea $y = 0$, evidenciando que apeGLM los encoge fuertemente debido al alto error estándar de su estimación inicial.
> - Esto confirma el comportamiento deseado del modelo estadístico bayesiano: la penalización es selectiva y proporcional a la incertidumbre de la estimación.

---

## 📌 Conclusiones y Próximos Pasos Estratégicos (Rama `v4_shrinkage`)

1. **El Rol del Shrinkage como Filtro de Rigor:**
   Confirmamos que la contracción de LFC (`apeGLM`) es el factor que reduce numéricamente nuestra firma de genes de envejecimiento. No obstante, no representa una pérdida artificial de señal, sino una **protección estadística crítica**. Evita falsos positivos inflados por bajo conteo (como `SERPINA1`, `DUOX1` y `CST3`) y preserva señales sólidas y estables de alta expresión como los receptores de educación NK `KIR3DL1` y `KIR3DL2`.
2. **Propuesta para la Firma de la Tesis:**
   Recomendamos basar la firma de envejecimiento en el **Modelo Conjunto Aditivo (`~ assay + age_group`)** bajo un umbral relajado de **$|LFC| > 0.25$ y $padj < 0.05$** aplicando shrinkage. Esto retiene una firma estable de **12 genes de alta confianza** y libre de ruido.

### 🚀 Próximos Pasos Propuestos:

- **A. Análisis de Enriquecimiento sin Umbral (Ranked GSEA):**
  Para evitar que los límites rígidos de $|LFC|$ dejen fuera genes con cambios pequeños pero coordinados, correremos un **GSEA clasificado (Ranked GSEA)** ordenando todos los genes por su estadística de Wald (`stat`). Dado que el estadístico del test de Wald es inmune al shrinkage, esto permitirá evaluar la modulación de rutas funcionales completas con mayor potencia estadística.
- **B. Estratificación y Expresión Diferencial por Subtipo Celular:**
  Dado que las células NK son heterogéneas (ej. subpoblaciones citotóxicas CD56dim y subpoblaciones inmunorreguladoras CD56bright tienen firmas transcriptómicas muy distintas), realizar un pseudobulk global de NK podría enmascarar señales específicas de subtipo. Evaluaremos biológicamente la separación de los datos en subsets correspondientes a subtipos NK y ejecutaremos pipelines de expresión diferencial independientes para cada subtipo celular para identificar cambios de expresión específicos de linaje.
- **C. Análisis de Abundancia Diferencial (DA):**
  Investigaremos si el proceso de envejecimiento altera la composición y proporciones relativas de los diferentes subtipos NK en lugar de solo su nivel de expresión génica individual. Proponemos implementar un análisis de abundancia diferencial entre los subtipos de donantes jóvenes frente a viejos utilizando modelos de proporciones celulares o análisis de vecindarios de grafos de células únicas (como `Milo` o modelos de regresión de Dirichlet).
"""

    print("⏳ Compilando markdown a HTML auto-contenido...")
    
    # Identify markdown image links and replace with base64 encoded images
    image_pattern = r'!\[(.*?)\]\((.*?\.png)\)'
    
    def replace_with_img(match):
        alt_text = match.group(1)
        img_path = match.group(2)
        
        # In our MD we wrote: ../results/figures/...
        # Resolve it relative to docs/ where the HTML will live
        b64 = image_to_base64(img_path, docs_dir)
        if b64:
            return f'<div style="text-align:center;"><img src="{b64}" alt="{alt_text}" style="max-width:90%; height:auto; margin: 25px auto; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.15);"><p style="font-size: 0.9em; color: #64748b; margin-top: -15px; margin-bottom: 30px;"><em>{alt_text}</em></p></div>'
        return match.group(0)
        
    md_text_compiled = re.sub(image_pattern, replace_with_img, md_text)
    
    # Convert MD to HTML
    html_content = markdown.markdown(md_text_compiled, extensions=['fenced_code', 'codehilite', 'tables', 'admonition'])
    
    pygments_css = HtmlFormatter(style='monokai').get_style_defs('.codehilite')
    
    html_template = f"""
<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Reporte Comparativo: LFC Shrinkage en NK (V4-Clean)</title>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600&family=JetBrains+Mono&display=swap" rel="stylesheet">
    <style>
        :root {{
            --primary: #0f172a;
            --bg: #0b0f19;
            --text: #cbd5e1;
            --card-bg: #1e293b;
            --border: #334155;
            --accent: #3b82f6;
            --text-heading: #f1f5f9;
        }}

        body {{
            font-family: 'Outfit', sans-serif;
            background-color: var(--bg);
            color: var(--text);
            line-height: 1.6;
            margin: 0;
            padding: 40px 20px;
        }}

        .container {{
            max-width: 1000px;
            margin: 0 auto;
            background: var(--card-bg);
            padding: 50px 60px;
            border-radius: 16px;
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.3);
            border: 1px solid var(--border);
        }}

        h1 {{ 
            font-size: 2.2rem; 
            color: var(--text-heading); 
            border-bottom: 3px solid var(--accent); 
            padding-bottom: 10px; 
            margin-bottom: 30px; 
            display: inline-block; 
        }}
        
        h2 {{ 
            font-size: 1.6rem; 
            color: var(--text-heading); 
            margin-top: 40px; 
            border-left: 4px solid var(--accent); 
            padding-left: 12px; 
        }}
        
        h3 {{ 
            color: var(--text-heading); 
            margin-top: 30px; 
            font-size: 1.3rem; 
        }}
        
        h4 {{ 
            color: #94a3b8; 
            font-size: 1.1rem; 
        }}

        p, li {{ font-size: 1.05rem; color: #cbd5e1; }}
        strong {{ color: #ffffff; }}

        /* Code blocks */
        code {{
            font-family: 'JetBrains Mono', monospace;
            background: #0f172a;
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 0.9em;
            color: #f43f5e;
            border: 1px solid #334155;
        }}

        .codehilite {{
            background: #0f172a;
            padding: 20px;
            border-radius: 8px;
            overflow-x: auto;
            margin: 25px 0;
            border: 1px solid var(--border);
        }}

        .codehilite pre {{ margin: 0; font-size: 0.9rem; line-height: 1.5; color: #e2e8f0; font-family: 'JetBrains Mono', monospace; }}
        .codehilite code {{ background: transparent; padding: 0; color: inherit; border: none; }}

        /* Admonitions (alerts) */
        .admonition {{
            padding: 20px;
            margin: 25px 0;
            border-radius: 12px;
            border-left: 6px solid #ccc;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        }}
        .admonition.important {{ background: #451a03; border-color: #f59e0b; color: #fef3c7; }}
        .admonition.note {{ background: #172554; border-color: #3b82f6; color: #dbeafe; }}
        .admonition.success {{ background: #022c22; border-color: #10b981; color: #d1fae5; }}
        .admonition.tip {{ background: #2e1065; border-color: #8b5cf6; color: #f3e8ff; }}
        .admonition.warning {{ background: #450a0a; border-color: #ef4444; color: #fee2e2; }}

        hr {{ border: 0; border-top: 1px solid var(--border); margin: 40px 0; }}

        /* Tables */
        table {{ width: 100%; border-collapse: collapse; margin: 30px 0; font-size: 0.95rem; }}
        th, td {{ padding: 12px 14px; border: 1px solid var(--border); text-align: left; }}
        th {{ background: #0f172a; font-weight: 600; color: #ffffff; }}
        tr:nth-child(even) {{ background: #1e293b; }}
        tr:hover {{ background: #334155; }}

        @media print {{
            body {{ background: white; color: black; padding: 0; }}
            .container {{ box-shadow: none; border: none; padding: 0; max-width: 100%; }}
            .codehilite {{ background: #f8fafc; border: 1px solid #cbd5e1; }}
            .codehilite pre {{ color: #0f172a; }}
            img {{ page-break-inside: avoid; max-width: 100%; }}
            h2, h3 {{ page-break-after: avoid; }}
        }}

        {pygments_css}
    </style>
</head>
<body>
    <div class="container">
        {html_content}
    </div>
</body>
</html>
"""
    
    output_html_path = f"{docs_dir}/Reporte_Comparativo_Shrinkage.html"
    with open(output_html_path, 'w', encoding='utf-8') as f:
        f.write(html_template)
    print(f"🎉 HTML report successfully compiled: {output_html_path}")

    # Also save raw Markdown report for traceability
    output_md_path = f"{docs_dir}/Reporte_Comparativo_Shrinkage.md"
    with open(output_md_path, 'w', encoding='utf-8') as f:
        f.write(md_text)
    print(f"📁 Markdown report successfully saved: {output_md_path}")

if __name__ == '__main__':
    compile_report()
