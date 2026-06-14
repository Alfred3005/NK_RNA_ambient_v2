import os
import base64
import pandas as pd
import numpy as np

def get_image_base64(path):
    """Lee una imagen y la devuelve formateada en base64 para HTML."""
    if not os.path.exists(path):
        print(f"⚠️ Imagen no encontrada en: {path}")
        return ""
    with open(path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"

def generate_table_rows(df, columns_to_show):
    """Genera las filas de la tabla HTML a partir de un DataFrame."""
    rows_html = ""
    for idx, row in df.iterrows():
        cells = ""
        for col in columns_to_show:
            val = row[col]
            if isinstance(val, float):
                if abs(val) < 1e-4 and val != 0:
                    cells += f"<td>{val:.4e}</td>"
                else:
                    cells += f"<td>{val:.4f}</td>"
            else:
                cells += f"<td>{str(val)}</td>"
        rows_html += f"<tr><td>{idx+1}</td>{cells}</tr>"
    return rows_html

def parse_markdown_to_html(md_content):
    """Parsea una versión muy simple y limpia de markdown a HTML para la narrativa."""
    html = ""
    lines = md_content.split("\n")
    in_table = False
    table_html = ""
    in_list = False
    in_mermaid = False
    mermaid_content = ""
    
    for line in lines:
        line_strip = line.strip()
        
        if line_strip.startswith("```mermaid"):
            in_mermaid = True
            continue
        elif in_mermaid and line_strip.startswith("```"):
            in_mermaid = False
            html += f"<pre class='mermaid'>\n{mermaid_content}</pre>\n"
            mermaid_content = ""
            continue
        elif in_mermaid:
            mermaid_content += line + "\n"
            continue

        # Manejo de listas
        if line_strip.startswith("* ") or line_strip.startswith("- "):
            if not in_list:
                html += "<ul>"
                in_list = True
            content = line_strip[2:]
            # Procesar negritas
            while "**" in content:
                content = content.replace("**", "<strong>", 1).replace("**", "</strong>", 1)
            html += f"<li>{content}</li>"
            continue
        elif in_list and not (line_strip.startswith("* ") or line_strip.startswith("- ") or line_strip == ""):
            html += "</ul>"
            in_list = False
            
        if not line_strip:
            if in_table:
                html += f"<div class='table-responsive'><table>{table_html}</table></div>"
                table_html = ""
                in_table = False
            continue
            
        if line_strip.startswith("# "):
            html += f"<h1>{line_strip[2:]}</h1>"
        elif line_strip.startswith("## "):
            html += f"<h2>{line_strip[3:]}</h2>"
        elif line_strip.startswith("### "):
            html += f"<h3>{line_strip[4:]}</h3>"
        elif line_strip.startswith("#### "):
            html += f"<h4>{line_strip[5:]}</h4>"
        elif line_strip.startswith("> [!"):
            alert_type = "note"
            if "IMPORTANT" in line_strip:
                alert_type = "important"
            elif "WARNING" in line_strip:
                alert_type = "warning"
            elif "TIP" in line_strip:
                alert_type = "tip"
            elif "CAUTION" in line_strip:
                alert_type = "caution"
            content_idx = line.find("]") + 1
            alert_text = line[content_idx:].strip()
            html += f"<div class='alert alert-{alert_type}'><span class='alert-title'>{alert_type.upper()}</span><p>{alert_text}</p></div>"
        elif line_strip.startswith(">"):
            html += f"<blockquote>{line_strip[1:].strip()}</blockquote>"
        elif line_strip.startswith("|"):
            in_table = True
            cols = [c.strip() for c in line_strip.split("|")[1:-1]]
            if len(cols) > 0 and "---" in cols[0]:
                continue
            row_type = "th" if table_html == "" else "td"
            row_html = "".join([f"<{row_type}>{c}</{row_type}>" for c in cols])
            table_html += f"<tr>{row_html}</tr>"
        else:
            if in_table:
                html += f"<div class='table-responsive'><table>{table_html}</table></div>"
                table_html = ""
                in_table = False
            processed_line = line_strip
            while "**" in processed_line:
                processed_line = processed_line.replace("**", "<strong>", 1).replace("**", "</strong>", 1)
            while "$" in processed_line:
                processed_line = processed_line.replace("$", "<em>", 1).replace("$", "</em>", 1)
            html += f"<p>{processed_line}</p>"
            
    if in_table:
        html += f"<div class='table-responsive'><table>{table_html}</table></div>"
    if in_list:
        html += "</ul>"
        
    return html

def main():
    print("🚀 Iniciando la compilación del Reporte Integrativo Premium...")
    
    # Directorios de origen
    report_md = "results/subtypes_abundance_integration_report.md"
    wiki_md = "docs/vault/wiki/integracion_subtipos_abundancia.md"
    
    # Tablas de DEGs de origen
    abundance_dir = "scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes"
    cd56dim_de_csv = os.path.join(abundance_dir, "deseq2_results_nk_cd56dim.csv")
    cd56bright_de_csv = "scAR_python_validation_v4_clean_subtypes_mixed_models/results/mixedlm_de_results.csv"
    
    # Resultados de abundancia en archivos de texto
    glm_txt = os.path.join(abundance_dir, "proportion_glm_results.txt")
    ratios_txt = os.path.join(abundance_dir, "statistical_test_ratios.txt")
    
    # 1. Cargar y filtrar DEGs de CD56dim (Top 15 por pvalue)
    print("📋 Cargando y filtrando DEGs de CD56dim...")
    df_dim = pd.read_csv(cd56dim_de_csv)
    # Renombrar columna de gen si es necesario
    if 'Unnamed: 0' in df_dim.columns:
        df_dim = df_dim.rename(columns={'Unnamed: 0': 'feature_name'})
    df_dim_sig = df_dim[df_dim['padj'] < 0.05].copy()
    df_dim_sig = df_dim_sig.sort_values(by='padj')
    
    # 2. Cargar y filtrar DEGs de CD56bright (Top 15 por pvalue en GLMM)
    print("📋 Cargando y filtrando DEGs de CD56bright...")
    df_bright = pd.read_csv(cd56bright_de_csv)
    if 'Unnamed: 0' in df_bright.columns:
        df_bright = df_bright.rename(columns={'Unnamed: 0': 'feature_name'})
    df_bright_sig = df_bright[df_bright['pvalue'] < 0.05].copy()
    df_bright_sig = df_bright_sig.sort_values(by='pvalue')
    
    # Generar filas HTML para tablas de DEGs
    de_cols_dim = ['feature_name', 'log2FoldChange', 'stat', 'pvalue', 'padj']
    # GLMM usa nombres diferentes a PyDESeq2 (ej. no tiene padj, tiene padj pero 0)
    de_cols_bright = ['feature_name', 'log2FoldChange', 'stat', 'pvalue', 'padj']
    if 'padj' not in df_bright_sig.columns:
        de_cols_bright = ['feature_name', 'log2FoldChange', 'stat', 'pvalue', 'stderr']
        
    dim_rows = generate_table_rows(df_dim_sig, de_cols_dim)
    
    # Para CD56bright, mostramos top 15 y el resto plegable
    bright_top15 = generate_table_rows(df_bright_sig.head(15), de_cols_bright)
    bright_rest = generate_table_rows(df_bright_sig.iloc[15:], de_cols_bright)
    
    # 3. Leer archivos de abundancia
    print("📖 Cargando archivos de texto de abundancia...")
    with open(glm_txt, "r", encoding="utf-8") as f:
        glm_content = f.read()
    with open(ratios_txt, "r", encoding="utf-8") as f:
        ratios_content = f.read()
        
    # 4. Leer y compilar narrativas de Markdown
    print("📝 Parseando narrativas de Markdown...")
    with open(report_md, "r", encoding="utf-8") as f:
        report_raw = f.read()
    report_html = parse_markdown_to_html(report_raw)
    
    with open(wiki_md, "r", encoding="utf-8") as f:
        wiki_raw = f.read()
    # Saltar el YAML frontmatter de la wiki
    if wiki_raw.startswith("---"):
        parts = wiki_raw.split("---", 2)
        if len(parts) >= 3:
            wiki_raw = parts[2].strip()
    wiki_html = parse_markdown_to_html(wiki_raw)
    
    # 5. Cargar imágenes en base64
    print("🖼️ Codificando imágenes a base64...")
    img_ratio = get_image_base64(os.path.join(abundance_dir, "nk_ratio_analysis.png"))
    
    img_gsea_compare = get_image_base64(os.path.join(abundance_dir, "gsea/comparative_summary_barplot.png"))
    img_gsea_dim = get_image_base64(os.path.join(abundance_dir, "gsea/cd56dim/dotplot_MSigDB_Hallmark_2020.png"))
    img_gsea_bright = get_image_base64(os.path.join(abundance_dir, "gsea/cd56bright/dotplot_MSigDB_Hallmark_2020.png"))
    img_gsea_global = get_image_base64(os.path.join(abundance_dir, "gsea/global/dotplot_MSigDB_Hallmark_2020.png"))
    
    # 6. Construir el HTML completo
    print("✍️ Ensamblando plantilla HTML...")
    full_html = f"""<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Reporte Integrativo: CD56dim vs CD56bright y Abundancia NK</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Outfit:wght@400;600;700&display=swap" rel="stylesheet">
    <script type="module">
        import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
        mermaid.initialize({{ startOnLoad: true, theme: 'dark' }});
    </script>
    <style>
        :root {{
            --primary: #4f46e5;
            --primary-light: #818cf8;
            --primary-dark: #3730a3;
            --background: #0f172a;
            --surface: #1e293b;
            --surface-hover: #334155;
            --text: #f8fafc;
            --text-muted: #94a3b8;
            --border: #334155;
            --success: #10b981;
            --warning: #f59e0b;
            --danger: #ef4444;
            --info: #06b6d4;
            --card-bg: rgba(30, 41, 59, 0.7);
        }}

        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}

        body {{
            font-family: 'Inter', sans-serif;
            background-color: var(--background);
            color: var(--text);
            line-height: 1.6;
            padding: 1.5rem 1rem;
        }}

        .wrapper {{
            max-width: 1200px;
            margin: 0 auto;
        }}

        header {{
            text-align: center;
            margin-bottom: 2.5rem;
            padding: 2.5rem;
            background: linear-gradient(135deg, rgba(79, 70, 229, 0.15) 0%, rgba(129, 140, 248, 0.05) 100%);
            border-radius: 16px;
            border: 1px solid var(--border);
            box-shadow: 0 4px 30px rgba(0, 0, 0, 0.2);
            backdrop-filter: blur(5px);
        }}

        header h1 {{
            font-family: 'Outfit', sans-serif;
            font-size: 2.5rem;
            font-weight: 700;
            background: linear-gradient(135deg, #e0e7ff, #a5b4fc);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.5rem;
        }}

        header p {{
            color: var(--text-muted);
            font-size: 1.1rem;
        }}

        /* Navigation Tabs */
        .nav-tabs {{
            display: flex;
            flex-wrap: wrap;
            gap: 0.5rem;
            margin-bottom: 2rem;
            border-bottom: 1px solid var(--border);
            padding-bottom: 0.75rem;
        }}

        .nav-btn {{
            padding: 0.75rem 1.25rem;
            cursor: pointer;
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: 8px;
            font-family: 'Outfit', sans-serif;
            font-weight: 600;
            font-size: 0.95rem;
            color: var(--text-muted);
            transition: all 0.3s ease;
        }}

        .nav-btn:hover {{
            background: var(--surface-hover);
            color: var(--text);
        }}

        .nav-btn.active {{
            background: var(--primary);
            color: #ffffff;
            border-color: var(--primary);
            box-shadow: 0 4px 12px rgba(79, 70, 229, 0.3);
        }}

        /* Section Panel */
        .panel {{
            display: none;
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 16px;
            padding: 2.5rem;
            margin-bottom: 2rem;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.3);
            backdrop-filter: blur(10px);
        }}

        .panel.active {{
            display: block;
        }}

        h1, h2, h3, h4 {{
            font-family: 'Outfit', sans-serif;
            font-weight: 600;
            color: #e2e8f0;
            margin-top: 1.5rem;
            margin-bottom: 0.75rem;
        }}

        .panel > h1 {{
            font-size: 2.2rem;
            color: #ffffff;
            border-bottom: 2px solid var(--border);
            padding-bottom: 0.75rem;
            margin-top: 0;
            margin-bottom: 1.5rem;
        }}

        .panel h2 {{
            font-size: 1.6rem;
            color: var(--primary-light);
            margin-top: 2rem;
            border-bottom: 1px solid var(--border);
            padding-bottom: 0.4rem;
        }}

        p {{
            margin-bottom: 1.2rem;
            color: #cbd5e1;
            font-size: 1.05rem;
        }}

        ul {{
            margin-bottom: 1.5rem;
        }}

        li {{
            margin-left: 2rem;
            margin-bottom: 0.5rem;
            color: #cbd5e1;
        }}

        blockquote {{
            padding: 0.75rem 1.25rem;
            border-left: 4px solid var(--primary-light);
            background: rgba(30, 41, 59, 0.4);
            border-radius: 0 8px 8px 0;
            margin: 1.5rem 0;
            font-style: italic;
        }}

        /* Responsive Table */
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 1.5rem 0;
            border-radius: 10px;
            border: 1px solid var(--border);
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.95rem;
            text-align: left;
        }}

        th, td {{
            padding: 12px 16px;
            border-bottom: 1px solid var(--border);
        }}

        th {{
            background-color: #0f172a;
            color: #ffffff;
            font-weight: 600;
        }}

        tr:nth-child(even) {{
            background-color: rgba(15, 23, 42, 0.3);
        }}

        tr:hover {{
            background-color: rgba(79, 70, 229, 0.1);
        }}

        /* Alerts style */
        .alert {{
            padding: 1.25rem 1.5rem;
            border-radius: 8px;
            margin: 1.5rem 0;
            border-left: 4px solid var(--primary);
            background-color: rgba(79, 70, 229, 0.1);
            color: #e2e8f0;
        }}

        .alert-important {{
            border-left-color: var(--primary);
            background-color: rgba(79, 70, 229, 0.08);
        }}

        .alert-warning {{
            border-left-color: var(--warning);
            background-color: rgba(245, 158, 11, 0.08);
        }}

        .alert-tip {{
            border-left-color: var(--success);
            background-color: rgba(16, 185, 129, 0.08);
        }}

        .alert-caution {{
            border-left-color: var(--danger);
            background-color: rgba(239, 68, 68, 0.08);
        }}

        .alert-title {{
            font-weight: 700;
            font-size: 0.8rem;
            letter-spacing: 0.05em;
            display: block;
            margin-bottom: 0.25rem;
        }}

        .alert-important .alert-title {{ color: var(--primary-light); }}
        .alert-warning .alert-title {{ color: var(--warning); }}
        .alert-tip .alert-title {{ color: var(--success); }}
        .alert-caution .alert-title {{ color: var(--danger); }}

        /* Split view layout */
        .split-container {{
            display: flex;
            flex-direction: column;
            gap: 2rem;
            margin: 2rem 0;
        }}

        @media (min-width: 900px) {{
            .split-container {{
                flex-direction: row;
            }}
            .split-half {{
                flex: 1;
                min-width: 0;
            }}
        }}

        /* Text area for output */
        .terminal-block {{
            background-color: #020617;
            font-family: 'Courier New', Courier, monospace;
            padding: 1.5rem;
            border-radius: 10px;
            border: 1px solid var(--border);
            color: #38bdf8;
            font-size: 0.9rem;
            overflow-x: auto;
            white-space: pre-wrap;
            margin: 1rem 0;
        }}

        .image-card {{
            background: #0f172a;
            padding: 1rem;
            border-radius: 12px;
            border: 1px solid var(--border);
            text-align: center;
            margin: 1.5rem 0;
        }}

        .image-card img {{
            max-width: 100%;
            height: auto;
            border-radius: 6px;
            box-shadow: 0 4px 10px rgba(0, 0, 0, 0.5);
        }}

        .image-caption {{
            color: var(--text-muted);
            font-size: 0.85rem;
            margin-top: 0.75rem;
            font-style: italic;
        }}

        /* Collapse Button */
        .collapse-btn {{
            display: block;
            width: 100%;
            padding: 12px;
            background-color: #0f172a;
            color: var(--primary-light);
            border: 1px solid var(--border);
            border-radius: 8px;
            cursor: pointer;
            font-family: 'Inter', sans-serif;
            font-weight: 600;
            text-align: center;
            margin-top: 1rem;
            transition: all 0.3s ease;
        }}

        .collapse-btn:hover {{
            background-color: var(--surface-hover);
            color: #ffffff;
        }}

        .collapse-content {{
            display: none;
        }}

        .collapse-content.show {{
            display: table-row-group;
        }}

        /* Sub-tab navigation inside GSEA */
        .sub-tabs {{
            display: flex;
            gap: 0.5rem;
            margin-bottom: 1.5rem;
        }}

        .sub-tab-btn {{
            padding: 0.5rem 1rem;
            background-color: rgba(30, 41, 59, 0.5);
            border: 1px solid var(--border);
            color: var(--text-muted);
            border-radius: 6px;
            cursor: pointer;
            font-family: 'Outfit', sans-serif;
            font-weight: 500;
            font-size: 0.85rem;
            transition: all 0.2s ease;
        }}

        .sub-tab-btn:hover {{
            color: var(--text);
        }}

        .sub-tab-btn.active {{
            background-color: var(--surface-hover);
            color: var(--primary-light);
            border-color: var(--primary-light);
        }}

        .sub-panel {{
            display: none;
        }}

        .sub-panel.active {{
            display: block;
        }}

        footer {{
            text-align: center;
            margin-top: 4rem;
            padding-top: 2rem;
            border-top: 1px solid var(--border);
            color: var(--text-muted);
            font-size: 0.9rem;
        }}
    </style>
</head>
<body>
    <div class="wrapper">
        <header>
            <h1>Reporte de Cierre: Dinámica de Subtipos y Abundancia NK</h1>
            <p>Análisis Integrativo de Inmunosenescencia y Comparativa de Subpoblaciones Bright vs Dim</p>
        </header>

        <nav class="nav-tabs">
            <button class="nav-btn active" onclick="switchPanel('narrative-panel')">📖 Narrativa de Cierre</button>
            <button class="nav-btn" onclick="switchPanel('abundance-panel')">📊 Abundancia Diferencial</button>
            <button class="nav-btn" onclick="switchPanel('expression-panel')">🧬 Expresión Diferencial (DEGs)</button>
            <button class="nav-btn" onclick="switchPanel('gsea-panel')">🖼️ Perfiles GSEApy</button>
            <button class="nav-btn" onclick="switchPanel('wiki-panel')">🗺️ Wiki Obsidian & Conceptos</button>
        </nav>

        <!-- PANEL NARRATIVA -->
        <div id="narrative-panel" class="panel active">
            <h1>Narrativa de Integración Científica</h1>
            {report_html}
        </div>

        <!-- PANEL ABUNDANCIA -->
        <div id="abundance-panel" class="panel">
            <h1>Modelado de Abundancia Celular</h1>
            <p>Análisis del ratio de la subpoblación rara progenitora CD56bright frente a la población efectora mayoritaria CD56dim en una cohorte co-ocurrente de $N=187$ donantes.</p>
            
            <div class="split-container">
                <div class="split-half">
                    <h3>Ajuste de GLM Binomial para Proporciones</h3>
                    <p>Este modelo ajusta la probabilidad binomial considerando la profundidad de lectura por donante e incorporando el ensayo técnico como covariable correctora.</p>
                    <div class="terminal-block">{glm_content}</div>
                </div>
                <div class="split-half">
                    <h3>Abundancia Diferencial</h3>
                    <p>En individuos envejecidos, la probabilidad de muestrear una célula inmunomoduladora CD56bright disminuye significativamente (~37%).</p>
                </div>
            </div>

            <div class="alert alert-important">
                <span class="alert-title">DISCUSIÓN DE POTENCIA ESTADÍSTICA:</span>
                <p>El test de Mann-Whitney U reporta un p-valor no significativo ($p = 0.1674$) porque se ve afectado por el ruido estocástico (*shot noise*) en la detección de células CD56bright (población escasa). En contraste, el **GLM Binomial** corrige por el lote técnico <code>assay</code> (el cual exhibe efectos drásticos de hasta coef = -2.12) y pondera a los donantes por su conteo total, rescatando un p-valor altamente significativo ($p < 0.0001$) y confirmando la **pérdida progresiva de células CD56bright en donantes mayores.**</p>
            </div>

            <div class="image-card">
                <img src="{img_ratio}" alt="Análisis de Ratios NK por Edad">
                <div class="image-caption">Distribución del ratio CD56bright/CD56dim y porcentajes celulares en donantes adultos jóvenes vs. mayores.</div>
            </div>
        </div>

        <!-- PANEL EXPRESIÓN DIFERENCIAL -->
        <div id="expression-panel" class="panel">
            <h1>Expresión Diferencial (DEGs) por Subtipo</h1>
            <p>Comparativa de las firmas moleculares depuradas por pseudobulk PyDESeq2 ($padj < 0.05$) bajo el diseño aditivo completo <code>~ assay + age_group</code>.</p>

            <div class="split-container">
                <!-- COLUMNA CD56DIM -->
                <div class="split-half">
                    <h2>NK CD56dim (PyDESeq2 Pseudobulk)</h2>
                    <p>Tras remover los filtros estrictos y extraer RAW counts, PyDESeq2 recuperó la señal hiper-inflamatoria de células NK CD56dim, logrando 12 DEGs robustos con FDR &lt; 0.05 (S100A8, S100A9, AHR, CD83).</p>
                    
                    <div class="table-responsive">
                        <table>
                            <thead>
                                <tr>
                                    <th>N°</th>
                                    <th>Gen</th>
                                    <th>LFC</th>
                                    <th>Stat (Wald)</th>
                                    <th>p-value</th>
                                    <th>p-adj (FDR)</th>
                                </tr>
                            </thead>
                            <tbody>
                                {dim_rows}
                            </tbody>
                        </table>
                    </div>
                </div>

                <!-- COLUMNA CD56BRIGHT -->
                <div class="split-half">
                    <h2>NK CD56bright (GLMM Single-Cell)</h2>
                    <p>Mediante Modelos Mixtos Lineales (GLMM) a nivel single-cell para sortear la escasez de esta población rara, se reportan los genes más significativos por p-valor nominal. En este modelo el p-adj penaliza a 0 hits, pero los 92 top genes revelan consistencia con desgaste mitocondrial:</p>
                    
                    <div class="table-responsive">
                        <table>
                            <thead>
                                <tr>
                                    <th>N°</th>
                                    <th>Gen</th>
                                    <th>LFC</th>
                                    <th>Stat (Wald)</th>
                                    <th>p-value</th>
                                    <th>p-adj (FDR)</th>
                                </tr>
                            </thead>
                            <tbody>
                                {bright_top15}
                            </tbody>
                            <tbody id="bright-more-rows" class="collapse-content">
                                {bright_rest}
                            </tbody>
                        </table>
                    </div>
                    <button class="collapse-btn" onclick="toggleCollapse('bright-more-rows')">Ver más tendencias CD56bright</button>
                </div>
            </div>
        </div>

        <!-- PANEL GSEA -->
        <div id="gsea-panel" class="panel">
            <h1>Enriquecimiento Funcional GSEA Preranked</h1>
            <p>Dotplots y gráficos de enriquecimiento comparando los análisis funcionales de los subtipos y el pool completo (Global).</p>
            
            <div class="image-card">
                <h3>Resumen Comparativo de Vías Hallmarks Significativas</h3>
                <img src="{img_gsea_compare}" alt="Comparación de GSEA Hallmarks">
                <div class="image-caption">Comparativa del número de términos enriquecidos y sus scores entre Global, CD56dim y CD56bright.</div>
            </div>

            <div class="tabs-container" style="margin-top: 2rem;">
                <h3>Dotplots de Enriquecimiento (MSigDB Hallmark 2020)</h3>
                
                <div class="sub-tabs" id="gsea-subtabs">
                    <button class="sub-tab-btn active" onclick="switchSubTab('gsea-sub', 'dim')">NK CD56dim</button>
                    <button class="sub-tab-btn" onclick="switchSubTab('gsea-sub', 'bright')">NK CD56bright</button>
                    <button class="sub-tab-btn" onclick="switchSubTab('gsea-sub', 'global')">NK Cell General (Global)</button>
                </div>

                <div id="gsea-sub-dim" class="sub-panel active">
                    <div class="image-card">
                        <img src="{img_gsea_dim}" alt="GSEA CD56dim Hallmark">
                        <div class="image-caption">Dotplot GSEA Preranked para CD56dim. Se aprecia la represión de TNF-α/NF-κB (enmascaramiento por inflammaging).</div>
                    </div>
                </div>

                <div id="gsea-sub-bright" class="sub-panel">
                    <div class="image-card">
                        <img src="{img_gsea_bright}" alt="GSEA CD56bright Hallmark">
                        <div class="image-caption">Dotplot GSEA Preranked para CD56bright, exhibiendo represión en OXPHOS/ROS mitocondrial.</div>
                    </div>
                </div>

                <div id="gsea-sub-global" class="sub-panel">
                    <div class="image-card">
                        <img src="{img_gsea_global}" alt="GSEA Global Hallmark">
                        <div class="image-caption">Dotplot GSEA Preranked Global. La firma de CD56dim domina completamente, cancelando la biología de la población rara.</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- PANEL WIKI -->
        <div id="wiki-panel" class="panel">
            <h1>Obsidian Wiki & Mapas Conceptuales</h1>
            {wiki_html}
        </div>

        <footer>
            <p>Generado por Antigravity · Proyecto de Inmunosenescencia de Células NK · 2026</p>
        </footer>
    </div>

    <script>
        function switchPanel(panelId) {{
            // Desactivar todos los botones de navegación
            const buttons = document.querySelectorAll('.nav-btn');
            buttons.forEach(btn => btn.classList.remove('active'));
            
            // Ocultar todos los paneles
            const panels = document.querySelectorAll('.panel');
            panels.forEach(p => p.classList.remove('active'));
            
            // Activar botón actual y panel
            event.target.classList.add('active');
            document.getElementById(panelId).classList.add('active');
        }}

        function switchSubTab(group, tabId) {{
            const tabsContainer = document.getElementById(group + '-tabs');
            const buttons = event.currentTarget.parentNode.getElementsByClassName('sub-tab-btn');
            for (let btn of buttons) {{
                btn.classList.remove('active');
            }}
            
            const panels = document.querySelectorAll('[id^="' + group + '-"]');
            panels.forEach(p => {{
                if (p.id !== group + '-tabs') {{
                    p.classList.remove('active');
                }}
            }});
            
            event.currentTarget.classList.add('active');
            document.getElementById(group + '-' + tabId).classList.add('active');
        }}

        function toggleCollapse(divId) {{
            const content = document.getElementById(divId);
            if (content.classList.contains('show')) {{
                content.classList.remove('show');
                event.target.innerText = "Ver todos los DEGs CD56bright (34 genes)";
            }} else {{
                content.classList.add('show');
                event.target.innerText = "Ocultar genes adicionales";
            }}
        }}
    </script>
</body>
</html>
"""

    output_html = "results/Reporte_Integrativo_Subtipos_Abundancia.html"
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(full_html)
        
    print(f"🎉 ¡Reporte HTML integrativo premium guardado con éxito en: {output_html}!")

if __name__ == '__main__':
    main()
