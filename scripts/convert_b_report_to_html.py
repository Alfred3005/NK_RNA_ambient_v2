import os
import re
import markdown
from pygments.formatters import HtmlFormatter

def convert_md_to_html():
    md_path = 'docs/reporte_diccionario_linaje_B.md'
    html_path = 'docs/reporte_diccionario_linaje_B.html'
    
    if not os.path.exists(md_path):
        print(f"Error: No se encontró el archivo Markdown en {md_path}")
        return
        
    print(f"Leyendo Markdown desde: {md_path}")
    with open(md_path, 'r', encoding='utf-8') as f:
        md_text = f.read()
        
    # Preprocesamiento de alertas estilo GitHub (> [!NOTE] / > [!IMPORTANT] / etc)
    # Buscaremos bloques de blockquotes y los convertiremos a divs de admonition correspondientes
    def replace_github_alerts(text):
        # Encontrar todas las citas que comienzan con > [!CLASE]
        # Regex para capturar todo el bloque de blockquote consecutivas
        lines = text.split('\n')
        new_lines = []
        in_blockquote = False
        blockquote_type = None
        blockquote_content = []
        
        for line in lines:
            match_alert = re.match(r'^>\s+\[!(NOTE|IMPORTANT|WARNING|TIP|CAUTION)\]', line, re.IGNORECASE)
            if match_alert:
                if in_blockquote:
                    # Cerrar bloque anterior
                    new_lines.append(f'<div class="admonition {blockquote_type.lower()}">')
                    # Procesar contenido interno del blockquote anterior como markdown
                    inner_md = markdown.markdown('\n'.join(blockquote_content), extensions=['fenced_code', 'tables'])
                    new_lines.append(inner_md)
                    new_lines.append('</div>')
                    blockquote_content = []
                in_blockquote = True
                blockquote_type = match_alert.group(1).upper()
                continue
                
            if in_blockquote and line.startswith('>'):
                # Quitar el prefijo de blockquote y agregar al contenido interno
                content_line = re.sub(r'^>\s?', '', line)
                blockquote_content.append(content_line)
            else:
                if in_blockquote:
                    # Cerrar bloque de alerta
                    new_lines.append(f'<div class="admonition {blockquote_type.lower()}">')
                    inner_md = markdown.markdown('\n'.join(blockquote_content), extensions=['fenced_code', 'tables'])
                    new_lines.append(inner_md)
                    new_lines.append('</div>')
                    in_blockquote = False
                    blockquote_type = None
                    blockquote_content = []
                new_lines.append(line)
                
        # Si quedó algún bloque abierto al final del archivo
        if in_blockquote:
            new_lines.append(f'<div class="admonition {blockquote_type.lower()}">')
            inner_md = markdown.markdown('\n'.join(blockquote_content), extensions=['fenced_code', 'tables'])
            new_lines.append(inner_md)
            new_lines.append('</div>')
            
        return '\n'.join(new_lines)
        
    # Reemplazar alertas de GitHub antes de pasarlo a markdown
    processed_md = replace_github_alerts(md_text)
    
    # Convertir el markdown a HTML
    html_content = markdown.markdown(processed_md, extensions=['fenced_code', 'codehilite', 'tables'])
    
    # Obtener estilos de Pygments para coloreado de código sintáctico
    try:
        pygments_css = HtmlFormatter(style='monokai').get_style_defs('.codehilite')
    except Exception:
        pygments_css = ""
        
    # Plantilla de diseño premium y responsiva con soporte especial para impresión (page-breaks, colores contrastantes)
    html_template = f"""<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Reporte de Pureza de Linaje: Exclusión de Células B</title>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
    <style>
        :root {{
            --primary: #1e3a8a;
            --primary-light: #eff6ff;
            --primary-border: #3b82f6;
            --bg: #f8fafc;
            --text: #1e293b;
            --text-muted: #64748b;
            --card-bg: #ffffff;
            --border: #e2e8f0;
            --accent: #2563eb;
            --shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
            --shadow-lg: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
        }}

        * {{
            box-sizing: border-box;
        }}

        body {{
            font-family: 'Outfit', sans-serif;
            background-color: var(--bg);
            color: var(--text);
            line-height: 1.6;
            margin: 0;
            padding: 40px 20px;
            -webkit-print-color-adjust: exact;
            print-color-adjust: exact;
        }}

        .container {{
            max-width: 900px;
            margin: 0 auto;
            background: var(--card-bg);
            padding: 50px 60px;
            border-radius: 16px;
            box-shadow: var(--shadow-lg);
            border: 1px solid var(--border);
        }}

        /* Encabezados */
        h1 {{
            font-size: 2.2rem;
            font-weight: 700;
            color: #0f172a;
            margin-top: 0;
            margin-bottom: 1.5rem;
            line-height: 1.25;
            border-bottom: 4px solid var(--accent);
            padding-bottom: 12px;
        }}

        h2 {{
            font-size: 1.6rem;
            font-weight: 600;
            color: #1e293b;
            margin-top: 2.5rem;
            margin-bottom: 1.2rem;
            border-left: 5px solid var(--accent);
            padding-left: 15px;
            padding-top: 2px;
            padding-bottom: 2px;
        }}

        h3 {{
            font-size: 1.25rem;
            font-weight: 600;
            color: #334155;
            margin-top: 1.8rem;
            margin-bottom: 0.8rem;
        }}

        p {{
            margin-top: 0;
            margin-bottom: 1.2rem;
            font-size: 1.05rem;
            color: #334155;
            text-align: justify;
        }}

        ul, ol {{
            margin-top: 0;
            margin-bottom: 1.2rem;
            padding-left: 24px;
        }}

        li {{
            margin-bottom: 0.6rem;
            font-size: 1.05rem;
            color: #334155;
        }}

        /* Código y pre */
        code {{
            font-family: 'JetBrains Mono', monospace;
            background: #f1f5f9;
            padding: 3px 6px;
            border-radius: 6px;
            font-size: 0.9em;
            color: #ef4444;
            word-break: break-word;
        }}

        pre code {{
            background: transparent;
            padding: 0;
            border-radius: 0;
            color: inherit;
            font-size: inherit;
        }}

        .codehilite {{
            background: #1e293b;
            padding: 20px 24px;
            border-radius: 12px;
            overflow-x: auto;
            margin: 25px 0;
            border: 1px solid #0f172a;
            box-shadow: inset 0 2px 4px rgba(0,0,0,0.2);
        }}

        .codehilite pre {{
            margin: 0;
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.9rem;
            line-height: 1.5;
            color: #f8fafc;
        }}

        /* Enlaces */
        a {{
            color: var(--accent);
            text-decoration: none;
            font-weight: 500;
            border-bottom: 1px dashed var(--accent);
            transition: all 0.2s ease;
        }}

        a:hover {{
            color: #1d4ed8;
            border-bottom-style: solid;
        }}

        /* Líneas horizontales */
        hr {{
            border: 0;
            border-top: 2px solid var(--border);
            margin: 40px 0;
        }}

        /* Tablas */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 30px 0;
            font-size: 0.95rem;
            box-shadow: var(--shadow);
            border-radius: 8px;
            overflow: hidden;
        }}

        th, td {{
            padding: 12px 16px;
            border: 1px solid var(--border);
            text-align: left;
        }}

        th {{
            background-color: #f1f5f9;
            color: #1e293b;
            font-weight: 600;
            font-size: 0.9rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}

        tr:nth-child(even) {{
            background-color: #f8fafc;
        }}

        /* Admonitions / Alertas */
        .admonition {{
            padding: 18px 24px;
            margin: 25px 0;
            border-radius: 10px;
            border-left: 6px solid #94a3b8;
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        }}

        .admonition p:last-child {{
            margin-bottom: 0;
        }}

        .admonition.note {{
            background-color: #eff6ff;
            border-color: #3b82f6;
            color: #1e40af;
        }}
        
        .admonition.note::before {{
            content: "💡 NOTA";
            font-weight: 700;
            display: block;
            margin-bottom: 8px;
            font-size: 0.9rem;
            letter-spacing: 0.05em;
        }}

        .admonition.important {{
            background-color: #fffbeb;
            border-color: #f59e0b;
            color: #92400e;
        }}
        
        .admonition.important::before {{
            content: "⚠️ IMPORTANTE";
            font-weight: 700;
            display: block;
            margin-bottom: 8px;
            font-size: 0.9rem;
            letter-spacing: 0.05em;
        }}

        .admonition.warning {{
            background-color: #fff5f5;
            border-color: #ef4444;
            color: #991b1b;
        }}
        
        .admonition.warning::before {{
            content: "🚨 ADVERTENCIA";
            font-weight: 700;
            display: block;
            margin-bottom: 8px;
            font-size: 0.9rem;
            letter-spacing: 0.05em;
        }}

        .admonition.tip {{
            background-color: #f5f3ff;
            border-color: #8b5cf6;
            color: #5b21b6;
        }}

        .admonition.tip::before {{
            content: "✨ CONSEJO";
            font-weight: 700;
            display: block;
            margin-bottom: 8px;
            font-size: 0.9rem;
            letter-spacing: 0.05em;
        }}

        /* Botón flotante para impresión rápida */
        .print-btn {{
            position: fixed;
            bottom: 30px;
            right: 30px;
            background-color: var(--accent);
            color: white;
            border: none;
            padding: 12px 24px;
            border-radius: 30px;
            font-family: 'Outfit', sans-serif;
            font-size: 1rem;
            font-weight: 500;
            cursor: pointer;
            box-shadow: 0 4px 10px rgba(37, 99, 235, 0.3);
            display: flex;
            align-items: center;
            gap: 8px;
            transition: all 0.2s ease;
            z-index: 100;
        }}

        .print-btn:hover {{
            background-color: #1d4ed8;
            transform: translateY(-2px);
            box-shadow: 0 6px 14px rgba(37, 99, 235, 0.4);
        }}

        /* Estilos de Impresión */
        @media print {{
            body {{
                background-color: #ffffff;
                color: #000000;
                padding: 0;
            }}

            .container {{
                box-shadow: none;
                border: none;
                padding: 0;
                max-width: 100%;
                width: 100%;
            }}

            .print-btn {{
                display: none; /* Ocultar botón al imprimir */
            }}

            h1, h2, h3, table, tr, img, .admonition {{
                page-break-inside: avoid;
            }}

            h2 {{
                border-left-width: 4px;
                background-color: transparent !important;
                padding-top: 0;
                padding-bottom: 0;
            }}

            .codehilite {{
                background-color: #f8fafc;
                border: 1px solid #cbd5e1;
                box-shadow: none;
            }}

            .codehilite pre {{
                color: #0f172a;
                white-space: pre-wrap;
                word-wrap: break-word;
            }}

            .admonition {{
                border-left-width: 4px;
                box-shadow: none;
            }}

            a {{
                color: #000000;
                border-bottom: none;
            }}
            
            /* Mostrar URL al lado del enlace en formato impreso */
            a::after {{
                content: " (" attr(href) ")";
                font-size: 0.85em;
                color: var(--text-muted);
            }}
            
            /* Exclusión de rutas locales del formateo para evitar URLs muy largas y feas en papel */
            a[href^="file://"]::after {{
                content: "";
            }}
        }}

        {pygments_css}
    </style>
</head>
<body>
    <button class="print-btn" onclick="window.print()">
        🖨️ Imprimir Reporte
    </button>
    <div class="container">
        {html_content}
    </div>
</body>
</html>
"""
    
    print(f"Escribiendo HTML final en: {html_path}")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_template)
    print("¡Proceso de conversión completado con éxito!")

if __name__ == "__main__":
    convert_md_to_html()
