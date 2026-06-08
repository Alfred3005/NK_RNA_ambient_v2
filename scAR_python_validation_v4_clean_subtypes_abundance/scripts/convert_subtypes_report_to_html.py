import markdown
import os
import re
import base64
from pygments.formatters import HtmlFormatter

def image_to_base64(path):
    # Limpiar la ruta
    path = path.replace('file:///', '').replace('\\', '/')
    
    # Manejar rutas relativas basadas en el directorio de la carpeta docs
    if not os.path.isabs(path):
        base_dir = os.path.dirname(os.path.abspath('scAR_python_validation_v4_clean_subtypes_abundance/docs/Reporte_Presentacion_Abundancia.md'))
        path = os.path.normpath(os.path.join(base_dir, path))
    
    if not os.path.exists(path):
        # Intentar ruta de WSL si está ejecutándose en Windows con WSL
        path_wsl = '/mnt/' + path[0].lower() + path[2:]
        if os.path.exists(path_wsl):
            path = path_wsl
        else:
            print(f"Advertencia: No se encontró la imagen en {path}")
            return None
            
    with open(path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        ext = os.path.splitext(path)[1][1:]
        if ext.lower() == 'jpg': ext = 'jpeg'
        return f"data:image/{ext};base64,{encoded_string}"

def convert_to_html():
    md_path = 'scAR_python_validation_v4_clean_subtypes_abundance/docs/Reporte_Completo_Subtipos_Abundancia.md'
    output_path = 'scAR_python_validation_v4_clean_subtypes_abundance/docs/Reporte_Completo_Subtipos_Abundancia.html'
    
    if not os.path.exists(md_path):
        print(f"Error: No se encontró el archivo markdown en {md_path}")
        return
        
    with open(md_path, 'r', encoding='utf-8') as f:
        md_text = f.read()
    
    # 1. Reemplazar imágenes de Markdown con imágenes base64 auto-contenidas
    image_pattern = r'!\[(.*?)\]\((.*?\.png)\)'
    
    def replace_with_img(match):
        alt_text = match.group(1)
        path = match.group(2)
        b64 = image_to_base64(path)
        if b64:
            return f'<div style="text-align:center; margin: 30px 0;"><img src="{b64}" alt="{alt_text}" style="max-width:100%; height:auto; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.1);"><p style="font-size: 0.9em; color: #64748b; margin-top: 8px;"><em>{alt_text}</em></p></div>'
        return match.group(0)
        
    md_text = re.sub(image_pattern, replace_with_img, md_text)
    
    # 2. Preprocesar alertas GFM (blockquotes con [!NOTE], [!WARNING], [!IMPORTANT])
    md_text = md_text.replace('> [!NOTE]', '> **NOTA:**')
    md_text = md_text.replace('> [!WARNING]', '> **ADVERTENCIA:**')
    md_text = md_text.replace('> [!IMPORTANT]', '> **IMPORTANTE:**')
    md_text = md_text.replace('> [!TIP]', '> **CONSEJO:**')
    md_text = md_text.replace('> [!CAUTION]', '> **ATENCIÓN:**')
    
    # 3. Convertir a HTML
    html_content = markdown.markdown(md_text, extensions=['fenced_code', 'codehilite', 'tables'])
    
    # CSS de Pygments para resaltar bloques de código si los hay
    pygments_css = HtmlFormatter(style='monokai').get_style_defs('.codehilite')
    
    # Plantilla HTML Premium con estilos optimizados para impresión PDF
    html_template = f"""<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Reporte Integrado: Composición, Expresión Diferencial y Enriquecimiento NK</title>
    <style>
        :root {{
            --primary: #0f172a;
            --bg: #f8fafc;
            --text: #334155;
            --card-bg: #ffffff;
            --border: #e2e8f0;
            --accent: #2563eb;
            --accent-light: #eff6ff;
            --warning: #eab308;
            --warning-light: #fefbeb;
        }}

        body {{
            font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
            background-color: var(--bg);
            color: var(--text);
            line-height: 1.6;
            margin: 0;
            padding: 40px 20px;
        }}

        .container {{
            max-width: 900px;
            margin: 0 auto;
            background: var(--card-bg);
            padding: 50px 60px;
            border-radius: 12px;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.05), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
        }}

        h1 {{
            font-size: 2.2rem;
            color: var(--primary);
            border-bottom: 2px solid var(--accent);
            padding-bottom: 12px;
            margin-top: 0;
            margin-bottom: 25px;
        }}

        h2 {{
            font-size: 1.5rem;
            color: #1e293b;
            margin-top: 40px;
            margin-bottom: 20px;
            border-left: 4px solid var(--accent);
            padding-left: 12px;
        }}

        h3 {{
            color: #475569;
            margin-top: 30px;
            font-size: 1.25rem;
        }}

        p, li {{
            font-size: 1.05rem;
            color: #475569;
        }}

        strong {{
            color: var(--primary);
        }}

        blockquote {{
            border-left: 4px solid var(--accent);
            background-color: var(--accent-light);
            padding: 15px 20px;
            margin: 25px 0;
            border-radius: 0 8px 8px 0;
        }}

        blockquote p {{
            margin: 0;
            font-style: normal;
        }}

        blockquote strong {{
            color: var(--accent);
        }}

        /* Estilo para advertencias específicas */
        blockquote:has(strong:contains('ADVERTENCIA')), 
        blockquote:has(strong:contains('ATENCIÓN')) {{
            border-left-color: #ef4444;
            background-color: #fef2f2;
        }}
        
        blockquote:has(strong:contains('ADVERTENCIA')) strong, 
        blockquote:has(strong:contains('ATENCIÓN')) strong {{
            color: #ef4444;
        }}

        /* Tablas elegantes */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 30px 0;
            font-size: 0.95rem;
        }}

        th {{
            background-color: #f1f5f9;
            color: var(--primary);
            text-align: left;
            padding: 12px 16px;
            font-weight: 600;
            border-bottom: 2px solid var(--border);
        }}

        td {{
            padding: 12px 16px;
            border-bottom: 1px solid var(--border);
            color: #475569;
        }}

        tr:hover {{
            background-color: #f8fafc;
        }}

        /* Código y bloques de comandos */
        code {{
            font-family: 'Consolas', 'Monaco', monospace;
            background: #f1f5f9;
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 0.9em;
            color: #e11d48;
        }}

        .codehilite {{
            background: #1e1e1e;
            padding: 20px;
            border-radius: 8px;
            overflow-x: auto;
            margin: 25px 0;
        }}

        .codehilite pre {{
            margin: 0;
            font-size: 0.9rem;
            line-height: 1.5;
            color: #d4d4d4;
        }}

        .codehilite code {{
            background: transparent;
            padding: 0;
            color: inherit;
        }}

        hr {{
            border: 0;
            border-top: 1px solid var(--border);
            margin: 40px 0;
        }}

        /* Estilos de Impresión PDF optimizados */
        @media print {{
            body {{
                background: white;
                padding: 0;
            }}
            .container {{
                box-shadow: none;
                border: none;
                padding: 0;
                max-width: 100%;
            }}
            .codehilite {{
                background: #f8fafc;
                border: 1px solid #cbd5e1;
                box-shadow: none;
            }}
            .codehilite pre {{
                color: var(--primary);
            }}
            img {{
                page-break-inside: avoid;
                max-width: 100%;
            }}
            h2, h3, table, blockquote {{
                page-break-inside: avoid;
            }}
            h2, h3 {{
                page-break-after: avoid;
            }}
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
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_template)
    
    print(f"✅ Reporte HTML de presentación generado exitosamente en: {output_path}")

if __name__ == "__main__":
    convert_to_html()
