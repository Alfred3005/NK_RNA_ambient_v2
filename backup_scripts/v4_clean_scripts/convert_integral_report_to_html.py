import markdown
import os
import re
import base64
from pygments.formatters import HtmlFormatter

def image_to_base64(path):
    # If path is relative in MD (like ../results/pydeseq2/figures/...), resolve it relative to docs/ folder
    if not os.path.isabs(path):
        base_dir = 'scAR_python_validation_v4_clean/docs'
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

def convert_to_html():
    md_path = 'scAR_python_validation_v4_clean/docs/Reporte_Integral_PyDESeq2_V4.md'
    output_path = 'scAR_python_validation_v4_clean/docs/Reporte_Integral_PyDESeq2_V4.html'
    
    print(f"⏳ Reading markdown: {md_path}")
    with open(md_path, 'r', encoding='utf-8') as f:
        md_text = f.read()
    
    # Identify markdown image links: ![alt](path)
    image_pattern = r'!\[(.*?)\]\((.*?\.png)\)'
    
    def replace_with_img(match):
        alt_text = match.group(1)
        path = match.group(2)
        b64 = image_to_base64(path)
        if b64:
            return f'<div style="text-align:center;"><img src="{b64}" alt="{alt_text}" style="max-width:90%; height:auto; margin: 25px auto; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.15);"><p style="font-size: 0.9em; color: #64748b; margin-top: -15px; margin-bottom: 30px;"><em>{alt_text}</em></p></div>'
        return match.group(0)
        
    md_text = re.sub(image_pattern, replace_with_img, md_text)
    
    # Convert to HTML with tables, code blocks, and admonitions
    html_content = markdown.markdown(md_text, extensions=['fenced_code', 'codehilite', 'tables', 'admonition'])
    
    # CSS for syntax highlighting
    pygments_css = HtmlFormatter(style='monokai').get_style_defs('.codehilite')
    
    # Premium template
    html_template = f"""
<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Reporte Integral: Expresión Diferencial NK (V4-Clean)</title>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600&family=JetBrains+Mono&display=swap" rel="stylesheet">
    <style>
        :root {{
            --primary: #0f172a;
            --bg: #f8fafc;
            --text: #334155;
            --card-bg: #ffffff;
            --border: #e2e8f0;
            --accent: #2563eb;
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
            max-width: 900px;
            margin: 0 auto;
            background: var(--card-bg);
            padding: 50px 60px;
            border-radius: 16px;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.05);
        }}

        h1 {{ font-size: 2.2rem; color: #0f172a; border-bottom: 3px solid var(--accent); padding-bottom: 10px; margin-bottom: 30px; display: inline-block; }}
        h2 {{ font-size: 1.6rem; color: #1e293b; margin-top: 40px; border-left: 4px solid var(--accent); padding-left: 12px; }}
        h3 {{ color: #475569; margin-top: 30px; font-size: 1.3rem; }}
        h4 {{ color: #64748b; font-size: 1.1rem; }}

        p, li {{ font-size: 1.05rem; }}
        strong {{ color: #0f172a; }}

        /* Code blocks */
        code {{
            font-family: 'JetBrains Mono', monospace;
            background: #f1f5f9;
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 0.9em;
            color: #ef4444;
        }}

        .codehilite {{
            background: #1e1e1e;
            padding: 20px;
            border-radius: 8px;
            overflow-x: auto;
            margin: 25px 0;
            box-shadow: inset 0 2px 4px rgba(0,0,0,0.1);
        }}

        .codehilite pre {{ margin: 0; font-size: 0.9rem; line-height: 1.5; color: #d4d4d4; font-family: 'JetBrains Mono', monospace; }}
        .codehilite code {{ background: transparent; padding: 0; color: inherit; }}

        /* Admonitions (alerts) */
        .admonition {{
            padding: 20px;
            margin: 25px 0;
            border-radius: 12px;
            border-left: 6px solid #ccc;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        .admonition.important {{ background: #fffbeb; border-color: #f59e0b; color: #92400e; }}
        .admonition.note {{ background: #eff6ff; border-color: #3b82f6; color: #1e40af; }}
        .admonition.success {{ background: #f0fdf4; border-color: #22c55e; color: #166534; }}
        .admonition.tip {{ background: #f5f3ff; border-color: #8b5cf6; color: #5b21b6; }}
        .admonition.warning {{ background: #fef2f2; border-color: #ef4444; color: #991b1b; }}

        hr {{ border: 0; border-top: 1px solid var(--border); margin: 40px 0; }}

        /* Tables */
        table {{ width: 100%; border-collapse: collapse; margin: 30px 0; font-size: 0.95rem; }}
        th, td {{ padding: 14px; border: 1px solid var(--border); text-align: left; }}
        th {{ background: #f1f5f9; font-weight: 600; color: #0f172a; }}

        @media print {{
            body {{ background: white; padding: 0; }}
            .container {{ box-shadow: none; border: none; padding: 0; max-width: 100%; }}
            .codehilite {{ background: #f8fafc; border: 1px solid #cbd5e1; box-shadow: none; }}
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
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_template)
    print(f"🎉 Self-contained HTML report successfully compiled: {output_path}")

if __name__ == "__main__":
    convert_to_html()
