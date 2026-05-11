import markdown
import os
import re
import base64
from pygments.formatters import HtmlFormatter

def image_to_base64(path):
    # Standardize path
    path = path.replace('file:///', '').replace('\\', '/')
    # If it's a relative path in the MD, we might need to prepend something.
    # But currently they are absolute C:/...
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
        return f"data:image/{ext};base64,{encoded_string}"

def convert_to_html():
    # Paths
    md_path = '/mnt/c/Users/PREDATOR/.gemini/antigravity/brain/6020230e-b638-44a6-829c-36118895bdfe/VISUAL_VALIDATION_AUDIT.md'
    output_path = 'scAR_python_validation_v4_clean/docs/VISUAL_VALIDATION_AUDIT.html'
    
    # Read Markdown
    with open(md_path, 'r', encoding='utf-8') as f:
        md_text = f.read()
    
    # 1. Pre-process: Identify image links and prepare for embedding
    # Format: [Label](file:///...)
    image_pattern = r'\[(.*?)\]\((file:///.*?\.png)\)'
    
    def replace_with_img(match):
        label = match.group(1)
        path = match.group(2)
        b64 = image_to_base64(path)
        if b64:
            return f'### {label}\n<img src="{b64}" alt="{label}" style="max-width:100%; height:auto; margin: 20px 0; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">'
        return match.group(0) # Keep link if conversion fails
        
    md_text = re.sub(image_pattern, replace_with_img, md_text)
    
    # Convert to HTML
    html_content = markdown.markdown(md_text, extensions=['fenced_code', 'codehilite', 'tables', 'admonition'])
    
    # Generate Pygments CSS
    pygments_css = HtmlFormatter(style='monokai').get_style_defs('.codehilite')
    
    # Premium HTML Boilerplate
    html_template = f"""
<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Auditoría de Validación Visual - Hito V4-Clean</title>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600&family=JetBrains+Mono&display=swap" rel="stylesheet">
    <style>
        :root {{
            --primary: #2563eb;
            --bg: #f8fafc;
            --text: #1e293b;
            --card-bg: #ffffff;
            --border: #e2e8f0;
            --accent: #3b82f6;
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
            padding: 50px;
            border-radius: 16px;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
        }}

        h1 {{ font-size: 2.5rem; font-weight: 600; color: #0f172a; margin-bottom: 0.5rem; border-bottom: 3px solid var(--primary); display: inline-block; }}
        h2 {{ font-size: 1.8rem; margin-top: 3rem; color: #334155; border-left: 5px solid var(--accent); padding-left: 15px; background: #f8fafc; padding-top: 10px; padding-bottom: 10px; border-radius: 0 8px 8px 0; }}
        h3 {{ color: #475569; margin-top: 2rem; font-size: 1.2rem; }}

        p, li {{ font-size: 1.05rem; color: #334155; }}

        /* Fix for code underlining */
        a code, a pre, code, pre {{
            text-decoration: none !important;
            display: inline-block;
        }}
        
        .codehilite, .codehilite pre, .codehilite span {{
            text-decoration: none !important;
        }}

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
            padding: 24px;
            border-radius: 12px;
            overflow-x: auto;
            margin: 25px 0;
            border: 1px solid #333;
        }}

        .codehilite pre {{ margin: 0; font-family: 'JetBrains Mono', monospace; font-size: 0.85rem; line-height: 1.5; }}

        blockquote {{
            margin: 25px 0;
            padding: 20px 30px;
            background: #fdf2f2;
            border-left: 6px solid #ef4444;
            border-radius: 12px;
            font-style: italic;
        }}
        
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

        hr {{ border: 0; border-top: 2px solid var(--border); margin: 50px 0; }}

        img {{
            max-width: 100%;
            display: block;
            margin: 30px auto;
            border-radius: 12px;
            box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04);
        }}

        a {{ color: var(--primary); text-decoration: none; font-weight: 500; border-bottom: 1px dashed var(--primary); }}
        a:hover {{ color: var(--accent); border-bottom-style: solid; }}

        table {{ width: 100%; border-collapse: collapse; margin: 30px 0; font-size: 0.95rem; }}
        th, td {{ padding: 14px; border: 1px solid var(--border); text-align: left; }}
        th {{ background: #f1f5f9; font-weight: 600; }}

        @media print {{
            body {{ background: white; padding: 0; }}
            .container {{ box-shadow: none; border: none; width: 100%; max-width: 100%; padding: 0; }}
            h2 {{ background: none; border-bottom: 1px solid var(--border); border-radius: 0; }}
            .codehilite {{ white-space: pre-wrap; word-wrap: break-word; background: #fff; border: 1px solid #ddd; color: #000; }}
            .codehilite pre {{ color: #000; }}
            img {{ page-break-inside: avoid; }}
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
    
    # Save HTML
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_template)
    
    print(f"Self-contained HTML Report generated at: {output_path}")

if __name__ == "__main__":
    convert_to_html()
