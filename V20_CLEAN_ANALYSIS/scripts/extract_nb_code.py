import json
import sys

nb_path = sys.argv[1]
try:
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
    print(f"--- Code for {nb_path} ---")
    for cell in nb.get('cells', []):
        if cell.get('cell_type') == 'code':
            print(''.join(cell.get('source', [])))
            print('\n# ---------------\n')
except Exception as e:
    print(f"Error: {e}")
