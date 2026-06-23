import requests
import json
import os
import urllib.parse
import time

citations = [
    "Iron metabolism in ferroptosis",
    "Ferroptosis: an iron-dependent form of nonapoptotic cell death",
    "Inflammaging: a new immune-metabolic viewpoint for age-related diseases",
    "A novel role for high-mobility group a proteins in cellular senescence and heterochromatin formation",
    "Mitoferrin is essential for erythroid iron assimilation",
    "Innate immunosenescence: effect of aging on cells and receptors of the innate immune system in humans",
    "S100A8/A9 in Inflammation",
    "Immunosenescence of human natural killer cells",
    "Human CD56bright NK cells: an update",
    "CD56bright natural killer (NK) cells: an important NK cell subset",
    "The Molecular Signatures Database (MSigDB) hallmark gene set collection",
    "A flux balance of glucose metabolism clarifies the requirements of the Warburg effect"
]

output_dir = r"C:\Users\PREDATOR\Documents\Antigravity_workspaces\NK_pipeline_RNA_ambient\referencias_bibliograficas"
os.makedirs(output_dir, exist_ok=True)

def search_epmc(title):
    query = f'TITLE:"{title}"'
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={urllib.parse.quote(query)}&format=json&resultType=lite"
    try:
        response = requests.get(url, timeout=10)
        data = response.json()
        results = data.get("resultList", {}).get("result", [])
        if results:
            return results[0]
    except Exception as e:
        print(f"Error searching {title}: {e}")
    return None

def download_pdf(pmcid, filename):
    url = f"https://europepmc.org/articles/{pmcid}?pdf=render"
    try:
        response = requests.get(url, stream=True, timeout=15)
        if response.status_code == 200:
            content_type = response.headers.get("Content-Type", "")
            if "application/pdf" in content_type:
                with open(filename, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                return True
            else:
                print(f"Not a PDF response for {pmcid}")
        else:
            print(f"HTTP {response.status_code} for {pmcid}")
    except Exception as e:
        print(f"Error downloading {pmcid}: {e}")
    return False

for title in citations:
    print(f"\nBuscando: {title}")
    result = search_epmc(title)
    if not result:
        print("  -> No encontrado en EuropePMC.")
        continue
    
    pmid = result.get("pmid", "unknown_pmid")
    pmcid = result.get("pmcid")
    is_oa = result.get("isOpenAccess") == "Y"
    has_pdf = result.get("hasPDF") == "Y"
    
    print(f"  -> PMID: {pmid}, PMCID: {pmcid}, OA: {is_oa}, PDF: {has_pdf}")
    
    if pmcid and has_pdf:
        # Create a safe filename
        safe_title = "".join([c if c.isalnum() else "_" for c in title[:50]])
        filename = os.path.join(output_dir, f"{safe_title}_{pmcid}.pdf")
        
        print(f"  -> Descargando PDF...")
        success = download_pdf(pmcid, filename)
        if success:
            print(f"  -> OK: {filename}")
        else:
            print(f"  -> FALLO en la descarga.")
    else:
        print("  -> Sin PMCID o PDF de acceso abierto disponible.")
    
    time.sleep(1) # Rate limit
