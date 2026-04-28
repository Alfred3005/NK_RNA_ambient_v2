import subprocess
import os
import time
import logging

# Configuración de logging
os.makedirs("scAR_python_validation/logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - 🚦 - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("scAR_python_validation/logs/00-orchestrator.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("Orchestrator")

# Rutas de ejecución (WSL path compatible)
VENV_PYTHON = "/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/V20_CLEAN_ANALYSIS/.v20_venv/bin/python3"
BASE_DIR = "/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation/scripts/"

SCRIPTS = [
    # "03-consolidate-scar.py", # Ya completado con éxito
    "04-purify-qc-lineage.py",
    "05-doublets-solo.py"
]

def run_script(script_name):
    script_path = os.path.join(BASE_DIR, script_name)
    logger.info(f"--- 🚀 EJECUTANDO: {script_name} ---")
    
    # Ejecutamos directamente ya que el orquestador ya corre en el entorno correcto
    cmd = f"{VENV_PYTHON} {script_path}"
    
    try:
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line.strip())
        
        process.wait()
        if process.returncode == 0:
            logger.info(f"--- ✅ COMPLETADO: {script_name} ---")
            return True
        else:
            logger.error(f"--- ❌ ERROR en {script_name} (Code: {process.returncode}) ---")
            return False
    except Exception as e:
        logger.error(f"Error lanzando {script_name}: {e}")
        return False

def main():
    logger.info("=== 🛡️ INICIANDO AUTOMATIZACIÓN TOTAL: PIPELINE PURE PYTHON V20 ===")
    
    start_time = time.time()
    
    for script in SCRIPTS:
        success = run_script(script)
        if not success:
            logger.error("🛑 Pipeline abortado por error en paso previo.")
            return

    end_time = time.time()
    duration = (end_time - start_time) / 3600
    logger.info(f"=== 🏁 PIPELINE FINALIZADO EXITOSAMENTE en {duration:.2f} horas ===")

if __name__ == "__main__":
    main()
