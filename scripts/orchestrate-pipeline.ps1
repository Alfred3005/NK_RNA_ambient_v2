# scripts/orchestrate-pipeline.ps1
# PHOENIX_PROTOCOL: Phase 02 Automation (PowerShell Version)

$InDir = "data/processed/segments"
$OutDir = "data/processed/corrected"
$LogFile = "logs/pipeline_20260415_resume.log"

if (-not (Test-Path $OutDir)) { New-Item -ItemType Directory -Path $OutDir | Out-Null }
if (-not (Test-Path "logs")) { New-Item -ItemType Directory -Path "logs" | Out-Null }

Write-Host "--- Starting PHOENIX Orchestrator: Phase 02 (Ambient Correction) ---"
"--- Starting PHOENIX Orchestrator: Phase 02 (Ambient Correction) ---" | Out-File -FilePath $LogFile -Append

$Files = Get-ChildItem -Path "$InDir\*.h5ad"
$Count = $Files.Count
$Current = 0

foreach ($f in $Files) {
    $Current++
    $BaseName = $f.BaseName
    $OutFile = Join-Path $OutDir "$BaseName.rds"
    
    if (Test-Path $OutFile) {
        $msg = "[$Current/$Count] Skipping $BaseName (Already exists)"
        Write-Host $msg
        $msg | Out-File -FilePath $LogFile -Append
        continue
    }
    
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $msg = "$timestamp - [$Current/$Count] Correcting $BaseName..."
    Write-Host $msg
    $msg | Out-File -FilePath $LogFile -Append
    
    # Execute Rscript
    & Rscript scripts/02-ambient-correction.R -i "$($f.FullName)" -o "$OutFile" 2>&1 | Out-File -FilePath $LogFile -Append
    
    if ($LASTEXITCODE -ne 0) {
        $err = "Error processing $BaseName (Exit code: $LASTEXITCODE)"
        Write-Host $err -ForegroundColor Red
        $err | Out-File -FilePath $LogFile -Append
    }
}

Write-Host "--- PHOENIX ORCHESTRATOR: Phase 02 Completed ---"
"--- PHOENIX ORCHESTRATOR: Phase 02 Completed ---" | Out-File -FilePath $LogFile -Append
