$ErrorActionPreference = "Stop"

$PythonExe = "C:\laragon\bin\python\python-3.13\python.exe"

if (-not (Test-Path $PythonExe)) {
    $PythonExe = "python"
}

& $PythonExe "Pathway Enrichment\pathway_enrichment_pipeline.py"
