$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

python "Using ML\01_build_patient_pseudobulk.py"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

python "Using ML\02_random_forest_knn_biomarkers.py"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

