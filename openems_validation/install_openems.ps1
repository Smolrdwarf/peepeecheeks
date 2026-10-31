# Installs openEMS (Windows x64 MSVC build) + a Python venv with its bindings.
#   powershell -ExecutionPolicy Bypass -File install_openems.ps1 [-Version v0.37.0] [-Prefix $HOME\tools]
param([string]$Version = "v0.37.0", [string]$Prefix = "$HOME\tools")
$ErrorActionPreference = "Stop"
$zipName = "openEMS_x64_$($Version)_msvc.zip"
$dl = Join-Path $Prefix "dl_openems"
New-Item -ItemType Directory -Force $dl | Out-Null
$zip = Join-Path $dl $zipName
if (-not (Test-Path $zip)) {
    Invoke-WebRequest "https://github.com/thliebig/openEMS-Project/releases/download/$Version/$zipName" -OutFile $zip
}
$dest = Join-Path $Prefix "openEMS"
if (-not (Test-Path (Join-Path $dest "openEMS.exe"))) {
    $tmp = Join-Path $dl "extract"
    Expand-Archive $zip -DestinationPath $tmp -Force
    Move-Item (Join-Path $tmp "openEMS") $dest
    Remove-Item $tmp -Recurse -Force
}
$venv = Join-Path $Prefix "openems-venv"
if (-not (Test-Path $venv)) { python -m venv $venv }
$py = Join-Path $venv "Scripts\python.exe"
& $py -m pip install -q -r (Join-Path $PSScriptRoot "requirements.txt")
& $py -m pip install -q --no-index --find-links (Join-Path $dest "python") openEMS
& $py -c "import os; os.environ.setdefault('CSXCAD_INSTALL_PATH', r'$dest'); import CSXCAD, openEMS; print('openEMS python bindings OK')"
Write-Host "Run the suite with:  $py run_all.py"
