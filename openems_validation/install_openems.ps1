# Installs openEMS (Windows x64 MSVC build), Elmer FEM, and a Python venv with the bindings.
#   powershell -ExecutionPolicy Bypass -File install_openems.ps1 [-Version v0.37.0] [-Prefix $HOME\tools]
param([string]$Version = "v0.37.0", [string]$Prefix = "$HOME\tools")
$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"   # progress bar makes Invoke-WebRequest very slow on PS 5.1
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
# Elmer FEM (native Windows build, no MPI). SourceForge is much faster than the CSC mirror.
$elmer = Join-Path $Prefix "elmer"
if (-not (Test-Path (Join-Path $elmer "bin\ElmerSolver.exe"))) {
    $ezip = Join-Path $dl "elmer.zip"
    Invoke-WebRequest "https://sourceforge.net/projects/elmerfem/files/WindowsBinaries/Win64Installer/ElmerFEM-nogui-nompi-Windows-AMD64.zip/download" -OutFile $ezip -UserAgent "Wget"
    $etmp = Join-Path $dl "elmer_extract"
    Expand-Archive $ezip -DestinationPath $etmp -Force
    Move-Item (Join-Path $etmp "ElmerFEM-nogui-nompi-Windows-AMD64") $elmer
    Remove-Item $etmp -Recurse -Force
}
Write-Host "Run the suite with:  $py run_all.py"
