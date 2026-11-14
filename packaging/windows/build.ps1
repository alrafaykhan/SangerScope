# Builds the Windows installer: a private Python runtime with SangerScope and its
# dependencies installed, packed with the reference sequences by Inno Setup.
# Run from the repository root on Windows:  powershell -File packaging\windows\build.ps1
# Requires Inno Setup 6 (iscc.exe on PATH or in its default folder).
param(
    [string]$PythonVersion = "3.12.8"
)
$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"

$root = (Resolve-Path "$PSScriptRoot\..\..").Path
$stage = Join-Path $root "build\windows\stage"
$python = Join-Path $stage "python"
if (Test-Path $stage) { Remove-Item $stage -Recurse -Force }
New-Item -ItemType Directory -Force $python | Out-Null

Write-Host "Downloading Python $PythonVersion (embeddable)"
$zip = Join-Path $stage "python.zip"
Invoke-WebRequest "https://www.python.org/ftp/python/$PythonVersion/python-$PythonVersion-embed-amd64.zip" -OutFile $zip
Expand-Archive $zip -DestinationPath $python
Remove-Item $zip

# let the embedded runtime use site-packages
$pth = Get-ChildItem $python -Filter "python*._pth" | Select-Object -First 1
(Get-Content $pth.FullName) -replace "^#\s*import site", "import site" | Set-Content $pth.FullName
Add-Content $pth.FullName "Lib\site-packages"

Write-Host "Installing pip and SangerScope"
$getpip = Join-Path $stage "get-pip.py"
Invoke-WebRequest "https://bootstrap.pypa.io/get-pip.py" -OutFile $getpip
& "$python\python.exe" $getpip --no-warn-script-location
if ($LASTEXITCODE -ne 0) { throw "get-pip failed" }
Remove-Item $getpip
# The embedded runtime cannot build packages from source (its ._pth file hides pip's build
# environment), so install the dependencies as wheels and copy the pure Python package.
$deps = & "$python\python.exe" -c "import tomllib; print('\n'.join(tomllib.load(open(r'$root\pyproject.toml','rb'))['project']['dependencies']))"
if ($LASTEXITCODE -ne 0) { throw "could not read the dependencies" }
& "$python\python.exe" -m pip install --no-warn-script-location --no-cache-dir --only-binary=:all: $deps
if ($LASTEXITCODE -ne 0) { throw "pip install failed" }
Copy-Item (Join-Path $root "sangerscope") (Join-Path $python "Lib\site-packages\sangerscope") -Recurse
Get-ChildItem (Join-Path $python "Lib\site-packages\sangerscope") -Filter "__pycache__" -Recurse -Directory | Remove-Item -Recurse -Force
& "$python\python.exe" -m sangerscope --version
if ($LASTEXITCODE -ne 0) { throw "SangerScope does not start in the embedded runtime" }

Copy-Item (Join-Path $PSScriptRoot "SangerScope.bat") $stage
$version = (& "$python\python.exe" -c "import sangerscope; print(sangerscope.__version__)").Trim()

$iscc = (Get-Command iscc.exe -ErrorAction SilentlyContinue).Source
if (-not $iscc) { $iscc = "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe" }
if (-not (Test-Path $iscc)) { throw "Inno Setup 6 was not found" }
Write-Host "Building the installer for version $version"
& $iscc "/DAppVersion=$version" "/DRepoRoot=$root" "/DStageDir=$stage" (Join-Path $PSScriptRoot "SangerScope.iss")
if ($LASTEXITCODE -ne 0) { throw "Inno Setup failed" }
Write-Host "Installer written to $root\build\windows"
