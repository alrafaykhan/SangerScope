# Installs SangerScope for the current user without the setup program.
#
# Paste this line into Windows PowerShell and press Enter:
#   irm https://raw.githubusercontent.com/alrafaykhan/SangerScope/master/packaging/windows/install.ps1 | iex
#
# It downloads the portable package of the latest release, unpacks it to
# %LOCALAPPDATA%\Programs\SangerScope and creates Start Menu and desktop shortcuts.
# No administrator rights are needed. To remove SangerScope, delete that folder and the shortcuts.
param(
    [string]$ZipPath = "",            # use a package that is already on this computer
    [string]$InstallDir = "$env:LOCALAPPDATA\Programs\SangerScope",
    [switch]$NoShortcuts
)
$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"

$temp = $null
if (-not $ZipPath) {
    Write-Host "Looking up the latest SangerScope release..."
    $release = Invoke-RestMethod "https://api.github.com/repos/alrafaykhan/SangerScope/releases/latest" -Headers @{ "User-Agent" = "sangerscope-install" }
    $asset = $release.assets | Where-Object { $_.name -like "*-portable.zip" } | Select-Object -First 1
    if (-not $asset) { throw "The latest release has no portable package." }
    $temp = Join-Path $env:TEMP $asset.name
    Write-Host "Downloading $($asset.name) ($([math]::Round($asset.size / 1MB)) MB). This can take a few minutes..."
    Invoke-WebRequest $asset.browser_download_url -OutFile $temp
    $ZipPath = $temp
}

Write-Host "Installing to $InstallDir"
if (Test-Path $InstallDir) {
    # replace the program files of an earlier version; downloaded references live elsewhere and are kept
    foreach ($item in "python", "reference", "references.zip", "SangerScope.bat", "SangerScope.ico") {
        $path = Join-Path $InstallDir $item
        if (Test-Path $path) { Remove-Item $path -Recurse -Force }
    }
}
New-Item -ItemType Directory -Force $InstallDir | Out-Null
Add-Type -AssemblyName System.IO.Compression.FileSystem
[System.IO.Compression.ZipFile]::ExtractToDirectory($ZipPath, $InstallDir, $true)
if ($temp) { Remove-Item $temp -Force }

$launcher = Join-Path $InstallDir "SangerScope.bat"
if (-not (Test-Path $launcher)) { throw "The package is incomplete: SangerScope.bat is missing." }

if (-not $NoShortcuts) {
    $shell = New-Object -ComObject WScript.Shell
    $targets = @(
        (Join-Path ([Environment]::GetFolderPath("Programs")) "SangerScope.lnk"),
        (Join-Path ([Environment]::GetFolderPath("Desktop")) "SangerScope.lnk")
    )
    foreach ($target in $targets) {
        $shortcut = $shell.CreateShortcut($target)
        $shortcut.TargetPath = $launcher
        $shortcut.WorkingDirectory = $InstallDir
        $shortcut.IconLocation = (Join-Path $InstallDir "SangerScope.ico")
        $shortcut.Description = "Sanger sequencing variant analysis"
        $shortcut.Save()
    }
}
Write-Host ""
Write-Host "SangerScope is installed. Start it from the desktop shortcut or the Start Menu."
