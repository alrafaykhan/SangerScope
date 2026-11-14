; Inno Setup script for the SangerScope Windows installer.
; Built by build.ps1, which passes AppVersion, RepoRoot and StageDir.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif

[Setup]
AppId={{6E0B5C0E-6B3F-4F0B-9B0A-5A2E6C1D7F31}
AppName=SangerScope
AppVersion={#AppVersion}
AppPublisher=SangerScope
AppPublisherURL=https://github.com/alrafaykhan/SangerScope
DefaultDirName={localappdata}\Programs\SangerScope
DefaultGroupName=SangerScope
DisableProgramGroupPage=yes
DisableDirPage=yes
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
LicenseFile={#RepoRoot}\LICENSE
OutputDir={#RepoRoot}\build\windows
OutputBaseFilename=SangerScope-{#AppVersion}-Setup
Compression=lzma2/max
SolidCompression=yes
LZMANumBlockThreads=4
WizardStyle=modern
SetupIconFile={#RepoRoot}\packaging\windows\SangerScope.ico
UninstallDisplayIcon={app}\SangerScope.ico
UninstallDisplayName=SangerScope

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"

[Files]
Source: "{#StageDir}\python\*"; DestDir: "{app}\python"; Flags: recursesubdirs ignoreversion
Source: "{#StageDir}\SangerScope.bat"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#RepoRoot}\packaging\windows\SangerScope.ico"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#StageDir}\references.zip"; DestDir: "{app}"; Flags: ignoreversion nocompression
Source: "{#RepoRoot}\LICENSE"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#RepoRoot}\NOTICE"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#RepoRoot}\Terms_of_use.txt"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#RepoRoot}\README.md"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\SangerScope"; Filename: "{app}\SangerScope.bat"; WorkingDir: "{app}"; IconFilename: "{app}\SangerScope.ico"; Comment: "Sanger sequencing variant analysis"
Name: "{group}\Uninstall SangerScope"; Filename: "{uninstallexe}"
Name: "{userdesktop}\SangerScope"; Filename: "{app}\SangerScope.bat"; WorkingDir: "{app}"; IconFilename: "{app}\SangerScope.ico"; Tasks: desktopicon

[Run]
Filename: "{app}\SangerScope.bat"; Description: "Start SangerScope now"; Flags: postinstall nowait skipifsilent

[InstallDelete]
; earlier versions installed the references as 53,000 separate files
Type: filesandordirs; Name: "{app}\reference"
Type: filesandordirs; Name: "{app}\python\Lib\site-packages\sangerscope"

[UninstallDelete]
Type: filesandordirs; Name: "{app}\python"
