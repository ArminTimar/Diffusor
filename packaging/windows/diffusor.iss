; Inno Setup script for the Windows installer. packaging/build.py compiles it with
; the version, the PyInstaller folder and the output name passed in as defines.
;
; The installer asks no administrator password: by default it installs for the
; current user into %LOCALAPPDATA%\Programs\Diffusor. The first page offers an
; installation for all users instead, which then asks for administrator rights.
; Installing a newer version over an older one replaces it in place.

#ifndef AppVersion
  #error AppVersion must be defined, e.g. /DAppVersion=0.5.0
#endif
#ifndef SourceDir
  #error SourceDir must be defined: the dist\Diffusor folder made by PyInstaller
#endif
#ifndef OutputDir
  #define OutputDir "..\..\dist"
#endif
#ifndef OutputName
  #define OutputName "Diffusor-" + AppVersion + "-Windows-Setup"
#endif

[Setup]
; never change AppId: Windows recognises earlier installations by it
AppId={{50BE6BC8-42ED-41CF-9CA6-D0CF41618FE5}
AppName=Diffusor
AppVersion={#AppVersion}
AppVerName=Diffusor {#AppVersion}
AppPublisher=Armin Timar
AppPublisherURL=https://github.com/ArminTimar/Diffusor
AppSupportURL=https://github.com/ArminTimar/Diffusor/issues
AppUpdatesURL=https://github.com/ArminTimar/Diffusor/releases
VersionInfoVersion={#AppVersion}
DefaultDirName={autopf}\Diffusor
DefaultGroupName=Diffusor
DisableProgramGroupPage=yes
DisableDirPage=auto
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
LicenseFile=..\..\LICENSE
SetupIconFile=..\..\diffusor\gui\icons\diffusor.ico
UninstallDisplayIcon={app}\Diffusor.exe
UninstallDisplayName=Diffusor
WizardStyle=modern
Compression=lzma2/max
SolidCompression=yes
CloseApplications=yes
OutputDir={#OutputDir}
OutputBaseFilename={#OutputName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[InstallDelete]
; the libraries of the previous version, so an update leaves none of them behind
Type: filesandordirs; Name: "{app}\_internal"

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
; the AppUserModelID matches the one Diffusor sets for its window, so a pinned
; taskbar button and the running window are one entry
Name: "{autoprograms}\Diffusor"; Filename: "{app}\Diffusor.exe"; AppUserModelID: "Diffusor.Diffusor"
Name: "{autodesktop}\Diffusor"; Filename: "{app}\Diffusor.exe"; AppUserModelID: "Diffusor.Diffusor"; Tasks: desktopicon

[Run]
Filename: "{app}\Diffusor.exe"; Description: "{cm:LaunchProgram,Diffusor}"; Flags: nowait postinstall skipifsilent
