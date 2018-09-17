Set oWS = WScript.CreateObject("WScript.Shell")
Set oFSO = CreateObject("Scripting.FileSystemObject")
Const Overwrite = True

matmodfit_dir = oWS.CurrentDirectory & "\"
startfolder = oWS.ExpandEnvironmentStrings("%USERPROFILE%\Start Menu\Programs\matmodfit\")
manual = matmodfit_dir & "matmodfit.pdf"



If Not oFSO.FolderExists(startfolder) Then
  oFSO.CreateFolder startfolder
End If

Set oLink = oWS.CreateShortcut(startfolder & "matmodfit manual.lnk")
oLink.TargetPath = manual
oLink.Save

Set oLink = oWS.CreateShortcut(startfolder & "matmodfit.lnk")
oLink.TargetPath = matmodfit_dir & "matmodfit.bat"
oLink.Description = "Material model fitting"
oLink.WorkingDirectory = matmodfit_dir
oLink.Save

x=msgbox("Installation completed", 0, "matmodfit installation")