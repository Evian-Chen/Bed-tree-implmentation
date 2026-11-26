# 直接用 g++ 編譯所有工具 (不需要 make / cmake)。輸出在 build\ 。
# 用法 (PowerShell, 在專案根目錄):  .\scripts\build.ps1
$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
New-Item -ItemType Directory -Force build | Out-Null
$flags = @("-std=c++14", "-O2", "-Wall", "-Wextra", "-Iinclude", "-static")
foreach ($t in @("ed_tool", "gram_tool", "zorder_tool", "bptree_tool", "bedtree_cli", "gen_data", "bench")) {
    Write-Host "building $t"
    & g++ @flags "src/tools/$t.cpp" "src/StringMapper.cpp" -o "build/$t.exe"
    if ($LASTEXITCODE -ne 0) { throw "build failed: $t" }
}
Write-Host "building test_all"
& g++ @flags "tests/test_all.cpp" "src/StringMapper.cpp" -o "build/test_all.exe"
if ($LASTEXITCODE -ne 0) { throw "build failed: test_all" }
Write-Host "done -> build\"
