param(
    [string]$BuildDirectory = "build-rc-serial"
)

$ErrorActionPreference = "Stop"
$Repository = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..\..")).Path
$Generator = Join-Path $Repository "$BuildDirectory\wcns_generate_release_cgns.exe"
$GridDirectory = Join-Path $PSScriptRoot "..\grid"

if (-not (Test-Path -LiteralPath $Generator)) {
    throw "grid generator not found: $Generator"
}
New-Item -ItemType Directory -Force -Path $GridDirectory | Out-Null

& $Generator cylinder-o `
    (Join-Path $GridDirectory "cylinder_o_32x20_r8.cgns") `
    32 20 4 1.0 8.0 2.5
if ($LASTEXITCODE -ne 0) { throw "coarse cylinder grid generation failed" }

& $Generator cylinder-o `
    (Join-Path $GridDirectory "cylinder_o_48x32_r8.cgns") `
    48 32 4 1.0 8.0 2.5
if ($LASTEXITCODE -ne 0) { throw "intermediate cylinder grid generation failed" }

& $Generator cylinder-o `
    (Join-Path $GridDirectory "cylinder_o_96x48_r10.cgns") `
    96 48 4 1.0 10.0 3.0
if ($LASTEXITCODE -ne 0) { throw "medium cylinder grid generation failed" }

Write-Host "Case07 grids generated in $GridDirectory"
