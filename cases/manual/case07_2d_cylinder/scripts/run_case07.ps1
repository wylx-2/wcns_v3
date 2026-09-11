param(
    [ValidateSet("re20", "re40", "re100", "re200", "mach5-euler", "mach5-robust")]
    [string]$Case,
    [string]$BuildDirectory = "build-rc-serial",
    [int]$Ranks = 1,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
if ($Ranks -lt 1) { throw "Ranks must be positive" }

$Repository = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..\..")).Path
$Executable = Join-Path $Repository "$BuildDirectory\wcns_run.exe"
$Configurations = @{
    "re20" = "cylinder_re20.wcns"
    "re40" = "cylinder_re40.wcns"
    "re100" = "cylinder_re100.wcns"
    "re200" = "cylinder_re200.wcns"
    "mach5-euler" = "cylinder_mach5_euler.wcns"
    "mach5-robust" = "cylinder_mach5_robust.wcns"
}
$Configuration = (Resolve-Path (
    Join-Path $PSScriptRoot "..\configs\$($Configurations[$Case])"
)).Path

if (-not (Test-Path -LiteralPath $Executable)) {
    throw "wcns_run not found: $Executable"
}
$Arguments = @("--config", $Configuration)
if ($DryRun) { $Arguments += "--dry-run" }

Push-Location $Repository
try {
    if ($Ranks -eq 1) {
        & $Executable @Arguments
    } else {
        $env:I_MPI_FABRICS = "shm"
        & mpiexec -n $Ranks $Executable @Arguments
    }
    if ($LASTEXITCODE -ne 0) {
        throw "Case07 run returned exit code $LASTEXITCODE"
    }
} finally {
    Pop-Location
}
