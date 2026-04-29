param(
    [string]$Python = "python",
    [switch]$RenderVisuals
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $repo

$configs = @(
    "configs/stage1dc_production_2000_lam4.yaml",
    "configs/stage1dc_production_4000_lam4.yaml",
    "configs/stage1dc_production_8000_lam4.yaml"
)

foreach ($cfg in $configs) {
    Write-Host "== Running $cfg =="
    & $Python -m acs.runner $cfg
    if ($LASTEXITCODE -ne 0) {
        throw "Simulation failed for $cfg"
    }

    if ($RenderVisuals) {
        $runName = [IO.Path]::GetFileNameWithoutExtension($cfg)
        $outDir = "results/$runName"
        if (Test-Path $outDir) {
            Write-Host "== Rendering visuals for $outDir =="
            & $Python scripts/render_run_visuals.py $outDir
            if ($LASTEXITCODE -ne 0) {
                throw "Visualization failed for $outDir"
            }
        }
    }
}

Write-Host "Stage 1d.c Lam4 cell-count sweep complete."
