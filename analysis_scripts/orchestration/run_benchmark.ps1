<#
.SYNOPSIS
    Canonical LineaLinkR 0.6.4 benchmark controller.
.DESCRIPTION
    Runs the formal benchmark chain in this order:
      1) competitors
      2) LineaLinkR
      3) exact full-cell fair comparison matrix
      4) benchmark metrics

    With no filters, the controller runs the 11 competitors represented in the
    12-algorithm independent fair comparison matrix. The complete-core contract
    is fixed before execution by canonical input size; the scale-extension panel
    preserves real missing or invalid states without entering marginal means.
    Dataset and algorithm filters are intended for targeted reruns against an
    otherwise complete canonical result tree.
.PARAMETER Datasets
    Optional comma-separated dataset names.
.PARAMETER Algos
    Optional comma-separated competitor names.
.PARAMETER Force
    Pass --force to the competitor runner.
.EXAMPLE
    .\run_benchmark.ps1
.EXAMPLE
    .\run_benchmark.ps1 -Datasets pancreas,dentategyrus -Algos slingshot,DPT
#>

param(
    [string]$Datasets = "",
    [string]$Algos = "",
    [switch]$Force
)

$ErrorActionPreference = "Stop"

$ProjectRoot = [System.IO.Path]::GetFullPath(
    (Join-Path -Path $PSScriptRoot -ChildPath "..\..")
)
$PipelineRoot = Join-Path -Path $ProjectRoot -ChildPath "pipeline"
$ResultsRoot = Join-Path -Path $ProjectRoot -ChildPath "results"
$CacheRoot = Join-Path -Path $ProjectRoot -ChildPath "cache"
$LogsRoot = Join-Path -Path $ProjectRoot -ChildPath "logs\pipeline\orchestration"
$CanonicalRLibrary = Join-Path -Path $ProjectRoot -ChildPath "rlib"

$Rscript = "D:\R\R-4.5.3\bin\Rscript.exe"
$Python = "E:\Python\Pythonbio\Scripts\python.exe"

$CompetitorScript = Join-Path -Path $PipelineRoot -ChildPath "manuscript\01_run_other_algos.R"
$LineaLinkRScript = Join-Path -Path $PipelineRoot -ChildPath "manuscript\02_run_linealinkr.R"
$FairMatrixScript = Join-Path -Path $PipelineRoot -ChildPath "qa\build_fair_comparison_matrix.R"
$EvaluationScript = Join-Path -Path $PipelineRoot -ChildPath "manuscript\03_evaluate_metrics.R"
$FairMatrix = Join-Path -Path $ResultsRoot -ChildPath "paper_tables\_registry\fair_comparison_matrix.csv"

$SelectedDatasets = $Datasets -replace "\s", ""

if (
    -not [System.IO.File]::Exists($Rscript) -or
    -not [System.IO.File]::Exists($Python) -or
    -not [System.IO.File]::Exists($CompetitorScript) -or
    -not [System.IO.File]::Exists($LineaLinkRScript) -or
    -not [System.IO.File]::Exists($FairMatrixScript) -or
    -not [System.IO.File]::Exists($EvaluationScript) -or
    -not [System.IO.Directory]::Exists($CanonicalRLibrary)
) {
    [Console]::Error.WriteLine("Required interpreter or benchmark stage is missing.")
    exit 1
}

[System.IO.Directory]::CreateDirectory($ResultsRoot) | Out-Null
[System.IO.Directory]::CreateDirectory($CacheRoot) | Out-Null
[System.IO.Directory]::CreateDirectory($LogsRoot) | Out-Null

$env:LLR064_WORK_ROOT = $ProjectRoot
$env:LLR_ASSETS = $ProjectRoot
$env:LLR_PIPELINE = $PipelineRoot
$env:LLR064_RLIB = $CanonicalRLibrary
$env:LLR_SMOKE_ROOT = ""
$env:LLR_BENCHMARK_TRACK = "annotation_free"
$env:LLR_BENCHMARK_ROOT_CELL = ""
$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"

$PackageContract = 'lib <- commandArgs(trailingOnly = TRUE)[1L]; .libPaths(c(lib, .libPaths())); stopifnot(as.character(packageVersion("LineaLinkR")) == "0.6.4"); cat("LINEALINKR_PACKAGE_CONTRACT_OK\n")'
& $Rscript --encoding=UTF-8 -e $PackageContract $CanonicalRLibrary
if ($LASTEXITCODE -ne 0) {
    [Console]::Error.WriteLine("LineaLinkR package contract failed with exit code $LASTEXITCODE")
    exit $LASTEXITCODE
}

& $Python -X utf8 -c 'import pandas; print("PYTHON_RUNTIME_CONTRACT_OK")'
if ($LASTEXITCODE -ne 0) {
    [Console]::Error.WriteLine("Python runtime contract failed with exit code $LASTEXITCODE")
    exit $LASTEXITCODE
}

$PanelCountContract = 'args<-commandArgs(trailingOnly=TRUE);zz<-textConnection("cfg_out","w");sink(zz);source(args[[1L]],encoding="UTF-8");sink();close(zz);cat(sprintf("PANEL_DATASET_COUNT=%d\nORDERING_DATASET_COUNT=%d\nCORE_DATASET_COUNT=%d\nEXTENSION_DATASET_COUNT=%d\nORDERING_DATASETS=%s\nCORE_DATASETS=%s\nEXTENSION_DATASETS=%s\nPSEUDOTIME_ALGORITHMS=%s\nCOMPETITOR_ALGORITHMS=%s\n",length(CONFIG$dataset_registry),length(BENCHMARK_ORDERING_DATASETS),length(BENCHMARK_COMPLETE_CORE_DATASETS),length(BENCHMARK_SCALE_EXTENSION_DATASETS),paste(BENCHMARK_ORDERING_DATASETS,collapse=","),paste(BENCHMARK_COMPLETE_CORE_DATASETS,collapse=","),paste(BENCHMARK_SCALE_EXTENSION_DATASETS,collapse=","),paste(BENCHMARK_PSEUDOTIME_ALGORITHMS,collapse=","),paste(CONFIG$algorithm_groups$ti_trajectory,collapse=",")))'
$PanelCountOutput = & $Rscript --encoding=UTF-8 -e $PanelCountContract (Join-Path $PipelineRoot "config\00_config.R")
if ($LASTEXITCODE -ne 0) {
    [Console]::Error.WriteLine("Dataset registry contract failed with exit code $LASTEXITCODE")
    exit $LASTEXITCODE
}
$PanelCountLine = @($PanelCountOutput | Where-Object { $_ -match '^PANEL_DATASET_COUNT=[0-9]+$' }) | Select-Object -Last 1
$OrderingCountLine = @($PanelCountOutput | Where-Object { $_ -match '^ORDERING_DATASET_COUNT=[0-9]+$' }) | Select-Object -Last 1
$CoreCountLine = @($PanelCountOutput | Where-Object { $_ -match '^CORE_DATASET_COUNT=[0-9]+$' }) | Select-Object -Last 1
$ExtensionCountLine = @($PanelCountOutput | Where-Object { $_ -match '^EXTENSION_DATASET_COUNT=[0-9]+$' }) | Select-Object -Last 1
$OrderingDatasetsLine = @($PanelCountOutput | Where-Object { $_ -match '^ORDERING_DATASETS=.+$' }) | Select-Object -Last 1
$CoreDatasetsLine = @($PanelCountOutput | Where-Object { $_ -match '^CORE_DATASETS=.+$' }) | Select-Object -Last 1
$ExtensionDatasetsLine = @($PanelCountOutput | Where-Object { $_ -match '^EXTENSION_DATASETS=.+$' }) | Select-Object -Last 1
$AlgorithmsLine = @($PanelCountOutput | Where-Object { $_ -match '^PSEUDOTIME_ALGORITHMS=.+$' }) | Select-Object -Last 1
$CompetitorsLine = @($PanelCountOutput | Where-Object { $_ -match '^COMPETITOR_ALGORITHMS=.+$' }) | Select-Object -Last 1
if (
    [string]::IsNullOrWhiteSpace($PanelCountLine) -or
    [string]::IsNullOrWhiteSpace($OrderingCountLine) -or
    [string]::IsNullOrWhiteSpace($CoreCountLine) -or
    [string]::IsNullOrWhiteSpace($ExtensionCountLine) -or
    [string]::IsNullOrWhiteSpace($OrderingDatasetsLine) -or
    [string]::IsNullOrWhiteSpace($CoreDatasetsLine) -or
    [string]::IsNullOrWhiteSpace($ExtensionDatasetsLine) -or
    [string]::IsNullOrWhiteSpace($AlgorithmsLine) -or
    [string]::IsNullOrWhiteSpace($CompetitorsLine)
) {
    [Console]::Error.WriteLine("Dataset registry contract did not emit the complete-core partition.")
    exit 1
}
$ExpectedPanelDatasetCount = [int]($PanelCountLine -replace '^PANEL_DATASET_COUNT=', '')
$ExpectedFairDatasetCount = [int]($OrderingCountLine -replace '^ORDERING_DATASET_COUNT=', '')
$ExpectedCoreDatasetCount = [int]($CoreCountLine -replace '^CORE_DATASET_COUNT=', '')
$ExpectedExtensionDatasetCount = [int]($ExtensionCountLine -replace '^EXTENSION_DATASET_COUNT=', '')
$OrderingDatasetCsv = $OrderingDatasetsLine -replace '^ORDERING_DATASETS=', ''
$CoreDatasetCsv = $CoreDatasetsLine -replace '^CORE_DATASETS=', ''
$ExtensionDatasetCsv = $ExtensionDatasetsLine -replace '^EXTENSION_DATASETS=', ''
$AlgorithmCsv = $AlgorithmsLine -replace '^PSEUDOTIME_ALGORITHMS=', ''
$CompetitorCsv = $CompetitorsLine -replace '^COMPETITOR_ALGORITHMS=', ''
$CanonicalAlgorithms = @($AlgorithmCsv -split ',')
$CanonicalCompetitors = @($CompetitorCsv -split ',')
$ExpectedAlgorithmCount = $CanonicalAlgorithms.Count
$SelectedAlgorithms = if ([string]::IsNullOrWhiteSpace($Algos)) {
    $CompetitorCsv
} else {
    $Algos -replace "\s", ""
}
$ExpectedPairCount = $ExpectedFairDatasetCount * $ExpectedAlgorithmCount
$ExpectedCorePairCount = $ExpectedCoreDatasetCount * $ExpectedAlgorithmCount
$ExpectedExtensionPairCount = $ExpectedExtensionDatasetCount * $ExpectedAlgorithmCount
if ([string]::IsNullOrWhiteSpace($SelectedDatasets)) {
    $SelectedDatasets = $OrderingDatasetCsv
}

$CompetitorArgs = @(
    "--encoding=UTF-8",
    $CompetitorScript,
    "--algos=$SelectedAlgorithms"
)
if (-not [string]::IsNullOrWhiteSpace($SelectedDatasets)) {
    $CompetitorArgs += "--datasets=$SelectedDatasets"
}
if ($Force) {
    $CompetitorArgs += "--force"
} else {
    $CompetitorArgs += "--reuse-existing"
}

Write-Host "[1/4] competitors"
$CompetitorProcess = Start-Process `
    -FilePath $Rscript `
    -ArgumentList $CompetitorArgs `
    -WorkingDirectory $ProjectRoot `
    -NoNewWindow `
    -RedirectStandardOutput (Join-Path $LogsRoot "competitors.stdout.log") `
    -RedirectStandardError (Join-Path $LogsRoot "competitors.stderr.log") `
    -PassThru `
    -Wait
if ($CompetitorProcess.ExitCode -ne 0) {
    [Console]::Error.WriteLine("competitors failed with exit code $($CompetitorProcess.ExitCode)")
    exit $CompetitorProcess.ExitCode
}

$LineaLinkRArgs = @("--encoding=UTF-8", $LineaLinkRScript)
if (-not [string]::IsNullOrWhiteSpace($SelectedDatasets)) {
    $LineaLinkRArgs += "--datasets=$SelectedDatasets"
}

Write-Host "[2/4] LineaLinkR"
$LineaLinkRProcess = Start-Process `
    -FilePath $Rscript `
    -ArgumentList $LineaLinkRArgs `
    -WorkingDirectory $ProjectRoot `
    -NoNewWindow `
    -RedirectStandardOutput (Join-Path $LogsRoot "linealinkr.stdout.log") `
    -RedirectStandardError (Join-Path $LogsRoot "linealinkr.stderr.log") `
    -PassThru `
    -Wait
if ($LineaLinkRProcess.ExitCode -ne 0) {
    [Console]::Error.WriteLine("LineaLinkR failed with exit code $($LineaLinkRProcess.ExitCode)")
    exit $LineaLinkRProcess.ExitCode
}

Write-Host "[3/4] fair comparison matrix"
$FairMatrixProcess = Start-Process `
    -FilePath $Rscript `
    -ArgumentList @("--encoding=UTF-8", $FairMatrixScript) `
    -WorkingDirectory $ProjectRoot `
    -NoNewWindow `
    -RedirectStandardOutput (Join-Path $LogsRoot "fair_matrix.stdout.log") `
    -RedirectStandardError (Join-Path $LogsRoot "fair_matrix.stderr.log") `
    -PassThru `
    -Wait
if ($FairMatrixProcess.ExitCode -ne 0) {
    [Console]::Error.WriteLine("fair comparison matrix failed with exit code $($FairMatrixProcess.ExitCode)")
    exit $FairMatrixProcess.ExitCode
}

$EvaluationArgs = @("--encoding=UTF-8", $EvaluationScript)

Write-Host "[4/4] benchmark metrics"
$EvaluationProcess = Start-Process `
    -FilePath $Rscript `
    -ArgumentList $EvaluationArgs `
    -WorkingDirectory $ProjectRoot `
    -NoNewWindow `
    -RedirectStandardOutput (Join-Path $LogsRoot "evaluation.stdout.log") `
    -RedirectStandardError (Join-Path $LogsRoot "evaluation.stderr.log") `
    -PassThru `
    -Wait
if ($EvaluationProcess.ExitCode -ne 0) {
    [Console]::Error.WriteLine("benchmark metrics failed with exit code $($EvaluationProcess.ExitCode)")
    exit $EvaluationProcess.ExitCode
}

$MatrixContract = @'
import pandas as pd
import sys
x = pd.read_csv(sys.argv[1], encoding="utf-8")
complete = x["output_status"].eq("complete_output")
applicable = x["applicable"].astype(str).str.lower().eq("true")
has_artifact = x["has_artifact"].astype(str).str.lower().eq("true")
expected_algorithms = pd.Index(sys.argv[11].split(",")).sort_values()
expected_dataset_count = int(sys.argv[2])
expected_core_dataset_count = int(sys.argv[3])
expected_extension_dataset_count = int(sys.argv[4])
expected_algorithm_count = int(sys.argv[5])
expected_pair_count = int(sys.argv[6])
expected_core_pair_count = int(sys.argv[7])
expected_extension_pair_count = int(sys.argv[8])
expected_core_datasets = set(sys.argv[9].split(","))
expected_extension_datasets = set(sys.argv[10].split(","))
allowed_statuses = {
    "complete_output", "missing_output", "invalid_output_schema",
    "stale_linealinkr_identity", "duplicate_cell_ids", "nonfinite_output",
    "unexpected_cell_ids", "partial_output", "constant_output"
}
required_columns = {
    "algorithm", "dataset", "dataset_role", "applicable",
    "applicability_reason", "has_artifact", "n_reference", "n_output",
    "n_aligned", "n_missing", "n_unexpected", "n_duplicate",
    "n_nonfinite", "nonzero_range", "output_status"
}
assert required_columns.issubset(x.columns)
assert x.shape[0] == expected_pair_count
assert x[["algorithm", "dataset"]].drop_duplicates().shape[0] == expected_pair_count
assert x["algorithm"].nunique() == expected_algorithm_count
assert x["dataset"].nunique() == expected_dataset_count
assert pd.Index(x["algorithm"].drop_duplicates()).sort_values().equals(expected_algorithms)
assert x.groupby("algorithm")["dataset"].nunique().eq(expected_dataset_count).all()
assert x.groupby("dataset")["algorithm"].nunique().eq(expected_algorithm_count).all()
assert applicable.all()
assert x["applicability_reason"].fillna("").astype(str).str.strip().eq("").all()
assert x["output_status"].isin(allowed_statuses).all()
core = x.loc[x["dataset_role"].eq("complete_core")]
extension = x.loc[x["dataset_role"].eq("scale_extension")]
assert core.shape[0] == expected_core_pair_count
assert extension.shape[0] == expected_extension_pair_count
assert core["dataset"].nunique() == expected_core_dataset_count
assert extension["dataset"].nunique() == expected_extension_dataset_count
assert set(core["dataset"].unique()) == expected_core_datasets
assert set(extension["dataset"].unique()) == expected_extension_datasets
assert core["output_status"].eq("complete_output").all()
assert x.loc[~complete, "dataset_role"].eq("scale_extension").all()
assert has_artifact.loc[complete].all()
assert x.loc[complete, "nonzero_range"].astype(str).str.lower().eq("true").all()
assert x.loc[complete, "n_reference"].gt(0).all()
assert x.loc[complete, "n_output"].eq(x.loc[complete, "n_reference"]).all()
assert x.loc[complete, "n_aligned"].eq(x.loc[complete, "n_reference"]).all()
assert x.loc[complete, ["n_missing", "n_unexpected", "n_duplicate", "n_nonfinite"]].eq(0).all().all()
missing = extension["output_status"].eq("missing_output")
assert (~has_artifact.loc[extension.index[missing]]).all()
assert extension.loc[missing, ["n_output", "n_aligned"]].eq(0).all().all()
assert extension.loc[missing, "n_missing"].eq(extension.loc[missing, "n_reference"]).all()
print(f"FAIR_MATRIX_CONTRACT_OK rows={expected_pair_count} core={expected_core_pair_count} extension={expected_extension_pair_count}")
'@

Write-Host "[contract] decoded fair comparison matrix"
& $Python -X utf8 -c $MatrixContract $FairMatrix $ExpectedFairDatasetCount $ExpectedCoreDatasetCount $ExpectedExtensionDatasetCount $ExpectedAlgorithmCount $ExpectedPairCount $ExpectedCorePairCount $ExpectedExtensionPairCount $CoreDatasetCsv $ExtensionDatasetCsv $AlgorithmCsv `
    1> (Join-Path $LogsRoot "fair_matrix_contract.stdout.log") `
    2> (Join-Path $LogsRoot "fair_matrix_contract.stderr.log")
if ($LASTEXITCODE -ne 0) {
    [Console]::Error.WriteLine("fair comparison matrix content contract failed with exit code $LASTEXITCODE")
    exit $LASTEXITCODE
}

Write-Host "BENCHMARK_CONTRACT_OK"
Write-Host "Results: $ResultsRoot"
Write-Host "Cache: $CacheRoot"
Write-Host "Logs: $LogsRoot"
exit 0
