$ErrorActionPreference = 'Stop'

$CONFIG = @{
    Python = 'E:\Python\Pythonbio\Scripts\python.exe'
    Rscript = 'D:\R\R-4.5.3\bin\Rscript.exe'
    PythonTask = 'F:\biotrade\0.64LLR\pipeline\manuscript\taskM_knn_precompute_one.py'
    RTask = 'F:\biotrade\0.64LLR\pipeline\manuscript\taskM_runtime_core_one.R'
    TempRoot = 'F:\biotrade\0.64LLR\cache\tmp\linealinkr_runtime_replicates'
    Output = 'F:\biotrade\0.64LLR\results\paper_tables\supplement\tables\runtime_benchmark\runtime_replicates_with_knn_memory.csv'
    Datasets = @('dentategyrus', 'pancreas', 'paul15')
    Replicates = @(1, 2, 3)
    PollMilliseconds = 250
}

$env:LLR064_WORK_ROOT = 'F:\biotrade\0.64LLR'
$env:LLR_ASSETS = 'F:\biotrade\0.64LLR'
$env:LLR064_RLIB = 'F:\biotrade\0.64LLR\rlib'
$env:R_LIBS_USER = 'F:/biotrade/0.64LLR/rlib;D:/R/win-library/4.5'

function Get-ProcessFamilyWorkingSet {
    param([System.Diagnostics.Process]$RootProcess)
    $RootProcess.Refresh()
    $children = Get-CimInstance Win32_Process -Filter "ParentProcessId=$($RootProcess.Id)"
    $childWorkingSet = ($children | Measure-Object -Property WorkingSetSize -Sum).Sum
    [int64]($RootProcess.WorkingSet64 + $childWorkingSet)
}

$resolvedTemp = [System.IO.Path]::GetFullPath($CONFIG.TempRoot)
$resolvedWorkspace = [System.IO.Path]::GetFullPath('F:\biotrade\0.64LLR\cache\tmp')
if (-not $resolvedTemp.StartsWith($resolvedWorkspace)) { throw 'TempRoot is outside workspace temp' }
New-Item -ItemType Directory -Path $CONFIG.TempRoot -Force | Out-Null
[System.IO.Directory]::CreateDirectory(
    [System.IO.Path]::GetDirectoryName($CONFIG.Output)
) | Out-Null
[System.IO.File]::WriteAllText(
    $CONFIG.Output,
    '',
    [System.Text.UTF8Encoding]::new($false)
)

$rows = @()
foreach ($dataset in $CONFIG.Datasets) {
    foreach ($replicate in $CONFIG.Replicates) {
        $replicateDir = Join-Path $CONFIG.TempRoot ($dataset + '_rep' + $replicate)
        New-Item -ItemType Directory -Path $replicateDir -Force | Out-Null
        $pythonStdout = Join-Path $replicateDir 'python_stdout.txt'
        $pythonStderr = Join-Path $replicateDir 'python_stderr.txt'
        $pythonTimer = [System.Diagnostics.Stopwatch]::StartNew()
        $pythonProcess = Start-Process -FilePath $CONFIG.Python -ArgumentList @('-X', 'utf8', $CONFIG.PythonTask, $dataset, $replicateDir) -PassThru -WindowStyle Hidden -RedirectStandardOutput $pythonStdout -RedirectStandardError $pythonStderr
        $pythonPeak = 0L
        while (-not $pythonProcess.HasExited) {
            $pythonPeak = [Math]::Max($pythonPeak, (Get-ProcessFamilyWorkingSet -RootProcess $pythonProcess))
            Start-Sleep -Milliseconds $CONFIG.PollMilliseconds
        }
        $pythonProcess.WaitForExit()
        $pythonTimer.Stop()
        if ($pythonProcess.ExitCode -ne 0) { throw [System.IO.File]::ReadAllText($pythonStderr, [System.Text.Encoding]::UTF8) }

        $rStdout = Join-Path $replicateDir 'r_stdout.txt'
        $rStderr = Join-Path $replicateDir 'r_stderr.txt'
        $rTimer = [System.Diagnostics.Stopwatch]::StartNew()
        $rProcess = Start-Process -FilePath $CONFIG.Rscript -ArgumentList @('--encoding=UTF-8', $CONFIG.RTask, $dataset, $replicateDir) -PassThru -WindowStyle Hidden -RedirectStandardOutput $rStdout -RedirectStandardError $rStderr
        $rPeak = 0L
        while (-not $rProcess.HasExited) {
            $rPeak = [Math]::Max($rPeak, (Get-ProcessFamilyWorkingSet -RootProcess $rProcess))
            Start-Sleep -Milliseconds $CONFIG.PollMilliseconds
        }
        $rProcess.WaitForExit()
        $rTimer.Stop()
        if ($rProcess.ExitCode -ne 0) { throw [System.IO.File]::ReadAllText($rStderr, [System.Text.Encoding]::UTF8) }

        $rLine = @([System.IO.File]::ReadAllLines($rStdout, [System.Text.Encoding]::UTF8) | Where-Object { $_ -match ('^' + [Regex]::Escape($dataset) + ',') })[-1]
        $rFields = $rLine.Split(',')
        $rows += [PSCustomObject]@{
            dataset = $dataset
            n_cells = [int]$rFields[1]
            replicate = $replicate
            knn_wall_sec = [Math]::Round($pythonTimer.Elapsed.TotalSeconds, 6)
            core_wall_sec = [Math]::Round($rTimer.Elapsed.TotalSeconds, 6)
            measured_total_wall_sec = [Math]::Round($pythonTimer.Elapsed.TotalSeconds + $rTimer.Elapsed.TotalSeconds, 6)
            python_peak_rss_mb = [Math]::Round($pythonPeak / 1MB, 3)
            r_peak_rss_mb = [Math]::Round($rPeak / 1MB, 3)
            peak_rss_mb = [Math]::Round([Math]::Max($pythonPeak, $rPeak) / 1MB, 3)
            graph_connect_sec = [double]$rFields[2]
            root_selection_sec = [double]$rFields[3]
            rust_pseudotime_sec = [double]$rFields[4]
            lineage_weights_sec = [double]$rFields[5]
            lr_static_sec = [double]$rFields[6]
            temporal_lr_core_sec = [double]$rFields[7]
            lineage_entropy_sec = [double]$rFields[8]
            bifurcation_sec = [double]$rFields[9]
            fork_contrast_sec = [double]$rFields[10]
            fork_evidence_sec = [double]$rFields[11]
            temporal_fork_sec = [double]$rFields[12]
            serialization_sec = [double]$rFields[13]
            trajectory_core_sec = [double]$rFields[14]
            llr_analysis_sec = [double]$rFields[15]
            package_version = $rFields[16]
            benchmark_protocol = $rFields[17]
            result_revision = $rFields[18]
            rds_source = $rFields[19]
            scope = 'Python normalization+HVG+PCA+kNN plus connected graph+expression-only state-aggregated attractor root+pseudotime+trunk-aware lineage+static/temporal LR+vectorized fork evidence+serialization; excludes H5AD extraction and optional UMAP visualization'
        }
        $csv = $rows | ConvertTo-Csv -NoTypeInformation
        [System.IO.File]::WriteAllLines(
            $CONFIG.Output,
            $csv,
            [System.Text.UTF8Encoding]::new($false)
        )
        Remove-Item -LiteralPath $replicateDir -Recurse -Force
    }
}

$rows | Format-Table dataset, replicate, measured_total_wall_sec, peak_rss_mb
