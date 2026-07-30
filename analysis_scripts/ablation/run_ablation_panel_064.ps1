param(
    [switch]$DryRun
)

$ErrorActionPreference = 'Stop'

$scriptPath = [System.IO.Path]::GetFullPath($MyInvocation.MyCommand.Path)
$projectRoot = [System.IO.Path]::GetFullPath((Join-Path (Split-Path $scriptPath -Parent) '..\..\..'))
$CONFIG = @{
    rscript = 'D:\R\R-4.5.3\bin\Rscript.exe'
    ablationScript = Join-Path $projectRoot 'pipeline\manuscript\ablation\run_linealinkr_ablation.R'
    panelAudit = Join-Path $projectRoot 'audit\panel_completion.tsv'
    cache = Join-Path $projectRoot 'cache\linealinkr_ablation'
    logs = Join-Path $projectRoot 'logs\ablation_064'
    localRlib = Join-Path $projectRoot 'rlib'
    minFreeMemoryGB = 10
    registeredSlots = 37
    formalSlots = 36
    excludedDatasets = @('GSE135893_IPF')
}
$CONFIG.rscriptEngine = Join-Path (Split-Path $CONFIG.rscript -Parent) 'x64\Rscript.exe'
$utf8NoBom = [System.Text.UTF8Encoding]::new($false)

[System.IO.Directory]::CreateDirectory($CONFIG.logs) | Out-Null
$env:LANG = ''
$env:LC_ALL = ''
$env:LC_CTYPE = ''
$env:R_LIBS_USER = ''
$env:LLR064_RLIB = $CONFIG.localRlib
$env:LLR_ASSETS = $projectRoot
$env:LLR_PIPELINE = Join-Path $projectRoot 'pipeline'

$registeredRows = @(
    Import-Csv -LiteralPath $CONFIG.panelAudit -Delimiter "`t" -Encoding UTF8
)
if ($registeredRows.Count -ne $CONFIG.registeredSlots) {
    throw 'The current LineaLinkR panel audit does not contain the 37 registered datasets.'
}
if (@($registeredRows | Where-Object { $_.dataset -in $CONFIG.excludedDatasets }).Count -ne $CONFIG.excludedDatasets.Count) {
    throw 'The formal ablation exclusion is absent from the registered panel audit.'
}
$panelRows = @(
    $registeredRows | Where-Object { $_.dataset -notin $CONFIG.excludedDatasets }
)
if ($panelRows.Count -ne $CONFIG.formalSlots) {
    throw 'The current LineaLinkR panel audit is empty.'
}
if (@($panelRows | Where-Object { $_.status -ne 'complete_current' }).Count -ne 0) {
    throw 'The current LineaLinkR panel audit contains missing or incompatible RDS files.'
}
$panel = @(
    $panelRows | Sort-Object { [double]$_.returned_cells }, dataset
)
$completed = [System.Collections.Generic.List[string]]::new()
$rdsContract = 'args<-commandArgs(TRUE);ds<-args[[1L]];csv<-args[[2L]];protocol<-args[[3L]];revision<-args[[4L]];module_csv<-args[[5L]];ordering_csv<-args[[6L]];modules<-args[7:10];nokmeans<-args[[11L]];graphs<-args[12:13];module_ids<-c(''no_lr_module'',''no_temporal_binning'',''no_bifurcation_module'',''no_lineage_specific_scoring'');expected_ids<-c(''full'',module_ids,''no_kmeans_refinement'',''bins_5'',''bins_10'',''bins_12'',''bins_20'',''graph_k_10'',''graph_k_30'');d<-data.table::fread(csv);m<-data.table::fread(module_csv)[dataset==ds];o<-data.table::fread(ordering_csv)[dataset==ds];graph_rows<-d[ablation_id%in%c(''graph_k_10'',''graph_k_30'')];stopifnot(nrow(d)==12L,data.table::uniqueN(d$ablation_id)==12L,setequal(as.character(d$ablation_id),expected_ids),all(d$dataset==ds),all(d$benchmark_protocol==protocol),all(d$result_revision==revision),nrow(m)==5L,data.table::uniqueN(m$ablation)==5L,setequal(as.character(m$ablation),c(''full'',module_ids)),all(m$benchmark_protocol==protocol),all(m$result_revision==revision),nrow(o)==12L,isTRUE(all.equal(o[order(ablation_id)],d[order(ablation_id)],check.attributes=FALSE)),d[ablation_id==''no_kmeans_refinement'',execution_status]==''pipeline_rerun_done'',nrow(graph_rows)==2L,all(graph_rows$execution_status==''pipeline_rerun_done''),all(as.logical(graph_rows$graph_k_rerun_done)));check_module<-function(p,id){x<-readRDS(p);stopifnot(x$dataset==ds,x$ablation_id==id,x$execution_status==''pipeline_rerun_done'',x$package_version==''0.6.4'',x$benchmark_protocol==protocol,x$result_revision==revision,length(x$pseudotime)==x$n_cells,all(is.finite(x$pseudotime)),all(x$pseudotime>=0&x$pseudotime<=1),nrow(x$lineage_weights)==x$n_cells,ncol(x$lineage_weights)==x$n_lineages,x$start_cell>=1L,x$start_cell<=x$n_cells)};invisible(mapply(check_module,modules,module_ids,SIMPLIFY=FALSE));lineage_free<-readRDS(modules[[4L]]);stopifnot(lineage_free$n_lineages==1L,ncol(lineage_free$lineage_weights)==1L,all(lineage_free$lineage_weights==1));z<-readRDS(nokmeans);stopifnot(z$dataset==ds,z$ablation_id==''no_kmeans_refinement'',z$execution_status==''pipeline_rerun_done'',z$package_version==''0.6.4'',z$benchmark_protocol==protocol,z$result_revision==revision,isTRUE(z$config_snapshot$reused_baseline_knn_pca),grepl(''^pre_refinement:.*:no_kmeans_refinement$'',z$root_cluster_source),length(z$pseudotime)==z$n_cells,length(z$root_clusters)==z$n_cells,length(z$lineage_clusters)==z$n_cells,all(is.finite(z$pseudotime)),all(z$pseudotime>=0&z$pseudotime<=1),nrow(z$root_selection)==4L,z$selected_root_potency_mode%in%z$root_selection$candidate);check_graph<-function(i){x<-readRDS(graphs[[i]]);stopifnot(x$dataset==ds,x$execution_status==''pipeline_rerun_done'',x$package_version==''0.6.4'',x$benchmark_protocol==protocol,x$result_revision==revision,all(c(''pseudotime'',''lineage_weights'',''root_selection'',''endpoint_selection'',''direction_diagnostics'',''root_clusters'',''lineage_clusters'')%in%names(x)),length(x$pseudotime)==nrow(x$lineage_weights),all(is.finite(x$pseudotime)),all(x$pseudotime>=0&x$pseudotime<=1),nrow(x$root_selection)==4L,x$root_strategy==''attractor_energy'',x$root_potency_mode==''detected_genes_ordering'',x$selected_root_potency_mode%in%x$root_selection$candidate)};invisible(lapply(1:2,check_graph));cat(''ABLATION_RDS_CONTRACT_OK\n'')'

function Test-AblationResult {
    param(
        [string]$Dataset,
        $PanelRow,
        [string]$ResultCsv,
        [string[]]$ModuleRds
    )

    $graphRoot = Join-Path $CONFIG.cache $Dataset
    $graphRds = @(
        Join-Path $graphRoot 'graph_k_10\graph_k_10_upstream.rds'
        Join-Path $graphRoot 'graph_k_30\graph_k_30_upstream.rds'
    )
    $noKmeansRds = Join-Path $graphRoot 'no_kmeans_refinement\no_kmeans_refinement_upstream.rds'
    $moduleAggregate = Join-Path $CONFIG.cache 'ablation_module_metrics.csv'
    $orderingAggregate = Join-Path $CONFIG.cache 'ablation_ordering_metrics.csv'
    $filesPresent = [System.IO.File]::Exists($ResultCsv) -and
        [System.IO.File]::Exists($moduleAggregate) -and
        [System.IO.File]::Exists($orderingAggregate) -and
        [System.IO.File]::Exists($noKmeansRds) -and
        ($ModuleRds | Where-Object { -not [System.IO.File]::Exists($_) }).Count -eq 0 -and
        ($graphRds | Where-Object { -not [System.IO.File]::Exists($_) }).Count -eq 0
    if (-not $filesPresent) { return $false }

    $resultRows = @(Import-Csv -LiteralPath $ResultCsv -Encoding UTF8)
    $tableValid = $resultRows.Count -eq 12 -and
        @($resultRows | Where-Object { $_.dataset -ne $Dataset }).Count -eq 0 -and
        @($resultRows | Where-Object {
            $_.benchmark_protocol -ne $PanelRow.expected_benchmark_protocol
        }).Count -eq 0 -and
        @($resultRows | Where-Object {
            $_.result_revision -ne $PanelRow.expected_result_revision
        }).Count -eq 0
    if (-not $tableValid) { return $false }

    $validationOutput = @(
        & $CONFIG.rscriptEngine --vanilla --encoding=UTF-8 -e $rdsContract `
            $Dataset $ResultCsv `
            $PanelRow.expected_benchmark_protocol `
            $PanelRow.expected_result_revision `
            $moduleAggregate $orderingAggregate `
            $ModuleRds[0] $ModuleRds[1] $ModuleRds[2] $ModuleRds[3] `
            $noKmeansRds `
            $graphRds[0] $graphRds[1] 2>&1
    )
    $validationExit = $LASTEXITCODE
    if (-not $DryRun) {
        $contractPath = Join-Path $CONFIG.logs ("{0}.contract.log" -f $Dataset)
        [System.IO.File]::WriteAllLines($contractPath, $validationOutput, $utf8NoBom)
    }
    $validationExit -eq 0 -and
        @($validationOutput | Where-Object {
            $_ -match 'ABLATION_RDS_CONTRACT_OK'
        }).Count -eq 1
}

foreach ($row in $panel) {
    $dataset = $row.dataset
    $resultCsv = Join-Path $CONFIG.cache ("ablation_{0}.csv" -f $dataset)
    $moduleRoot = Join-Path $CONFIG.cache $dataset
    $moduleRds = @(
        Join-Path $moduleRoot 'no_lr_module\linealinkr_complete.rds'
        Join-Path $moduleRoot 'no_temporal_binning\linealinkr_complete.rds'
        Join-Path $moduleRoot 'no_bifurcation_module\linealinkr_complete.rds'
        Join-Path $moduleRoot 'no_lineage_specific_scoring\linealinkr_complete.rds'
    )
    $stderr = Join-Path $CONFIG.logs ("{0}.stderr.log" -f $dataset)
    $resultValid = Test-AblationResult $dataset $row $resultCsv $moduleRds
    if ($resultValid) { $completed.Add($dataset) }
}

$summary = [ordered]@{
    formal_slots = $panel.Count
    completed = $completed.Count
    pending = $panel.Count - $completed.Count
}
if ($DryRun) {
    $summary | ConvertTo-Json -Compress
    exit 0
}

$statusPath = Join-Path $CONFIG.logs 'panel_status.tsv'
if (-not [System.IO.File]::Exists($statusPath)) {
    [System.IO.File]::WriteAllText($statusPath, "timestamp`tdataset`tstate`texit_code`tstdout`tstderr`n", $utf8NoBom)
}

$pendingRows = @($panel | Where-Object { -not $completed.Contains($_.dataset) })
foreach ($row in $pendingRows) {
    $dataset = $row.dataset
    $freeMemoryGB = [math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB, 2)
    while ($freeMemoryGB -lt $CONFIG.minFreeMemoryGB) {
        Start-Sleep -Seconds 30
        $freeMemoryGB = [math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB, 2)
    }

    $stdout = Join-Path $CONFIG.logs ("{0}.stdout.log" -f $dataset)
    $stderr = Join-Path $CONFIG.logs ("{0}.stderr.log" -f $dataset)
    $processArgs = @{
        FilePath = $CONFIG.rscriptEngine
        ArgumentList = @('--encoding=UTF-8', $CONFIG.ablationScript, $dataset)
        WorkingDirectory = $projectRoot
        RedirectStandardOutput = $stdout
        RedirectStandardError = $stderr
        WindowStyle = 'Hidden'
        Wait = $true
        PassThru = $true
    }
    $process = Start-Process @processArgs

    $resultCsv = Join-Path $CONFIG.cache ("ablation_{0}.csv" -f $dataset)
    $moduleRoot = Join-Path $CONFIG.cache $dataset
    $moduleRds = @(
        Join-Path $moduleRoot 'no_lr_module\linealinkr_complete.rds'
        Join-Path $moduleRoot 'no_temporal_binning\linealinkr_complete.rds'
        Join-Path $moduleRoot 'no_bifurcation_module\linealinkr_complete.rds'
        Join-Path $moduleRoot 'no_lineage_specific_scoring\linealinkr_complete.rds'
    )
    $resultValid = $process.ExitCode -eq 0 -and
        (Test-AblationResult $dataset $row $resultCsv $moduleRds)
    if ($resultValid) { $completed.Add($dataset) }
    $state = if ($resultValid) { 'completed' } else { 'failed' }
    $statusLine = "{0}`t{1}`t{2}`t{3}`t{4}`t{5}`n" -f ([DateTime]::UtcNow.ToString('o')), $dataset, $state, $process.ExitCode, $stdout, $stderr
    [System.IO.File]::AppendAllText($statusPath, $statusLine, $utf8NoBom)
}

$finalCompleted = $completed.Count
$finalSummary = [ordered]@{
    formal_slots = $panel.Count
    completed_files = $finalCompleted
    finished_at_utc = [DateTime]::UtcNow.ToString('o')
}
[System.IO.File]::WriteAllText(
    (Join-Path $CONFIG.logs 'panel_runner_summary.json'),
    ($finalSummary | ConvertTo-Json),
    $utf8NoBom
)
