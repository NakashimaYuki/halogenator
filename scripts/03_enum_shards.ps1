# ==============================================================================
# Parallel Flavonoid Shard Enumeration (Windows PowerShell)
# ==============================================================================
#
# Purpose: Run halogenation enumeration on all shard files in parallel
#
# Features:
# - Controlled parallelism using PowerShell jobs
# - Individual log files per shard
# - Progress tracking with progress bar
# - Error handling and resume capability
#
# Usage:
#   powershell -ExecutionPolicy Bypass -File scripts/03_enum_shards.ps1 -MaxParallel 4
#
# Parameters:
#   -MaxParallel: Number of parallel jobs (default: 4)
#   -RdkitThreads: RDKit threads per job (default: 4)
#
# ==============================================================================

param(
    [int]$MaxParallel = 4,
    [int]$RdkitThreads = 4,
    [switch]$Resume
)

# ==============================================================================
# Configuration
# ==============================================================================

$ShardDir = "data\work\shards"
$OutputRoot = "data\output\cnpd_flav_k2"
$Config = "configs\flavonoids_k2_prod.yaml"
$LogDir = "data\output\cnpd_flav_k2_logs"

# ==============================================================================
# Display Configuration
# ==============================================================================

Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host "Flavonoid Shard Enumeration - Parallel Execution (Windows)" -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host ""

# ==============================================================================
# Validation
# ==============================================================================

# Check shard directory
if (-not (Test-Path $ShardDir)) {
    Write-Host "ERROR: Shard directory not found: $ShardDir" -ForegroundColor Red
    Write-Host "Please run scripts/02_make_shards.py first" -ForegroundColor Yellow
    exit 1
}

# Check config file
if (-not (Test-Path $Config)) {
    Write-Host "ERROR: Configuration file not found: $Config" -ForegroundColor Red
    exit 1
}

# Get shard files
$ShardFiles = Get-ChildItem -Path $ShardDir -Filter "flav_shard_*.smi"
$NumShards = $ShardFiles.Count

if ($NumShards -eq 0) {
    Write-Host "ERROR: No shard files found in $ShardDir" -ForegroundColor Red
    exit 1
}

Write-Host "Configuration:"
Write-Host "  Shard directory: $ShardDir"
Write-Host "  Number of shards: $NumShards"
Write-Host "  Output root: $OutputRoot"
Write-Host "  Config file: $Config"
Write-Host "  Max parallel jobs: $MaxParallel"
Write-Host "  RDKit threads per job: $RdkitThreads"
Write-Host ""

# Create output directories
New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

# ==============================================================================
# Progress Tracking
# ==============================================================================

$ProgressFile = Join-Path $LogDir "progress.txt"
$CompletedFile = Join-Path $LogDir "completed_shards.txt"
$FailedFile = Join-Path $LogDir "failed_shards.txt"

# Initialize progress files
New-Item -ItemType File -Force -Path $ProgressFile | Out-Null

if (-not (Test-Path $CompletedFile)) {
    New-Item -ItemType File -Force -Path $CompletedFile | Out-Null
}

if (-not (Test-Path $FailedFile)) {
    New-Item -ItemType File -Force -Path $FailedFile | Out-Null
}

# Load completed shards (for resume capability)
$CompletedShards = @{}
if (Test-Path $CompletedFile) {
    Get-Content $CompletedFile | ForEach-Object {
        if ($_ -ne "") {
            $CompletedShards[$_] = $true
        }
    }
}

Write-Host "Resume mode: $(if ($Resume) { 'ENABLED' } else { 'DISABLED' })"
Write-Host "Already completed: $($CompletedShards.Count) shards"
Write-Host ""

# ==============================================================================
# Enumeration Script Block
# ==============================================================================

$EnumerationScriptBlock = {
    param($SmiFile, $OutputRoot, $Config, $LogDir, $RdkitThreads)

    $ShardName = [System.IO.Path]::GetFileNameWithoutExtension($SmiFile)
    $OutDir = Join-Path $OutputRoot $ShardName
    $LogFile = Join-Path $LogDir "$ShardName.log"

    # Create output directory
    New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

    # Build command
    $cmd = "python"
    $args = @(
        "-m", "halogenator.cli", "enum",
        "-c", $Config,
        "--rdkit-threads", $RdkitThreads,
        "--outdir", $OutDir,
        "-i", $SmiFile
    )

    # Run enumeration
    try {
        $process = Start-Process -FilePath $cmd -ArgumentList $args `
            -RedirectStandardOutput $LogFile `
            -RedirectStandardError "$LogFile.err" `
            -NoNewWindow -Wait -PassThru

        $exitCode = $process.ExitCode

        # Combine stdout and stderr
        if (Test-Path "$LogFile.err") {
            Get-Content "$LogFile.err" | Add-Content $LogFile
            Remove-Item "$LogFile.err"
        }

        return @{
            ShardName = $ShardName
            Success = ($exitCode -eq 0)
            LogFile = $LogFile
        }
    }
    catch {
        return @{
            ShardName = $ShardName
            Success = $false
            LogFile = $LogFile
            Error = $_.Exception.Message
        }
    }
}

# ==============================================================================
# Parallel Execution
# ==============================================================================

Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host "Starting parallel enumeration..." -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host ""

$StartTime = Get-Date
$Jobs = @()
$ProcessedCount = 0
$SuccessCount = 0
$FailCount = 0

# Filter shards to process
$ShardsToProcess = $ShardFiles | Where-Object {
    $shardName = [System.IO.Path]::GetFileNameWithoutExtension($_.Name)
    -not $CompletedShards.ContainsKey($shardName) -or -not $Resume
}

$TotalToProcess = $ShardsToProcess.Count
Write-Host "Shards to process: $TotalToProcess (skipping $($NumShards - $TotalToProcess) already completed)"
Write-Host ""

foreach ($SmiFile in $ShardsToProcess) {
    $ShardName = [System.IO.Path]::GetFileNameWithoutExtension($SmiFile.Name)

    # Wait if max jobs reached
    while ((Get-Job -State Running).Count -ge $MaxParallel) {
        Start-Sleep -Milliseconds 500

        # Check completed jobs
        $CompletedJobs = Get-Job -State Completed
        foreach ($job in $CompletedJobs) {
            $result = Receive-Job -Job $job
            Remove-Job -Job $job

            $ProcessedCount++

            if ($result.Success) {
                $SuccessCount++
                Write-Host "[DONE] $($result.ShardName) ($ProcessedCount/$TotalToProcess)" -ForegroundColor Green

                # Log completion
                Add-Content -Path $CompletedFile -Value $result.ShardName
                $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
                Add-Content -Path $ProgressFile -Value "$timestamp SUCCESS $($result.ShardName)"
            }
            else {
                $FailCount++
                Write-Host "[FAIL] $($result.ShardName) ($ProcessedCount/$TotalToProcess)" -ForegroundColor Red
                Write-Host "       See: $($result.LogFile)" -ForegroundColor Yellow

                # Log failure
                Add-Content -Path $FailedFile -Value $result.ShardName
                $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
                Add-Content -Path $ProgressFile -Value "$timestamp FAILED $($result.ShardName)"
            }

            # Update progress bar
            $percentComplete = ($ProcessedCount / $TotalToProcess) * 100
            Write-Progress -Activity "Enumerating shards" `
                -Status "$ProcessedCount/$TotalToProcess complete (Success: $SuccessCount, Failed: $FailCount)" `
                -PercentComplete $percentComplete
        }
    }

    # Start new job
    Write-Host "[START] $ShardName" -ForegroundColor Cyan
    $job = Start-Job -ScriptBlock $EnumerationScriptBlock `
        -ArgumentList $SmiFile.FullName, $OutputRoot, $Config, $LogDir, $RdkitThreads
    $Jobs += $job
}

# Wait for remaining jobs
Write-Host ""
Write-Host "Waiting for remaining jobs to complete..." -ForegroundColor Yellow

while ((Get-Job -State Running).Count -gt 0) {
    Start-Sleep -Milliseconds 500

    # Process completed jobs
    $CompletedJobs = Get-Job -State Completed
    foreach ($job in $CompletedJobs) {
        $result = Receive-Job -Job $job
        Remove-Job -Job $job

        $ProcessedCount++

        if ($result.Success) {
            $SuccessCount++
            Write-Host "[DONE] $($result.ShardName) ($ProcessedCount/$TotalToProcess)" -ForegroundColor Green

            Add-Content -Path $CompletedFile -Value $result.ShardName
            $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
            Add-Content -Path $ProgressFile -Value "$timestamp SUCCESS $($result.ShardName)"
        }
        else {
            $FailCount++
            Write-Host "[FAIL] $($result.ShardName) ($ProcessedCount/$TotalToProcess)" -ForegroundColor Red
            Write-Host "       See: $($result.LogFile)" -ForegroundColor Yellow

            Add-Content -Path $FailedFile -Value $result.ShardName
            $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
            Add-Content -Path $ProgressFile -Value "$timestamp FAILED $($result.ShardName)"
        }

        $percentComplete = ($ProcessedCount / $TotalToProcess) * 100
        Write-Progress -Activity "Enumerating shards" `
            -Status "$ProcessedCount/$TotalToProcess complete (Success: $SuccessCount, Failed: $FailCount)" `
            -PercentComplete $percentComplete
    }
}

# Clean up any remaining jobs
Get-Job | Remove-Job -Force

Write-Progress -Activity "Enumerating shards" -Completed

$EndTime = Get-Date
$Elapsed = $EndTime - $StartTime

# ==============================================================================
# Summary
# ==============================================================================

Write-Host ""
Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host "Enumeration Complete" -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host ""

Write-Host "Results:"
Write-Host "  Total shards: $NumShards"
Write-Host "  Processed: $ProcessedCount"
Write-Host "  Completed: $SuccessCount" -ForegroundColor $(if ($SuccessCount -eq $ProcessedCount) { "Green" } else { "Yellow" })
Write-Host "  Failed: $FailCount" -ForegroundColor $(if ($FailCount -eq 0) { "Green" } else { "Red" })
Write-Host ""
Write-Host "Elapsed time: $($Elapsed.ToString('hh\:mm\:ss'))"
Write-Host ""

if ($FailCount -gt 0) {
    Write-Host "Failed shards:" -ForegroundColor Red
    Get-Content $FailedFile
    Write-Host ""
    Write-Host "Review log files in $LogDir for error details" -ForegroundColor Yellow
    exit 1
}

Write-Host "All shards completed successfully!" -ForegroundColor Green
Write-Host ""
Write-Host "Next step: Merge shards and run QC" -ForegroundColor Cyan
Write-Host "  python scripts/04_merge_and_qc.py" -ForegroundColor White
Write-Host ""
