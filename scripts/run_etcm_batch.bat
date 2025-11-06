@echo off
REM ETCM 2000 Flavonoids Batch Processing Script (Windows)
REM Complete workflow from Excel ingestion to halogenated product generation

setlocal enabledelayedexpansion

REM Configuration - ADJUST THESE FOR YOUR SETUP
set ETCM_INPUT=flavonoid_smiles_data_ETCM.xlsx
set K2_CONFIG=configs\etcm_flavonoids.yaml
set K3_CONFIG=configs\etcm_flavonoids_k3.yaml
set OUTPUT_BASE=out\etcm2000
set PARENTS_DIR=out\etcm2000\parents
set P1_K2_DIR=out\etcm2000\p1_k2
set P1_K3_DIR=out\etcm2000\p1_k3
set AUDIT_DIR=out\etcm2000\audit
set SAMPLE_SIZE=0

REM QA conflict tolerance settings
set K2_TOL=1
set K3_TOL=2

REM Create directories
if not exist "%OUTPUT_BASE%" mkdir "%OUTPUT_BASE%"
if not exist "%PARENTS_DIR%" mkdir "%PARENTS_DIR%"
if not exist "%P1_K2_DIR%" mkdir "%P1_K2_DIR%"
if not exist "%P1_K3_DIR%" mkdir "%P1_K3_DIR%"
if not exist "%AUDIT_DIR%" mkdir "%AUDIT_DIR%"

echo [INFO] Starting ETCM 2000 flavonoids batch processing workflow...

REM Check pandas/pyarrow availability for statistics
set PANDAS_AVAILABLE=0
python "scripts\batch_stats.py" --check-deps >nul 2>&1
if %errorlevel% equ 0 (
    set PANDAS_AVAILABLE=1
    echo [INFO] Pandas/PyArrow available for statistics
) else (
    echo [WARN] Pandas/PyArrow not available - statistics will be simplified
)

REM Check prerequisites
if not exist "%ETCM_INPUT%" (
    echo [ERROR] ETCM input file not found: %ETCM_INPUT%
    echo Please ensure the ETCM Excel file is in the project root directory.
    exit /b 1
)

if not exist "%K2_CONFIG%" (
    echo [ERROR] k=2 configuration not found: %K2_CONFIG%
    exit /b 1
)

REM Step 1: ETCM Data Ingestion
echo [INFO] Step 1: Converting ETCM Excel to SMILES...

set SAMPLE_ARGS=
if %SAMPLE_SIZE% gtr 0 (
    set SAMPLE_ARGS=--sample %SAMPLE_SIZE%
    echo [WARN] Running in SAMPLE MODE with %SAMPLE_SIZE% molecules
)

python -m halogenator.cli etcm-ingest ^
    -i "%ETCM_INPUT%" ^
    -o "%PARENTS_DIR%\parents.smi" ^
    --out-meta "%PARENTS_DIR%\parents_meta.csv" ^
    %SAMPLE_ARGS%

if %errorlevel% neq 0 (
    echo [ERROR] ETCM ingestion failed
    exit /b 1
)

REM Count molecules
for /f %%i in ('type "%PARENTS_DIR%\parents.smi" ^| find /c /v ""') do set MOLECULE_COUNT=%%i
echo [INFO] Prepared %MOLECULE_COUNT% flavonoid parent molecules

REM Step 2: k=2 Halogenation Enumeration
echo [INFO] Step 2: Running k=2 halogenation enumeration...
echo [INFO] Using parent molecules from: %PARENTS_DIR%\parents.smi

python -m halogenator.cli enum ^
    -c "%K2_CONFIG%" ^
    --parents "%PARENTS_DIR%\parents.smi" ^
    --subset flavonoids ^
    --outdir "%P1_K2_DIR%" ^
    --qa-completion-mode distribute ^
    --qa-conflict-tolerance %K2_TOL% ^
    --qa-max-warnings 500

if %errorlevel% neq 0 (
    echo [ERROR] k=2 enumeration failed
    exit /b 1
)

REM Step 3: Quality checks and reporting
echo [INFO] Step 3: Performing quality checks...

if not exist "%P1_K2_DIR%\products_k2.parquet" (
    echo [ERROR] Expected output file not found: %P1_K2_DIR%\products_k2.parquet
    exit /b 1
)

if not exist "%P1_K2_DIR%\qa_summary.json" (
    echo [WARN] QA summary not found, enumeration may have issues
)

REM Create standardized file names for k=2 products
echo [INFO] Creating standardized product file names...
copy "%P1_K2_DIR%\products_k2.parquet" "%P1_K2_DIR%\products.parquet" >nul

REM Generate CSV version using Python
if %PANDAS_AVAILABLE% equ 1 (
    python -c "import pandas as pd; df = pd.read_parquet('%P1_K2_DIR%/products_k2.parquet'); df.to_csv('%P1_K2_DIR%/products.csv', index=False); print('Created %P1_K2_DIR%/products.csv')"
) else (
    echo [WARN] Pandas not available - skipping CSV generation
)

REM Basic statistics
python "scripts\batch_stats.py" --products "%P1_K2_DIR%\products_k2.parquet" --qa "%P1_K2_DIR%\qa_summary.json" --label "k=2"

REM Step 4: Generate audit sample
echo [INFO] Step 4: Generating audit sample...
python "scripts\batch_stats.py" --products "%P1_K2_DIR%\products_k2.parquet" --audit "%AUDIT_DIR%" --audit-size 100

REM Step 5: Optional k=3 enumeration
echo [INFO] Step 5: Checking k=3 enumeration option...

REM Check for non-interactive mode via environment variable
if not defined RUN_K3 (
    echo Do you want to run k=3 enumeration (tri-substituted products)? Warning: This can be very slow! (y/N):
    set /p RUN_K3="Enter choice: "
    if not defined RUN_K3 set RUN_K3=n
) else (
    echo [INFO] RUN_K3 environment variable set to: %RUN_K3%
)

if /i "%RUN_K3%"=="y" (
    echo [INFO] Running k=3 enumeration (tri-substituted products)...
    echo [INFO] Using parent molecules from: %PARENTS_DIR%\parents.smi

    if not exist "%K3_CONFIG%" (
        echo [ERROR] k=3 configuration not found: %K3_CONFIG%
        echo [WARN] k=3 enumeration failed but continuing with k=2 results
        goto :skip_k3
    )

    python -m halogenator.cli enum ^
        -c "%K3_CONFIG%" ^
        --parents "%PARENTS_DIR%\parents.smi" ^
        --subset flavonoids ^
        --outdir "%P1_K3_DIR%" ^
        --qa-completion-mode distribute ^
        --qa-conflict-tolerance %K3_TOL% ^
        --qa-max-warnings 500

    if %errorlevel% equ 0 (
        echo [INFO] k=3 enumeration completed successfully

        if exist "%P1_K3_DIR%\products_k3.parquet" (
            REM Create standardized file names for k=3 products
            echo [INFO] Creating standardized k=3 product file names...
            copy "%P1_K3_DIR%\products_k3.parquet" "%P1_K3_DIR%\products.parquet" >nul

            REM Generate CSV version using Python
            if %PANDAS_AVAILABLE% equ 1 (
                python -c "import pandas as pd; df = pd.read_parquet('%P1_K3_DIR%/products_k3.parquet'); df.to_csv('%P1_K3_DIR%/products.csv', index=False); print('Created %P1_K3_DIR%/products.csv')"
            ) else (
                echo [WARN] Pandas not available - skipping k=3 CSV generation
            )

            echo [INFO] k=3 statistics:
            python "scripts\batch_stats.py" --products "%P1_K3_DIR%\products_k3.parquet" --qa "%P1_K3_DIR%\qa_summary.json" --label "k=3"
        )
    ) else (
        echo [WARN] k=3 enumeration failed but continuing with k=2 results
    )
) else (
    echo [INFO] Skipping k=3 enumeration (user choice or default)
)

:skip_k3

REM Final Summary
echo.
echo ==================================================================
echo             ETCM 2000 FLAVONOIDS BATCH COMPLETE
echo ==================================================================
echo.
echo Processing Summary:
echo   Input file: %ETCM_INPUT%
echo   Parent molecules: %MOLECULE_COUNT%
echo   k=2 configuration: %K2_CONFIG%

if %SAMPLE_SIZE% gtr 0 (
    echo   MODE: SAMPLE RUN (%SAMPLE_SIZE% molecules)
    echo   To run full batch: set SAMPLE_SIZE=0 and rerun
) else (
    echo   MODE: FULL BATCH
)

echo.
echo Generated Outputs:
echo   Parents:       %PARENTS_DIR%\parents.smi
echo   Parents meta:  %PARENTS_DIR%\parents_meta.csv
echo   k=2 Products:  %P1_K2_DIR%\products_k2.parquet
echo   k=2 CSV:       %P1_K2_DIR%\summary_k2.csv
echo   k=2 QA:        %P1_K2_DIR%\qa_summary.json
if /i "%RUN_K3%"=="y" if exist "%P1_K3_DIR%\products_k3.parquet" (
    echo   k=3 Products:  %P1_K3_DIR%\products_k3.parquet
    echo   k=3 CSV:       %P1_K3_DIR%\summary_k3.csv
    echo   k=3 QA:        %P1_K3_DIR%\qa_summary.json
)
echo   Audit sample:  %AUDIT_DIR%\sample_100.csv
echo.
echo Next Steps:
echo 1. Review QA summaries for data quality issues
echo 2. Load products into your analysis pipeline
echo 3. Consider ADMET/property prediction on the generated library
echo.
echo [INFO] Batch processing workflow complete!

pause