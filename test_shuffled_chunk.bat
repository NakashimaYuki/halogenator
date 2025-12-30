@echo off
REM Test run a single shuffled chunk to validate processing time

echo ================================================================================
echo TESTING SHUFFLED CHUNK
echo ================================================================================
echo.
echo This will test Chunk 8 (lowest complexity: 163.2)
echo Expected time: ~6-7 hours based on complexity analysis
echo.
echo Starting at: %TIME%
echo.

python scripts/08_transform_library_v2.py apply ^
  --input data/output/transforms/polyphenol-2X_SHUFFLED/chunk_008_input.parquet ^
  --outdir data/output/transforms/SHUFFLED_TEST_CHUNK8 ^
  --xf-config configs/transforms.yaml ^
  --xf-name FG_PHENOL_OH__OH__TO__OMe ^
  --workers 16 ^
  --max-in-flight 6 ^
  --batch-size 50000 ^
  --use-bloom-filter ^
  --bloom-expected-items 5000000

echo.
echo Completed at: %TIME%
echo.
echo Check results in: data/output/transforms/SHUFFLED_TEST_CHUNK8/
pause
