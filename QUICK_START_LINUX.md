# Quick Start on Linux Server

## Step 1: Pull Latest Changes

```bash
cd /path/to/halogenator
git pull origin main
```

If the environment files are not in main yet, pull from the feature branch:

```bash
git pull origin fix/pr2-contract-and-sugar-accept
```

## Step 2: Setup Environment

```bash
# Create conda environment
conda env create -f environment.yml

# Activate it
conda activate halogenator

# Verify
python -c "from rdkit import Chem; import pyarrow; print('Ready!')"
```

## Step 3: Check Data Location

The shuffled chunks should be at:
```bash
ls data/output/transforms/polyphenol-2X_SHUFFLED_500K/
```

You should see:
- 28 chunk files (chunk_000 to chunk_027)
- Each ~94-95 MB
- Plus chunks_metadata.json

## Step 4: Test Run One Chunk

```bash
# Test chunk 0 (should take ~3-4 hours)
python scripts/08_transform_library_v2.py apply \
  --input data/output/transforms/polyphenol-2X_SHUFFLED_500K/chunk_000_input.parquet \
  --outdir data/output/transforms/TEST_CHUNK_000 \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 16 \
  --max-in-flight 6 \
  --batch-size 50000 \
  --use-bloom-filter \
  --bloom-expected-items 5000000
```

## Step 5: Monitor Progress

```bash
# Check log file
tail -f data/output/transforms/TEST_CHUNK_000/transform.log

# Check memory usage
free -h

# Check CPU usage
htop
```

## Expected Results

- Processing time: 3-4 hours per chunk
- Memory usage: Should stay under 70%
- All 28 chunks should complete without timeout
- Total time: ~100 hours for all chunks

## Troubleshooting

See `LINUX_SETUP.md` for detailed troubleshooting steps.

## Next Steps

After successful test run:
1. Set up batch processing for all 28 chunks
2. Monitor first few chunks to confirm timing
3. Let it run!

The shuffle-based rechunking has reduced variance by 94.7%, so all chunks should have similar processing times.
