# CI Improvements and Pre-Push Workflow

This document outlines the comprehensive CI layering and pre-push validation setup for the halogenator project.

## CI Pipeline Layers

### Layer 1: Smoke Tests (Fast Feedback - 2-5 minutes)
**File**: `.github/workflows/smoke-tests.yml`

- **Basic Smoke Tests**: Core imports and CLI without RDKit dependencies  
- **RDKit Smoke Tests**: Essential RDKit-dependent functionality
- **Performance Smoke Tests**: Benchmark and comparison functionality
- **Integration Smoke Tests**: P0/P1 workflow validation

### Layer 2: Unit & Integration Tests (Medium - 8-15 minutes)  
**File**: `.github/workflows/ci.yml`

- **Unit Tests**: Fast unit test suite with comprehensive coverage
- **Integration Tests**: Cross-component integration validation
- **P0 Demo Workflow**: Full end-to-end P0 demonstration with artifacts

### Layer 3: P1 Extended Tests (Longer - 10-30 minutes)
**File**: `.github/workflows/p1-tests.yml`  

- **P1 Smoke Tests**: Extended P1 functionality validation
- **Performance Benchmarks**: k=2 baseline with 10 flavonoids
- **Constraint Validation**: Complex constraint handling

## Pre-Push Validation

### Automated Pre-Commit Hooks
**File**: `.pre-commit-config.yaml`

```bash
# Install pre-commit hooks
pip install pre-commit
pre-commit install
```

Hooks include:
- **ASCII Compliance**: Ensures all files use ASCII encoding
- **Python Import Check**: Validates core module imports  
- **Schema Validation**: Tests schema record creation

### Pre-Push Script
**File**: `scripts/pre-push.sh`

```bash
# Run comprehensive pre-push validation
./scripts/pre-push.sh
```

The pre-push script performs 5 essential checks:

1. **ASCII Compliance** - Validates encoding standards
2. **Basic Imports** - Tests core module imports
3. **Schema Validation** - Validates data schemas  
4. **Core Unit Tests** - Runs critical unit tests
5. **P0 Smoke Test** - Single molecule enumeration test

### Git Hook Integration (Optional)

To automatically run pre-push validation:

```bash
# Install as git pre-push hook
ln -sf ../../scripts/pre-push.sh .git/hooks/pre-push
chmod +x .git/hooks/pre-push
```

## CI Strategy Benefits

### 1. Fast Feedback Loop
- Smoke tests provide immediate feedback (2-5 minutes)
- Developers get quick validation on basic functionality
- Reduced wait time for common issues

### 2. Layered Validation
- Progressive complexity: smoke -> unit -> integration -> P1
- Each layer catches different types of issues
- Resource-efficient execution

### 3. Local Validation
- Pre-commit hooks catch issues before commit
- Pre-push script prevents broken pushes
- Reduces CI queue contention

### 4. Comprehensive Coverage
- Multiple environments (with/without RDKit, different Python versions)
- Both unit and integration testing
- Performance validation with real datasets

## Usage Guidelines

### For Developers
1. **Before committing**: Pre-commit hooks run automatically
2. **Before pushing**: Run `./scripts/pre-push.sh` manually
3. **Monitor CI**: Check all layers pass before merging

### For CI Maintenance  
1. **Keep smoke tests fast**: < 5 minutes total
2. **Maintain layer separation**: Don't duplicate tests across layers
3. **Update pre-push script**: Keep aligned with smoke tests

### For Releases
1. **All layers must pass**: No exceptions for releases
2. **P1 benchmarks**: Must meet performance targets
3. **Artifact validation**: Verify P0 demo outputs

## Troubleshooting

### Common Pre-Push Failures
- **ASCII compliance**: Use `scripts/check_ascii.sh` to identify files
- **Import errors**: Check for missing dependencies or circular imports
- **Schema validation**: Verify schema definitions match usage
- **Unit test failures**: Run `python -m pytest tests/` for details

### CI Layer Failures  
- **Smoke test failures**: Usually import/dependency issues
- **Unit test failures**: Logic or regression issues
- **Integration failures**: Cross-component compatibility
- **P1 test failures**: Performance or constraint violations

### Performance Targets
- **P0 operations**: < 1 second per molecule (k=1)
- **P1 baseline**: 10-15 minutes for 10 flavonoids (k=2) 
- **Memory usage**: <= 8GB peak for P1 baseline
- **CPU utilization**: Efficient on 4-core systems