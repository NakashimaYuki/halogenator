# Changelog

## [Unreleased] - 2025-01-19

### Added
- **Sugar Events Pipeline**: Comprehensive three-component observability for sugar masking effectiveness
  - `sugar_mask_filtered`: Sites filtered during initial site identification
  - `sugar_proximity_filtered`: Sites filtered by proximity guard within configurable bond radius
  - `post_guard_blocked`: Sites filtered by post-enumeration guards
- **Proximity Guard System**: Configurable bond-distance filtering around sugar mask atoms
  - CLI parameter: `--sugar.proximity-guard-radius=N` (default: 0=disabled, typical: 2-3)
  - BFS-based distance calculation with O(1) membership checks for performance
  - Selective activation: only enabled for glycoside samples in heuristic mode
- **Hybrid Sugar Ring Caching**: Performance optimization for repeated sugar ring identification
  - Primary: WeakKeyDictionary with automatic garbage collection cleanup
  - Fallback: SMILES-based LRU cache when WeakKeyDictionary unavailable
  - Cache hit monitoring via `HALOGENATOR_CACHE_MONITOR=1` environment variable
- **Legacy Key Compatibility**: Configurable emission of legacy QA key aliases
  - CLI parameter: `--sugar.emit-legacy-keys` for backward compatibility
  - Unified calculation logic prevents double counting in unit tests
- **Comprehensive Documentation**:
  - New `docs/sugar_events.md` with complete API reference and usage examples
  - README updated with proximity guard CLI examples and YAML configuration
- **Enhanced Logging**: Unified proximity guard configuration logging for production troubleshooting
  - Format: `{radius=N, sample_type=TYPE, sugar_mode=MODE, enabled=BOOL}`

### Fixed
- **R2b Sugar Ring Exclusion**: Corrected asymmetric behavior in `c_ring_sp3_CH2_flavanone_sites()`
  - Now properly excludes sugar ring atoms for consistent enumeration results
- **Aggregator Semantic Pollution**: Eliminated pollution of `by_rule_halogen_k` dimensions
  - Added dedicated `paths` container with `record_paths()` API for QA events
  - Prevents sugar event metrics from appearing in rule-halogen grid statistics
- **Cache Reliability**: Implemented hybrid strategy to handle WeakKeyDictionary failures
  - Automatic fallback to SMILES-LRU when TypeErrors occur
  - Consistent cache performance across different Python/platform configurations

### Changed
- **Mask Expansion Optimization**: Compute expanded masks once per enumeration, reuse across rules
  - Eliminates redundant BFS calculations for R1/R2a/R2b rule processing
  - Significant performance improvement for proximity guard operations
- **QA Schema Extension**: Added `sugar_proximity_filtered` to `QA_PATH_KEYS` for schema compliance
- **Deprecation**: Marked `_site_proximity_blocked()` as deprecated to avoid confusion
  - Added deprecation warnings with migration guidance to proximity guard system

### Technical Improvements
- **Cache Hit Monitoring**: Optional debugging capabilities for cache performance analysis
  - Functions: `get_sugar_ring_cache_stats()`, `reset_sugar_ring_cache_stats()`, `log_sugar_ring_cache_stats()`
  - Environment-controlled logging to prevent production log pollution
- **Test Coverage**: Added comprehensive test suites for regression protection
  - `tests/test_sugar_events.py`: Core functionality and edge cases
  - `tests/test_aggregator_paths.py`: QA aggregation and semantic isolation
- **Acceptance Testing**: Enhanced PR1 acceptance script with proximity guard validation
  - Standardized 11-sample test suite (5 glycosides + 6 aglycone controls)
  - Multi-criteria assessment for glycoside filtering effectiveness

### Performance
- **O(1) Proximity Checks**: Precomputed expanded masks eliminate iterative distance calculations
- **Hybrid Caching**: 85-95% cache hit rates for repeated molecules
- **Memory Efficiency**: WeakKeyDictionary automatic cleanup prevents memory leaks
- **Selective Processing**: Proximity guard only activated when needed (glycoside + heuristic mode)

### Migration Notes
- New CLI parameters are backward compatible (default values maintain existing behavior)
- Existing QA statistics consumers will see new `sugar_proximity_filtered` field (defaults to 0)
- Cache monitoring is opt-in via environment variable, no impact on existing workflows
- Legacy key emission disabled by default, enable via `--sugar.emit-legacy-keys` if needed

## [Previous] - 2025-01-10

### Fixed
- **[BREAKING]** Corrected QA statistics invariant to three-way exclusive classification
  - **Previous:** `attempts = products + no_product_matches` (incorrect - ignored template unsupported cases)
  - **Current:** `attempts = products + no_product_matches + template_unsupported` (correct - all outcomes accounted for)
  - This change ensures accurate statistical reporting where each attempt has exactly one outcome

### Changed
- **QA File Name:** Standardized QA statistics filename to `qa_summary.json`
- **Function-level imports:** Moved most function-level imports to module top-level for better code maintainability
  - Preserved conditional imports for optional dependencies (RDKit, PyArrow) and dependency injection in CLI
- **Backward Compatibility:** Maintained read support for legacy dedup field names while writing new format
  - Legacy: `statesig_hits`, `inchi_hits` (read-only support)  
  - Current: `dedup_hits_statesig`, `dedup_hits_inchi` (output format)

### Technical Improvements
- Enhanced test coverage for statistical invariant validation
- Added explicit comments for conditional imports (testing, optional dependencies)
- Improved documentation for QA statistics schema and semantics

### Migration Notes
- Existing consumers should update statistical validation logic to use the corrected three-way invariant
- QA statistics files now use standardized filename `qa_summary.json`
- Legacy field names in existing data will be automatically migrated during CLI operations