# Contributing Guidelines

## Code Style and Standards

### ASCII-only Requirement

ASCII-only applies to source code and repository scripts, not to user interface text or commit messages.

This means:
- All Python source files (`.py`) must contain only ASCII characters
- Shell scripts in the `scripts/` directory must contain only ASCII characters  
- Configuration files and repository metadata should use ASCII encoding

However, the following are exempt from ASCII-only requirements:
- User-facing text and messages displayed to users
- Commit messages and PR descriptions
- Documentation that includes scientific notation or chemical formulas
- Data files containing chemical structures or molecular representations