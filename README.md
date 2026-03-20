# DixonRes_Win

Standalone Windows build package for DixonRes.

## Prerequisites

- `x86_64-w64-mingw32-gcc` installed on your system

## Quick Start

Simply run `make` in this directory to build everything:

```bash
make
```

## Build Outputs

After a successful build, you'll get:

| Location | Files |
|----------|-------|
| Root directory | `dixon.exe`, `dixon_win_gui.exe` |
| `bin/` | `dixon_cli_real.exe` |
| `dll/` | `libdixon-1.dll` and all runtime DLLs |
| `lib/` | `libdixon-1.a`, `libdixon-1.dll.a` |

## Architecture

### Launcher Mechanism

- **`dixon.exe`** is a thin launcher that:
  1. Adds the `dll/` directory to `PATH`
  2. Launches `bin/dixon_cli_real.exe`

This design ensures that DLL dependencies in `dll/` are properly resolved when users double-click `dixon.exe`. Windows doesn't automatically search subdirectories for DLLs before process startup, so this launcher layer is necessary for a seamless user experience.

- **`dixon_win_gui.exe`** internally calls `dixon.exe` from the same directory, inheriting the same DLL search path behavior.

### Using libdixon-1.dll in Your Own Applications

If you want to create your own Windows applications that use `libdixon-1.dll`:
- Add the `dixon_win/dll` directory to your `PATH` before running your application
- Or explicitly set the DLL search directory in your host program

## Included Components

This directory already contains:

- Project sources: `dixon.c`, `include/`, `src/`, `win_gui/`
- PML headers and libraries: `third_party/pml/`
- MinGW includes for Windows: `third_party/mingw/include/`
- Windows import libraries for FLINT/GMP/MPFR/OpenBLAS: `third_party/flint/`
- Runtime DLL sources: `runtime/`

## Available Make Commands

| Command | Description |
|---------|-------------|
| `make` | Build everything |
| `make clean` | Remove all build artifacts |
| `make show-config` | Display build configuration |
| `make package-layout` | Show output directory structure |

## Directory Structure

```
dixon_win/
├── dixon.exe                 # Launcher (sets DLL path)
├── dixon_win_gui.exe         # Windows GUI application
├── bin/
│   └── dixon_cli_real.exe    # Actual CLI program
├── dll/
│   └── *.dll                 # All runtime DLLs
├── lib/
│   ├── libdixon-1.a          # Static library
│   └── libdixon-1.dll.a      # DLL import library
├── include/                  # Header files
├── src/                      # Source files
├── win_gui/                  # GUI source files
├── third_party/              # Third-party dependencies
│   ├── pml/                  # PML headers and libs
│   ├── mingw/                # MinGW includes
│   └── flint/                # FLINT import libs
└── runtime/                  # Runtime DLL sources
```

## License

DixonRes is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.
```
