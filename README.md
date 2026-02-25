# Dixon Resultant & Polynomial System Solver - Windows Version

A Windows GUI application for computing Dixon resultants and solving polynomial systems over finite fields using FLINT library.

## Features

- **Modern Windows GUI**: User-friendly tabbed interface with progress indicators
- **Three Computation Modes**:
  - Polynomial System Solver (n×n systems)
  - Dixon Resultant (basic variable elimination)
  - Dixon with Ideal Reduction (triangular ideal reduction)
- **Finite Field Support**: Prime fields and extension fields (F_p, F_{p^k})
- **File I/O**: Load input from files and save results
- **Multi-threaded**: Non-blocking computation with progress feedback
- **DLL Architecture**: Separated math library (dixon.dll) and GUI

## System Requirements

- Windows 7 or later (64-bit)
- MinGW-w64 cross-compiler (for building)
- FLINT, GMP, MPFR libraries (DLLs must be accessible)

## Pre-built Binaries

If you have pre-built binaries, you need:
- `dixon_gui.exe` - Main GUI application
- `dixon.dll` - Math computation library
- `dixon.exe` - Command-line version (optional)
- FLINT/GMP/MPFR DLLs in the same directory or in `dll/` subdirectory
- You may still need libopenblas.dll, which is not included here

## Building from Source

### Build Process

Run the provided build script:
```cmd
build_gui.bat
```

This will:
1. Compile all math source files into object files
2. Link them into `dixon.dll` (math library)
3. Compile the GUI application
4. Link `dixon_gui.exe` with static dixon library
5. Build `dixon.exe` (command-line version)

### Build Output

```
dixon.dll         - Math computation DLL (contains all FLINT operations)
dixon.lib         - Import library for the DLL
dixon_gui.exe     - Windows GUI application
dixon.exe         - Command-line interface (in dll/ directory)
build/*.o         - Object files (preserved for incremental builds)
```

## Directory Structure

```
.
├── dixon_gui.exe           # GUI application
├── dixon.exe               # CLI application (in dll/)
├── dixon.dll               # Math library
├── dixon.lib               # Import library
├── src/                    # Source files
│   ├── dixon_wrapper.c     # DLL wrapper functions
│   ├── dixon_flint.c
│   ├── polynomial_system_solver.c
│   └── ...
├── include/                # Header files
│   └── ...
├── lib/                    # External libraries (FLINT, GMP, MPFR)
├── local/include/          # External headers
├── dll/                    # Runtime DLLs directory
│   └── ...
└── build/                  # Build artifacts
    └── *.o
```

## Using the GUI

### 1. Polynomial System Solver Tab

Solve n×n polynomial systems (n equations in n variables).

**Example:**
- Polynomials: `x + y - 3, 2*x - y`
- Field size: `257`
- Click "Compute"

The solver automatically detects all variables and finds all solutions.

### 2. Dixon Resultant Tab

Eliminate variables using Dixon resultant method.

**Example:**
- Polynomials: `x + y + z, x*y + y*z + z*x, x*y*z + 1`
- Variables to eliminate: `x, y`
- Field size: `257`

**Important**: Number of variables to eliminate must be (number of equations - 1).

### 3. Dixon with Ideal Reduction Tab

Compute resultants with triangular ideal reduction.

**Example:**
- Polynomials: `a1^2 + a2^2 + a3^2 + a4^2 - 10, a4^3 - a1 - a2*a3 - 5`
- Variables to eliminate: `a4`
- Ideal generators: `a2^3 = 2*a1 + 1, a3^3 = a1*a2 + 3, a4^3 = a1 + a2*a3 + 5`
- Field size: `257`

### Field Size Formats

All formats are supported:
- `257` - Prime field F_257
- `256` - Automatically detected as F_{2^8}
- `2^8` - Extension field F_{2^8}
- `3^5` - Extension field F_{3^5}

### Loading Files

Click "Load" to load input files. File format:
```
257
x + y - 3
2*x - y
```

For Dixon mode, add a third line with variables to eliminate.

### Saving Results

Click "Save" after computation to export results to a text file.

## Command-Line Interface

The package also includes `dixon.exe` for command-line usage:

```cmd
# Polynomial solver
dixon.exe --solve "x + y - 3, 2*x - y" 257

# Basic Dixon
dixon.exe "x + y + z, x*y + y*z + z*x, x*y*z + 1" "x, y" 257

# File input
dixon.exe input.dat
```

See the main README for complete CLI documentation.

## DLL Architecture

The application uses a clean separation:

1. **dixon.dll**: Contains all FLINT math operations
   - No Windows GUI code
   - Pure computational functions
   - Can be used by other applications

2. **dixon_gui.exe**: Windows GUI application
   - Statically links dixon.lib
   - No FLINT headers in GUI code
   - Calls wrapped functions from the DLL

This architecture allows:
- Independent updates to math library
- Reuse of math library in other projects
- Clean separation of concerns
- Easier debugging

## Troubleshooting

### DLL Not Found Error

If you see "dixon.dll not found":
1. Ensure `dixon.dll` is in the same directory as `dixon_gui.exe`
2. Check that FLINT/GMP/MPFR DLLs are accessible
3. Try placing all DLLs in the `dll/` subdirectory

The application will search:
- Current directory
- `dll/` subdirectory
- System PATH

### Function Import Error

If you see "Failed to get function addresses":
1. Rebuild `dixon.dll` with `build_gui.bat`
2. Ensure the DLL export definitions are correct
3. Check for version mismatch between DLL and GUI

### Computation Hangs

If computation appears frozen:
- The progress bar shows activity
- Large systems may take significant time
- You can close the application (with confirmation)

### Field Size Invalid

If field size validation fails:
- Use supported formats: `257`, `2^8`, `256`
- Prime must be less than 2^63
- For extension fields, base must be prime

## Performance Notes

- Extension fields are slower than prime fields
- Parallel computation is enabled (OpenMP)
- Progress bar shows activity during computation
- Results appear in the text area when complete

## Building Tips

### Incremental Builds

Object files are preserved in `build/`. To rebuild only modified files:
1. Manually compile changed `.c` files
2. Relink the DLL or executable

### Static vs Dynamic Linking

The build script creates:
- `dixon.dll` - Dynamic library (shared)
- `dixon.lib` - Import library (for static linking)
- `dixon_gui.exe` - GUI with statically linked dixon.lib

### Debug Build

For debugging, modify `build_gui.bat`:
```batch
REM Add -g flag for debug symbols
x86_64-w64-mingw32-gcc -c -g src\file.c ...
```

## Additional Resources

- Main documentation: See README.md
- FLINT documentation: https://flintlib.org/

## License

DixonRes is distributed under the GNU General Public License version 3.0. See the file COPYING.
