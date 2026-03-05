@echo off
setlocal enabledelayedexpansion
echo Building Dixon Math DLL and GUI separately...

REM Create build directory
if not exist build mkdir build

REM Clean previous files
if exist dixon.dll del dixon.dll
if exist dixon.lib del dixon.lib
if exist dixon_gui.exe del dixon_gui.exe
if exist build\*.o del build\*.o

REM 1. Compile all math-related source files into object files
echo [1/4] Compiling all math source files...
set MATH_SOURCES=dixon_wrapper.c dixon_complexity.c dixon_flint.c dixon_interface_flint.c dixon_test.c dixon_with_ideal_reduction.c fq_mat_det.c fq_mpoly_mat_det.c fq_multivariate_interpolation.c fq_mvpoly.c fq_nmod_roots.c fq_poly_mat_det.c fq_sparse_interpolation.c fq_unified_interface.c gf2128_mpoly.c gf28_mpoly.c gf2n_field.c gf2n_poly.c polynomial_system_solver.c resultant_with_ideal_reduction.c unified_mpoly_det.c unified_mpoly_interface.c unified_mpoly_resultant.c

set COMPILE_ERROR=0
set OBJECT_FILES=

for %%f in (%MATH_SOURCES%) do (
    if exist src\%%f (
        echo Compiling src\%%f...
        x86_64-w64-mingw32-gcc -c src\%%f -o build\%%~nf.o ^
            -O3 -march=native -static-libgcc -static-libstdc++ ^
            -I./include -I./local/include -DDLL_EXPORT ^
            -fPIC -Wno-format -Wno-unused-function -Wno-unused-variable
        
        if !ERRORLEVEL! neq 0 (
            echo ERROR: Failed to compile src\%%f
            set COMPILE_ERROR=1
        ) else (
            echo OK src\%%f compiled successfully
            set OBJECT_FILES=!OBJECT_FILES! build\%%~nf.o
        )
    ) else (
        echo Warning: src\%%f not found, skipping...
    )
)

if %COMPILE_ERROR% neq 0 (
    echo Some files failed to compile. Stopping.
    pause
    exit /b 1
)

REM Check if any object files were generated
if "!OBJECT_FILES!"=="" (
    echo ERROR: No object files were created
    pause
    exit /b 1
)

REM 2. Link all object files into a DLL
echo [2/4] Creating dixon.dll...
x86_64-w64-mingw32-gcc -shared -o dixon.dll !OBJECT_FILES! ^
    -O3 -march=native -static-libgcc -static-libstdc++ ^
    -L./lib -lflint -lmpfr -lgmp -lm -lpthread -lstdc++ ^
    -Wl,--enable-runtime-pseudo-reloc ^
    -Wl,--out-implib,dixon.lib

if !ERRORLEVEL! neq 0 (
    echo ERROR: Failed to create DLL
    pause
    exit /b 1
)

echo OK dixon.dll created successfully

REM 3. Compile GUI program (no math library dependency)
echo [3/4] Compiling dixon_gui.c (GUI only)...
x86_64-w64-mingw32-gcc -c dixon_gui.c -o build\dixon_gui.o ^
    -O3 -march=native -static-libgcc -static-libstdc++ ^
    -I./include -I./local/include ^
    -mwindows

if !ERRORLEVEL! neq 0 (
    echo ERROR: Failed to compile dixon_gui.c
    pause
    exit /b 1
)

echo OK dixon_gui.c compiled successfully

REM 4. Link GUI program (statically link dixon library)
echo [4/4] Creating dixon_gui.exe...
x86_64-w64-mingw32-gcc -o dixon_gui.exe build\dixon_gui.o ^
    -O3 -march=native -static-libgcc -static-libstdc++ ^
    -mwindows -L./lib -L./ ^
    -lflint -lmpfr -lgmp -lm -lpthread -lstdc++ ^
    -lcomctl32 -lcomdlg32 -luser32 -lgdi32 -lshlwapi -ldixon

if !ERRORLEVEL! neq 0 (
    echo ERROR: Failed to link GUI
    pause
    exit /b 1
)
x86_64-w64-mingw32-gcc -o ./dll/dixon.exe dixon.c -O3 -march=native -static-libgcc -static-libstdc++  -L./lib -L./ -I./include -I./local/include -lflint -lmpfr -lgmp -lm -lpthread -lstdc++  -lcomctl32 -lcomdlg32 -luser32 -lgdi32 -lshlwapi -ldixon
echo OK dixon_gui.exe created successfully

echo.
echo ================================================================
echo BUILD SUCCESSFUL!
echo ================================================================
echo Created files:
echo - dixon.dll       (All FLINT math operations)
echo - dixon.lib       (Import library for DLL)
echo - dixon_gui.exe   (GUI with dixon.lib statically linked)
echo.
echo Directory structure:
echo - Source files (*.c): ./src/
echo - GUI file: ./dixon_gui.c
echo - Headers (dixon): ./include/
echo - Headers (flint/gmp): ./local/include/
echo - Libraries: ./lib/
echo - Object files: ./build/ (preserved)
echo - Output: current directory
echo.
echo Architecture:
echo 1. dixon.dll contains all FLINT math operations
echo 2. dixon_gui.exe has dixon.lib statically linked
echo 3. Object files preserved in build/ for incremental builds
echo.
echo Run dixon_gui.exe to test
echo ================================================================
pause
