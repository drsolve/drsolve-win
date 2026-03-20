# DixonRes-Windows

Standalone Windows build package for DixonRes, including:

- `dixon.exe`: user-facing CLI launcher
- `dixon_win_gui.exe`: Windows GUI frontend
- `dll/libdixon-1.dll`: reusable Dixon library DLL

## Requirements

- `x86_64-w64-mingw32-gcc`
- GNU `make`

The third-party Windows headers, import libraries, and runtime DLLs used by this package are already included in this directory.

## Build

Run:

```bash
make
```

Optional helper targets:

```bash
make clean
make show-config
make package-layout
```

## Output Layout

After a successful build:

| Location | Description |
| --- | --- |
| `dixon.exe` | CLI launcher shown to end users |
| `dixon_win_gui.exe` | GUI application |
| `bin/dixon_cli_real.exe` | Internal CLI executable |
| `dll/libdixon-1.dll` | Dixon shared library |
| `dll/*.dll` | Runtime DLL dependencies |
| `lib/libdixon-1.a` | Static archive |
| `lib/libdixon-1.dll.a` | Import library for `libdixon-1.dll` |

## How The Windows Packaging Works

`dixon.exe` is a small launcher. It resolves the local `dll/` directory, prepends that directory to `PATH`, and then starts `bin/dixon_cli_real.exe`.

This layout keeps the top-level directory clean while still allowing:

- double-click execution of `dixon.exe`
- GUI execution through `dixon_win_gui.exe`
- reuse of `libdixon-1.dll` by other Windows programs

`dixon_cli_real.exe` links against `libdixon-1.dll` through `lib/libdixon-1.dll.a`, so the solver code is shared through the DLL instead of being duplicated into the CLI executable.

## Rebuilding The Windows Version Later

If you want to rebuild the Windows package in the future:

1. Keep this directory structure unchanged.
2. Make sure the bundled `third_party/` and `runtime/` directories are present.
3. Run `make clean`.
4. Run `make`.

That produces a fresh `dixon.exe`, `dixon_win_gui.exe`, `bin/dixon_cli_real.exe`, `dll/libdixon-1.dll`, and the import/static libraries under `lib/`.

## Using The DLL In Another Program

If another Windows application wants to use `libdixon-1.dll`:

- link against `lib/libdixon-1.dll.a`
- ensure `dll/libdixon-1.dll` and the other runtime DLLs are available at runtime
- either add the `dll/` directory to `PATH` or set the DLL search path in the host program

## Repository Contents

This package contains:

- project sources in `dixon.c`, `include/`, `src/`, and `win_gui/`
- bundled Windows-facing third-party headers and link libraries in `third_party/`
- runtime DLL files in `runtime/`

## License

DixonRes is distributed under the GNU General Public License v2.0 or later. See `COPYING` for details.
