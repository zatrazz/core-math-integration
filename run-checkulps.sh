#!/bin/bash
#
# Wrapper around the checkulps tool.
#
#  * Builds the project on demand (into a caller-supplied directory or a
#    throw-away temporary one created per invocation).
#  * Runs checkulps either against the system libc or against a given glibc
#    build directory.  In the latter case the tool is launched with the glibc
#    build's own dynamic loader and a --library-path that also resolves the
#    libstdc++ and libgomp used by checkulps itself.
#  * Optionally checks the floating-point exceptions (-e) and/or errno (-E)
#    semantics on top of the plain ULP check.

set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage () {
  cat <<EOF
Usage: ${0##*/} [OPTIONS] [FUNCTION...]

Run checkulps for the given FUNCTIONs (names of <desc-dir>/<name>.json files).
With no FUNCTION, every description file found in the description directory is
checked.

Options:
  -b, --build-dir DIR   Build (or reuse a build) in DIR.  When omitted a
                        temporary directory is created for this invocation and
                        removed on exit (see --keep).
  -g, --glibc-dir DIR   glibc build directory to test.  checkulps is run under
                        that build's loader and --library-path.  When omitted
                        the system libc is tested.
  -d, --desc-dir DIR    Directory with the JSON description files
                        (default: ${SCRIPT_DIR}/src/checkulps/description).
  -o, --output-dir DIR  Write per-function results to DIR/<name>.out instead of
                        streaming everything to stdout.
  -e, --exceptions      Also check raised floating-point exceptions.
  -E, --errno           Also check errno (EDOM/ERANGE).
  -f, --full            Shorthand for -e -E.
  -j, --jobs N          Parallel build jobs (default: nproc).
  -k, --keep            Do not delete an auto-created temporary build dir.
      --clang           Configure the build with clang/clang++.
  -h, --help            This help.

Anything after a literal -- is forwarded verbatim to checkulps.
EOF
  exit "${1:-0}"
}

die () { echo "${0##*/}: error: $*" >&2; exit 1; }

# --- argument parsing -------------------------------------------------------

BUILD_DIR=""
GLIBC_DIR=""
DESC_DIR="${SCRIPT_DIR}/src/checkulps/description"
OUTPUT_DIR=""
CHECK_EXC=0
CHECK_ERRNO=0
JOBS="$(nproc 2>/dev/null || echo 1)"
KEEP=0
USE_CLANG=0
declare -a FUNCS=()
declare -a EXTRA=()

while [ $# -gt 0 ]; do
  case "$1" in
    -b|--build-dir)  BUILD_DIR="${2:?}"; shift 2 ;;
    -g|--glibc-dir)  GLIBC_DIR="${2:?}"; shift 2 ;;
    -d|--desc-dir)   DESC_DIR="${2:?}"; shift 2 ;;
    -o|--output-dir) OUTPUT_DIR="${2:?}"; shift 2 ;;
    -e|--exceptions) CHECK_EXC=1; shift ;;
    -E|--errno)      CHECK_ERRNO=1; shift ;;
    -f|--full)       CHECK_EXC=1; CHECK_ERRNO=1; shift ;;
    -j|--jobs)       JOBS="${2:?}"; shift 2 ;;
    -k|--keep)       KEEP=1; shift ;;
    --clang)         USE_CLANG=1; shift ;;
    -h|--help)       usage 0 ;;
    --)              shift; EXTRA=("$@"); break ;;
    -*)              die "unknown option '$1' (see --help)" ;;
    *)               FUNCS+=("$1"); shift ;;
  esac
done

[ -d "$DESC_DIR" ] || die "description directory not found: $DESC_DIR"

# --- build ------------------------------------------------------------------

TMP_BUILD=""
cleanup () { [ -n "$TMP_BUILD" ] && [ "$KEEP" -eq 0 ] && rm -rf "$TMP_BUILD"; }
trap cleanup EXIT

if [ -z "$BUILD_DIR" ]; then
  TMP_BUILD="$(mktemp -d "${TMPDIR:-/tmp}/checkulps-build.XXXXXX")"
  BUILD_DIR="$TMP_BUILD"
fi

CHECKULPS="$BUILD_DIR/src/checkulps/checkulps"

if [ ! -x "$CHECKULPS" ]; then
  echo ">>> building checkulps in $BUILD_DIR" >&2
  declare -a CMAKE_ARGS=(-S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release)
  if [ "$USE_CLANG" -eq 1 ]; then
    CMAKE_ARGS+=(-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++)
  fi
  cmake "${CMAKE_ARGS[@]}" >&2 || die "cmake configuration failed"
  cmake --build "$BUILD_DIR" --target checkulps -j "$JOBS" >&2 \
    || die "build failed"
fi
[ -x "$CHECKULPS" ] || die "checkulps binary not found after build: $CHECKULPS"

# --- glibc loader wrapper ---------------------------------------------------

declare -a LAUNCH=()
if [ -n "$GLIBC_DIR" ]; then
  G="$(cd "$GLIBC_DIR" 2>/dev/null && pwd)" || die "glibc dir not found: $GLIBC_DIR"

  # Locate the built dynamic loader (ld-linux-<arch>.so.<n>).
  LOADER=""
  for f in "$G"/elf/ld-linux*.so.* "$G"/elf/ld-*.so.[0-9]; do
    [ -f "$f" ] || continue
    case "$f" in *.dyn|*.dynsym|*.jmprel|*.note|*.phdr) continue ;; esac
    LOADER="$f"; break
  done
  [ -n "$LOADER" ] || die "no dynamic loader found under $G/elf"

  # glibc library components (only the directories that actually exist).
  LIBPATH="$G"
  for sub in math mathvec elf dlfcn nss nis rt resolv support misc debug nptl; do
    [ -d "$G/$sub" ] && LIBPATH="$LIBPATH:$G/$sub"
  done

  # checkulps is a C++/OpenMP program: make sure its libstdc++ and libgomp
  # resolve while running under the fresh glibc.
  CXX_BIN="${CXX:-g++}"
  [ "$USE_CLANG" -eq 1 ] && CXX_BIN="${CXX:-clang++}"
  for lib in libstdc++.so.6 libgomp.so.1; do
    d=$(dirname "$("$CXX_BIN" -print-file-name="$lib" 2>/dev/null)")
    case ":$LIBPATH:" in *":$d:"*) ;; *) [ -d "$d" ] && LIBPATH="$LIBPATH:$d" ;; esac
  done

  LAUNCH=("$LOADER" --library-path "$LIBPATH")
  echo ">>> glibc loader : $LOADER" >&2
  echo ">>> library path : $LIBPATH" >&2
fi

# --- checkulps flags --------------------------------------------------------

declare -a CHECK_FLAGS=()
[ "$CHECK_EXC" -eq 1 ]   && CHECK_FLAGS+=(-e)
[ "$CHECK_ERRNO" -eq 1 ] && CHECK_FLAGS+=(-E)

# --- function list ----------------------------------------------------------

if [ "${#FUNCS[@]}" -eq 0 ]; then
  for j in "$DESC_DIR"/*.json; do
    [ -e "$j" ] || die "no description files in $DESC_DIR"
    FUNCS+=("$(basename "$j" .json)")
  done
fi

[ -n "$OUTPUT_DIR" ] && mkdir -p "$OUTPUT_DIR"

# --- run --------------------------------------------------------------------

rc=0
for func in "${FUNCS[@]}"; do
  desc="$DESC_DIR/$func.json"
  [ -f "$desc" ] || { echo "${0##*/}: skip '$func': no $desc" >&2; rc=1; continue; }

  if [ -n "$OUTPUT_DIR" ]; then
    echo ">>> $func" >&2
    "${LAUNCH[@]}" "$CHECKULPS" -d "$desc" "${CHECK_FLAGS[@]}" "${EXTRA[@]}" \
      > "$OUTPUT_DIR/$func.out" || rc=1
  else
    echo "==================== $func ===================="
    "${LAUNCH[@]}" "$CHECKULPS" -d "$desc" "${CHECK_FLAGS[@]}" "${EXTRA[@]}" || rc=1
  fi
done

exit "$rc"
