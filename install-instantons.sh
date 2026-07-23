#!/usr/bin/env bash

# Install MAMBO in a configured Herwig 7.3.0 tree and build the
# HerwigQCDInstantons contrib plugin.

set -Eeuo pipefail

PROGRAM_NAME=${0##*/}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)

HERWIG_BUILD=
JOBS=
BACKUP_USED=0
BACKUP_ROOT=

report_failed_install() {
  local status=$?
  if (( status != 0 && BACKUP_USED )); then
    printf '%s: backups from the incomplete installation are in %s\n' \
      "$PROGRAM_NAME" "$BACKUP_ROOT" >&2
  fi
}
trap report_failed_install EXIT

usage() {
  cat <<EOF
Usage:
  $PROGRAM_NAME [options] HERWIG_SOURCE HERWIG_INSTALL

Install the QCD-instanton model into an existing, configured Herwig 7.3.0
build. HERWIG_INSTALL must be the prefix recorded when that build was
configured.

Options:
  -b, --build-dir DIR  Herwig build directory (default: HERWIG_SOURCE)
  -j, --jobs N         Parallel build jobs (default: detected CPU count)
  -h, --help           Show this help

Examples:
  ./$PROGRAM_NAME /path/to/Herwig-7.3.0 /path/to/herwig

  ./$PROGRAM_NAME --build-dir /path/to/build -j 8 \\
    /path/to/Herwig-7.3.0 /path/to/herwig
EOF
}

log() {
  printf '==> %s\n' "$*"
}

die() {
  printf '%s: error: %s\n' "$PROGRAM_NAME" "$*" >&2
  exit 1
}

require_argument() {
  local option=$1
  local count=$2
  (( count >= 2 )) || die "$option requires an argument"
}

canonical_directory() {
  local directory=$1
  local description=$2
  [[ -d "$directory" ]] || die "$description is not a directory: $directory"
  (cd "$directory" && pwd -P)
}

make_value() {
  local makefile=$1
  local variable=$2
  awk -v wanted="$variable" '
    $0 ~ "^[[:space:]]*" wanted "[[:space:]]*=" {
      line = $0
      sub("^[[:space:]]*" wanted "[[:space:]]*=[[:space:]]*", "", line)
      print line
      exit
    }
  ' "$makefile"
}

detect_jobs() {
  local detected=

  if command -v getconf >/dev/null 2>&1; then
    detected=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
  fi
  if [[ ! "$detected" =~ ^[1-9][0-9]*$ ]] &&
     command -v sysctl >/dev/null 2>&1; then
    detected=$(sysctl -n hw.ncpu 2>/dev/null || true)
  fi
  if [[ ! "$detected" =~ ^[1-9][0-9]*$ ]]; then
    detected=1
  fi

  printf '%s\n' "$detected"
}

backup_file() {
  local source_file=$1
  local backup_name=$2
  local backup_file_path=$BACKUP_ROOT/$backup_name

  [[ -e "$source_file" ]] || return 0
  [[ ! -L "$source_file" ]] ||
    die "refusing to replace symbolic link: $source_file"

  if [[ ! -e "$backup_file_path" ]]; then
    mkdir -p "$(dirname "$backup_file_path")"
    cp -p "$source_file" "$backup_file_path"
  fi
  BACKUP_USED=1
}

copy_managed_file() {
  local source_file=$1
  local destination_file=$2
  local backup_name=$3

  [[ -f "$source_file" ]] || die "required repository file is missing: $source_file"
  [[ ! -d "$destination_file" ]] ||
    die "expected a file but found a directory: $destination_file"
  [[ ! -L "$destination_file" ]] ||
    die "refusing to replace symbolic link: $destination_file"

  if [[ -f "$destination_file" ]] &&
     cmp -s "$source_file" "$destination_file"; then
    log "Already current: $destination_file"
    return 0
  fi

  backup_file "$destination_file" "$backup_name"
  mkdir -p "$(dirname "$destination_file")"
  cp "$source_file" "$destination_file"
  chmod 644 "$destination_file"
  log "Installed source: $destination_file"
}

patch_make_list() {
  local makefile=$1
  local variable=$2
  local entry=$3
  local backup_name=$4
  local temporary

  [[ -f "$makefile" ]] || die "required makefile is missing: $makefile"
  [[ ! -L "$makefile" ]] || die "refusing to edit symbolic link: $makefile"

  temporary=$(mktemp "${TMPDIR:-/tmp}/herwig-instantons.XXXXXX")
  if ! awk -v wanted="$variable" -v entry="$entry" '
    BEGIN {
      in_list = 0
      saw_list = 0
      completed = 0
      found = 0
    }

    !in_list && $0 ~ "^[[:space:]]*" wanted "[[:space:]]*=" {
      saw_list = 1
      in_list = 1
      found = 0
      if ($0 !~ /\\[[:space:]]*$/) {
        exit 7
      }
      print
      next
    }

    in_list {
      line = $0
      content = $0
      sub(/^[[:space:]]*/, "", content)
      sub(/[[:space:]]*\\[[:space:]]*$/, "", content)
      sub(/[[:space:]]*$/, "", content)
      if (content == entry) {
        found = 1
      }

      if ($0 !~ /\\[[:space:]]*$/) {
        if (!found) {
          sub(/[[:space:]]*$/, "", line)
          print line " \\"
          print entry
        } else {
          print line
        }
        in_list = 0
        completed = 1
        next
      }

      print
      next
    }

    {
      print
    }

    END {
      if (!saw_list || !completed) {
        exit 8
      }
    }
  ' "$makefile" >"$temporary"; then
    rm -f "$temporary"
    die "could not update $variable safely in $makefile"
  fi

  if cmp -s "$makefile" "$temporary"; then
    rm -f "$temporary"
    log "Already listed: $entry in $makefile"
    return 0
  fi

  backup_file "$makefile" "$backup_name"
  cp "$temporary" "$makefile"
  chmod 644 "$makefile"
  rm -f "$temporary"
  log "Added $entry to $variable in $makefile"
}

positionals=()
while (( $# > 0 )); do
  case $1 in
    -b|--build-dir)
      require_argument "$1" "$#"
      HERWIG_BUILD=$2
      shift 2
      ;;
    -j|--jobs)
      require_argument "$1" "$#"
      JOBS=$2
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while (( $# > 0 )); do
        positionals[${#positionals[@]}]=$1
        shift
      done
      ;;
    -*)
      die "unknown option: $1"
      ;;
    *)
      positionals[${#positionals[@]}]=$1
      shift
      ;;
  esac
done

(( ${#positionals[@]} == 2 )) || {
  usage >&2
  exit 2
}

HERWIG_SOURCE=$(canonical_directory "${positionals[0]}" "Herwig source")
HERWIG_INSTALL=$(canonical_directory "${positionals[1]}" "Herwig installation")
if [[ -z "$HERWIG_BUILD" ]]; then
  HERWIG_BUILD=$HERWIG_SOURCE
else
  HERWIG_BUILD=$(canonical_directory "$HERWIG_BUILD" "Herwig build")
fi

case "$HERWIG_SOURCE$HERWIG_BUILD$HERWIG_INSTALL" in
  *[[:space:]]*)
    die "Herwig's contrib build machinery does not support whitespace in paths"
    ;;
esac

if [[ -z "$JOBS" ]]; then
  JOBS=$(detect_jobs)
fi
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] ||
  die "--jobs must be a positive integer"

SOURCE_PHASESPACE=$HERWIG_SOURCE/MatrixElement/Matchbox/Phasespace
BUILD_PHASESPACE=$HERWIG_BUILD/MatrixElement/Matchbox/Phasespace
BUILD_MAKEFILE=$HERWIG_BUILD/Makefile
CONTRIB_ROOT=$HERWIG_BUILD/Contrib
CONTRIB_SETUP=$CONTRIB_ROOT/make_makefiles.sh
PLUGIN_BUILD=$CONTRIB_ROOT/HerwigQCDInstantons

[[ -f "$HERWIG_SOURCE/configure.ac" ]] ||
  die "not a Herwig source tree: $HERWIG_SOURCE"
grep -Fq 'AC_INIT([Herwig],[7.3.0]' "$HERWIG_SOURCE/configure.ac" ||
  die "the source tree is not Herwig 7.3.0"
[[ -f "$BUILD_MAKEFILE" && -f "$HERWIG_BUILD/config.status" ]] ||
  die "the Herwig build is not configured: $HERWIG_BUILD"
[[ -f "$BUILD_PHASESPACE/Makefile" ]] ||
  die "the configured Matchbox phase-space build directory is missing"
[[ -f "$CONTRIB_SETUP" ]] ||
  die "the configured contrib setup script is missing: $CONTRIB_SETUP"

BUILD_VERSION=$(make_value "$BUILD_MAKEFILE" PACKAGE_VERSION)
[[ "$BUILD_VERSION" == 7.3.0 ]] ||
  die "the configured build is Herwig ${BUILD_VERSION:-unknown}, not 7.3.0"

CONFIGURED_CXX=$(make_value "$BUILD_MAKEFILE" CXX)
[[ -n "$CONFIGURED_CXX" ]] ||
  die "could not read the configured C++ compiler from $BUILD_MAKEFILE"
CXX_COMMAND=${CONFIGURED_CXX%%[[:space:]]*}

TOP_SOURCE=$(make_value "$BUILD_MAKEFILE" top_srcdir)
[[ -n "$TOP_SOURCE" ]] || die "could not read top_srcdir from $BUILD_MAKEFILE"
case $TOP_SOURCE in
  /*) CONFIGURED_SOURCE=$TOP_SOURCE ;;
  *)  CONFIGURED_SOURCE=$HERWIG_BUILD/$TOP_SOURCE ;;
esac
CONFIGURED_SOURCE=$(canonical_directory "$CONFIGURED_SOURCE" \
  "configured Herwig source")
[[ "$CONFIGURED_SOURCE" == "$HERWIG_SOURCE" ]] ||
  die "build directory was configured from $CONFIGURED_SOURCE, not $HERWIG_SOURCE"

CONFIGURED_PREFIX=$(make_value "$BUILD_MAKEFILE" prefix)
[[ -n "$CONFIGURED_PREFIX" ]] ||
  die "could not read the installation prefix from $BUILD_MAKEFILE"
CONFIGURED_PREFIX=$(canonical_directory "$CONFIGURED_PREFIX" \
  "configured Herwig installation")
[[ "$CONFIGURED_PREFIX" == "$HERWIG_INSTALL" ]] ||
  die "build prefix is $CONFIGURED_PREFIX, not $HERWIG_INSTALL"

for required_command in awk cmp cp make mktemp; do
  command -v "$required_command" >/dev/null 2>&1 ||
    die "required command is not available: $required_command"
done

export PATH="$HERWIG_INSTALL/bin${PATH:+:$PATH}"
if [[ -r "$HERWIG_INSTALL/bin/activate" ]]; then
  log "Activating $HERWIG_INSTALL"
  set +u
  # shellcheck disable=SC1090
  source "$HERWIG_INSTALL/bin/activate"
  set -u
fi
export PATH="$HERWIG_INSTALL/bin${PATH:+:$PATH}"
export LD_LIBRARY_PATH="$HERWIG_INSTALL/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="$HERWIG_INSTALL/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"

if [[ "$CXX_COMMAND" == */* ]]; then
  [[ -x "$CXX_COMMAND" ]] ||
    die "configured C++ compiler is unavailable: $CXX_COMMAND; reconfigure Herwig or restore that compiler"
else
  command -v "$CXX_COMMAND" >/dev/null 2>&1 ||
    die "configured C++ compiler is unavailable: $CXX_COMMAND; reconfigure Herwig or restore that compiler"
fi

[[ -x "$HERWIG_INSTALL/bin/Herwig" ]] ||
  die "Herwig executable is missing from $HERWIG_INSTALL/bin"
HERWIG_VERSION_OUTPUT=$("$HERWIG_INSTALL/bin/Herwig" --version 2>&1) ||
  die "the installed Herwig executable could not be run"
INSTALLED_VERSION=$(printf '%s\n' "$HERWIG_VERSION_OUTPUT" |
  awk '$1 == "Herwig" { print $2; exit }')
[[ "$INSTALLED_VERSION" == 7.3.0 ]] ||
  die "the installation contains Herwig ${INSTALLED_VERSION:-unknown}, not 7.3.0"

command -v gsl-config >/dev/null 2>&1 ||
  die "gsl-config is required to link Instantons.so"

BACKUP_ROOT=$HERWIG_SOURCE/.herwig-qcd-instantons-backup/$(date -u +%Y%m%dT%H%M%SZ)-$$

log "Installing MAMBO sources"
copy_managed_file \
  "$SCRIPT_DIR/Phasespace/MamboPhasespace.h" \
  "$SOURCE_PHASESPACE/MamboPhasespace.h" \
  "source/MatrixElement/Matchbox/Phasespace/MamboPhasespace.h"
copy_managed_file \
  "$SCRIPT_DIR/Phasespace/MamboPhasespace.cc" \
  "$SOURCE_PHASESPACE/MamboPhasespace.cc" \
  "source/MatrixElement/Matchbox/Phasespace/MamboPhasespace.cc"

log "Updating Herwig phase-space source lists"
patch_make_list \
  "$SOURCE_PHASESPACE/Makefile.am" ALL_H_FILES MamboPhasespace.h \
  "source/MatrixElement/Matchbox/Phasespace/Makefile.am"
patch_make_list \
  "$SOURCE_PHASESPACE/Makefile.am" ALL_CC_FILES MamboPhasespace.cc \
  "source/MatrixElement/Matchbox/Phasespace/Makefile.am"
patch_make_list \
  "$SOURCE_PHASESPACE/Makefile.in" ALL_H_FILES MamboPhasespace.h \
  "source/MatrixElement/Matchbox/Phasespace/Makefile.in"
patch_make_list \
  "$SOURCE_PHASESPACE/Makefile.in" ALL_CC_FILES MamboPhasespace.cc \
  "source/MatrixElement/Matchbox/Phasespace/Makefile.in"
patch_make_list \
  "$BUILD_PHASESPACE/Makefile" ALL_H_FILES MamboPhasespace.h \
  "build/MatrixElement/Matchbox/Phasespace/Makefile"
patch_make_list \
  "$BUILD_PHASESPACE/Makefile" ALL_CC_FILES MamboPhasespace.cc \
  "build/MatrixElement/Matchbox/Phasespace/Makefile"

log "Building MAMBO and the Herwig core with $JOBS job(s)"
make -C "$BUILD_PHASESPACE" -j"$JOBS"
make -C "$HERWIG_BUILD" -j"$JOBS"

log "Installing the rebuilt Herwig core into $HERWIG_INSTALL"
make -C "$HERWIG_BUILD" install

log "Preparing the instanton contrib plugin"
mkdir -p "$PLUGIN_BUILD"
copy_managed_file \
  "$SCRIPT_DIR/MEInstanton.cc" "$PLUGIN_BUILD/MEInstanton.cc" \
  "contrib/HerwigQCDInstantons/MEInstanton.cc"
copy_managed_file \
  "$SCRIPT_DIR/MEInstanton.h" "$PLUGIN_BUILD/MEInstanton.h" \
  "contrib/HerwigQCDInstantons/MEInstanton.h"
copy_managed_file \
  "$SCRIPT_DIR/Makefile.in" "$PLUGIN_BUILD/Makefile.in" \
  "contrib/HerwigQCDInstantons/Makefile.in"

(
  cd "$CONTRIB_ROOT"
  bash "$CONTRIB_SETUP"
)

[[ -f "$PLUGIN_BUILD/Makefile" ]] ||
  die "Herwig did not generate the contrib plugin Makefile"
PLUGIN_PREFIX=$(make_value "$PLUGIN_BUILD/Makefile" HERWIGINSTALL)
[[ "$PLUGIN_PREFIX" == "$HERWIG_INSTALL" ]] ||
  die "generated plugin prefix is ${PLUGIN_PREFIX:-unknown}, not $HERWIG_INSTALL"

log "Building Instantons.so"
make -C "$PLUGIN_BUILD" clean
make -C "$PLUGIN_BUILD" -j"$JOBS"

INSTALLED_PLUGIN=$HERWIG_INSTALL/lib/Herwig/Instantons.so
if [[ -f "$INSTALLED_PLUGIN" ]] &&
   cmp -s "$PLUGIN_BUILD/Instantons.so" "$INSTALLED_PLUGIN"; then
  log "Already current: $INSTALLED_PLUGIN"
else
  backup_file "$INSTALLED_PLUGIN" "install/lib/Herwig/Instantons.so"
  make -C "$PLUGIN_BUILD" HERWIGINSTALL="$HERWIG_INSTALL" install
fi

INSTALLED_MAMBO_HEADER=$HERWIG_INSTALL/include/Herwig/MatrixElement/Matchbox/Phasespace/MamboPhasespace.h
[[ -s "$INSTALLED_MAMBO_HEADER" ]] ||
  die "MAMBO header was not installed: $INSTALLED_MAMBO_HEADER"
cmp -s "$SCRIPT_DIR/Phasespace/MamboPhasespace.h" "$INSTALLED_MAMBO_HEADER" ||
  die "installed MAMBO header does not match this repository"
[[ -s "$BUILD_PHASESPACE/Phasespace__all.cc" ]] ||
  die "the combined Matchbox phase-space source was not generated"
grep -Fq 'describeHerwigMamboPhasespace' \
  "$BUILD_PHASESPACE/Phasespace__all.cc" ||
  die "the rebuilt Matchbox source does not contain MAMBO"
[[ -s "$PLUGIN_BUILD/Instantons.so" && -s "$INSTALLED_PLUGIN" ]] ||
  die "Instantons.so was not built and installed"
cmp -s "$PLUGIN_BUILD/Instantons.so" "$INSTALLED_PLUGIN" ||
  die "installed Instantons.so does not match the new build"

printf '\nInstalled Herwig QCD Instantons successfully.\n'
printf '  Herwig:    %s\n' "$HERWIG_INSTALL"
printf '  MAMBO:     %s\n' "$INSTALLED_MAMBO_HEADER"
printf '  Plugin:    %s\n' "$INSTALLED_PLUGIN"
if (( BACKUP_USED )); then
  printf '  Backups:   %s\n' "$BACKUP_ROOT"
fi
if [[ -r "$HERWIG_INSTALL/bin/activate" ]]; then
  printf '\nActivate this installation in new shells with:\n'
  printf '  source %q\n' "$HERWIG_INSTALL/bin/activate"
fi
