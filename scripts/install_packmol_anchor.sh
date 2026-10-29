#!/usr/bin/env bash
# Local compiler and separately versioned geometry generator; no system install.
set -euo pipefail
project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
tools_dir=/workspace/md-tools
source_dir="$tools_dir/packmol-source"
source_commit=14e50c65fa9120b58e1ba33ad482c7c7260f72b2
test -f "$tools_dir/apt/apt.conf" || { echo 'Run the base project setup first.' >&2; exit 1; }
if test ! -d "$source_dir"; then
    git clone https://github.com/m3g/packmol.git "$source_dir"
    git -C "$source_dir" checkout --detach "$source_commit"
fi
test "$(git -C "$source_dir" rev-parse HEAD)" = "$source_commit" || {
    echo 'Existing Packmol source differs; preserve it and inspect before building.' >&2; exit 1;
}
APT_CONFIG="$tools_dir/apt/apt.conf" apt-get -o Debug::NoLocking=true --download-only --no-install-recommends -y install \
    gfortran-14=14.2.0-19 gfortran-14-x86-64-linux-gnu=14.2.0-19 libgfortran-14-dev=14.2.0-19
for package in "$tools_dir"/apt/archives/*gfortran*.deb; do
    dpkg-deb -x "$package" "$tools_dir/root"
done
cd "$project_dir"
source scripts/activate.sh
if test ! -f "$tools_dir/packmol-anchor-pbc-v1/packmol"; then
    python scripts/build_packmol_anchor_variant.py
fi
python scripts/validate_packmol_anchor.py
