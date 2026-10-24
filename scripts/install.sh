#!/usr/bin/env bash
# Debian 13 CPU-only setup; uses signed Debian packages without system writes.
set -euo pipefail
project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export MD_TOOLS_DIR="${MD_TOOLS_DIR:-/workspace/md-tools}"
mkdir -p "$MD_TOOLS_DIR/apt/lists/partial" "$MD_TOOLS_DIR/apt/archives/partial" "$MD_TOOLS_DIR/root" "$MD_TOOLS_DIR/bin" "$MD_TOOLS_DIR/apt/empty"
cat > "$MD_TOOLS_DIR/apt/sources.list" <<'APT'
deb [signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] https://deb.debian.org/debian trixie main
APT
cat > "$MD_TOOLS_DIR/apt/apt.conf" <<APT
Dir::Etc::Parts "$MD_TOOLS_DIR/apt/empty";
Dir::Etc::main "-";
Dir::Etc::sourcelist "$MD_TOOLS_DIR/apt/sources.list";
Dir::Etc::sourceparts "$MD_TOOLS_DIR/apt/empty";
Dir::State::lists "$MD_TOOLS_DIR/apt/lists";
Dir::State::status "/var/lib/dpkg/status";
Dir::Cache::archives "$MD_TOOLS_DIR/apt/archives";
APT::Sandbox::User "$(id -un)";
APT::Get::List-Cleanup "false";
APT
APT_CONFIG="$MD_TOOLS_DIR/apt/apt.conf" /usr/bin/apt-get update
# NoLocking is used only for downloading into our private cache; no dpkg installation.
APT_CONFIG="$MD_TOOLS_DIR/apt/apt.conf" /usr/bin/apt-get -o Debug::NoLocking=true --download-only --no-install-recommends -y install gromacs=2025.2-1
for package in "$MD_TOOLS_DIR"/apt/archives/*.deb; do
    dpkg-deb -x "$package" "$MD_TOOLS_DIR/root"
done
cat > "$MD_TOOLS_DIR/bin/gmx" <<'GMX'
#!/usr/bin/env bash
set -euo pipefail
tools_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export GMXLIB="$tools_dir/root/usr/share/gromacs/top"
export LD_LIBRARY_PATH="$tools_dir/root/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$tools_dir/root/usr/bin/gmx" "$@"
GMX
chmod +x "$MD_TOOLS_DIR/bin/gmx"
python3 -m venv "$MD_TOOLS_DIR/chem-venv"
"$MD_TOOLS_DIR/chem-venv/bin/python" -m pip install --disable-pip-version-check -r "$project_dir/config/requirements-chem.lock"
"$MD_TOOLS_DIR/bin/gmx" --version
"$MD_TOOLS_DIR/chem-venv/bin/python" -c 'from rdkit import rdBase; print("RDKit", rdBase.rdkitVersion)'
bash "$project_dir/scripts/smoke.sh"
