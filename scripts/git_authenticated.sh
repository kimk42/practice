#!/usr/bin/env bash
# Use this worktree's private credentials without changing other worktrees.
set -euo pipefail
task_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export GH_CONFIG_DIR="$task_root/.local-auth/github"
unset GH_TOKEN GITHUB_TOKEN
exec git -c credential.helper= -c 'credential.helper=!gh auth git-credential' "$@"
