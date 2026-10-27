# Repository and calculation backup

The MD work uses branch `md/cbaa-surface` in `kimk42/practice`. The GitHub repository is public. Only tracked code,model inputs,configuration,derived numerical reports and their history belong on that branch. Original documents,downloaded articles,raw structures and runs/ remain excluded. Never upload a private recovery archive to this public repository.

`scripts/backup_md_workspace.py` creates a recovery directory under `/workspace/md-backups`,outside the worktree. It preserves a compressed worktree including ignored calculation data/private sources,a self-contained Git bundle,the uncommitted worktree patch and checkpoint/file hashes. It checks archive decompression,Git bundle integrity and unchanged checkpoint sizes/timestamps during creation. It refuses to run while a GROMACS/SQM/Antechamber calculation is active;it does not stop or restart calculations. Existing pending tracked deletions are preserved in the patch rather than silently committed/restored.

This archive is still in the same cloud environment. It is not an independent external backup until it has been transferred to a user-selected private destination and verified there. Installed shared tools are not archived;recreate them using the project's setup instructions. The archive contains a RESTORE.txt with the restore order,including application of the worktree patch before extracting the working files.

To make a new local snapshot when calculations are idle:

```bash
source scripts/activate.sh
python scripts/backup_md_workspace.py
```

The backup destination's latest-backup.json records its actual verified location. Backup status is separate from MD scientific validation and separate from GitHub branch synchronization.
