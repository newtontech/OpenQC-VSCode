# Agent Development & Git Workflow Standards

This document outlines the standards and best practices for agents working on this project, with a focus on the Git workflow.

## Git Workflow & `git worktree`

To maintain a clean and efficient development environment, we will adopt a `git worktree`-based workflow. This allows for parallel development on different features and fixes without needing to stash changes or switch branches in the main project directory.

### Branching Strategy

- **`feat/<feature-name>`:** For developing new features.
- **`fix/<bug-name>`:** For fixing bugs.
- **`docs/<doc-name>`:** For updating documentation.

### PR Documentation Requirements

When implementing a PR, keep user-facing claims synchronized with the actual
implementation:

- If the PR adds, removes, or materially changes supported software, commands,
  LSP behavior, workflows, screenshots, setup steps, or user-visible features,
  update `README.md` in the same PR.
- If `README.md` already accurately describes the changed behavior, state that
  explicitly in the PR summary instead of leaving the documentation check
  implicit.
- Avoid README-only claims about broad support unless the code, tests, or
  compatibility matrix in the same PR proves the claim.

### `git worktree` Usage

When starting a new task (feature, fix, or docs), create a new worktree. The worktree will be located in the parent directory of the main project clone to keep it isolated and clearly named.

**Example: Starting a new feature**

```bash
# 1. Ensure your main branch is up to date
git checkout master
git pull origin master

# 2. Create a new worktree for the feature branch
# This creates a new directory ../OpenQC-VSCode-feat-new-ui
# and a new branch feat/new-ui based on the current master
git worktree add -b feat/new-ui ../OpenQC-VSCode-feat-new-ui master

# 3. Navigate to the new worktree to start development
cd ../OpenQC-VSCode-feat-new-ui
```

**Example: Cleaning up after a merge**

After your feature or fix branch has been merged into `master`, you should clean up the worktree.

```bash
# 1. Remove the worktree directory (ensure you are NOT inside it)
rm -rf ../OpenQC-VSCode-feat-new-ui

# 2. Prune the worktree metadata from git
git worktree prune

# 3. (Optional) Delete the local branch, as it's been merged
git branch -d feat/new-ui
```
