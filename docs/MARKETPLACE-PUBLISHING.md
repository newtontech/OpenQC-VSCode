# Marketplace Publishing Guide

OpenQC's release identity is `newtontech.openqc@0.0.1`, displayed as **OpenQC - DFT/MD/Quantum Chemistry Suite**. The supported release path publishes one verified VSIX to the VS Code Marketplace and then attaches the same bytes, checksum, and SBOM to the GitHub Release.

## Repository setup

- Node is fixed to major 22 through `.nvmrc`, `package.json`, and GitHub Actions.
- npm is declared as `npm@10.9.8` and dependencies are installed with `npm ci`.
- `@vscode/vsce` is a dev dependency pinned exactly to `3.9.2`; do not use a global or floating `vsce` for release work.
- Configure a protected GitHub environment named `marketplace` and store a publisher-scoped `VSCE_PAT` secret in that environment. The token needs **Marketplace > Manage** permission for publisher `newtontech`.

## Build and verify locally

From a clean checkout:

```bash
nvm use
npm ci
make format lint typecheck test check
npm run check:release
```

The release gate produces ignored local artifacts:

- `vsix/openqc-0.0.1.vsix`
- `vsix/openqc-0.0.1.vsix.sha256`
- `vsix/openqc-0.0.1.sbom.json`

`verify:vsix` reads the packaged archive, confirms `newtontech.openqc@0.0.1`, enforces size/file-count limits, and rejects development files, nested VSIX files, credentials, caches, tests, scripts, docs, and the retired `core/` package.

## Create a release tag

Tagging is an explicit release action. Do it only after the release PR is merged and `origin/master` is green:

```bash
git checkout master
git pull --ff-only origin master
npm ci
npm run check:release
npm run release:tag
git push origin v0.0.1
```

`release:tag` fetches `origin/master`, refuses a dirty tree, refuses a commit that is not exactly the current remote master, rejects an existing tag, and derives `v0.0.1` from `package.json`. `release:tag:push` combines the final two commands, but the separate push is easier to review.

## Automated release stages

1. **Build and verify** checks out full tag history, confirms the tag points to the clean current `origin/master`, installs the lockfile on Node 22, and runs `npm run check:release`.
2. **Publish Marketplace** downloads those exact artifacts, verifies SHA-256, waits for the protected `marketplace` environment, and publishes with the pinned local `vsce` plus `VSCE_PAT`.
3. **Finalize GitHub Release** runs only after Marketplace publication succeeds, re-verifies SHA-256, and creates the GitHub Release with the VSIX, checksum, and CycloneDX SBOM.

The workflow does not publish Open VSX and does not rewrite `CHANGELOG.md` after tagging.

## Failure boundaries

| Failure | Result |
|---------|--------|
| Dirty tree, stale master, or mismatched tag | Release stops before building |
| Quality, package, identity, hygiene, checksum, or SBOM failure | No external publication |
| Missing/unauthorized `VSCE_PAT` or environment rejection | Marketplace is unchanged and no GitHub Release is created |
| GitHub Release finalization failure | Marketplace may already contain the version; rerun only the failed job after checking attached artifact hashes |

Never print, commit, or place the PAT in `.env`; `.env` and common key formats are ignored by Git and rejected from the VSIX.
