# OpenQC 0.0.1 Release Operator Guide

Normal contributors only need to run:

```bash
nvm use
npm ci
make format lint typecheck test check
npm run check:release
```

This validates and packages `newtontech.openqc@0.0.1`; it does not publish anything.

Repository administrators must configure a protected `marketplace` GitHub environment with required reviewers and a `VSCE_PAT` secret authorized for the `newtontech` publisher. Keep the token in GitHub secrets only.

After the release PR is merged and CI is green, a release operator may create and push `v0.0.1` from the updated `master` branch:

```bash
git checkout master
git pull --ff-only origin master
npm ci
npm run check:release
npm run release:tag
git push origin v0.0.1
```

The tag workflow builds first, pauses for Marketplace approval, publishes the verified VSIX, and creates the GitHub Release last. If any guard or verification fails, fix the source and use a new version; do not move or force-push a published release tag.

See [Marketplace Publishing](MARKETPLACE-PUBLISHING.md) for artifact details and failure handling.
