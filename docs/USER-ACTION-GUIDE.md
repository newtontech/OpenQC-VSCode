# User Action Guide - Publishing OpenQC-VSCode to Marketplace

## Quick Start (5 Steps to Publish)

### 1. Install vsce Tool
```bash
npm install -g @vscode/vsce
```

### 2. Create Azure DevOps PAT
1. Go to https://dev.azure.com/
2. Click "User Settings" > "Personal Access Tokens"
3. Create new token with:
   - Name: OpenQC-VSCode Publishing
   - Scopes: Marketplace > Manage
4. Copy the token (save it securely!)

### 3. Login to vsce
```bash
cd /home/yhm/desktop/code/OpenQC-VSCode
vsce login newtontech
# Paste your PAT when prompted
```

### 4. Publish Extension
```bash
vsce publish
```

### 5. Verify
Go to https://marketplace.visualstudio.com/ and search for "OpenQC-VSCode"

## What's Been Completed

### ✅ Functional Implementation
- **viewResults**: Complete results viewing with job information and output display
- **exportData**: Export completed job data to JSON/CSV format
- **Helper Functions**: Results HTML generation with dark theme styling

### ✅ Test Coverage
- **327 tests passing** with 95.98% statement coverage
- New tests for DefinitionProvider
- Enhanced tests for StructureViewer and DataPlotter
- All 7 quantum chemistry packages tested

### ✅ Documentation
- CHANGELOG.md updated for v2.0.0
- Marketplace publishing guide created
- Completion summary with all details

### ✅ Package Ready
- openqc-vscode-2.0.0.vsix created (41.84 MB)
- All required fields in package.json verified
- Extension icon present (128x128 pixels)

## Pre-Publication Checklist

Before publishing, verify:

- [ ] Azure DevOps PAT created with "Marketplace > Manage" scope
- [ ] vsce tool installed globally
- [ ] Extension packaged successfully (openqc-vscode-2.0.0.vsix exists)
- [ ] All tests passing (run `npm run test:unit`)
- [ ] TypeScript compiles without errors (run `npm run compile`)

## After Publishing

1. **Create GitHub Release**:
   ```bash
   git tag -a v2.0.0 -m "Release v2.0.0"
   git push origin v2.0.0
   ```

2. **Add Marketplace Badge to README**:
   ```markdown
   [![Version](https://img.shields.io/visual-studio-marketplace/v/newtontech.openqc-vscode)](https://marketplace.visualstudio.com/items?itemName=newtontech.openqc-vscode)
   ```

3. **Create Screenshots** for Marketplace listing:
   - Molecular structure visualization
   - Data plotting interface
   - Input file preview
   - Sidebar panels (Molecules and Jobs)

## Troubleshooting

### "Authentication failed"
- Verify PAT is correct
- Check PAT has "Marketplace > Manage" scope
- Try `vsce logout` then `vsce login` again

### "Validation failed"
- Run `vsce publish --dry-run` first
- Check icon.png is 128x128 pixels
- Verify all required fields in package.json

### "Version already exists"
- Increment version in package.json
- Run `npm run compile`
- Package again: `npx vsce package`

## Support

For detailed information, see:
- `/home/yhm/desktop/code/OpenQC-VSCode/docs/MARKETPLACE-PUBLISHING.md` - Full publishing guide
- `/home/yhm/desktop/code/OpenQC-VSCode/docs/ISSUE-2-COMPLETION-SUMMARY.md` - Complete task summary

## Ready to Publish?

Run these commands in order:

```bash
# 1. Navigate to project directory
cd /home/yhm/desktop/code/OpenQC-VSCode

# 2. Login to vsce (replace with your PAT)
vsce login newtontech

# 3. Publish!
vsce publish
```

That's it! Your extension will be live on the marketplace within minutes.
