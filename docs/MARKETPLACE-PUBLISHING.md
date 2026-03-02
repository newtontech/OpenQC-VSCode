# Marketplace Publishing Guide

This guide walks through publishing OpenQC-VSCode to the Visual Studio Code Marketplace.

## Prerequisites

1. **Azure DevOps Account**: You need a Personal Access Token (PAT) for publishing
2. **vsce Tool**: The Visual Studio Code Extension Manager

### Setting up Azure DevOps PAT

1. Go to [Azure DevOps](https://dev.azure.com/)
2. Sign in or create an account
3. Click on "User Settings" (top right) > "Personal Access Tokens"
4. Create a new token with:
   - **Organization**: All accessible organizations
   - **Scopes**: Marketplace > Manage
5. Copy the token (you won't see it again!)

### Installing vsce

```bash
npm install -g @vscode/vsce
```

Or install locally:

```bash
npm install --save-dev @vscode/vsce
```

## Pre-Publication Checklist

### 1. Verify package.json

Run this checklist:

```bash
# Verify all required fields
cat package.json | grep -E "(name|displayName|description|version|publisher|repository|bugs|homepage|license|icon)"
```

Required fields:
- ✅ `name`: openqc-vscode
- ✅ `displayName`: OpenQC-VSCode
- ✅ `description`: Clear and concise
- ✅ `version`: Follows semver (current: 2.0.0)
- ✅ `publisher`: newtontech
- ✅ `repository`: Git repository URL
- ✅ `bugs`: Issue tracker URL
- ✅ `homepage`: Project homepage
- ✅ `license`: MIT
- ✅ `icon`: icon.png (128x128 pixels)

### 2. Verify Icon

```bash
file icon.png
# Should show: PNG image data, 128 x 128
```

If you need to resize:

```bash
# Using ImageMagick
convert icon.svg -resize 128x128 icon.png
```

### 3. Build and Test

```bash
# Compile TypeScript
npm run compile

# Run tests
npm run test:unit

# Check for linting errors
npm run lint

# Package the extension
npx vsce package
```

### 4. Test the Extension Locally

1. Install the packaged extension:
   ```bash
   code --install-extension openqc-vscode-2.0.0.vsix
   ```

2. Test core functionality:
   - Open a quantum chemistry file (POSCAR, input.com, etc.)
   - Verify syntax highlighting
   - Test molecular visualization
   - Test data plotting
   - Test sidebar panels
   - Test LSP management

## Publishing Process

### Step 1: Create Publisher (First Time Only)

If you haven't created a publisher yet:

```bash
vsce create-publisher your-publisher-name
```

For OpenQC-VSCode, the publisher is `newtontech`.

### Step 2: Login to vsce

```bash
vsce login your-publisher-name
```

When prompted, paste your Azure DevOps PAT.

### Step 3: Publish

```bash
# Publish the extension
vsce publish

# Or publish a specific version
vsce publish patch
vsce publish minor
vsce publish major
```

### Step 4: Verify

1. Go to [Visual Studio Code Marketplace](https://marketplace.visualstudio.com/)
2. Search for "OpenQC-VSCode"
3. Verify the extension appears correctly
4. Check the extension page for proper display

## Post-Publication

### 1. Create a Release on GitHub

```bash
# Create a git tag
git tag -a v2.0.0 -m "Release v2.0.0 - Universal Quantum Chemistry Platform"
git push origin v2.0.0
```

Then create a GitHub release with:
- Version number
- Release notes from CHANGELOG.md
- Screenshots
- Installation instructions

### 2. Update Documentation

- Verify README.md is up to date
- Check that all links work
- Update any version-specific information

### 3. Announce

- Post on relevant forums (Computational Chemistry, VS Code extensions)
- Update project website
- Send to interested users

## Troubleshooting

### Authentication Failed

```
Error: Authentication failed. Please check your PAT.
```

Solution:
1. Verify your PAT is correct
2. Check that the PAT has "Marketplace > Manage" scope
3. Try logging in again: `vsce logout` then `vsce login`

### Package Validation Failed

```
Error: Validation failed
```

Solution:
1. Run `vsce publish --dry-run` to check for issues
2. Verify all required fields in package.json
3. Check that icon.png exists and is 128x128 pixels
4. Verify README.md doesn't have broken links

### Version Conflict

```
Error: Version x.y.z already exists
```

Solution:
1. Increment the version in package.json
2. Run `npm run compile` to update compiled files
3. Package again: `npx vsce package`

## Maintenance

### Updating the Extension

1. Make changes
2. Update version in package.json
3. Update CHANGELOG.md
4. Run tests: `npm run test:unit`
5. Package: `npx vsce package`
6. Test locally
7. Publish: `vsce publish`

### Versioning Strategy

- **Major version (X.0.0)**: Breaking changes, major features
- **Minor version (0.X.0)**: New features, backward compatible
- **Patch version (0.0.X)**: Bug fixes, small improvements

## Resources

- [vsce Documentation](https://github.com/microsoft/vscode-vsce)
- [VS Code Extension API](https://code.visualstudio.com/api)
- [Marketplace Publishing Guidelines](https://code.visualstudio.com/api/working-with-extensions/publishing-extension)
