# Remote LSP Smoke Test

Use this checklist after packaging or installing the OpenQC extension locally. The extension is configured with `extensionKind: ["workspace", "ui"]`, so LSP processes prefer the workspace extension host. In Remote SSH, WSL, Dev Containers, and Codespaces, `openqc.lsp.<language>.command` and deprecated `openqc.lsp.<language>.path` resolve in the workspace-side environment, not on the local desktop.

## Local

1. Open the repository in a local VS Code window.
2. Install one language server on the local machine or set `openqc.lsp.gaussian.command` to a known local executable.
3. Open `tests/fixtures/gaussian/water_optimization.com`.
4. Run `OpenQC: Start Language Server`.
5. Check the OpenQC output channel for `local extension host` and the resolved command.
6. If startup fails, confirm the error names the local extension host and the `openqc.lsp.gaussian.command` setting.

## Remote SSH

1. Connect with the VS Code Remote SSH extension.
2. Open a workspace folder on the SSH host.
3. Install the selected language server on the SSH host or set `openqc.lsp.<language>.command` to a path that exists on the SSH host.
4. Open a matching OpenQC fixture from the remote workspace.
5. Run `OpenQC: Start Language Server`.
6. Check the OpenQC output channel for `remote extension host (Remote SSH)` and confirm the command is the remote host command.
7. Temporarily set the command to a missing executable and confirm the startup error says commands and paths resolve in the remote workspace environment.

## WSL

1. Open the repository through `WSL: Open Folder in WSL`.
2. Install the selected language server inside the WSL distribution or configure `openqc.lsp.<language>.command` to a WSL path.
3. Open a matching fixture and start the LSP.
4. Check the output channel for `remote extension host (WSL)`.
5. Confirm missing executable errors refer to the WSL workspace-side environment.

## Dev Containers

1. Reopen the repository in a dev container.
2. Install the selected language server in the container image or set `openqc.lsp.<language>.command` to a container path.
3. Open a matching fixture and start the LSP.
4. Check the output channel for `remote extension host (Dev Container)`.
5. Confirm missing executable errors refer to the container workspace-side environment.

## Codespaces

1. Open the repository in GitHub Codespaces.
2. Install the selected language server in the Codespace or configure `openqc.lsp.<language>.command` to a Codespace path.
3. Open a matching fixture and start the LSP.
4. Check the output channel for `remote extension host (Codespaces)`.
5. Confirm missing executable errors refer to the Codespace workspace-side environment.

## Webview Resource Audit

- `MoleculeViewerWebview`, `DataPlotter`, and `StructureViewer` load bundled JavaScript through `webview.asWebviewUri(...)` with `localResourceRoots` scoped to extension resources.
- Inline setup, AI, sidebar, and results webviews do not load local extension files.
- Extension persistence uses VS Code `workspaceState` or `globalState`; no browser `localStorage` or `sessionStorage` use was found.
