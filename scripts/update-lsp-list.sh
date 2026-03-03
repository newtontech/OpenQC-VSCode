#!/bin/bash
#
# Update LSP List Script
# Fetches LSP repositories from OpenQuantumChemistry GitHub organization
#

set -e

GITHUB_API_URL="https://api.github.com/orgs/OpenQuantumChemistry/repos"

echo "Fetching LSP repositories from OpenQuantumChemistry..."

response=$(curl -s -H "Accept: application/vnd.github.v3+json" \
  -H "User-Agent: OpenQC-VSCode-UpdateScript" \
  "\${GITHUB_API_URL}?per_page=100&sort=updated" 2>/dev/null || true)

if [ -z "\$response" ]; then
  echo "Error: Failed to fetch repositories"
  exit 1
fi

echo "LSP Repositories:"
echo "\$response" | grep -o '"name":"[^"]*-lsp[^"]*"' | sed 's/"name":"//g; s/"\$//g'
