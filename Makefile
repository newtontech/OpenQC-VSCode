.DEFAULT_GOAL := help

.PHONY: help install build format format-write lint typecheck test check release-check package clean
.PHONY: worktree-list cleanup-merged

help: ## Show available targets
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "%-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install the locked Node dependencies
	npm ci

build: ## Compile the extension
	npm run compile

format: ## Check formatting without changing the tree
	npm run format:check

format-write: ## Apply Prettier formatting
	npm run format:write

lint: ## Run ESLint
	npm run lint

typecheck: ## Run the TypeScript type checker
	npm run typecheck

test: ## Run the maintained Jest test suite
	npm run test:all

check: ## Run the canonical pull-request quality gate
	npm run ci:pr

release-check: ## Run quality gates and build verified release assets
	npm run check:release

package: ## Build and verify the VSIX without publishing it
	npm run package:vsix
	npm run verify:vsix
	npm run release:assets

clean: ## Remove generated outputs
	rm -rf out coverage output
	rm -f vsix/*.vsix vsix/*.sha256 vsix/*.sbom.json

worktree-list: ## List repository worktrees
	git worktree list

cleanup-merged: ## Remove worktrees for merged issue branches
	bash scripts/cleanup_merged_worktrees.sh
