# OpenQC-VSCode Makefile
# 常用开发命令和清理命令

.PHONY: help install install-dev build test lint format clean clean-all
.PHONY: worktree-create worktree-clean worktree-list
.PHONY: python-install python-test python-lint python-format
.PHONY: vscode-compile vscode-test vscode-package vscode-publish docs

# 默认目标
.DEFAULT_GOAL := help

# 颜色定义
BLUE := \033[36m
GREEN := \033[32m
YELLOW := \033[33m
RED := \033[31m
RESET := \033[0m

##@ 帮助

help: ## 显示本帮助信息
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "$(BLUE)%-20s$(RESET) %s\n", $$1, $$2}' $(MAKEFILE_LIST)

##@ 安装与依赖

install: ## 安装所有依赖（Node.js + Python）
	@echo "$(GREEN)Installing Node.js dependencies...$(RESET)"
	npm install
	@echo "$(GREEN)Installing Python dependencies...$(RESET)"
	cd core && pip install -e ".[dev]"

install-dev: install ## 安装开发依赖（包含 husky 钩子）
	npm run prepare

python-install: ## 仅安装 Python 依赖
	cd core && pip install -e ".[dev]"

##@ 编译与构建

build: vscode-compile ## 构建项目（编译 TypeScript）

vscode-compile: ## 编译 VS Code 扩展 TypeScript
	@echo "$(GREEN)Compiling TypeScript...$(RESET)"
	npm run compile

watch: ## 监听 TypeScript 文件变化并自动编译
	npm run watch

python-build: ## 构建 Python 包
	cd core && python -m build

vscode-package: build ## 打包 VS Code 扩展
	@echo "$(GREEN)Packaging VS Code extension...$(RESET)"
	vsce package

##@ 测试

test: test-unit test-integration ## 运行所有测试

test-unit: ## 运行单元测试
	@echo "$(GREEN)Running unit tests...$(RESET)"
	npm run test:unit
	cd core && pytest tests/unit -v

test-integration: ## 运行集成测试
	@echo "$(GREEN)Running integration tests...$(RESET)"
	npm run test:integration
	cd core && pytest tests/integration -v

test-coverage: ## 运行测试并生成覆盖率报告
	@echo "$(GREEN)Running tests with coverage...$(RESET)"
	npm run test:coverage
	cd core && pytest --cov=openqc --cov-report=html --cov-report=term

test-e2e: build ## 运行端到端测试
	@echo "$(GREEN)Running E2E tests...$(RESET)"
	npm run test:e2e

python-test: ## 仅运行 Python 测试
	cd core && pytest -v

##@ 代码检查与格式化

lint: ## 运行所有 lint 检查
	@echo "$(GREEN)Linting TypeScript...$(RESET)"
	npm run lint
	@echo "$(GREEN)Linting Python...$(RESET)"
	cd core && ruff check .
	cd core && mypy openqc

format: format-write ## 格式化所有代码（等同于 format-write）

format-check: ## 检查代码格式
	@echo "$(GREEN)Checking TypeScript format...$(RESET)"
	npm run format:check
	@echo "$(GREEN)Checking Python format...$(RESET)"
	cd core && black --check .
	cd core && isort --check-only .

format-write: ## 格式化代码并保存
	@echo "$(GREEN)Formatting TypeScript...$(RESET)"
	npm run format:write
	@echo "$(GREEN)Formatting Python...$(RESET)"
	cd core && black .
	cd core && isort .

python-lint: ## 仅运行 Python lint
	cd core && ruff check .
	cd core && mypy openqc

python-format: ## 仅格式化 Python 代码
	cd core && black .
	cd core && isort .

##@ 文档

docs-dev: ## 启动文档开发服务器
	npm run docs:dev

docs-build: ## 构建文档
	npm run docs:build

docs-preview: ## 预览构建的文档
	npm run docs:preview

##@ Git Worktree 工作流（AGENTS.md）

WORKTREE_PREFIX := ../OpenQC-VSCode-

worktree-list: ## 列出所有 worktree
	@echo "$(BLUE)Current worktrees:$(RESET)"
	@git worktree list

worktree-create: ## 创建新的 worktree（用法: make worktree-create TYPE=feat NAME=new-feature）
	@if [ -z "$(TYPE)" ] || [ -z "$(NAME)" ]; then \
		echo "$(RED)Error: TYPE and NAME are required$(RESET)"; \
		echo "用法: make worktree-create TYPE=feat NAME=new-feature"; \
		echo "TYPE: feat, fix, docs"; \
		exit 1; \
	fi
	@echo "$(GREEN)Creating $(TYPE)/$(NAME) worktree...$(RESET)"
	git checkout master
	git pull origin master
	git worktree add -b $(TYPE)/$(NAME) $(WORKTREE_PREFIX)$(TYPE)-$(NAME) master
	@echo "$(GREEN)Worktree created at: $(WORKTREE_PREFIX)$(TYPE)-$(NAME)$(RESET)"
	@echo "$(YELLOW)进入目录: cd $(WORKTREE_PREFIX)$(TYPE)-$(NAME)$(RESET)"

worktree-clean: ## 清理已合并的 worktree（用法: make worktree-clean TYPE=feat NAME=new-feature）
	@if [ -z "$(TYPE)" ] || [ -z "$(NAME)" ]; then \
		echo "$(RED)Error: TYPE and NAME are required$(RESET)"; \
		echo "用法: make worktree-clean TYPE=feat NAME=new-feature"; \
		exit 1; \
	fi
	@echo "$(YELLOW)Removing worktree $(WORKTREE_PREFIX)$(TYPE)-$(NAME)...$(RESET)"
	@rm -rf $(WORKTREE_PREFIX)$(TYPE)-$(NAME)
	git worktree prune
	@git branch -d $(TYPE)/$(NAME) 2>/dev/null || echo "$(YELLOW)Branch $(TYPE)/$(NAME) already deleted or not found$(RESET)"
	@echo "$(GREEN)Worktree cleaned$(RESET)"

##@ 清理

clean: ## 清理构建输出和缓存
	@echo "$(YELLOW)Cleaning build outputs...$(RESET)"
	rm -rf out/
	rm -rf coverage/
	rm -f .coverage
	rm -f openqc-vscode-*.vsix
	cd core && rm -rf .coverage htmlcov/ .pytest_cache/
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	@echo "$(GREEN)Clean completed$(RESET)"

clean-all: clean ## 清理所有生成文件和依赖（完全重置）
	@echo "$(RED)Cleaning all generated files and dependencies...$(RESET)"
	rm -rf node_modules/
	rm -rf core/openqc.egg-info/
	rm -rf core/build/
	rm -rf core/dist/
	rm -rf .ruff_cache/
	cd core && pip uninstall -y openqc || true
	@echo "$(GREEN)Full clean completed$(RESET)"

clean-vscode: ## 仅清理 VS Code 相关输出
	rm -rf out/
	rm -rf coverage/
	rm -f openqc-vscode-*.vsix

clean-python: ## 仅清理 Python 相关缓存
	cd core && rm -rf .coverage htmlcov/ .pytest_cache/ build/ dist/
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true

##@ 发布

vscode-publish: ## 发布 VS Code 扩展（需要令牌）
	vsce publish

##@ 快捷命令

dev: watch ## 启动开发模式（监听文件变化）

check: lint typecheck format-check test ## 运行完整检查（lint + typecheck + 格式检查 + 测试）

pr-ready: format lint test ## 准备 PR（格式化 + lint + 测试）
	@echo "$(GREEN)PR is ready!$(RESET)"

# repo-governance-kit:standard-targets-v1
.PHONY: install format lint typecheck test check cleanup-merged

typecheck:
	bash scripts/typecheck.sh

cleanup-merged:
	bash scripts/cleanup_merged_worktrees.sh
