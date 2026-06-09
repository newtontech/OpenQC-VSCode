GitHub Ops GitHub Fix CI merge current PRs with all CI green, and use parallel subagents to do gh issues in different .worktrees in different branches then PR and merge if all tests pass and CI green

1:13 AM

I’ll use github-ops first for the repo/PR inventory and gh-fix-ci for CI status interpretation. I’m staying read-only while we’re in planning mode: I can inspect GitHub, checks, issues, branches, and worktrees, then produce a decision-complete execution plan.





Explored 1 search, ran 3 commands

The local repo is on feat/97-crossplatform-commands and already has uncommitted changes, so any execution plan needs to preserve that work and use separate worktrees for issue work. GitHub auth is present, so I’m checking open PRs and issue queues read-only next.





Ran 4 commands

Open PRs split cleanly: #131, #128, and #127 are mergeable with all checks green; #126 has all checks green but GitHub reports mergeability as UNKNOWN; #130 and #129 have real failed checks. I’m checking the repo test surface and the active branch PR/issue link so the plan can define exact gates before any merge or subagent work.





Ran 4 commands

The active branch appears to be work for issue #97, but it has no PR yet. The repo’s PR gate is npm run ci:pr; that’s the local verification command I’ll use in the plan for every issue worktree before creating PRs.





Ran 4 commands

There’s a dependency chain in the issues: #97 is the foundation, then #98 and #100, then #101/#102. That means “parallel subagents” should be used only for independent branches once prerequisites are merged, otherwise they will conflict or stack on unmerged code.



I’ve got enough repo truth to avoid guessing on mechanics. The one high-impact product decision left is which issue slice the parallel agents should consume first, because the open queue includes both an LSP dependency chain and separate viewer/scientific-platform work.





Asked 2 questions

I’ll lock the plan around the LSP chain and squash merges. One more quick check is for local agent/worktree configuration so the plan matches this repo’s ECC setup instead of inventing a process.





Explored 1 file, 1 search, 1 list

Worked for 4m 41s

This checkout does not currently contain .codex/ or .agents/skills/, despite the AGENTS instructions describing them. The plan will therefore use plain git worktree plus GitHub CLI, and treat Codex multi-agent roles as optional only if configured by the execution surface.



Plan

OpenQC-VSCode PR Merge and LSP Issue Worktree Plan

Summary

First merge only open PRs that are non-draft, mergeable, and have every required GitHub check green.

Use squash merge and delete branches where GitHub allows it.

Then work the LSP issue chain, starting with the existing local feat/97-crossplatform-commands work for issue #97; only parallelize follow-up issues once their dependencies are merged.

Preserve the current dirty worktree and never mix unrelated issue work into it.

Key Changes

Merge candidates now: #131, #128, #127; recheck #126 because CI is green but GitHub currently reports mergeability as UNKNOWN.

Do not merge #130 or #129; both have failing checks and need separate CI investigation.

Finish #97 in the current branch, covering the command resolver, deprecated path aliases, command/args/env config, cross-platform executable checks, and removal of the fixed restart delay.

After #97 is merged, create isolated worktrees from updated master:

../OpenQC-VSCode-fix-98-language-feature-mode on fix/98-language-feature-mode

../OpenQC-VSCode-fix-100-remote-lsp on fix/100-remote-lsp

../OpenQC-VSCode-feat-101-lsp-status-ux on feat/101-lsp-status-ux, after #100

../OpenQC-VSCode-test-102-fake-lsp on test/102-fake-lsp, after #98

Execution Flow

Before each PR merge: run gh pr view <n> --json mergeable,isDraft,statusCheckRollup and gh pr checks <n>.

If all checks pass and mergeability is MERGEABLE, run gh pr merge <n> --squash --delete-branch.

If mergeability is UNKNOWN, refresh/requery; if it remains unknown or blocked, skip and report it rather than forcing.

For each issue branch: create a fresh worktree from updated origin/master, implement only that issue, run local gates, push, open PR, wait for CI, and squash merge only if all checks are green.

Use parallel subagents only for independent branches; do not run #101 before #100, or #102 before #98.

Test Plan

Local gate for every issue branch: npm run ci:pr.

Targeted checks:

#97: npm run compile, npm run typecheck, npm run test:unit -- LSPManager, npm run test:unit -- commandResolver

#98: npm run test:unit -- extension, plus provider/mode tests

#100: compile/typecheck plus manual remote-safety doc review

#101: npm run test:unit -- LSPManager and command/status tests

#102: npm run test:integration

After each PR opens: wait for GitHub Actions and merge only if CI is green.

Assumptions

The default branch is master.

Merge policy is squash merge.

The first issue batch is the LSP chain.

The current uncommitted branch is intended for issue #97.

This checkout lacks .codex/ and .agents/skills/, so execution should use plain worktrees and GitHub CLI unless a separate Codex multi-agent surface is available.
