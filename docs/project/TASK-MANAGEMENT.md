# Task Management Guide

> How we organize, track, and complete tasks in OpenQC-VSCode

## Overview

We use GitHub's built-in project management tools combined with a structured workflow to ensure efficient task tracking and completion. This guide covers our task management philosophy, tools, and processes.

---

## GitHub Projects Board

### Board Structure

We use a **Kanban-style** board with the following columns:

| Column | Purpose | WIP Limit |
|--------|---------|-----------|
| **Backlog** | Future tasks, not yet prioritized | ∞ |
| **Ready** | Prioritized tasks ready to start | 10 |
| **In Progress** | Currently being worked on | 3 per person |
| **Review** | PR submitted, awaiting review | 5 |
| **Testing** | Feature complete, needs testing | 5 |
| **Done** | Completed and merged | ∞ |

### Workflow Stages

```
Backlog → Ready → In Progress → Review → Testing → Done
   ↓                  ↓           ↓         ↓
   └──────────────────┴───────────┴─────────┘
              (Issues can return to previous stages)
```

---

## Issue Management

### Issue Types

We use **labels** to categorize issues:

| Label | Description | Color |
|-------|-------------|-------|
| `type/bug` | Something isn't working | 🔴 Red |
| `type/feature` | New feature request | 🟢 Green |
| `type/enhancement` | Improvement to existing feature | 🔵 Blue |
| `type/documentation` | Documentation improvements | 🟡 Yellow |
| `type/test` | Testing-related tasks | 🟣 Purple |
| `type/refactor` | Code refactoring | 🟠 Orange |

### Priority Labels

| Label | Description | Response Time |
|-------|-------------|---------------|
| `priority/critical` | Blocks release or critical bug | < 4 hours |
| `priority/high` | Important for current milestone | < 24 hours |
| `priority/medium` | Normal priority | < 1 week |
| `priority/low` | Nice to have | Backlog |

### Area Labels

| Label | Description |
|-------|-------------|
| `area/parser` | File parsing logic |
| `area/visualization` | 3D rendering and UI |
| `area/conversion` | Format conversion |
| `area/ai` | AI/ML features |
| `area/testing` | Test infrastructure |
| `area/docs` | Documentation |
| `area/ci` | CI/CD pipeline |

---

## Issue Templates

### Bug Report

```markdown
---
name: Bug Report
about: Report a bug in OpenQC-VSCode
title: '[BUG] '
labels: 'type/bug'
assignees: ''
---

## Description
Clear description of the bug.

## Steps to Reproduce
1. Open file '...'
2. Run command '....'
3. See error

## Expected Behavior
What should happen.

## Actual Behavior
What actually happens.

## Screenshots
If applicable, add screenshots.

## Environment
- OS: [e.g., Windows 10, macOS 13, Ubuntu 22.04]
- VSCode Version: [e.g., 1.85.0]
- Extension Version: [e.g., 0.1.0]
- File Format: [e.g., VASP POSCAR]

## Additional Context
Any other relevant information.

## Sample File
If possible, attach a sample file that triggers the bug.
```

### Feature Request

```markdown
---
name: Feature Request
about: Suggest a new feature for OpenQC-VSCode
title: '[FEATURE] '
labels: 'type/feature'
assignees: ''
---

## Problem Statement
What problem does this feature solve?

## Proposed Solution
Describe the feature you'd like.

## Use Case
How would this feature be used?

## Alternatives Considered
Other solutions you've considered.

## Additional Context
Any other relevant information or screenshots.

## Would you be willing to submit a PR?
[ ] Yes, I'd like to contribute this feature
```

### Documentation Improvement

```markdown
---
name: Documentation
about: Improve or add documentation
title: '[DOCS] '
labels: 'type/documentation'
assignees: ''
---

## What needs to be documented?
[ ] New feature
[ ] API reference
[ ] Tutorial/Guide
[ ] Example
[ ] Other: ___

## Current State
What's missing or incorrect?

## Proposed Changes
What should be added or changed?

## Location
Which file(s) need updates?
```

---

## Task Workflow

### 1. Creating a Task

1. **Choose the right template**
   - Use issue templates for consistency
   - Fill in all required fields

2. **Add appropriate labels**
   - Type (bug/feature/etc.)
   - Priority (critical/high/etc.)
   - Area (parser/visualization/etc.)

3. **Assign to a milestone** (if applicable)
   - Link to project phase or release

4. **Add to project board**
   - New issues start in **Backlog**

### 2. Starting Work

1. **Move issue to "Ready"**
   - Ensure task is well-defined
   - Dependencies are resolved

2. **Self-assign the issue**
   - Click "Assign yourself"

3. **Create a branch**
   ```bash
   git checkout -b feature/issue-123-vasp-parser
   ```

4. **Move to "In Progress"**
   - Update project board

### 3. Development

1. **Follow TDD workflow**
   - Write tests first
   - Implement feature
   - Ensure all tests pass

2. **Commit frequently**
   ```bash
   git commit -m "feat: add INCUT parameter validation (#123)"
   ```

3. **Reference issue in commits**
   - Use `#123` in commit messages

### 4. Pull Request

1. **Create PR**
   - Reference the issue: `Closes #123`
   - Fill in PR template
   - Request reviewers

2. **Move to "Review"**
   - Update project board

3. **Address feedback**
   - Make requested changes
   - Push new commits

4. **Ensure CI passes**
   - All tests green
   - Coverage maintained

### 5. Completion

1. **Merge PR**
   - Squash and merge (clean history)
   - Delete branch

2. **Move to "Testing"**
   - Manual testing if needed

3. **Close issue**
   - Automated by PR merge
   - Move to "Done"

---

## Project Milestones

### Current Milestone: v0.1.0 - Foundation

**Goal**: Basic parsing and syntax highlighting

**Issues**:
- [ ] #1: Setup VSCode extension scaffold
- [ ] #2: Configure build system
- [ ] #3: Setup Jest testing framework
- [ ] #4: Create INCAR parser
- [ ] #5: Create POSCAR parser
- [ ] #6: Add syntax highlighting for VASP

**Due Date**: Week 4

---

## Sprint Planning

### Sprint Cadence
- **Duration**: 2 weeks
- **Planning**: Every other Monday
- **Review/Demo**: Every other Friday

### Sprint Process

1. **Planning Meeting** (1 hour)
   - Review backlog
   - Prioritize tasks
   - Estimate effort
   - Assign to sprint

2. **Daily Standups** (async in GitHub)
   - What I did yesterday
   - What I'll do today
   - Blockers

3. **Sprint Review** (30 min)
   - Demo completed features
   - Gather feedback
   - Update backlog

4. **Retrospective** (30 min)
   - What went well
   - What could improve
   - Action items

---

## Estimation

### Story Points

We use Fibonacci sequence for estimation:

| Points | Effort | Example |
|--------|--------|---------|
| 1 | Trivial | Fix typo in docs |
| 2 | Small | Add new parameter to parser |
| 3 | Medium | Implement new command |
| 5 | Large | Add new file format support |
| 8 | Very Large | Integrate external tool |
| 13 | Epic | Complete feature phase |

### Velocity Tracking

- Track story points completed per sprint
- Use for future sprint planning
- Adjust estimates based on actuals

---

## Definition of Done

A task is **Done** when:

- [ ] Code implemented
- [ ] Unit tests written (≥ 80% coverage)
- [ ] Integration tests (if applicable)
- [ ] Documentation updated
- [ ] PR reviewed and approved
- [ ] CI/CD passes
- [ ] Manual testing completed
- [ ] Merged to main branch
- [ ] Issue closed

---

## Automation

### GitHub Actions

We automate the following:

1. **Issue Created**
   - Add to project board (Backlog)
   - Apply default labels based on title
   - Notify team via Slack

2. **PR Opened**
   - Move linked issue to "Review"
   - Run tests
   - Check coverage
   - Request review from code owners

3. **PR Merged**
   - Close linked issue
   - Move to "Done"
   - Delete branch
   - Update changelog

### GitHub Bots

- **@dependabot**: Dependency updates
- **@codecov**: Coverage reports
- **@stale**: Close inactive issues

---

## Metrics & Reporting

### Key Metrics

1. **Velocity**: Story points per sprint
2. **Cycle Time**: Time from "In Progress" to "Done"
3. **Lead Time**: Time from issue creation to completion
4. **Bug Rate**: Bugs per release
5. **Coverage**: Test coverage percentage

### Dashboards

- **GitHub Projects**: Task status
- **GitHub Insights**: Contribution stats
- **Codecov**: Coverage trends
- **Custom Dashboard**: Combined metrics

---

## Best Practices

### Do's ✅

- Keep issues small and focused
- Write clear, actionable descriptions
- Update issue status regularly
- Link related issues
- Use templates
- Estimate tasks before starting
- Document decisions in issues

### Don'ts ❌

- Create vague or oversized issues
- Ignore stale issues (close or reprioritize)
- Work on multiple high-priority tasks simultaneously
- Skip the review process
- Merge without tests
- Forget to update project board

---

## Tools & Integrations

### Primary Tools
- **GitHub Issues**: Task tracking
- **GitHub Projects**: Kanban board
- **GitHub Actions**: Automation
- **GitHub Discussions**: Q&A and ideas

### External Integrations
- **Slack**: Notifications
- **CodeCov**: Coverage tracking
- **Dependabot**: Dependency updates

---

## Getting Help

- **Documentation**: Check `docs/` folder
- **GitHub Discussions**: Ask questions
- **Slack**: #openqc-dev channel
- **Office Hours**: Weekly sync (Fridays 3pm)

---

## Templates

### Issue Checklist

```markdown
- [ ] Issue created with template
- [ ] Labels applied
- [ ] Milestone assigned
- [ ] Added to project board
- [ ] Estimated story points
- [ ] Dependencies identified
```

### PR Checklist

```markdown
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] Issue linked (Closes #123)
- [ ] CI passes
- [ ] Coverage maintained
- [ ] Reviewers assigned
- [ ] Ready for review
```

---

## Resources

- [GitHub Projects Documentation](https://docs.github.com/en/issues/planning-and-tracking-with-projects)
- [GitHub Actions](https://docs.github.com/en/actions)
- [Kanban Guide](https://www.atlassian.com/agile/kanban)

---

## Questions?

- Open an issue with label `question`
- Ask in GitHub Discussions
- Join our Slack channel

---

**Last Updated**: 2026-02-28
**Version**: 1.0.0
