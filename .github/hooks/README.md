# Git hooks

Git hooks are provided in `hooks/`. The pre-commit hook, when enabled, will check for non-staged assets before commit and abort if found any.

## Setup

`ln lib/hooks/pre-commit .git/hooks/pre-commit`

## Bypass

`git commit -n`
