# Contributing to HARD

## Workflow (fork model)

We work in a **fork model**. Do not push branches to the main repo.

1. Fork the project into your own namespace.
2. Create a branch in your fork and commit your work there.
3. Open a Merge Request against the `main` branch of the main repo.
4. **Assign the group** to your MR so reviewers can pick it up.
5. Every MR needs **at least one review** before it can be merged.

## Project structure

Keep physics and applications separated:

- **`modules/`** — the core physics tasks, reusable across applications.
- **`apps/`** — one directory per application (e.g. `hydro`, `radhydro`), wiring
  modules together.

**Avoid code duplication.** Duplicated code is hard to maintain and
error prone. Factor shared apps logic into `actions/` instead of copying it between
apps.

## Code style

- Use **snake_case** for all names (variables, functions, files).
- Header files use `.hh`, source files use `.cc`.
- Guard headers with `#ifndef <NAME>_HH` / `#define <NAME>_HH`.
- Namespaces mirror the module (e.g. `spec::io`).

## Formatting

CI rejects unformatted code, so format **before** committing — don't add commits
just to fix formatting afterward.

We use **clang-format 20** (`20...<21`). Check your version:

```bash
clang-format --version
```

Run the formatter from a build directory configured with `ENABLE_FORMAT=ON`:

```bash
make format
```

Commit formatting together with the related work, not as a separate commit —
this keeps the history readable and makes reviewing easier.

## Citing papers

When using papers for code development (e.g., adapting a recent algorithm) you
can cite them in the code directly, adding a link to the `doi` database or to
the preprint itself if not yet published. Please, also cite the paper in place
so that other people can find the work even if the link breaks.

Optionally, the paper can be added to the code documentation if it does not
apply to only an easily identifiable part of the code.
