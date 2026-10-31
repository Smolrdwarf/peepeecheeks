# peepeecheeks
for acacia hackerthon

## Notes for agents working in this repo

Several agents may be working here at the same time. To avoid stepping on each other:

- **Do not commit work directly to `main`.** Use your own branch and merge or PR when you're done. Only coordination docs like this one go straight to `main`.
- **Add yourself to the table below** (commit + push to `main`) before you start, and update your row when your status changes.
- `git pull --rebase` before you push. Never force-push to someone else's branch.
- Keep your work in its own top-level folder (e.g. `hfss/`) so merges stay conflict-free.

### Active work

| Agent / owner | Branch | Folder | Task | Status |
|---|---|---|---|---|
| Claude (Opus 5.5) | `hfss-re` | `hfss/` | Reverse engineering Ansys HFSS (file formats, solver behaviour, scripting interfaces) | In progress |
| Claude (Opus 5.5) #2 | `openems-validation` | `openems_validation/` | Install openEMS, build benchmark cases, compare against Ansys HFSS reference results | In progress |

### Environment / tooling (on the host machine)

- RE tools: `~/tools/re/ghidra_12.2_PUBLIC` (Ghidra), `~/tools/re/recaf` (Recaf, a Java bytecode editor)
- Ansys HFSS / AEDT is **not** installed on this machine as of 2026-10-31
- openEMS: being installed to `~/tools/openEMS` (see `openems_validation/README.md`)
- `gh` CLI is not installed. Use plain `git` against `origin` (github.com/Smolrdwarf/peepeecheeks)
