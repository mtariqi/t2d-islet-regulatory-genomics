#!/usr/bin/env python3
"""Update the README hypothesis table: H3 and H5 status, new H6 row, and the note below it.

Run from the repository root:
    python3 tools/update_readme_status.py           # preview
    python3 tools/update_readme_status.py --apply   # write README.md
Safe to rerun: does nothing if the update is already present.
"""
import pathlib, sys

H3_ADD = " Disease × chemistry interaction test planned."
H5_STATUS = "Replication pre-registered (tag `prereg-validation-v1`); not yet run. Functional experiments planned."
H6_CELLS = ["H6",
            "Robust T2D-associated genes converge on coherent biological pathways.",
            "Ranked fgsea (Reactome, GO BP) across three models; redundancy clustering; leading-edge analysis.",
            "Evaluated; excitability and potassium-channel genes higher in T2D; a 29-gene proteasome–ubiquitin "
            "module lower in T2D (small, chemistry-sensitive effect)."]
OLD_NOTE = "H4 and H5 remain untested."
NEW_NOTE = ("H6 is exploratory; the proteasome–ubiquitin module was selected from discovery data "
            "and is tested prospectively under H5. H4 remains untested.")

def cells(line):
    parts = line.strip().split("|")
    return [c.strip() for c in parts[1:-1]]

def row(c):
    return "| " + " | ".join(c) + " |"

def main():
    apply = "--apply" in sys.argv
    p = pathlib.Path("README.md")
    if not p.exists():
        sys.exit("README.md not found; run from the repository root.")
    lines = p.read_text(encoding="utf-8").split("\n")
    if any("prereg-validation-v1" in l for l in lines):
        print("README already updated; nothing to do."); return
    idx = {k: next((i for i, l in enumerate(lines) if l.strip().startswith(f"| {k} ")
                    or l.strip().startswith(f"|{k}|")), None) for k in ("H3", "H5")}
    for k, i in idx.items():
        if i is None:
            sys.exit(f"Could not find the {k} row (a line starting with '| {k} |'). README unchanged.")
    ncol = len(cells(lines[idx["H5"]]))
    if ncol != 4:
        sys.exit(f"Expected 4 columns in the hypothesis table, found {ncol}. README unchanged.")

    c3 = cells(lines[idx["H3"]]); old3 = c3[-1]; c3[-1] = c3[-1].rstrip() + H3_ADD
    c5 = cells(lines[idx["H5"]]); old5 = c5[-1]; c5[-1] = H5_STATUS
    lines[idx["H3"]] = row(c3)
    lines[idx["H5"]] = row(c5)
    lines.insert(idx["H5"] + 1, row(H6_CELLS))

    text = "\n".join(lines)
    note_found = OLD_NOTE in text
    text = text.replace(OLD_NOTE, NEW_NOTE, 1)

    print(f"H3 status: '{old3}'\n        -> '{c3[-1]}'")
    print(f"H5 status: '{old5}'\n        -> '{c5[-1]}'")
    print("H6 row: added after H5")
    print("Note under table: " + ("updated" if note_found else
          "NOT FOUND (edit by hand: replace 'H4 and H5 remain untested.')"))
    if apply:
        p.write_text(text, encoding="utf-8"); print("\nREADME.md written.")
    else:
        print("\nPreview only. Run with --apply to write README.md.")

if __name__ == "__main__":
    main()
