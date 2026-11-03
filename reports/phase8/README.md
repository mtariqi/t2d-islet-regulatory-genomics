# Phase 8 deliverables

- `Phase8_Results_Report.docx`: manuscript-oriented report based on the user's completed Phase 8 terminal results.
- `create_phase8_figures.R`: generates three publication figure pairs (PDF + PNG) from the **existing local repository TSVs**. No figures have been calculated from original TSVs here because those files were not uploaded.
- `Independent_Validation_Plan.md`: dataset shortlist, prespecified hypothesis, donor-aware replication protocol and criteria.

## Run figures

Copy `create_phase8_figures.R` to the root of `~/t2d-islet-regulatory-genomics`, then execute:

```bash
cd ~/t2d-islet-regulatory-genomics
Rscript create_phase8_figures.R 2>&1 | tee reports/phase8/publication_figures.log
```

Requires `ggplot2`; expected results are under `figures/phase8/publication/`. Inspect the generated figures, particularly the gene order and nonredundant pathway labels. The script fails rather than silently guessing if essential files or columns are missing.
