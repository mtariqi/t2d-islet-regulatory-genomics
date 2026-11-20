# Deviations from pre-specified plans

| Date | Change | Reason |
|---|---|---|
| 2026-11-19 | v1 full-pathway proteasome comparator not run | The plan did not name a specific Reactome pathway; choosing one after data access would add an undisclosed analytic choice. The frozen 29-gene module remains the v1 primary test. |
| 2026-11-19 | GSE153855 annotation parsing: strip CR characters before reading | Final annotation line uses CR separators and was split into three rows; detected by the row-count check before any results were produced. Affects one non-beta cell only. |
