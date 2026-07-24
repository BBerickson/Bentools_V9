# BenTools v9 — Smoke Test

Two layers of regression protection during the refactor:

## 1. Automated data-layer fingerprint (fast, runnable headless)

Exercises the perf-critical data functions in `R_scripts/functions.R` on `test_files/`
and compares a numeric fingerprint against a saved golden copy.

```bash
Rscript tests/regression_harness.R          # print current fingerprint
Rscript tests/regression_harness.R --save   # capture golden (do once, on a known-good tree)
Rscript tests/regression_harness.R --check  # PASS/FAIL vs golden (exit 1 on mismatch)
```

Run `--check` after every data-layer change (Phases 5 & 6). If a change is *intended*
to alter numeric output, re-run `--save` and note why in the commit.

Golden values on the pre-refactor baseline (branch `refactor`, first commit):
`load`: 4 sets, 66000 rows, 610 genes, 80 bins, score_sum 6886.13; `active`: 640 rows / 2 genes;
`apply_math_mean`: value_sum 28.282; `find_clusters`: 101 genes / 100 merges; `find_groups`: 101 genes.

## 2. Manual UI smoke path (the acceptance gate for reactivity/UI changes)

The automated harness does not cover Shiny reactivity or the UI, so run this by hand
after Phases 2, 3, 4, and 7. Launch:

```r
shiny::runApp()
```

Then walk the path — each step must succeed without console errors:

1. **Load Data**: Browse → select all of `test_files/A_test.matrix.gz`, `B_`, `C_`, `D_test.matrix.gz`.
   Progress bar completes; the loaded-files table populates; disabled tabs light up.
2. **Plot**: open the Plot tab → the lines-and-labels popup appears → **SET and Plot** → a line plot renders.
   - Toggle samples in the dropdown; **Update Plot** re-renders.
   - Change a cosmetic option (color, y-range, label) — plot updates (should be instant after Phase 4).
   - Try the smoothing / log2 / abs / AUC sub-panels and the Violin sub-panel.
3. **Compare Lists**: load `test_files/test_genelist.txt` as a gene list; it appears as a new list.
4. **Filter Tool**: pick a bin range + sample → **filter sum** → preview plot + gene count; new list appears back on Plot.
5. **Ratio Tool**: pick numerator/denominator samples → run → ratio plot renders.
6. **Cluster Tools**: pick a sample + cluster count → run → cluster plots + table; change cluster number re-renders.
7. **Groups Tools**: same shape as Cluster — run and confirm plots/table.
8. **CDF Tools**: run → ECDF plot renders.
9. **Norm data** and **Group data**: run each once; confirm no errors and outputs appear.
10. **Data Table View**: table renders; row filtering works.

A step "passes" if it produces the expected plot/table with no error modal or console error.
