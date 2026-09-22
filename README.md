# BMAT-ID

Analysis code for the site-resolved transcriptomics of **human bone marrow
adipocytes**: how much of the transcriptional difference between anatomical
sites is a property of the adipocyte itself, and how much is the composition of
the tissue it was isolated from.

Every number, figure and table of the accompanying manuscript traces to a script
below. Three scripts render a figure from a released, machine-readable table
rather than from numbers typed into the source (`51_*`, `52_*`, `53_*`), so the
chain raw data → result table → figure is unbroken and can be re-run end to end
from public data. The released figures and result tables ship alongside the code.

- **Bulk data analysed:** [GSE291355](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291355) — 24 libraries, 4 donors × 3 anatomical sites × 2 culture states (freshly isolated *in vivo* and donor-matched *in vitro* differentiated adipocytes), fully paired.
- **Primary analysis:** threshold-free variance decomposition (`variancePartition` + `lme4`) of the anatomical-site signal, and its two-engine sensitivity (`limma-voom`, `DESeq2`).
- **Composition:** hand-defined marker axes and signatures derived from an independent single-cell reference, used to ask whether the site term survives adjustment.
- **Headline result:** the site signal is large *in vivo*, is quantitatively abolished *in vitro*, and *in vivo* is not separable from cell composition — the composition covariates absorb it.

---

## Repository layout

```
BMAT-ID/
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── data/
│   ├── README.md          # what to download and where to put it
│   ├── raw/               # GSE291355 counts go here (not tracked)
│   └── processed/         # written at run time (not tracked)
├── scripts/               # 41 numbered analysis scripts
├── results/
│   ├── figures/           # released figures, PNG + PDF
│   └── tables/            # released result tables
└── logs/                  # written at run time (not tracked)
```

`results/` in this repository is the **released** output: the figures as they
appear in the manuscript and the machine-readable tables behind them. Running
the pipeline rewrites those two directories in place, so a fresh run reproduces
the shipped files; a few scripts additionally write diagnostic tables that are
not part of the released set (listed under *Notes and limitations*).

---

## Running the pipeline

### Paths — nothing to edit

Every script finds the project root from its own location (the parent of
`scripts/`), for R and Python alike:

```bash
cd BMAT-ID          && Rscript scripts/21_varpart_main_models.R   # works
cd BMAT-ID/scripts  && Rscript 21_varpart_main_models.R           # works
Rscript /abs/path/to/BMAT-ID/scripts/21_varpart_main_models.R     # works
```

Set `BMAT_ID_DIR` to override, if you keep the scripts somewhere else:

```bash
export BMAT_ID_DIR=/path/to/BMAT-ID
```

`logs/`, `results/tables/` and `results/figures/` are created automatically on
first run.

### External inputs

Nothing is bundled. `data/README.md` describes the two downloads, the exact URLs,
and the environment variables (`BMAT_ID_DIR` for the project root,
`BMAT_REF_RDS` and `BMAT_REF_FEATURES` for the reference atlas) the scripts read.
Only `30_*` and `31_*` need the reference variables; `32_*` reads the output of
`31_*` and runs without them.

### Requirements

**R 4.6.1** — the variance decomposition is the primary analysis and its
rounding is what the manuscript reports:

| package | version |
|---|---|
| limma | 3.68.x |
| variancePartition | 1.42.0 |
| lme4 | 2.0.6 |
| DESeq2 | 1.52.0 |
| BiocParallel | 1.46.0 |
| Seurat | 5.5.1 (reference processing only) |
| Matrix | 1.7.5 |

**Python 3.13** — `pip install -r requirements.txt` (numpy, scipy, pandas,
matplotlib).

Two steps reach out to the internet: `08_*` (MyGene.info, symbol → Ensembl) and
`26_*`/`29_*`/`40_*`/`41_*`/`42_*`/`50_*` (Enrichr, and the same MyGene mapping
cached by `50_*`).

### Run order

The numbering is the run order, with one exception: the single step added after
the numbering was frozen carries a letter suffix, and `07b_*` runs between `07_*`
and `08_*`. Within a stage the scripts are independent unless noted.

#### Stage 1 — data, diagnostics and engine comparison

| script | what it does | key output |
|---|---|---|
| `01_load_and_annotate.R` | loads the counts, applies the CPM filter, fits limma-voom with quantile normalisation | `results/tables/BMATID_limma_allgenes.csv` |
| `02_site_identity_quantification.py` | how much of the site-identity signal the culture removes, and which site it drifts towards | `logs/quant_out.txt` |
| `03_site_identity_symbols.py` | top site-identity genes → symbols, grouped EPI / META / SCAT-high | `logs/symbols_out.txt` |
| `04_sequencing_depth_diagnostics.R` | library sizes, within-group correlation, p-value distributions, limma-trend robustness | `logs/diag_out.txt` |
| `05_partial_F_validation.py` | checks the partial-F implementation against limma — **run after `01_*`** | `logs/diag2_out.txt` |
| `06_immunoglobulin_origin_qc.py` | are the immunoglobulin genes a property of the adipocyte or of contaminating haematopoietic cells? | `logs/qc_out.txt` |

#### Stage 2 — annotation, composition axes and the residual schemes

| script | what it does | key output |
|---|---|---|
| `07_primary_annotation.R` | p-value annotation of the primary samples, CD45 co-variation | `BMATID_primary_cd45_annotated.csv`, `BMATID_residual_site_genes.tsv` |
| `07b_residual_scheme_tables.R` | the two CD45-decontaminated residual-site gene schemes, from the raw counts (caliber A: filter then fit; caliber B: fit then filter) | `BMATID_residual494_schemeA.tsv`, `BMATID_residual1248_schemeB.tsv` |
| `08_axis_marker_ensembl_map.py` | marker sets → Ensembl ids via MyGene.info (**network**) | `BMATID_marker_ensg.tsv` |
| `09_composition_axes.R` | axis scores and per-gene axis correlations — **run after `07_*`, `07b_*` and `08_*`** | `BMATID_axis_scores.rds`, `BMATID_axis_correlations.rds` |
| `10_axis_definition_robustness.R` | re-defines the axes (drop immunoglobulins, drop collagens, MEPE/DMP1/PHEX only) | `logs/axis_robust_out.txt` |

#### Stage 3 — the two engines

| script | what it does | key output |
|---|---|---|
| `11_deseq2_site_effect.R` | DESeq2 LRT `~ donor + localization` vs `~ donor` | `BMATID_deseq2_diff_residual_genes.csv` |
| `12_deseq2_effect_sizes.R` | pairwise Wald tests with `lfcShrink(type = "normal")` | `BMATID_deseq2_primary_shrunk.csv`, `BMATID_deseq2_diff_shrunk.csv` |
| `13_deseq2_composition_attribution.R` | refits the LRT after removing axis-co-varying genes | `BMATID_deseq2_axis_attribution.csv` |
| `14_deseq2_outlier_robustness.R` | outlier and low-count robustness of the *in vitro* residual genes | `logs/deseq2_robust_out.txt` |
| `15_deseq2_residual_symbols.py` | symbols for the *in vitro* residual genes | `logs/deseq2_residual_symbols.txt` |
| `16_limma_threshold_caliber.R` | re-derives the counts under the original publication's threshold (FDR < 0.05, \|log2FC\| > 2); reads the raw counts only | `logs/caliber_out.txt` |
| `17_limma_matched_caliber.py` | matched-caliber recount from the same starting point | `logs/caliber_head.txt` |
| `18_contamination_attribution_ladder.py` | CD45 co-variation fingerprint, and site-dependent counts with and without it | `logs/clean_out.txt` |
| `19_contamination_filter_limma.R` | the same removal ladder in the limma arm, with the equal-number random control (fixed seed) | `logs/clean_out2.txt`, `results/tables/BMATID_ladder_data.json` |
| `20_limma_axis_removal_sensitivity.R` | limma sensitivity under axis removal | `logs/sens_out.txt` |

#### Stage 4 — variance decomposition (the primary analysis)

| script | what it does | key output |
|---|---|---|
| `21_varpart_main_models.R` | the three mixed models `~ (1\|donor) + (1\|localization)`: *in vivo*, *in vitro*, and *in vivo* + composition covariates; plus the VST sensitivity. Accepts `smoke [n]` for a fast trial run. **≈ 50 min per model at full size** | `BMATID_varpart_invivo.csv`, `BMATID_varpart_invitro.csv`, `BMATID_varpart_invivo_comp.csv`, `BMATID_varpart_inputs.rds` |
| `22_varpart_collinearity.py` | canonical correlations and η² between site and the axis scores | `logs/varpart_collinearity_out.txt` |
| `23_varpart_permutation_2x2.R` | 2 × 2 control: with and without composition covariates, true versus shuffled site labels | `logs/varpart_m6b_out.txt` |
| `24_varpart_low_depth_sensitivity.R` | drops the single low-depth library and refits | `logs/varpart_m5b_out.txt` |
| `25_varpart_surviving_genes.py` | the strong-site genes that survive composition adjustment | `BMATID_varpart_survivors.csv` |
| `26_varpart_residual_genes.py` | residual site genes after adjustment, against the permutation null | `BMATID_varpart_residual.csv` |
| `27_varpart_diagnostics.py` | diagnostics for the bimodal variance allocation | `logs/varpart_diag_out.txt` |
| `28_varpart_summary.py` | the headline quantities, machine-readable | `BMATID_varpart_summary.json` |
| `29_varpart_gene_set_enrichment.py` | Enrichr on the survivor gene set (**network**) | `BMATID_varpart_enrichment_survivors233.csv` |

#### Stage 5 — reference-anchored signatures

Requires `BMAT_REF_RDS` and `BMAT_REF_FEATURES` (see `data/README.md`).

| script | what it does | key output |
|---|---|---|
| `30_reference_signature_export.R` | pseudo-bulk profiles per cluster, cluster composition, marker scoring on the reference atlas | `BMATID_ref_pseudobulk_clusters.csv` + 3 more |
| `31_reference_anchored_scores.py` | data-driven signatures per super-type and the relative scores of each bulk library; also attempts the absolute (NNLS) deconvolution | `BMATID_deconv_L1_scores.csv`, `BMATID_ref_supertype_specificity.csv`, `BMATID_deconv_proportions.csv` |
| `32_varpart_reference_anchored.R` | the mixed model with reference-anchored covariates instead of the hand-built axes. **≈ 55 min** | `BMATID_varpart_M7_reference.csv` |

#### Stage 6 — modules and enrichment

| script | what it does | key output |
|---|---|---|
| `40_functional_modules.py` | *in vitro*/*in vivo* effect-size retention per functional module — **run after `01_*`** | `logs/modules_out.txt`, `results/tables/BMATID_module_retention.csv` |
| `41_enrichr_residual_genes.py` | Enrichr on the residual sets (**network**) | `BMATID_enrichment_residual494.csv`, `BMATID_residual494_symbols.json` |
| `42_enrichr_sensitivity_genes.py` | Enrichr on the alternative-caliber residual set (**network**) | `BMATID_enrichment_residual1248.csv` |

#### Stage 7 — figures and panel data

| script | what it does | key output |
|---|---|---|
| `50_fig1_panel_data.py` | maps every gene to a symbol and merges the four variance tables (**network** for the mapping, cached) — **run after `07_*`, `21_*`, `25_*`** | `BMATID_varpart_allgenes_merged.csv`, `BMATID_allgenes_symbols.csv` |
| `51_fig1_site_composition.py` | **Figure 1** (a–d): variance composition across the four models, per-gene ECDF, positive-control panel, surviving genes | `BMATID_Fig1_site_composition.{png,pdf}` |
| `52_fig2_attribution_ladder.py` | **Figure 2**: the attribution ladder and its specificity control; reads `BMATID_ladder_data.json` | `BMATID_Fig2_attribution_ladder.{png,pdf}` |
| `53_fig3_module_retention.py` | **Figure 3**: what survives adipogenic culture, by module; reads `BMATID_module_retention.csv` | `BMATID_Fig3_module_retention.{png,pdf}` |
| `54_figS1_reference_scores.py` | **Figure S1**: per-library reference-anchored scores | `BMATID_FigS1_reference_scores.{png,pdf}` |

### Reproducing only the figures

The released figures in `results/figures/` were rendered from the released
tables, so the four figure scripts (`51`–`54`) can be re-run on their own
without fitting a single mixed model — every input they need is already in
`results/tables/`. `52_*` needs `BMATID_ladder_data.json` (shipped) and `53_*`
needs `BMATID_module_retention.csv` (shipped); both fail with a clear message if
the file is missing.

---

## Expected values

If you re-run the pipeline, these are the quantities to compare against
(`28_varpart_summary.py` writes them as json). They are the numbers quoted in
the manuscript.

| quantity | value | produced by |
|---|---|---|
| genes retained (CPM ≥ 1 in ≥ 4 samples) | 20,073 | `01_*` |
| genes in the fitted universe (20,073 minus 32 all-zero *in vivo* genes) | 20,041 | `21_*` → `28_*` |
| mean site-attributable variance, *in vivo* | 17.1% | `21_*` → `28_*` |
| median, *in vivo* | 6.3% | `21_*` → `28_*` |
| mean, *in vitro* | 3.6% | `21_*` → `28_*` |
| mean, *in vivo* + composition covariates | 4.2% | `21_*` → `28_*` |
| genes whose site term is exactly zero after adjustment | 77.3% (15,501 of 20,041) | `21_*` → `28_*` |
| the same quantity at a 1e-9 tolerance (Figure 1 legend only) | 91.7% | `21_*` → `28_*` |
| genes above 50% *in vivo* → surviving above 50% | 2,059 → 233 | `21_*` → `25_*` |
| VST sensitivity of the *in vivo* estimate | 18.2% | `21_*` |
| composition terms, *in vivo* + axes | haematopoietic 20.2%, bone 12.7% | `21_*` |
| canonical correlation, site vs axes | 0.887 (haematopoietic), 0.944 (bone) | `22_*` |
| permutation null *without* composition covariates: observed / null median | 17.3% vs 3.1% | `23_*` |
| permutation null *with* composition covariates: observed / null range (median) | 4.0% inside 2.4–14.0% (median 5.0%) | `23_*` |
| contamination ladder (site genes left; limma arm) | 2,712 → 1,222 / 494 / 207 / 128 | `19_*` |
| equal-number random removal (same arm) | 2,320 / 1,777 | `19_*` |
| candidate signature genes, immune super-type | 1,165 → top 200 by specificity | `31_*` |
| correlation, immune signature vs CD45 | +0.972 | `31_*` |
| reference-anchored model: site / exact-zero | 7.2% / 66.8% | `32_*` |
| module retention (*in vitro*/*in vivo* effect size) | adipogenic core 62.9%, haematopoietic/plasma-cell 3.5%, bone/mineralisation 10.6%, MSC/stemness 12.3%, haematopoietic niche 13.8%, endothelial 96.4% | `40_*` |

The ladder in `19_*` uses `set.seed(42)`, so the equal-number random control is
reproducible; `20,041` is the universe the mixed models are fitted on, whereas
`20,073` is what survives the CPM filter.

---

## Figures

| file | content |
|---|---|
| `BMATID_Fig1_site_composition` | (a) variance composition across four mixed models; (b) ECDF of per-gene site variance, with the exact-zero fraction (and the 1e-9 comparison) marked; (c) positive control — axis-patterning transcription factors keep their site variance, haematopoietic/bone/adipogenic markers lose it; (d) survival of the strongest site genes |
| `BMATID_Fig2_attribution_ladder` | site-dependent gene count as a function of how many CD45-co-varying genes are removed, with the equal-number random control |
| `BMATID_Fig3_module_retention` | *in vitro*/*in vivo* effect-size retention by functional module, and the per-gene *in vivo* versus *in vitro* F comparison |
| `BMATID_FigS1_reference_scores` | reference-anchored score of every library, and the donor-matched *in vivo* versus *in vitro* comparison |

Each figure ships as PNG (400 dpi for Figure 1, 300 dpi for the rest) and
vector PDF. The PNGs are byte-reproducible; the PDFs differ between runs only in
the `/CreationDate` field that the matplotlib PDF backend writes into every file.

## Result tables

`results/tables/` holds the released tables, in five groups:

- **Variance decomposition** — `BMATID_varpart_invivo.csv`, `BMATID_varpart_invitro.csv`, `BMATID_varpart_invivo_comp.csv` (per-gene variance share of every term, the three main models), `BMATID_varpart_M7_reference.csv` (reference-anchored model), `BMATID_varpart_summary.json` (headline quantities), `BMATID_varpart_inputs.rds`.
- **Survivors and residuals** — `BMATID_varpart_survivors.csv`, `BMATID_varpart_residual.csv`, `BMATID_residual494_schemeA.tsv`, `BMATID_residual1248_schemeB.tsv`, `BMATID_residual_site_genes.tsv`, plus the Enrichr outputs.
- **Engines** — `BMATID_limma_allgenes.csv`, `BMATID_deseq2_primary_shrunk.csv`, `BMATID_deseq2_diff_shrunk.csv`, `BMATID_deseq2_axis_attribution.csv`, `BMATID_deseq2_diff_residual_genes.csv`.
- **Composition** — `BMATID_axis_scores.rds`, `BMATID_axis_correlations.rds`, `BMATID_marker_ensg.tsv`, `BMATID_primary_cd45_annotated.csv`, `BMATID_allgenes_symbols.csv`.
- **Reference atlas** — `BMATID_ref_pseudobulk_clusters.csv`, `BMATID_ref_cluster_by_group.csv`, `BMATID_ref_cluster_markerscores.csv`, `BMATID_ref_axismarkers_log2cpm.csv`, `BMATID_ref_supertype_specificity.csv`, `BMATID_deconv_L1_scores.csv`, `BMATID_deconv_proportions.csv`.
- **Figure inputs** — `BMATID_ladder_data.json` (Figure 2), `BMATID_module_retention.csv` (Figure 3).

A few scripts report their result as a plain-text table in `logs/` rather than
as a csv (the diagnostics of `02_*`, `03_*`, `04_*`, `06_*`, `10_*`, `14_*`,
`16_*`, `18_*`, `20_*`, `22_*`, `23_*`, `24_*`, `27_*`, `40_*`); the numbers
quoted from those analyses are printed there.

---

## Notes and limitations

These are carried over from the manuscript and matter for anyone re-using the
code.

1. **Absolute deconvolution is not identifiable.** `31_*` also fits a
   non-negative least-squares mixture (written to
   `BMATID_deconv_proportions.csv`), but the reference is 3′-biased 10x data
   while the bulk libraries are full-length; the residuals are large
   (0.53–0.97). Only the relative scores are reported and only they should be
   used.
2. **The reference atlas cannot anchor the bone compartment.** GSE169396 is
   overwhelmingly immune (15,508 of 18,205 cells; 19 of 26 clusters), with a
   single osteo-chondral cluster and no pure MSC/stromal cluster. The skeletal
   signature is therefore a weak anchor, and this is reported as a limitation.
3. **Gene-level variance shares are not interpretable** once the composition
   covariates are in the model: the covariates are strongly collinear with
   site *in vivo* (canonical correlation 0.887 and 0.944), so the site term and
   the composition terms cannot be separated per gene. Only the aggregate site
   term is reported, and survivors are given as categories rather than as a
   ranking.
4. **Bordering on the released set.** Three scripts write further tables that are
   *not* shipped, because they are per-contrast or alternative-caliber variants
   rather than the figures and tables of the manuscript: `12_*` writes one
   `..._primary_shrunk_<contrast>.csv` and one `..._diff_shrunk_<contrast>.csv`
   per pairwise contrast, `26_*` writes `BMATID_varpart_residual_enrichment.csv`
   for the alternative residual set, and `29_*` writes
   `BMATID_varpart_enrichment_residual1591.csv` when run on the second gene set.
   `18_*` additionally caches `data/processed/noncontam_top.npy` (not tracked).
5. **The scripts log in Chinese.** The external documents (`README.md`,
   `data/README.md`, `CITATION.cff`) are in English, but most of the per-script
   log text in `logs/` is written in Chinese — it is the working record of the
   analysis rather than documentation.
6. **Script numbers are not contiguous.** The analysis grew in stages, so
   `scripts/` runs 01–32, 40–42 and 50–54 (33–39 and 43–49 were never used).
   The gaps carry no meaning beyond the order in which the stages were added;
   the one step added after the numbering was frozen is `07b_*`.

## Citation

If you use this code, please cite the manuscript and this repository
(see `CITATION.cff`).

## License

MIT — see `LICENSE`. You are free to use, modify and redistribute the code
and tables, including commercially, provided the copyright notice and the
license text are kept. The work is provided as is; no warranty of any kind.
Citation of the manuscript and of this repository is requested, not required.

## Contact

Qingqing Li — liqingqing_1999@163.com
Jiangsu Province Hospital (The First Affiliated Hospital with Nanjing Medical University)
