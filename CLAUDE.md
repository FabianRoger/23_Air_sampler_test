# CLAUDE.md — Project Context for AI-Assisted Development

## Quick Summary

This is an R-based analysis repository for a scientific study comparing 8 air sampling devices for airborne eDNA metabarcoding across 5 organism groups (Bacteria, Fungi, Plants, Insects, Vertebrates), sampled at 8 sites in Switzerland. The repository lives at https://github.com/FabianRoger/23_Air_sampler_test.

## Workflow Recommendations

- **Use a new chat for each discrete task** (e.g., one chat per script harmonization, one for figure work, one for manuscript edits). This keeps context manageable.
- **Use Sonnet for coding tasks** (writing/editing scripts, fixing bugs, reformatting).
- **Use Opus for planning tasks** (deciding analytical strategy, reviewing contamination approaches, discussing manuscript structure).
- **Read this file at the start of each new chat** to restore context.

---

## Study Design

### Samplers (8 total + Drone)

| Sampler | Type | Flow rate (L/min) | Durations |
|---------|------|-------------------|-----------|
| Kärcher | active | 3000 | 30 min, 5 hours |
| Coriolis | active | 300 | 30 min, 5 hours |
| Sass | active | 300 | 30 min, 5 hours |
| Hepa | active | 60 | 30 min, 5 hours |
| Burkhart | active | 16.5 | 30 min, 5 hours |
| Electrostatic | active | 10 | 30 min, 5 hours |
| MWAC | passive | — | 2 weeks |
| WSL_Filter | passive | — | 2 weeks |
| Drone | active (aerial) | 300 | 10 min (1 site only) |

### Sites (8)

Dübendorf, Talent, Coruz, Haslibach, Leugene, Gäbelbach, Weierbach, Riedgraben (greater Bern/Zürich area, Switzerland).

### Taxonomic Groups and Primers

| Group | Marker | Primer pair | Reference DB |
|-------|--------|-------------|-------------|
| Bacteria | 16S rRNA | 515F / 806R | RDP trainset 19 |
| Fungi | ITS2 | fITS7 / ITS4 | UNITE + GBIF |
| Plants | ITS2 | ITS-S2F / ITS4 | PLANiTS + GBIF |
| Insects | COI | mlCOIintF / jgHCO2198 | BOLD + MIDORI + GBIF |
| Vertebrates | 16S large ribosomal | 16Smam1 / 16Smam2 | MIDORI + GBIF |

Plants are also detected from the Fungi ITS2 primer (cross-amplification); results from both primers are merged and compared.

### Controls

- **Lab controls**: 9 extraction blanks + 3 PCR negative controls (per primer)
- **Field controls**: 2 per sampler type (deployed at Haslibach and Weierbach)
- 152 total samples in `Meta_data.tsv` (columns: Plate, row, col, Sample_names, Sites, Sampler, duration, Type, Index)
- `Type` values: active, passive, field_control, extr_control, PCR_control

---

## Repository Structure

```
├── Scripts/
│   ├── Prep_Metadata.Rmd              # Source → Meta_data.tsv
│   ├── Prepare_RefDB.Rmd              # Build ref DBs (DO NOT RUN — resource heavy)
│   ├── Demultiplexing_by_Primer.Rmd   # Demux raw reads
│   ├── Seq_Analysis_*.Rmd             # Per-group bioinformatics (DO NOT RUN)
│   ├── Seq_analysis_helper.R          # Shared functions
│   ├── Data_Analysis_Fungi.Rmd        # ★ TEMPLATE for all Data_Analysis scripts
│   ├── Data_Analysis_Bacteria.Rmd     # Harmonized (2025-02)
│   ├── Data_Analysis_Insects.Rmd      # Needs harmonization
│   ├── Data_Analysis_Vertebrates.Rmd  # Needs harmonization
│   ├── Data_Analysis_Plants.Rmd       # Needs harmonization (has unique sections)
│   └── Joined_Figures.Rmd             # Cross-group richness figures
├── Data/
│   ├── Meta_data.tsv                  # THE metadata file (all scripts must use this)
│   ├── sampler_comp_meta.xlsx         # Source workbook
│   ├── Barcodes_4_plates.xlsx         # Barcode plate layout
│   ├── DNA_conc_air_sampler_comp.txt  # DNA concentrations
│   ├── bbduk_stats/                   # Demultiplexing statistics
│   └── {Group}/                       # Per-group data folders (Fungi, Bacteria, etc.)
│       ├── demultiplexed/             # Raw demultiplexed reads (not in repo)
│       ├── filtered/                  # Quality-filtered reads (not in repo)
│       ├── in_files/                  # Small lookup files + large ref DBs (ignored)
│       ├── out_final/                 # Final processed data (what Data_Analysis reads)
│       │   ├── {Group}_ASV_glom.txt   # OTU count table (note: Bacteria = "ASW")
│       │   ├── {Group}_taxa_*_glom.txt # Taxonomy (80 = SINTAX, blast = BLAST)
│       │   ├── {Group}_Seq_glom       # Representative sequences (FASTA, no extension)
│       │   └── {Group}_tree.nwk       # Phylogenetic tree (placeholder, not yet built)
│       └── out_intermediate/          # Intermediate outputs & exports for cross-group figs
│           ├── {Group}_richness.tsv
│           ├── {Group}_contamination.tsv
│           ├── {Group}_completeness.tsv
│           ├── {Group}_nestedness_turnover.tsv
│           └── {Group}_rank_correlation.tsv
├── GENERAL_ANALYSIS_CODING_STYLE.md   # Coding style — MUST FOLLOW
├── README.md
├── .gitignore
├── .claudeignore
└── 23_Air_sampler_test.Rproj
```

### Data File Naming Quirks

- Bacteria ASV file is `Bacteria_ASW_glom.txt` (typo: ASW not ASV)
- Vertebrates taxonomy file is `Vertebrates_taxa_blast_glom.txt` (blast, not 80)
- Plants has extra files: `Plants_ASV_glom_fungi.txt` and `Plants_taxa_blast_glom_fungi.txt` (from fungi primer cross-amplification)
- Seq_glom files have no file extension (they are FASTA format)

---

## Coding Style (from GENERAL_ANALYSIS_CODING_STYLE.md)

- **Format**: RMarkdown (.Rmd), html_notebook output
- **Style**: tidyverse, pipes (`%>%`), ggplot2
- **Structure**: One chunk per task, one displayed output per chunk
- **Markdown**: Bullet-point style explanations before each chunk
- **Paths**: Use `here()` for all file paths
- **Outputs**: TSV for tabular data
- **Contamination**: Explicit reporting of what's removed and why
- **End**: Session info at end of script

---

## Contamination Assessment Strategy (Agreed Upon)

All Data_Analysis scripts must implement this standardized contamination approach:

### 1. Decontam (prevalence-based)

```r
# Field controls set to NA (excluded from decontam — they contain real species)
isContaminant(matrix, neg = blanks$is_neg, method = "prevalence")
```

### 2. Lab Controls (extraction blanks + PCR negatives)

- Compute mean and max reads per OTU across all lab controls
- Compute ratio = mean(lab control reads) / mean(sample reads)
- Visualize as heatmap
- For ALL OTUs found in lab controls: subtract `2 × max(lab control reads)` from ALL samples globally

### 3. Field Controls

- Compute ratio = mean(field control reads) / mean(sample reads) per OTU
- Only flag OTUs where ratio > 1 (more abundant in controls than samples)
- For flagged OTUs: subtract `2 × max(field control reads)` from ALL samples globally
- OTUs with ratio ≤ 1: keep untouched (likely environmental signal that happened to appear in controls)

### 4. Post-cleaning

- Set any negative values to 0
- Remove empty OTUs (colSums = 0)
- Remove control samples from dataset
- Remove samples with < 1000 reads

### Key Implementation (from Fungi template)

```r
Remove_reads_lab <- Lab_only %>% filter(name == "max_val") %>% select(seq, value)
Remove_reads_field <- Field_only %>%
  group_by(seq) %>% filter(min(ratio) > 1) %>% ungroup() %>%
  filter(name == "max_val") %>% select(seq, value)
Remove_reads <- bind_rows(Remove_reads_lab, Remove_reads_field) %>%
  group_by(seq) %>% filter(value == max(value)) %>%
  mutate(seq_remove = ceiling(value * 2))
```

---

## Data_Analysis Script Template — REVISED (Feb 2026)

The analysis plan was revised to focus on four core questions per species group:
1. **How many species?** — Richness by sampler × duration, formal model
2. **Does more air = more species?** — Richness vs volume, 30min vs 5h completeness
3. **Do samplers agree on what's there?** — Rank-abundance correlation + nestedness/turnover
4. **Do we see real site differences through all samplers?** — Richness correlation across sites

### Diversity Metrics Strategy

- **Primary: Observed richness (S)** — justified by rarefaction curves reaching saturation (sequencing depth is not the bottleneck; sampler capture is)
- **Secondary: Faith's PD** — phylogenetic diversity from FastTree NJ trees (built in Seq_Analysis scripts). Shown alongside S in main richness figure
- **Sensitivity (supplement): iNEXT coverage-based q=0** — confirms sampler rankings hold under coverage standardization
- **No Hill q=1/q=2** — read abundance is unreliable in metabarcoding (PCR bias, copy number variation), especially for multi-copy markers (ITS, 16S)

### Phylogenetic Trees

FastTree NJ trees from representative sequences. One tree per marker, built in the respective `Seq_Analysis_*.Rmd` script (code to be added but NOT run — resource heavy). Tree files stored as `Data/{Group}/out_final/{Group}_tree.nwk`. Data_Analysis scripts load trees with `ape::read.tree()`.

### Section Order (Template)

Each Data_Analysis script follows this structure:

1. **Setup**: libraries (including `betapart`, `ape`, `picante`, `iNEXT`), `primer <- "{Group}"`, load Meta_data.tsv
2. **Load data**: ASV table, taxonomy, sequences (Biostrings), DNA concentrations, phylogenetic tree (placeholder — commented out until trees are built)
3. **Fragment length distribution**: sequence length histogram, reads vs length
4. **Rarefaction (pre-cleaning)**: by major phyla/groups — demonstrates sequencing saturation, justifies use of observed S
5. **Clean data**:
   - Index hopping check (unused indices)
   - Remove known bad samples if any
   - Decontam (prevalence-based, field controls as NA)
   - Visualize control contamination ratios
   - Lab controls: heatmap of contaminant OTUs
   - Field controls: by-sampler contamination patterns (report contamination levels per sampler — this is itself a result). Export contamination summary to `Data/{Group}/out_intermediate/{Group}_contamination.tsv`
   - Exclude contaminants (2× max subtraction)
   - Remove empty OTUs, remove controls, check empty samples
6. **Subset for target group**: filter to target kingdom/phylum
7. **Prevalence**: histogram
8. **Species richness (S)**: observed S by sampler × duration × site. Faceted by Sampler, colored by Site, with duration on x-axis. Include GLM: `S ~ Sampler * duration + Sites` (or GLMM with `(1|Sites)` if appropriate). Report estimated richness and CIs. Include Burkhart comparison table (fold-differences)
9. **Faith's PD**: same plot structure as S, using `picante::pd()`. Shown as parallel panel to S figure. (Placeholder until trees are built)
10. **Richness vs air volume**: scatter of S vs log(air_vol) for active samplers + linear model fit. Report formal test
11. **Sampling completeness**: fraction of pooled 5h community captured per sampler × duration. Pool all samplers at each site for 5h as "reference community," compute fraction detected by each sampler × duration combination. Directly answers "what fraction of the detectable community do you capture in 30 min vs 5 hours?"
12. **Richness correlation across sites**: corrplot of S across sampler × duration combinations. Shows whether samplers agree on which sites are species-rich
13. **Richness by taxonomic subgroup** (group-specific): repeat richness analysis for major subdivisions to test whether sampler effects differ by subgroup. See group-specific details below
14. **Nestedness vs turnover (betapart)**: pairwise Sørensen decomposition across sampler pairs, averaged across sites. Shows whether low-richness samplers detect a subset of what high-richness samplers detect (nestedness) or genuinely different species (turnover). Display as matrix or grouped bar chart. Replaces the old NMDS + PERMANOVA section
15. **Rank-abundance correlation**: Spearman correlation on relative abundance (or rank) across all shared OTUs between sampler × duration pairs. Display as corrplot. Shows whether samplers agree on which species are dominant. Replaces the old similarity ranking section
16. **Venn/Euler**: one key comparison — 4 main active samplers (Kärcher, Coriolis, Sass, Hepa) at 5 hours
17. **iNEXT sensitivity (supplement)**: coverage-based rarefaction/extrapolation for q=0. Confirm sampler ranking matches observed S
18. **Group-specific extras** (short — one figure max): see below
19. **Export results for cross-group figures**: write TSV files to `Data/{Group}/out_intermediate/` for use in `Joined_Figures.Rmd`. See export spec below
20. **Session info**

### Sections REMOVED from Old Template

These are intentionally dropped — do NOT re-add:
- NMDS ordination (replaced by nestedness/turnover + rank-abundance correlation)
- PERMANOVA + betadisper (confounds richness with composition for binary Bray-Curtis)
- specaccum curves (replaced by pooled-reference completeness)
- Old similarity ranking (replaced by rank-abundance correlation)
- Full red-list sub-pipeline for fungi (condensed to one-figure extra)
- Aquatic fungi, plant pathogens sections

### Group-Specific Sections (condensed)

- **Fungi**: Richness by subgroup = Ascomycota vs Basidiomycota (+ Agaricomycetes specifically). Tests whether sampler effects differ by spore type (unicellular Asco vs fruiting-body Basidio). Red list: one summary (% of Swiss macrofungi checklist detected, by sampler) — no full sub-pipeline
- **Plants**: Merge plant + fungi ITS2 data (phyloseq tax_glom). Vegetation survey comparison (InfoFlora) — kept but condensed. Invasive species: simple detection list
- **Insects**: Genus-level bubble plot (ggh4x::facet_nested). Primer specificity: brief check
- **Vertebrates**: Species bubble plot (ggh4x::facet_nested). Homo sapiens exclusion before analysis
- **Bacteria**: No kingdom-level filter (just phylum != NA). Subgroup richness by major phyla (Proteobacteria, Actinobacteria, Firmicutes, etc.)

### Exports for Cross-Group Figures

Each Data_Analysis script must export TSV files to `Data/{Group}/out_intermediate/` so that `Joined_Figures.Rmd` can produce multi-panel figures across all 5 species groups. All exports include a `taxon` column (e.g., "Fungi", "Bacteria") for binding across groups.

| File | Contents | Key columns |
|------|----------|-------------|
| `{Group}_richness.tsv` | Observed S (and Faith's PD when available) per sample | taxon, Index, Sites, Sampler, duration, S, PD |
| `{Group}_contamination.tsv` | Contamination summary per sampler: n OTUs in field controls, n with ratio > 1, n with ratio > 10, total reads in controls | taxon, Sampler, n_OTUs_field, n_ratio_gt1, n_ratio_gt10, total_reads_field, n_OTUs_lab, n_ratio_gt1_lab |
| `{Group}_completeness.tsv` | Fraction of pooled 5h reference community detected | taxon, Sites, Sampler, duration, fraction_detected, S_sampler, S_reference |
| `{Group}_nestedness_turnover.tsv` | Betapart decomposition per sampler pair | taxon, Sampler1, Sampler2, duration, beta_sor, beta_sim (turnover), beta_sne (nestedness) |
| `{Group}_rank_correlation.tsv` | Rank-abundance Spearman rho per sampler pair | taxon, Sampler1, Sampler2, rho, p_value |

Export using `write_tsv(df, here("Data", primer, "out_intermediate", filename))` where `primer` is the group name (e.g., "Fungi").

**Contamination as a result**: The field control contamination levels differ substantially across samplers and are themselves a finding (some samplers are inherently "dirtier" than others). The contamination export must include per-sampler summaries so the manuscript can report contamination patterns across all 5 species groups in a single table or figure.

### Variable Naming Convention

| Variable | Description |
|----------|-------------|
| `{Group}_ASV` | ASV tibble (Index + OTU columns) |
| `{Group}_tax` | Taxonomy tibble (seq + rank columns) |
| `{Group}_seq` | DNAStringSet of representative sequences |
| `{Group}_tree` | phylo object (ape) — phylogenetic tree |
| `{Group}_ASV_clean` | ASV tibble after target-group filtering |
| `{Group}_ASV_clean_long` | Long-format of clean ASV tibble |
| `S_df` | Species richness data frame (observed S + Faith's PD) |
| `Cor_mat` | Richness correlation matrix |
| `beta_df` | Betapart nestedness/turnover decomposition |
| `rank_cor` | Rank-abundance correlation matrix across samplers |
| `cont_ratio` | Control contamination ratio table |
| `Lab_only` / `Field_only` | Subsets of cont_ratio |
| `completeness_df` | Sampling completeness (fraction of pooled reference) |

---

## Harmonization Status

| Script | Status | Key Issues |
|--------|--------|-----------|
| Data_Analysis_Fungi.Rmd | **REWRITE IN PROGRESS** | Template script — being rewritten to match revised analysis plan (Feb 2026). Old version had NMDS/PERMANOVA, full red-list sub-pipeline, aquatic fungi — all to be replaced |
| Data_Analysis_Bacteria.Rmd | **TODO** | Previously harmonized (Feb 2025) but needs rewrite to match new template. Was: Google Sheets metadata, wrong paths |
| Data_Analysis_Insects.Rmd | **TODO** | Google Sheets metadata, needs full rewrite to new template |
| Data_Analysis_Vertebrates.Rmd | **TODO** | Google Sheets metadata, wrong data paths, very incomplete — needs full rewrite |
| Data_Analysis_Plants.Rmd | **TODO** | xlsx metadata (not Meta_data.tsv), needs full rewrite. Has unique: merge ITS data, vegetation survey, invasives |

### Implementation Order

1. **Fungi first** — rewrite as template following revised section order
2. **Then harmonize remaining scripts** one at a time, modelling on Fungi template
3. **Trees**: Add FastTree code blocks to Seq_Analysis scripts (NOT run) after Fungi template is stable

### Common Issues to Fix in All Scripts

1. **Metadata**: Use `read_tsv(here("Data", "Meta_data.tsv"))` everywhere
2. **Data paths**: Use `here("Data", primer, "out_final", ...)` consistently
3. **Contamination**: Full decontam + lab/field control + 2× max subtraction pipeline
4. **Libraries**: Must include: forcats, purrr, corrplot, Biostrings, decontam, eulerr, betapart, ape, picante, iNEXT
5. **New sections to add**: GLM for richness, Faith's PD (placeholder), sampling completeness (pooled reference), betapart nestedness/turnover, rank-abundance correlation, iNEXT sensitivity
6. **Sections to remove**: NMDS, PERMANOVA, betadisper, specaccum, old similarity ranking
7. **Remove**: References to `here("Documents", ...)` (doesn't exist in repo)

---

## Scripts NOT to Modify (Only Cosmetic Cleanup)

- `Seq_Analysis_*.Rmd` — resource-heavy bioinformatics (DADA2, BLAST), DO NOT RUN
- `Prepare_RefDB.Rmd` — downloads large databases, DO NOT RUN
- `Demultiplexing_by_Primer.Rmd` — raw read processing
- `Seq_analysis_helper.R` — shared functions, stable

---

## Git & GitHub

- **Remote**: https://github.com/FabianRoger/23_Air_sampler_test
- **Branch**: main
- **Credentials setup**:
  ```bash
  git config credential.helper store
  echo "https://FabianRoger:<YOUR_TOKEN>@github.com" > ~/.git-credentials
  ```
- **Index workaround**: If you get `index.lock` errors, use `GIT_INDEX_FILE=/tmp/git_alt_index` for all git commands
- **Total tracked**: ~238 files, ~30 MB (large files excluded via .gitignore)

---

## R Session Requirements

Key packages (all available on CRAN/Bioconductor):

```
# Core tidyverse + plotting
dplyr, tidyr, purrr, forcats, readr, readxl, here, ggplot2, ggh4x

# Ecology / diversity
vegan, betapart, iNEXT, picante

# Phylogenetics
ape, Biostrings

# Contamination
decontam, microDecon

# Visualization
corrplot, eulerr

# Other
phyloseq, rgbif, data.table
```

Plants additionally uses: rvest, sf, leaflet (for vegetation survey comparison).
