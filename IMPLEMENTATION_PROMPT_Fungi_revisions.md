# Implementation Prompt: Revise Data_Analysis_Fungi.Rmd

Read `CLAUDE.md` first, then `GENERAL_ANALYSIS_CODING_STYLE.md`, then the current `Scripts/Data_Analysis_Fungi.Rmd` in full. Apply all changes below. Do not reorder sections that aren't mentioned — only modify the specific sections listed. Preserve the overall script structure and style.

---

## Change 1: Add libraries (Setup chunk, ~line 10)

Add `glmmTMB` and `DHARMa` to the library block:

```r
library(glmmTMB)
library(DHARMa)
```

---

## Change 2: Condense decontam (Section "Check controls > Decontam", ~lines 228–286)

**Replace** the current decontam implementation (which includes a threshold sweep, a plot, and extraction of contaminant sequences) with a single condensed chunk:

- Keep the `Blanks` definition chunk exactly as is (lines 232–244).
- Replace everything from the `isContaminant` call through the `contam_seqs` extraction (lines 248–286) with ONE chunk that:
  1. Runs `isContaminant(as.matrix(Fungi_ASV[,-1]), neg = Blanks$is_neg, method = "prevalence")`
  2. Counts contaminants at p < 0.1: `sum(contamdf.prev$p < 0.1, na.rm = TRUE)`
  3. Reports the result in a markdown sentence after the chunk (e.g., "Decontam (prevalence method, p < 0.1) identifies 0 contaminants. This is expected for air eDNA where contaminants are also environmental organisms, making the prevalence signal weak.")
- **Remove**: The threshold sweep (`p_thresholds`, `contam_counts`, the ggplot of threshold vs n_contaminants). Remove the `contam_seqs` extraction chunk.

---

## Change 3: Replace lab control heatmap with taxonomy commentary (~lines 323–374)

**Remove** the heatmap chunk (the chunk that filters to lab control sequences, joins taxonomy, creates `spec_name`, and plots with `geom_tile` + `scale_fill_viridis_c`).

**Replace** with a taxonomy summary chunk:

- Join `Lab_only` (filtered to `name == "mean_val"`) to `Fungi_tax` by `seq`.
- Group by phylum (or class/order where phylum is too broad) and summarise: n OTUs, n with ratio > 1, total reads.
- Display as a simple table (print the tibble).
- Add markdown commentary after the chunk identifying expected contaminant categories:
  - Molds (Aspergillus, Penicillium, Cladosporium) — omnipresent lab/air contaminants
  - Skin-associated fungi (Malassezia) — classical human contaminants
  - Soil fungi (Mortierella) — may indicate sample handling contamination
  - Adapt the specific genera to whatever is actually found in the data.

Keep the summary statistics chunks before this (counting unique sequences, counting ratio > 1) — those are fine.

---

## Change 4: Add taxonomy commentary for field controls (~after line 470)

After the existing field control summary chunks (which show contamination by ratio categories per sampler — keep those), add a new chunk:

- Join field control OTUs (those with ratio > 1) to `Fungi_tax`.
- Group by phylum/class and sampler, summarise n OTUs.
- Display as a summary table.
- Add markdown commentary noting which contaminant types dominate in which samplers and whether they are expected (same logic as lab controls).

Keep all the existing field control visualization chunks (the jitter plot of reads by sampler, the n OTUs by sampler plot, the ratio category table). These are compact and show contamination as a result.

---

## Change 5: Replace GLM with GLMM + diagnostics (Section "Formal statistical model", ~lines 778–803)

**Replace** the current `glm(S ~ Sampler * duration + Sites, family = poisson)` with:

Chunk 1 — Fit GLMM:
```r
richness_glmm <- glmmTMB(S ~ Sampler * duration + (1|Sites),
                          data = S_df_model,
                          family = poisson(link = "log"))
summary(richness_glmm)
```

Chunk 2 — Check overdispersion:
```r
sim_res <- simulateResiduals(richness_glmm, n = 1000)
testDispersion(sim_res)
plot(sim_res)
```

Chunk 3 — Refit with negative binomial if overdispersed (conditional):
```r
# If overdispersion is detected (p < 0.05), refit:
richness_glmm_nb <- glmmTMB(S ~ Sampler * duration + (1|Sites),
                              data = S_df_model,
                              family = nbinom2)
summary(richness_glmm_nb)

# Re-check residuals
sim_res_nb <- simulateResiduals(richness_glmm_nb, n = 1000)
plot(sim_res_nb)
```

Add markdown between chunks explaining the rationale: "Sites are a blocking factor (random sample of locations), not of direct interest. A random intercept is more appropriate than fixed effects and saves degrees of freedom. We check for overdispersion — common in metabarcoding richness data — and refit with negative binomial if needed."

Chunk 4 — Extract predictions: Update the prediction chunk to use whichever model is final (the nb model if overdispersed, otherwise the Poisson). Use `predict()` from glmmTMB — note that `se.fit` works differently; use `predict(model, type = "conditional", se.fit = TRUE)`.

Also update the subgroup GLMs (Ascomycota, Basidiomycota, ~lines 1113–1122) to use the same GLMM approach: `glmmTMB(S ~ Sampler * duration + (1|Sites), family = poisson)` with a brief overdispersion note. These don't need the full diagnostic treatment — just note that the same approach is used.

---

## Change 6: Fix sampling completeness (Section "Sampling completeness", ~lines 926–1049)

### 6a. Simplify the fraction_detected computation (~lines 950–968)

The current code computes `S_shared` which is always identical to `S_sampler` (since the reference is the union of all samplers at the site). Simplify:

```r
completeness_df <-
  Fungi_ASV_clean_long %>%
  left_join(Meta) %>%
  filter(Sampler != "Drone") %>%
  filter(reads > 0) %>%
  group_by(Sites, Sampler, duration) %>%
  summarise(S_sampler = n_distinct(seq), .groups = "drop") %>%
  left_join(ref_community_richness) %>%
  mutate(
    fraction_detected = S_sampler / S_reference,
    taxon = "Fungi"
  )
```

Drop the `S_shared` column entirely.

### 6b. Replace the estimateD block (~lines 986–1009)

The current code pools across sites incorrectly and uses `base = "size"`. Replace with per-site estimation at 95% coverage:

```r
# Pool reads across all samplers per site, then estimate richness at 95% coverage
site_abundance_list <-
  Fungi_ASV_clean_long %>%
  left_join(Meta) %>%
  filter(Sampler != "Drone") %>%
  group_by(Sites, seq) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  filter(reads > 0) %>%
  group_by(Sites) %>%
  group_split()

# Convert each site to a named abundance vector
names(site_abundance_list) <- map_chr(site_abundance_list, ~unique(.x$Sites))

site_vectors <- map(site_abundance_list, ~{
  v <- .x$reads
  names(v) <- .x$seq
  v
})

inext_estimates <-
  map_df(names(site_vectors), function(site) {
    est <- estimateD(site_vectors[[site]],
                     datatype = "abundance",
                     base = "coverage",
                     level = 0.95) %>%
      filter(Order.q == 0)
    tibble(
      Sites = site,
      S_estimated = est$qD,
      SE = est$qD.LCL  # or extract s.e. as needed
    )
  })
```

Note: Check the exact column names returned by `estimateD` with `base = "coverage"` — they may be `qD`, `qD.LCL`, `qD.UCL` rather than `Observed`, `Estimator`, `s.e.`. Adapt column selection accordingly.

### 6c. Update fraction_of_estimated computation (~lines 1014–1023)

Remove the `group_by() %>% summarise(S_sampler = mean(S_sampler))` — there should be no averaging, since each row in `completeness_df` is already one sampler × duration × site observation:

```r
completeness_df <-
  completeness_df %>%
  left_join(inext_estimates, by = "Sites") %>%
  mutate(fraction_of_estimated = S_sampler / S_estimated)
```

### 6d. Update the export (~lines 1041–1049)

Export `completeness_df` directly — it now contains all needed columns: taxon, Sites, Sampler, duration, S_sampler, S_reference, fraction_detected, S_estimated, fraction_of_estimated. No separate `completeness_export` join needed.

```r
write_tsv(completeness_df, here("Data", "Fungi", "out_intermediate", "Fungi_completeness.tsv"))
```

---

## Change 7: Fix betapart row alignment (Section "Nestedness vs turnover", ~lines 1157–1248)

The current code builds `comm_mat` from `Fungi_ASV_clean` and `sample_meta` from a separate `Meta` filter. The loop indexes both with the same `i, j` indices, which assumes identical row order. This is fragile.

**Fix**: Build `sample_meta` directly from `Fungi_ASV_clean$Index` to guarantee alignment:

```r
sample_meta <-
  tibble(Index = Fungi_ASV_clean$Index) %>%
  left_join(Meta) %>%
  filter(Sampler != "Drone") %>%
  mutate(Sampler_dur = paste(Sampler, duration))
```

Then also filter `comm_mat` to only include non-Drone samples (currently `comm_mat` includes all samples in `Fungi_ASV_clean`, which may include Drone). Make sure `comm_mat` and `sample_meta` are built from the same filtered set:

```r
non_drone_idx <-
  Meta %>%
  filter(Index %in% Fungi_ASV_clean$Index) %>%
  filter(Sampler != "Drone") %>%
  pull(Index)

comm_mat <-
  Fungi_ASV_clean %>%
  filter(Index %in% non_drone_idx) %>%
  data.frame() %>%
  `rownames<-`(.$Index) %>%
  select(-Index) %>%
  as.matrix()

sample_meta <-
  tibble(Index = rownames(comm_mat)) %>%
  left_join(Meta) %>%
  mutate(Sampler_dur = paste(Sampler, duration))
```

This guarantees `sample_meta[i,]` corresponds to `comm_mat[i,]`.

### Fix export columns

The current export embeds duration in `Sampler1`/`Sampler2` (e.g., "Kärcher 5 hours"). Instead, store sampler name and duration separately:

In the loop, replace:
```r
Sampler1 = meta1$Sampler_dur,
Sampler2 = meta2$Sampler_dur,
```
with:
```r
Sampler1 = meta1$Sampler,
Sampler2 = meta2$Sampler,
duration1 = meta1$duration,
duration2 = meta2$duration,
```

Update `beta_summary` grouping and export accordingly.

---

## Change 8: Fix rank-abundance correlation row alignment (~lines 1251–1344)

Apply the same row-alignment fix as Change 7. The `rel_abund_mat` must be built from the same filtered index set as `sample_meta`. Reuse the `non_drone_idx` and `sample_meta` from the betapart section (they should already exist in the environment).

```r
rel_abund_mat <-
  Fungi_ASV_clean %>%
  filter(Index %in% non_drone_idx) %>%
  data.frame() %>%
  `rownames<-`(.$Index) %>%
  select(-Index) %>%
  as.matrix()

rel_abund_mat <- sweep(rel_abund_mat, 1, rowSums(rel_abund_mat), "/")
```

Also fix the export columns: store Sampler1, Sampler2, duration1, duration2 separately (same as betapart).

---

## Change 9: Replace iNEXT sensitivity section (~lines 1370–1408)

**Remove** the entire current implementation (which uses `iNEXT()` with pooled-across-sites data and `ggiNEXT(type = 1)`).

**Replace** with `estimateD`-based per-sample coverage-standardized richness:

```r
# Prepare abundance data per sample (sampler × duration × site)
sample_abundance <-
  Fungi_ASV_clean_long %>%
  left_join(Meta) %>%
  filter(Sampler != "Drone") %>%
  filter(reads > 0) %>%
  group_by(Index) %>%
  group_split()

# Estimate richness at 95% coverage for each sample
inext_sensitivity <-
  map_df(sample_abundance, function(samp) {
    idx <- unique(samp$Index)
    abund_vec <- samp$reads

    est <- tryCatch(
      estimateD(abund_vec,
                datatype = "abundance",
                base = "coverage",
                level = 0.95) %>%
        filter(Order.q == 0),
      error = function(e) NULL
    )

    if(is.null(est)) return(NULL)

    tibble(
      Index = idx,
      S_coverage95 = est$qD
    )
  }) %>%
  left_join(Meta)
```

Then plot:
```r
inext_sensitivity %>%
  filter(!is.na(S_coverage95)) %>%
  left_join(S_df %>% select(Index, S)) %>%
  pivot_longer(c(S, S_coverage95), names_to = "metric", values_to = "richness") %>%
  mutate(metric = ifelse(metric == "S", "Observed", "Coverage-standardized (95%)")) %>%
  ggplot(aes(y = richness, x = duration, colour = Sites)) +
  facet_grid(metric ~ Sampler, scales = "free_y") +
  geom_line(aes(group = paste(Sites, metric)), linewidth = 0.2, colour = "grey") +
  geom_point() +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "Duration", y = "Species richness")
```

Add markdown: "Coverage-standardized richness at 95% coverage confirms that sampler rankings are robust to differences in sequencing depth. [Describe whether rankings match or not based on the figure.]"

Remove the `ggiNEXT` chunk and the `$AsyEst` extraction chunk.

---

## Change 10: Fix final export section (~lines 1476–1505)

The `if(!file.exists(...))` guards are unnecessary — just write unconditionally. The exports already happened inline earlier in the script. But if keeping a final "ensure all exports exist" section, remove the conditional logic (the files should be written during the analysis, not checked at the end). Alternatively, just remove this entire section since all exports are written inline.

---

## General reminders

- Follow `GENERAL_ANALYSIS_CODING_STYLE.md`: one chunk per task, one output per chunk, bullet-point markdown before each chunk.
- Use `here()` for all paths.
- Preserve all existing markdown commentary that isn't explicitly replaced.
- Do NOT add or remove sections beyond what's specified above.
- Do NOT change the Venn diagram, red list, richness correlation, or richness vs volume sections.
- Test that the `estimateD` column names match what iNEXT 3.x returns — they may be `qD`, `qD.LCL`, `qD.UCL` for coverage-based estimation rather than `Observed`, `Estimator`, `s.e.`.
