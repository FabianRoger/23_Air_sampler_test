# General Analysis Coding Style Guide

## Purpose
- Write analysis scripts that are easy for scientists and reviewers to read, run, and audit.
- Prioritize scientific clarity, reproducibility, and transparent analytical decisions.
- Treat scripts as manuscript companions, not production software.

## Core Principles
- Keep workflows linear and analysis-focused.
- Prefer readability over cleverness.
- Prefer simplicity over abstraction.
- Explain analytical choices and interpretation, not coding tricks.
- Keep intermediate outputs visible and saved when useful.

## Recommended Script Format
- Preferred format: `RMarkdown` (`.Rmd`) for analysis narratives.
- Use short sections that follow analysis logic from input to result.
- Place clear markdown before each code chunk.
- Write markdown in bullet-point style.
- Keep each code chunk focused on one task.
- Keep each chunk to one displayed output (one table or one figure).

## Data Workflow
- Load all required data near the top of the script.
- Read files directly from repository paths.
- Use stable, human-readable object names.
- Do not over-engineer defensive checks for known schema issues.
- Fix broken metadata/data structure upstream instead of adding many in-script workarounds.

## Coding Conventions
- Prefer tidyverse workflows.
- Prefer pipes (`%>%`) for sequential analysis steps.
- Use `ggplot2` for graphics.
- Keep transformations explicit and easy to follow.
- Avoid unnecessary helper functions for one-time analyses.
- Avoid dense one-liners if they reduce interpretability.

## Explanation Style
- Explain what is being tested or summarized.
- Explain why that step is needed for interpretation.
- Report key outcomes in markdown bullets with inline values when helpful.
- Prefer markdown reporting over `cat()`-based prose chunks.

## Quality Control and Reporting
- Include early QC sections (for example read depth, controls, detection patterns).
- Output manuscript-reportable summary statistics (means, medians, counts, thresholds).
- Explicitly report flagged samples/outliers and criteria used.
- Save QC summary tables to file for transparent reporting.

## Contamination Assessment
- Use explicit contamination checks relevant to the project (for example blanks, known contaminants, host/human sequences).
- Report what contaminants are removed and why.
- Report contaminant taxonomy when available.
- Report where contaminants occur (number of samples, controls, and sample IDs).
- Report before/after filtering effects, including residual signal in controls.

## Output Management
- Save intermediate tabular results as TSV by default.
- Save final figures with clear, descriptive names.
- Use consistent output folders (for example `Data/derived/`, `Figures/...`).
- Ensure outputs are directly reusable in manuscript writing.

## Editing Existing Scripts
- Preserve analytical intent while improving structure.
- Reorder code to match analysis flow when needed.
- Remove redundancy and dead code.
- Rename ambiguous variables to clear domain-specific names.
- Keep style consistent throughout the script.

## Reproducibility Practices
- Keep all required analysis inputs inside the repository when possible.
- End scripts with session/version information.
- Avoid hidden manual steps.
- Make reruns deterministic where practical.

## Checklist Before Finalizing
- Data load is at top and complete.
- Section order follows scientific logic.
- Every chunk has one clear purpose.
- Every chunk has at most one displayed output.
- Markdown explains analytical purpose and interpretation.
- QC and contamination results are explicitly reported.
- Intermediate and final outputs are written with clear names.
- Script runs from repository inputs without ad hoc editing.
