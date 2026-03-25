# Feature Engineering and Statistical Analysis

This document describes the feature engineering pipeline and statistical methods used for sleep architectural instability analysis in the AD continuum. For the full analytical workflow, see the figure in the main repository.

---

## Table of Contents

1. [Subject Selection and Baseline Definition](#1-subject-selection-and-baseline-definition)
2. [Hypnogram Parsing and Cycle Detection](#2-hypnogram-parsing-and-cycle-detection)
3. [Insomnia and Sleep-Onset Features](#3-insomnia-and-sleep-onset-features)
4. [Temporal Deviation Metrics (Stage Duration)](#4-temporal-deviation-metrics-stage-duration)
5. [Transition Metrics](#5-transition-metrics)
6. [Stage Composition Metrics](#6-stage-composition-metrics)
7. [Stage Stability and Fragmentation via Run-Length Encoding](#7-stage-stability-and-fragmentation-via-run-length-encoding)
8. [Weighted Sleep-Quality Indices (Cycle-Level)](#8-weighted-sleep-quality-indices-cycle-level)
9. [Autonomic Function (RMSSD) Features](#9-autonomic-function-rmssd-features)
10. [Moving-Window Analysis](#10-moving-window-analysis)

---

## 1. Subject Selection and Baseline Definition

To establish reference distributions for normalisation, baseline-qualified nights were identified among participants labelled as Healthy control that met the following predefined criteria:

- **3–7 detected sleep cycles** per night (`n_distinct[cycle_id]`)
- **Awake epochs**, if present, restricted to the very first or very last 5-minute epoch of the night (or absent)

Stage-specific means (μ) and standard deviations (σ) were computed from these baseline-qualified Healthy control nights and used as reference values for z-score normalisation and deviation-based feature construction.

> **Note:** To prevent reference leakage, all baseline-qualified participants were excluded from subsequent feature construction and inferential analyses.

---

## 2. Hypnogram Parsing and Cycle Detection

Wearable-derived hypnogram data were recorded at a **5-minute epoch resolution**. To summarise sleep-stage dynamics, hypnograms were decomposed using **run-length encoding (RLE)** within each `email × date × cycle_id` segment. Within each cycle, the following were derived from the RLE representation:

- Stage durations (minutes)
- Transition events
- Stage composition (proportions)

---

## 3. Insomnia and Sleep-Onset Features

Sleep-onset-related features were constructed using two components:

| Component | Definition |
|---|---|
| **Awake initial** | Duration (minutes) of the initial consecutive Awake run at the start of cycle 1 (derived via RLE; 0 if the first run was not Awake) |
| **SOL** | Sleep onset latency reported by the wearable device, converted to minutes (SOL / 60) |

An insomnia score was calculated as the sum of SOL (minutes) and Awake initial (minutes). To capture shared variance between these two components:

1. SOL and Awake initial were z-standardised
2. PCA was performed
3. The first principal component (**Insomnia PC1**) was retained as an integrated insomnia-related feature

---

## 4. Temporal Deviation Metrics (Stage Duration)

Stage-specific durations (minutes) within each sleep cycle were standardised using z-scores derived from the healthy sleep baseline.

### Z-score computation

$$Z_{\text{stage}} = \frac{\text{Duration}_{\text{obs, stage}} - \mu_{\text{Healthy, stage}}}{\sigma_{\text{Healthy, stage}}^{*}}$$

where:

$$\sigma_{\text{Healthy, stage}}^{*} = \begin{cases} \sigma_{\text{Healthy, stage}} & \text{if } \sigma_{\text{Healthy, stage}} > 0 \\ \varepsilon & \text{if } \sigma_{\text{Healthy, stage}} = 0 \end{cases}$$

and **ε = 0.1** was used to ensure numerical stability.

### Cycle-level temporal deviation

For each sleep cycle, the absolute z-scores across all sleep stages were summed:

$$D_{\text{cycle}} = \sum_{\text{stage}=1}^{4} |Z_{\text{stage}}|$$

At the night level:

- **Mean and median** of D_cycle → overall temporal deviation
- **Standard deviation** of D_cycle → within-night variability across cycles

---

## 5. Transition Metrics

Sleep-stage transitions were defined as changes in sleep stage between **adjacent 5-minute epochs** within a cycle.

### Feature construction

For each cycle:
1. All `previous stage → current stage` transition types were extracted
2. Frequency of each transition was normalised by the total number of transitions within that cycle to obtain **transition-type proportions**

Cycle-level transition proportions were aggregated to the night level by computing the **mean proportion** for each transition type.

### Feature selection for transitions

Transition-specific linear mixed-effects models were fitted comparing Healthy control vs Dementia (random intercept for participant). Transition types with **p < 0.01** were retained and converted to wide-format nightly features.

---

## 6. Stage Composition Metrics

Within each sleep cycle, the **proportion of time** spent in each sleep stage was calculated as:

$$\text{Proportion}_{\text{stage}} = \frac{\text{Number of 5-min epochs in stage}}{\text{Total epochs in cycle}}$$

Stage-specific proportions were standardised using healthy baseline reference distributions.

### Compositional deviation score

For each cycle, a compositional deviation score was defined as the **sum of absolute stage-wise proportion z-scores**. This was aggregated to the night level using:

- Mean and median summaries of cycle-level compositional deviation
- Nightly mean stage-wise proportion z-scores (averaged across cycles)

> **Exclusion criterion:** Nights with one or fewer detected sleep cycles were excluded due to insufficient structural information.

---

## 7. Stage Stability and Fragmentation via Run-Length Encoding

RLE was applied to characterise sleep-stage stability within individual sleep cycles by quantifying how long a given sleep stage was maintained without interruption. Consecutive epochs assigned to the same sleep stage were grouped into **runs**, each representing a continuous period in that stage.

### Derived cycle-level metrics

| Metric | Definition |
|---|---|
| **Mean run length** | Average number of consecutive 5-minute epochs within the same sleep stage |
| **Maximum run length** | Longest continuous run observed within a cycle |
| **Dominant stable stage** | Sleep stage corresponding to the maximum run length |

These metrics quantify the extent to which sleep stages are stably maintained versus fragmented within each sleep cycle.

### Deviation from baseline and nightly aggregation

Stage proportion metrics at the cycle level were compared with healthy baseline distributions to compute standardised deviation scores. Cycle-level compositional deviation was aggregated to the night level using mean and median summaries.

> **Exclusion criterion:** Nights with one or fewer detected sleep cycles were excluded from further analyses.

---

## 8. Weighted Sleep-Quality Indices (Cycle-Level)

Several weighted, cycle-level composite indices were constructed to explore whether sleep architecture unfolded with appropriate timing and minimal disruption. These indices were designed based on physiological considerations, with auxiliary weights assigned to:

- **Deep sleep** in early cycles
- **Light sleep** in intermediate cycles
- **REM sleep** in late cycles
- Penalty for excessive Awake dominance

### Component indices

| Index | Definition |
|---|---|
| **Stability index** | Average duration of uninterrupted sleep-stage runs within a cycle (higher = more stable) |
| **Fragmentation index** | Number of within-cycle transitions normalised by number of runs: `frag_idx = transitions / runs` (where `runs = transitions + 1`) |
| **Deep-sleep dominance index** | Longest continuous run weighted by its dominant sleep stage |

### Composite score

For each index, healthy baseline distributions were used to derive z-scores. These components were combined into an exploratory composite score:

$$\text{run\_frag\_score} = -Z_{\text{frag}} + Z_{\text{stable}} + Z_{\text{deep}}$$

Cycle-level composite indices were summarised at the night level using mean, SD, and CV.

> **Note:** These weighted composite indices were used exclusively during the feature engineering and exploratory phase and were **not included in the final inferential analyses**. Nights with one or fewer detected sleep cycles were excluded.

---

## 9. Autonomic Function (RMSSD) Features

Stage-specific RMSSD reference means and standard deviations were computed from **baseline-qualified Healthy control nights**. For non-baseline participants, stage-specific mean RMSSD values were transformed into **z-scores** relative to the baseline reference, capturing relative parasympathetic activity during sleep.

### RMSSD-based recovery scoring

To enhance physiological interpretability of RMSSD deviations, a discrete weighting scheme was defined:

| Score | Condition | Interpretation |
|---|---|---|
| **−2** (impaired recovery) | Mean RMSSD < baseline mean − 1 SD | Insufficient parasympathetic activation; persistent physiological stress |
| **+1** (normal recovery) | Mean RMSSD within baseline range | Typical nocturnal autonomic recovery |
| **+2** (high resilience) | Mean RMSSD > baseline mean + 1 SD | Enhanced parasympathetic activity; autonomic resilience |

An **integrated autonomic regulation metric** was generated by multiplying the RMSSD z-score by the corresponding recovery score, thereby amplifying clinically meaningful deviations in autonomic control.

Cycle-level autonomic metrics were summarised at the night level using mean, SD, and CV. In addition to these derived measures, conventional sleep and autonomic variables provided directly by the wearable device were included in the feature set.

---

## 10. Moving-Window Analysis

To mitigate day-to-day variability and measurement noise in longitudinal sleep data, a **right-aligned 7-day rolling window** was applied within each participant after ordering nights by date.

### Aggregation

For each nightly feature, the following rolling statistics were computed:

- Median
- Mean
- SD
- CV (SD / mean; set to `NA` when |mean| ≤ 1 × 10⁻⁸)

### Retention criterion

A windowed observation was retained only when **at least 60% of nights were non-missing** (≥ 5 of 7) for every feature within that window.

These rolling-window summaries constituted the **final feature set (124 features)**.

---

## Summary: Feature Domains

| Domain | Source | # Features (before aggregation) |
|---|---|---|
| Stage duration (z-score) | Hypnogram + HC baseline | Per stage × cycle |
| Stage composition | Hypnogram | Per stage × cycle |
| Stability / fragmentation | RLE | 3 metrics × cycle |
| Transitions | Hypnogram | Selected transition types |
| RMSSD autonomic | RMSSD + HC baseline | Per stage × cycle |
| Composite index | Derived | 1 per cycle |
| Device features | Oura Ring output | 8 variables |

**Total: 33 source features → 124 derived features** (after moving-window aggregation with mean, median, SD, CV)
