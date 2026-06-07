# Can you smell the difference? Identifying fungal discriminating volatiles using Random Forest classification

Fungi communicate chemically through volatile organic compounds (VOCs). This repo asks: can VOC profiles alone reveal what kind of fungus is emitting them — and which specific compounds carry the discriminatory signal?

📄 **[View the full analysis](fungal_voc_random_forest.md)**

---

## Structure

The analysis has two parts that build on each other:

**Part 1 — Replication** replicates the Random Forest analysis from
[Katariya et al. (2017)](https://doi.org/10.1007/s10886-017-0902-4) using
the published supplementary data from that paper. In the original study, we
showed that fungus-farming termites discriminate their crop fungus
(*Termitomyces*) from a parasitic weed (*Pseudoxylaria*) using VOC profiles
alone — without visual contact. Here we reproduce that analysis to confirm
that the key discriminating compounds are stably selected across 100
`varSelRF` runs.

**Part 2 — Extension** applies the same analytical framework to a broader
public dataset from [Guo et al. (2021)](https://doi.org/10.1038/s42003-021-02198-8),
asking whether VOC profiles discriminate **Ascomycota** from
**Basidiomycota** across 39 fungal species. This tests whether the
chemical discrimination signal observed at the species level in Part 1
generalises to the phylum level.

---

## Data

| Part | Source | Data |
|------|--------|------|
| Part 1 | Katariya et al. (2017). *J. Chem. Ecol.* 43:986–995 | Supplementary Table S1 (41 VOCs, 16 replicates) — embedded in script |
| Part 2 | Guo et al. (2021). *Commun. Biol.* 4:673 | Supplementary Table S4 ([osf.io/bva2q](https://osf.io/bva2q), CC-BY 4.0) |

> Part 1 data is entered directly in the script — no download needed.
> Part 2 requires downloading `TableS4_GCMSdata.xlsx` from the link above.

---

## Methods

| Step | Tool | Purpose |
|------|------|---------|
| Variable selection | `varSelRF` × 100 runs | Identifies the minimal stable set of discriminating compounds |
| Classification | `randomForest` on selected variables | Obtains OOB error, permutation importance, proximity matrix |
| Importance metric | Mean Decrease in Accuracy | Measures how much accuracy drops when each compound is randomly shuffled |
| Ordination | MDS on RF proximity matrix | Non-parametric visualisation of similarity in VOC space |
| Parallelisation | `doParallel` + `foreach` | Distributes 100 `varSelRF` runs across CPU cores |

---

## Repo structure

```
├── fungal_voc_random_forest.Rmd   # source — edit and knit here
├── fungal_voc_random_forest.md    # rendered output — view on GitHub
└── figures/                       # auto-generated plots
    ├── fig-varselrf-p1-1.png      # Part 1: varSelRF stability
    ├── fig-importance-p1-1.png    # Part 1: feature importance
    ├── fig-mds-p1-1.png           # Part 1: MDS proximity plot
    ├── fig-varselrf-p2-1.png      # Part 2: varSelRF stability
    ├── fig-importance-p2-1.png    # Part 2: feature importance
    └── fig-mds-p2-1.png           # Part 2: MDS proximity plot
```

---

## How to reproduce

1. Download `TableS4_GCMSdata.xlsx` from [osf.io/bva2q](https://osf.io/bva2q)
   and place it in the same directory as the `.Rmd` file.
   Part 1 data is embedded directly in the script.

2. Install required packages:
```r
install.packages(c("tidyverse", "readxl", "randomForest",
                   "varSelRF", "janitor", "parallel",
                   "doParallel", "foreach"))
```

3. Open `fungal_voc_random_forest.Rmd` in RStudio and click **Knit**

> **Note:** Both `varSelRF` chunks are parallelised across available CPU
> cores and cached after the first run (`cache=TRUE`). Expect 10–20 minutes
> on first knit depending on your machine.

---

## Author

**Lakshya Katariya**
[LinkedIn](https://be.linkedin.com/in/lakshyakatariya) ·
[Google Scholar](https://scholar.google.com/citations?user=iEpb8ycAAAAJ&hl=en) ·
[ORCID](https://orcid.org/0000-0001-5667-8281)
