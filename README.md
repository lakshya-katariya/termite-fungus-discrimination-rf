# Fungal VOC profiles discriminate between phyla: a Random Forest approach

Can a fungus's scent reveal what kind of organism it is?

This repo applies Random Forest classification to gas chromatography–mass
spectrometry (GC-MS) volatile organic compound (VOC) profiles from 43 fungal
species to ask whether VOC emissions systematically differ between
**Ascomycota** and **Basidiomycota** — and which specific compounds drive that
discrimination.

📄 **[View the full analysis](fungal_voc_random_forest.md)**

---

## Scientific background

This analysis replicates the analytical approach of:

> Katariya L, Ramesh PB, Gopalappa T et al. (2017). *Fungus-Farming Termites
> Selectively Bury Weedy Fungi that Smell Different from Crop Fungi.*
> Journal of Chemical Ecology 43:986–995.
> [https://doi.org/10.1007/s10886-017-0902-4](https://doi.org/10.1007/s10886-017-0902-4)

In that study, we showed that fungus-farming termites discriminate their
crop fungus (*Termitomyces*) from a parasitic weed (*Pseudoxylaria*) using
VOC profiles alone — without visual contact. Random Forest classification on
GC-MS chemical profiles identified the candidate volatile compounds driving
that discrimination.

Here I apply the same framework to a broader public dataset to ask whether
phylogenetic identity (phylum) is similarly encoded in fungal VOC chemistry.

---

## Data

- **Source:** Guo et al. (2021). *Volatile organic compound patterns predict
  fungal trophic mode and lifestyle.* Communications Biology 4:673.
  [https://doi.org/10.1038/s42003-021-02198-8](https://doi.org/10.1038/s42003-021-02198-8)
- **Data:** Supplementary Table S4 (GC-MS emission intensities) and Table S1
  (species metadata including phylum), available at
  [https://osf.io/bva2q](https://osf.io/bva2q) under CC-BY 4.0.
- **Species:** 39 species retained after excluding Zygomycota
  (26 Ascomycota, 13 Basidiomycota)
- **Features:** Named volatile compounds only (m/z-only entries excluded)

---

## Methods

| Step | Approach | Rationale |
|------|----------|-----------|
| Variable selection | `varSelRF` repeated 100 times | Identifies minimal stable discriminating compound set across runs |
| Classification | `randomForest` on selected variables | Obtains OOB error, feature importance, proximity matrix |
| Importance metric | Mean Decrease in Accuracy | Measures genuine compound contribution, not spurious correlation |
| Visualisation | MDS on RF proximity matrix | Non-parametric species similarity in VOC space |

---

## Repo structure

```
├── fungal_voc_random_forest.Rmd   # source — edit and knit here
├── fungal_voc_random_forest.md    # rendered output — view on GitHub
└── figures/                       # auto-generated plots
    ├── fig-varselrf-1.png         # varSelRF stability plot
    ├── fig-importance-1.png       # feature importance plot
    └── fig-mds-1.png              # MDS proximity plot
```

---

## How to reproduce

1. Download `TableS4_GCMSdata.xlsx` from [osf.io/bva2q](https://osf.io/bva2q)
   and place it in the same directory as the `.Rmd` file
2. Install required packages:
```r
install.packages(c("tidyverse", "readxl", "randomForest",
                   "varSelRF", "janitor", "parallel",
                   "doParallel", "foreach"))
```
3. Open `fungal_voc_random_forest.Rmd` in RStudio and click **Knit**

> **Note:** The `varSelRF` chunk is parallelised across available CPU cores
> and cached after the first run (`cache=TRUE`). Expect 5–15 minutes on first
> knit depending on your machine.

---

## Author

**Lakshya Katariya**
[LinkedIn](https://be.linkedin.com/in/lakshyakatariya) ·
[Google Scholar](https://scholar.google.com/citations?user=iEpb8ycAAAAJ&hl=en) ·
[ORCID](https://orcid.org/0000-0001-5667-8281)
