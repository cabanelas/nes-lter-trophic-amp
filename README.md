# NES-LTER Trophic Amplification
## Pelagic Community Structure: Interannual Variability and Long-Term Change in Pelagic Community Structure Across a Latitudinal Gradient
### LTER Pelagic Synthesis Working Group

**Author:** Alexandra C. Cabanelas  
**Project:** LTER Synthesis Working Group: Pelagic Community Structure  
**LTER Site:** Northeast U.S. Shelf Long-Term Ecological Research (NES-LTER)  
**Created:** August 2025 | **Last updated:** March 2026

---

## Working Group Context

This repository is part of a cross-site **LTER Pelagic Synthesis Working Group**: *Interannual variability and long-term change in pelagic community structure across a latitudinal gradient*. The working group brings together four marine LTER sites spanning a wide latitudinal gradient to compare how pelagic communities respond to environmental variability at annual and longer time scales:

| Site | Description |
|------|-------------|
| **NES** - Northeast U.S. Shelf | Rapidly warming temperate shelf in the northwest Atlantic Ocean. Supports productive fisheries and dynamic planktonic food webs shaped by strong seasonal forcing. |
| **CCE** - California Current Ecosystem | Eastern boundary upwelling system off the U.S. West Coast. High spatial variability in ocean conditions supports diverse assemblages across trophic levels. |
| **NGA** - Northern Gulf of Alaska | Subarctic coastal system with strong freshwater influence from glacial runoff. Characterized by a large spring bloom and diapausing zooplankton, transitioning to a smaller-bodied summer community. |
| **PAL** - Palmer, Antarctica | Rapidly warming polar system west of the Antarctic Peninsula. Sea ice seasonality is a dominant driver of food web structure across trophic levels. |

The working group is pursuing three interconnected projects:

1. **[Normalized Biomass Size Spectra (NBSS)](https://github.com/cabanelas/nes-lter-zp-sizespectra.git)**
2. **Trophic Amplification** = *this repository*
3. **[Double Integration Hypothesis](https://github.com/cabanelas/doubleintegration.git)**

This repository contains the **NES-LTER pipeline** for Project 2. Equivalent pipelines for the other projects are maintained in separate repositories.

---

## Overview

This repository tests the **trophic amplification hypothesis**: that climate-driven biomass variability is amplified at higher trophic levels and at lower latitudes (Kwiatkowski et al. 2019; Lotze et al. 2019). We quantify interannual variability in biomass at four trophic levels — phytoplankton (not on this repo), zooplankton, ichthyoplankton, and forage fish — across three regions of the Northeast U.S. Shelf (Gulf of Maine, Georges Bank, Mid-Atlantic Bight) and compare the standard deviation of log-transformed annual biomass means across trophic levels and time periods.

All trophic levels follow a common processing pipeline:
1. **Log-transform**
2. **Annual mean** across stations per region and season
3. **5-year running mean**
4. **Standard deviation** of the annual log-mean time series

---

## Data Sources

| Trophic Level | Metric | Dataset | Source |
|---|---|---|---|
| Zooplankton | Displacement volume (mL mm<sup>-2</sup>) | NOAA EcoMon Plankton Survey | [NCEI: 0187513](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513) |
| Ichthyoplankton | Larval abundance (ind. 10 m<sup>-2</sup>) | NOAA EcoMon Plankton Survey | [NCEI: 0187513](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513) |
| Forage fish | Biomass CPUE (kg tow<sup>-1</sup>) | NOAA NEFSC Bottom Trawl Survey | [InPort: 22560 (Fall)](https://inport.nmfs.noaa.gov/inport/item/22560), [22561 (Spring)](https://inport.nmfs.noaa.gov/inport/item/22561) |

**Ichthyoplankton and Forage fish species included:** Atlantic herring (*Clupea harengus*), Atlantic mackerel (*Scomber scombrus*), Atlantic butterfish (*Peprilus spp.*), sand lance (2 *Ammodytes* spp.)

**Auxiliary files:**

| File | Source |
|---|---|
| `NEFSC_conversion_factors.csv` | [OceanAdapt GitHub](https://github.com/pinskylab/OceanAdapt/blob/master/data_raw/NEFSC_conversion_factors.csv) |
| `Miller2010_Bigelow_calibration_factors.csv` | [Miller et al. 2010](https://repository.library.noaa.gov/view/noaa/3726), NEFSC Ref Doc 10-05, Tables 56–58 |
| `EcomonStrata_v4.shp`, `EcomonStrata_v4b.shp` | [NOAA EcoMon strata shapefiles](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513) |

---

## Repository Structure

```
nes-lter-trophic-amp/
├── scripts/
│   ├── trophic_amp_zooplankton.R       # log transform → running mean → SD → outputs for zooplankton
│   ├── trophic_amp_ichthyoplankton.R   # log transform → running mean → SD → outputs for ichthyoplankton
│   ├── trophic_amp_foragefish.R        # vessel/gear corrections; log transform → running mean → SD → outputs for forage fishes 
├── data/
│   ├── raw/                            # raw, publicly available data (not tracked by git)
│   ├── output/                         # processed CSVs               (not tracked by git)
├── .gitignore
└── README.md
```

---

## Regions

Regions are assigned spatially using EcoMon strata shapefiles. Southern New England (SNE) is merged into the Mid-Atlantic Bight.

| Region | Code |
|---|---|
| Gulf of Maine | GOM |
| Georges Bank | GB |
| Mid-Atlantic Bight (+ SNE) | MAB |

---

## Forage Fish Vessel & Gear Corrections

The NEFSC bottom trawl survey spans multiple vessels and gear configurations. Catches are standardized to the **Albatross IV** baseline using the following corrections applied in order:

| Step | Correction | Applied when | Direction |
|---|---|---|---|
| 1 | Door Conversion Factor (DCF) | year < 1985, all vessels | multiply |
| 2 | Gear Conversion Factor (GCF) | spring surveys, 1973–1981 | divide |
| 3 | Vessel Conversion Factor (VCF) | vessel == Delaware II ("DE") | multiply |
| 4 | Bigelow/Pisces calibration (ρ_W) | vessel == HB or PC (post-2008) | divide |

Sources: `NEFSC_conversion_factors.csv` (steps 1–3); Miller et al. 2010 Tables 56–58 (step 4).  
Season-specific ρ_W values used where available (herring, butterfish); combined estimate used otherwise (mackerel, sand lance). No ρ_W correction applied to *Ammodytes americanus* — no published value.

---

## How to Run

1. Clone this repository
2. Place all raw data files in `raw/` (see Data Sources above)
3. Open R and set your working directory to the project root, or open the `.Rproj` file
4. Run scripts (can be run in any order):

```r
source("R/trophic_amp_zooplankton.R")
source("R/trophic_amp_ichthyoplankton.R")
source("R/trophic_amp_foragefish.R")
```

---

## Methods Summary

### Log transformation
Station-level data is log-transformed as:

```
log10(x + min_nonzero / 2)
```

where `min_nonzero` is the smallest positive observed value within each region × season × trophic level group.

### Annual means
Log-transformed values are averaged across all stations within a region × season × year group.

### 5-year running mean
A right-aligned, index-aware running mean (k = 5) is applied using `runner::mean_run()`. The window uses all available years up to k, so early years in a time series receive shorter windows rather than being dropped.

### Standard deviation
SD is computed across the full annual log-mean time series within each region × season group. This is the primary metric for the trophic amplification test: if variability increases at higher trophic levels, SD should be larger for forage fish than for phytoplankton.

---

## Dependencies

R version used: 4.5.2 (2025-10-31 ucrt)

This project uses `renv` for package version management. To restore the exact environment:

```r
install.packages("renv")
renv::restore()
```

> Note: `renv::restore()` will automatically install all required packages (tidyverse, here, runner, janitor, data.table, sf) at the correct versions.  

---

## Citation

If you use this code, please cite:
<!-- Add citation here -->

---

## License

<!-- Add license here -->
