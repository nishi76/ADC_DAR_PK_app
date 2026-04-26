# ADC In Vivo DAR Monitoring & PK Tracker

![R](https://img.shields.io/badge/R-4.2%2B-blue)
![Shiny](https://img.shields.io/badge/Shiny-App-green)
![License](https://img.shields.io/badge/Use-Research%20Only-orange)

An interactive **R Shiny** dashboard for integrating **Total ADC**, **Total Antibody**, and **Free Payload** concentration-time data from in vivo ADC studies.

The app brings these analytes into a single workflow to calculate **average DAR over time**, assess **payload release kinetics**, run **non-compartmental analysis (NCA)**, and generate **report-ready outputs** for study review and submission support.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Who It Is For](#who-it-is-for)
- [Installation](#installation)
- [Usage](#usage)
- [Input Format](#input-format)
- [Analysis Workflow](#analysis-workflow)
- [Output](#output)
- [Limitations](#limitations)
- [Project Structure](#project-structure)
- [Citation](#citation)
- [Research Use Only](#research-use-only)

## Overview

ADC biotransformation in vivo is often monitored through separate assays for Total ADC, Total Antibody, and Free Payload. In many workflows, these results are combined manually in Excel to understand deconjugation, linker cleavage, payload loss, and DAR stability over time.

This app automates that process by combining the three analytes in one dashboard and turning them into a reproducible analysis workflow.

## Features

- Import multi-analyte PK data from a CSV file.
- Calculate average DAR over time from Total ADC and Total Antibody.
- Estimate DAR stability trends and decay over time.
- Quantify Free Payload kinetics.
- Run NCA for all three analytes using PKNCA.
- Generate publication-quality plots, tables, PDF reports, and Excel exports.
- Provide an interactive dashboard for rapid review and comparison.

## Who It Is For

This tool is designed for teams working with ADC PK and biotransformation data, including:

- Bioanalytical scientists developing or interpreting ADC assays.
- DMPK and pharmacokinetics scientists analyzing ADC disposition.
- Translational and pharmacology teams linking exposure to payload release.
- Study directors preparing integrated study summaries.
- CROs and internal teams that currently reconcile results manually in Excel.

## Installation

### Requirements

- R 4.2 or later.
- The following R packages:
  - `shiny`
  - `shinydashboard`
  - `shinyjs`
  - `shinycssloaders`
  - `PKNCA` (>= 0.11)
  - `dplyr`
  - `tidyr`
  - `ggplot2`
  - `plotly`
  - `DT`
  - `openxlsx`
  - `rmarkdown`
  - `knitr`
  - `scales`

### PDF export requirements

To generate PDF reports, you also need a LaTeX installation:

- Linux: `sudo apt-get install texlive-latex-base texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended lmodern`
- macOS: install MacTeX
- Windows: install MiKTeX

### Run the app

From R or RStudio:

```r
shiny::runApp("path/to/adc_pk_app/")
```

Or open `app.R` in RStudio and click **Run App**.

## Usage

1. Launch the app.
2. Upload a CSV file with time-matched ADC data.
3. Confirm the concentration unit used in the study.
4. Review plots, DAR trajectory, and NCA outputs.
5. Export the PDF report or Excel workbook.

## Input Format

Upload a CSV file with the following columns:

| Column | Required | Description |
|---|---:|---|
| `Time_h` | Yes | Nominal sampling time in hours |
| `Total_ADC_mean` | Yes | Mean Total ADC concentration |
| `Total_Ab_mean` | Yes | Mean Total Antibody concentration |
| `FreePayload_mean` | Yes | Mean Free Payload concentration |
| `Dose_mg_kg` | Yes | Administered dose in mg/kg |
| `Total_ADC_SD` | Optional | SD for Total ADC, used for error bars |
| `Total_Ab_SD` | Optional | SD for Total Antibody |
| `FreePayload_SD` | Optional | SD for Free Payload |
| `Species` | Optional | Species label shown in reports |
| `Study_ID` | Optional | Study identifier |

All concentration columns must use the same unit, selected in the app sidebar. Supported units include `ug/mL`, `ng/mL`, `nM`, and `pM`.

A blank template can be downloaded from the app’s Upload tab.

## Analysis Workflow

### DAR calculation

The app estimates average DAR over time using:

```text
avg_DAR(t) = [Total_ADC(t) / Total_Ab(t)] x Nominal_DAR
```

This assumes both analytes are measured in compatible mass-based units and come from the same antibody backbone. The app also fits an exponential decay model to summarize DAR stability over time.

### NCA

NCA is performed using PKNCA. The default interval is from the first observed timepoint to the last observed timepoint, and the app reports parameters such as AUClast, Cmax, Tmax, half-life, lambda-z, and R-squared.

For IV bolus datasets, `tmax = 0` is expected when the first measured timepoint is also the peak observed value relative to the interval start. AUCinf, CL, and Vz are not computed by default unless a T=0 row is included in the input data.

## Output

### Interactive dashboard tabs

The app includes seven main tabs:

1. Upload & Configure.
2. PK Curves.
3. DAR Trajectory.
4. Payload Kinetics.
5. PK/PD Overlay.
6. NCA Parameters.
7. Export.

### PDF report

The PDF report includes:

- Study metadata.
- Publication-quality figures.
- DAR stability metrics.
- NCA parameter tables for all analytes.
- A methods section with software versions and PKNCA citation.

### Excel workbook

The Excel export includes:

- Raw Data.
- DAR Trajectory.
- NCA - Total ADC.
- NCA - Total Antibody.
- NCA - Free Payload.
- Biotransformation Summary.

## Limitations

- The app is designed for mean concentration profiles, not individual-subject PK analysis.
- DAR calculations assume mass-compatible assays; unit conversion may be needed if assay formats differ.
- AUCinf, CL, and Vz are not calculated by default unless the dataset includes a T=0 row.
- Multi-dose and steady-state PK are not currently supported.
- The current brand color palette may not be optimal for all color-vision conditions.

## Project Structure

```text
adc_pk_app/
├── app.R
├── report_template.Rmd
├── sample_data.csv
└── README.md
```

## Citation

PK analysis is powered by the PKNCA R package:

Denney WS, Duvvuri S, Buckeridge C. Simple, Automatic Noncompartmental Analysis: The PKNCA R Package. *J Pharmacokinet Pharmacodyn.* 2015. https://doi.org/10.1007/s10928-015-9432-2

## Research Use Only

This software is for research use only. It has not been validated for regulatory submission without independent verification against a qualified PK software system such as Phoenix WinNonlin or SAS.
