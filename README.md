# ADC In Vivo DAR Monitoring & PK Tracker

An interactive R Shiny dashboard for integrating **Total ADC**, **Total Antibody**, and **Free Payload** concentration-time data from in vivo ADC studies.

The app combines these analytes into a single workflow to calculate **average DAR over time**, assess **payload release kinetics**, run **non-compartmental analysis (NCA)**, and generate **report-ready outputs** for study review and submission support [file:1].

## What this app does

ADC biotransformation in vivo is often tracked through separate assays for Total ADC, Total Antibody, and Free Payload. In many workflows, these results are then merged manually in Excel to understand deconjugation, linker cleavage, payload loss, and DAR stability over time [file:1].

This app automates that workflow by:

- Importing multi-analyte PK data from a CSV file [file:1].
- Calculating average DAR over time from Total ADC and Total Antibody measurements [file:1].
- Estimating DAR stability trends and decay over time [file:1].
- Quantifying Free Payload kinetics [file:1].
- Running NCA for all three analytes using the PKNCA R package [file:1].
- Producing publication-quality plots, tables, PDF reports, and Excel exports [file:1].

## Who this is for

This tool is designed for teams working with ADC PK and biotransformation data, including:

- Bioanalytical scientists developing or interpreting ADC assays [file:1].
- DMPK and pharmacokinetics scientists analyzing ADC disposition [file:1].
- Translational and pharmacology teams linking exposure to payload release [file:1].
- Study directors preparing integrated study summaries [file:1].
- CROs and internal teams that currently reconcile results manually in Excel [file:1].

## Key features

- Interactive multi-tab Shiny dashboard [file:1].
- DAR over time visualization with stability fitting [file:1].
- Free payload kinetics plots with Cmax and Tmax annotation [file:1].
- PK/PD overlay views combining DAR and payload behavior [file:1].
- NCA tables for Total ADC, Total Antibody, and Free Payload [file:1].
- Export to PDF report and Excel workbook [file:1].

## Installation

### Requirements

- R 4.2 or later [file:1].
- The following R packages:
  - shiny
  - shinydashboard
  - shinyjs
  - shinycssloaders
  - PKNCA (>= 0.11)
  - dplyr
  - tidyr
  - ggplot2
  - plotly
  - DT
  - openxlsx
  - rmarkdown
  - knitr
  - scales [file:1]

### PDF export requirements

To generate PDF reports, you also need a LaTeX installation [file:1]:

- Linux: `sudo apt-get install texlive-latex-base texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended lmodern`
- macOS: install MacTeX
- Windows: install MiKTeX [file:1]

### Run the app

From R or RStudio:

```r
shiny::runApp("path/to/adc_pk_app/")
```

Or open `app.R` in RStudio and click **Run App** [file:1].

## Usage

1. Launch the app.
2. Upload a CSV file with time-matched ADC data.
3. Confirm the concentration unit used in the study.
4. Review plots, DAR trajectory, and NCA outputs.
5. Export the PDF report or Excel workbook [file:1].

## Input format

Upload a CSV file with the following columns:

| Column | Required | Description |
|---|---:|---|
| `Time_h` | Yes | Nominal sampling time in hours [file:1]. |
| `Total_ADC_mean` | Yes | Mean Total ADC concentration [file:1]. |
| `Total_Ab_mean` | Yes | Mean Total Antibody concentration [file:1]. |
| `FreePayload_mean` | Yes | Mean Free Payload concentration [file:1]. |
| `Dose_mg_kg` | Yes | Administered dose in mg/kg [file:1]. |
| `Total_ADC_SD` | Optional | SD for Total ADC, used for error bars [file:1]. |
| `Total_Ab_SD` | Optional | SD for Total Antibody [file:1]. |
| `FreePayload_SD` | Optional | SD for Free Payload [file:1]. |
| `Species` | Optional | Species label shown in reports [file:1]. |
| `Study_ID` | Optional | Study identifier [file:1]. |

All concentration columns must use the same unit, selected in the app sidebar. Supported units include `ug/mL`, `ng/mL`, `nM`, and `pM` [file:1].

A blank template can be downloaded from the app’s Upload tab [file:1].

## Analysis workflow

### DAR calculation

The app estimates average DAR over time using:

```text
avg_DAR(t) = [Total_ADC(t) / Total_Ab(t)] x Nominal_DAR
```

This approach assumes both analytes are measured in compatible mass-based units and come from the same antibody backbone [file:1]. The app also fits an exponential decay model to summarize DAR stability over time [file:1].

### NCA

NCA is performed using PKNCA [file:1]. The default interval is from the first observed timepoint to the last observed timepoint, and the app reports parameters such as AUClast, Cmax, Tmax, half-life, lambda-z, and R-squared [file:1].

For IV bolus datasets, `tmax = 0` is expected when the first measured timepoint is also the peak observed value relative to the interval start [file:1]. AUCinf, CL, and Vz are not computed by default unless a T=0 row is included in the input data [file:1].

## Output

### Interactive dashboard tabs

The app includes seven main tabs [file:1]:

1. Upload & Configure.
2. PK Curves.
3. DAR Trajectory.
4. Payload Kinetics.
5. PK/PD Overlay.
6. NCA Parameters.
7. Export.

### PDF report

The PDF report includes [file:1]:
- Study metadata.
- Publication-quality figures.
- DAR stability metrics.
- NCA parameter tables for all analytes.
- A methods section with software versions and PKNCA citation.

### Excel workbook

The Excel export includes [file:1]:
- Raw Data.
- DAR Trajectory.
- NCA - Total ADC.
- NCA - Total Antibody.
- NCA - Free Payload.
- Biotransformation Summary.

## Limitations

- The app is designed for mean concentration profiles, not individual-subject PK analysis [file:1].
- DAR calculations assume mass-compatible assays; unit conversion may be needed if assay formats differ [file:1].
- AUCinf, CL, and Vz are not calculated by default unless the dataset includes a T=0 row [file:1].
- Multi-dose and steady-state PK are not currently supported [file:1].
- The current brand color palette may not be optimal for all color-vision conditions [file:1].

## Files

```text
adc_pk_app/
├── app.R
├── report_template.Rmd
├── sample_data.csv
└── README.md
```

## Citation

PK analysis is powered by the PKNCA R package [file:1]:

Denney WS, Duvvuri S, Buckeridge C. Simple, Automatic Noncompartmental Analysis: The PKNCA R Package. *J Pharmacokinet Pharmacodyn.* 2015. https://doi.org/10.1007/s10928-015-9432-2

## Research use only

This software is for research use only. It has not been validated for regulatory submission without independent verification against a qualified PK software system such as Phoenix WinNonlin or SAS [file:1].
