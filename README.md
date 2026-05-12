# Projecting climate change impacts on health: A tutorial integrating the latest climate and demographic scenarios

Fully reproducible code for the illustrative example presented in:  
**[ADD CITATION AND DOI]**

This repository illustrates the methodology presented in the paper through a projection study of heat-related mortality in London under a middle-of-the-road scenario (SSP2-4.5).

## What this repository includes

The code:

- Provides all files necessary to run the projection study and demonstrates **how to download all open-access data** from the original sources.  
  **Note:** Downloading and processing climate model data can be substantially slower than the rest of the workflow.

- Shows how to fit the **epidemiological models**.

- Processes raw **RCP climate projections** to derive regional temperature time series and applies bias correction using the temperature observations used in the epidemiological analysis.

- Calibrates raw **SSP demographic projections** to match the spatial and temporal resolution of the study.

- Computes **climate- and demographic-driven health impacts** while accounting for uncertainty from both epidemiological and climate models.

- Summarises health impact projections for specific time windows (e.g., decades) and aggregated populations, reporting point estimates and 95% empirical confidence intervals (CIs) that combine all sources of uncertainty.

- **Generates all figures** presented in the paper.


