# SpatialDengueRefit
 
A spatially explicit, patch-based stochastic model of dengue transmission, extended to evaluate *Wolbachia* replacement and reactive case-area targeted interventions (CATI) in Rio de Janeiro, Brazil.
 
This repository builds on the **SpatialDengue** package developed by Oliver Brady and colleagues and described in:
 
> Brady OJ et al. (2021). Case-area targeted interventions (CATI) for reactive dengue control: Modelling effectiveness of vector control and prophylactic drugs in Singapore. *PLoS Neglected Tropical Diseases* 15(8): e0009562. https://doi.org/10.1371/journal.pntd.0009562
 
The underlying model (human SIR and mosquito SEI compartments in 500 m–1 km patches, daily stochastic time steps, and exponential, gravity or radiation human movement) is unchanged. The additions below are optional: unless they are used, `DEN.spatial()` behaves as in the original package.
 
---
 
## Changes from the original SpatialDengue package
 
### Model code (`R/`)
 
| Change | Where | Description |
|---|---|---|
| Wolbachia transmission reduction | `DEN.spatial()` | New optional `paramsList` element `betaEnv_reduction`: a proportional reduction in each patch's transmission coefficient, either a single value or one value per patch. |
| Area-based CATI coverage | `DEN.spatial()` | `vectreat` can include `Cov_index`, `Cov_edge` and `Cov_corner`: the share of the index case's patch and of its edge and corner neighbours covered by a circle of the CATI radius. Without these columns, CATI behaves as in the original package. |
| `weekdates` renamed to `timestep` | `DEN.spatial()` | The model can be used with case data at any time-step resolution (e.g. weeks or months). The first element selects the seed rows of `fitdat` (`fitdat$Week == timestep[1]`). Calls that pass this argument by position are unaffected. |
| Safe binomial draws | `DEN.spatial()`, `Multinom.R` | `rb_safe()` guards binomial draws against invalid sizes or probabilities. |
| Parallel ensembles | `DEN.spatial.ensemble()` | Simulations run in parallel with `parallel::mclapply()`. |
 
Other functions (`DEN.spatial.seasonal()`, `DEN.spatial.RCT()`, `DEN.spatial.UCT()`, `make.eval.set()`) are unchanged and still assume weekly time steps.
 
### Data (`data/`)
 
- **`finalWeights`** now contains the Rio de Janeiro fit: round-3 ABC-SMC parameter sets for the radiation movement model, weighted 1 − d/max(d) following Brady et al. (2021). It replaces the Singapore parameter sets shipped with the original package; `DEN.spatial()` samples one parameter set per simulation from it.
- **Rio de Janeiro inputs** (`*_rio.rda`): dengue case and *Wolbachia* data (`wmel_rio`), model patches (`unipix_rio`), distances (`pixdistmat_rio`), population (`riopop_rio`), immunity landscape (`stim_vec_rio`), seed infections (`fitdat_rio`, `pastdat_rio`), time steps (`timestep_rio`) and the daily simulation calendar (`sim_calendar_rio`).
`wmel_rio` is derived from the publicly available data of Ribeiro dos Santos et al. (2022), *The Lancet Infectious Diseases* 22(11): 1587–1595 (https://github.com/pdgcam/wMel_Rio_dengue_chik).
 
---
 
## Rio de Janeiro analysis
 
The two analysis files support the manuscript:
 
> Tiley K et al. *[Title TBC]*. *[Journal TBC]*.
>
> **Abstract:** *[TBC]*
 
| File | Purpose |
|---|---|
| `Rio_refit_workflow.Rmd` | Prepares the Rio de Janeiro inputs and refits the model to monthly dengue case data (November 2011–October 2012) by ABC-SMC for three human movement models, then selects the best-fitting model and runs posterior predictive checks. Writes the `*_rio.rda` inputs and `finalWeights` to `data/`. |
| `install_and_run_simulation_Rio.R` | Calibrates the *Wolbachia* transmission reduction against published field estimates and runs the intervention scenarios: (RQ1) observed vs homogeneous *Wolbachia* coverage, (RQ2) indirect protection of untreated patches, (RQ3) heterogeneous vs homogeneous coverage at matched mean prevalence, and (RQ4) combined *Wolbachia* and CATI. |
| `smc_checkpoints/smc_*_fullfit.rds` | Saved ABC-SMC fits for the exponential, gravity and radiation movement models, so the Rmd can be knitted without refitting. |
 
### Running the analysis
 
1. Clone this repository and open `SpatialDengueRefit.Rproj` in RStudio, so that the working directory is the repository root.
2. Install the packages loaded at the top of each file (`library()` calls), plus `devtools`.
3. Both files load the package with `devtools::load_all()`, so they use the version of the model in this repository. No separate installation is needed.
4. **Refitting workflow (optional):** knit `Rio_refit_workflow.Rmd`. The three fitting chunks are not evaluated when knitting; the saved fits in `smc_checkpoints/` are loaded instead. To refit from scratch, set those chunks to `eval = TRUE` (this takes many hours per movement model).
5. **Scenario analysis:** run `install_and_run_simulation_Rio.R` from top to bottom. It uses the inputs and `finalWeights` in `data/`, so it can be run without first running the Rmd. Outputs (`.rds` summaries and figures) are saved to the repository root with the date prefix set in `OUT`. A full run takes many hours so it's recommended to do a small test run on your machine to better estimate total duration.
### Notes
 
- Each scenario is simulated 100 times, and simulation *n* of every scenario uses the same random seed as simulation *n* of the no-intervention scenario (same posterior draw, transmission landscape and initial infections), so that effectiveness is calculated from paired simulations.
- Simulations run in parallel on 8 cores on Mac using `parallel::mclapply()`, which relies on forking. On Windows it runs on a single core; reduce the number of cores if memory is limited.
---
 
## Licence
 
MIT.
 
