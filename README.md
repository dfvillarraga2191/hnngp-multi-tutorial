# HNNGP tutorial (multinomial mode choice with rstan)

Tutorial materials for estimating a **hierarchical nearest neighbor Gaussian process (HNNGP)** for three-alternative mode choice using [`rstan`](https://mc-stan.org/users/interfaces/rstan.html). This is the **IIA contrast** specification: car and transit utilities are nearest-neighbor Gaussian processes measured against a non-motorized utility fixed at zero. This repository accompanies Villarraga and Daziano (2025).

The binary tutorial is published separately at <https://github.com/dfvillarraga2191/hnngp-tutorial>.

## Citation

If you use this tutorial, please cite the paper:

> Villarraga, D. F., & Daziano, R. A. (2025). Hierarchical Nearest Neighbor Gaussian Process models for discrete choice: Mode choice in New York City. *Transportation Research Part B: Methodological*, 191, 103132. <https://doi.org/10.1016/j.trb.2024.103132>

A machine-readable citation is in [`CITATION.cff`](CITATION.cff). The code is released under the [MIT License](LICENSE).

Repository: <https://github.com/dfvillarraga2191/hnngp-multi-tutorial>

## Contents

| Path | Description |
|------|-------------|
| `HNNGP_multi_tutorial.Rmd` | Tutorial source (R Markdown) |
| `docs/HNNGP_multi_tutorial.html` | Self-contained HTML (for a website or GitHub Pages) |
| `docs/tutorial.css` | Stylesheet used by the HTML |
| `data/mode_multi.csv` | Example three-alternative mode-choice data |
| `stan/nngp_multinomial.stan` | IIA multinomial Stan model |
| `R/NNMatrix.R` | Nearest-neighbor matrix helpers (`spNNGP`) |
| `R/ess_rhats.R` | ESS (`sns`) and $\widehat{R}$ (`asbio`) diagnostics |
| `scripts/fit_hnngp.R` | End-to-end estimation script |

## Requirements

- R ≥ 4.x with a working C++ toolchain (for Stan)
- Packages: `rstan`, `spNNGP`, `tidyverse`, `HDInterval`, `ggridges`, `loo`, `sns`, `asbio`, `rmarkdown`

```r
install.packages(c(
  "rstan", "spNNGP", "tidyverse", "HDInterval",
  "ggridges", "loo", "sns", "asbio", "rmarkdown"
))
```

## Quick start

From this directory:

```r
source("scripts/fit_hnngp.R")
```

For a faster trial run, uncomment the subsample block inside `scripts/fit_hnngp.R` (or in the Rmd).

`R/ess_rhats.R` assumes 4 chains, 3000 iterations, and 1000 warmup draws. Those are the settings in `scripts/fit_hnngp.R`. Change them together if you edit the `stan()` call.

## Rebuild the HTML

```r
rmarkdown::render(
  "HNNGP_multi_tutorial.Rmd",
  output_dir = "docs",
  output_file = "HNNGP_multi_tutorial.html"
)
```

The HTML is **self-contained** (`self_contained: true`), so you can upload `docs/HNNGP_multi_tutorial.html` alone, link it from a site, or serve the `docs/` folder with GitHub Pages.

## Website embedding

- **Link / download:** host `docs/HNNGP_multi_tutorial.html` and link to it.
- **iframe:** `<iframe src="path/to/HNNGP_multi_tutorial.html" width="100%" height="900" style="border:0"></iframe>`
- **GitHub Pages:** push this folder and set Pages to `/docs`. The site root opens the tutorial.

## Diagnostics

After fitting, with `fit` in the workspace:

```r
source("R/ess_rhats.R")
```

This prints and returns ESS (summed across chains via `sns::ess`) and Gelman–Rubin $\widehat{R}$ (`asbio::R.hat`) for $\beta$, $\gamma$, the alternative-specific constants, $\sigma^2$, and $\phi$.

## License / data

Example data are an extract for this tutorial and still contain trip origin and destination coordinates. Household, person, and census-tract identifiers are not included. Adapt paths and priors to your own study as needed.
