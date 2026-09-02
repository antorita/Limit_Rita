# Merged CYGNO Bayesian fit

This folder combines the Samuele command-line workflow with the Rita BAT model idea.
It fits the same toy-file format as Samuele while using the requested priors:

- background: Poisson prior with mean equal to the second value in the toy file;
- signal: uniform prior between zero and the observed total number of events.

## Input format

A toy file contains:

1. signal metadata count (`ns`);
2. background metadata count (`nb`);
3. one observed count per bin.

`background.txt` and `signal.txt` contain one template value per line. Their lengths must equal the number of toy bin values.

## Build

Run from this directory, with ROOT and BAT installed:

```bash
cd /Users/rita/Documents/limit/Limit/Merged
make clean
make
```

The Makefile uses `root-config` and links BAT from `/usr/local/lib`.

## Run

Run from `Merged`, because result paths are relative:

```bash
./runfit \
  ../Samuele/toyMC_txt/1_10/0.txt \
  ../Samuele/background.txt \
  ../Samuele/signal.txt \
  0 1_10 NID
```

The executable requires exactly six arguments after `runfit`:

```text
DATA BKG SIG RUN_INDEX REDUCTION_EXPOSURE FLAG
```

For example, `1_10` is split into reduction `1` and exposure `10`.

## Output

Results are written under:

```text
results_Merged/run_1_10_0_800000/
```

The directory contains BAT Markov-chain ROOT files, including `mcmc_neutrinofit_SB.root` and the background-only chain, plus the BAT model-comparison output. The Bayes-factor summary is appended to `pippo_newAR_10_NID.txt` in the current directory.

## Statistical model

For bin `i`, the expected count is:

```text
lambda_i = nb * bkg_i / sum(bkg) + ns * sig_i / sum(sig)
```

The likelihood is the product of Poisson probabilities over bins. The S+B model has parameters `nb` and `ns`; the background-only model has only `nb`.

The models are passed to BAT's MCMC marginalization and then integrated to compute the Bayes factor `B(S+B, B)`.

## Compare the models

The two MCMC ROOT files can be compared with:

```bash
python3 plot.py \
  results_Merged/run_1_10_0_800000/mcmc_neutrinofit_SB.root \
  results_Merged/run_1_10_0_800000/mcmc_neutrinofit_B.root
```

The script reads `neutrinofit_SB_mcmc` from the S+B file and `neutrinofit_B_mcmc` from the background-only file. It overlays their `nb` posterior distributions in `nb_model_comparison.pdf`. The `ns` posterior exists only in the S+B model and is saved as `ns.pdf`.

This plot compares parameter posteriors. The formal comparison between the two hypotheses is the Bayes factor printed by `runfit` and stored in `pippo_newAR_<exposure>_<flag>.txt`.
