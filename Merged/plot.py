import sys

import matplotlib.pyplot as plt
import numpy as np
import uproot

signal_background_filename = sys.argv[1] if len(sys.argv) > 1 else (
    "results_Merged/run_1_10_0_800000/mcmc_neutrinofit_SB.root"
)
background_filename = sys.argv[2] if len(sys.argv) > 2 else (
    "results_Merged/run_1_10_0_800000/mcmc_neutrinofit_B.root"
)

def read_chain(filename, tree_name):
    data = uproot.open(filename)
    if tree_name not in data:
        raise KeyError(
            f"Tree '{tree_name}' not found in {filename}. "
            "The parameters tree is not the MCMC tree."
        )
    tree = data[tree_name]

    required_branches = {"Phase", "nb"}
    missing_branches = required_branches.difference(tree.keys())
    if missing_branches:
        raise KeyError(
            f"Missing branches in '{tree_name}': "
            f"{', '.join(sorted(missing_branches))}"
        )
    arrays = tree.arrays(["Phase", "nb"], library="np")
    mask = arrays["Phase"] > 0
    if not np.any(mask):
        raise RuntimeError(f"The MCMC tree in {filename} has no Phase > 0 entries.")
    return arrays["nb"][mask], tree

nb_sb, sb_tree = read_chain(signal_background_filename, "neutrinofit_SB_mcmc")
nb_b, b_tree = read_chain(background_filename, "neutrinofit_B_mcmc")

if "ns" not in sb_tree.keys():
    raise KeyError("Branch 'ns' not found in the signal-plus-background MCMC tree.")
sb_arrays = sb_tree.arrays(["Phase", "ns"], library="np")
ns = sb_arrays["ns"][sb_arrays["Phase"] > 0]

print("90% posterior upper limit on ns (S+B):", np.quantile(ns, 0.90))
print("90% posterior upper quantile on nb (S+B):", np.quantile(nb_sb, 0.90))
print("90% posterior upper quantile on nb (B):", np.quantile(nb_b, 0.90))

plt.figure()
plt.hist(nb_sb, bins=50, density=True, alpha=0.6, label="S+B")
plt.hist(nb_b, bins=50, density=True, alpha=0.6, label="B only")
plt.xlabel("nb")
plt.ylabel("Posterior density")
plt.legend()
plt.tight_layout()
plt.savefig("nb_model_comparison.pdf")

plt.figure()
plt.hist(ns, bins=50, density=True, alpha=0.8, label="S+B")
plt.xlabel("ns")
plt.ylabel("Posterior density")
plt.legend()
plt.tight_layout()
plt.savefig("ns.pdf")
