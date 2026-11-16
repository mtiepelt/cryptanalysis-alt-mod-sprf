# <CryptanalysisAltModSPRF>

Artifact for the paper **"Cryptanalysis of the Alternating-Moduli Strong PRF"** by Mehdi Tibouchi, Marcel Tiepelt,
IACR Transactions on Symmetric Cryptology 2026, Issue 04.

The code implements the collision attack (CollisionSearch, Fig. 1 of the paper) on the
(2,3) strong PRF of [Def. 9, https://eprint.iacr.org/2024/582.pdf]. It contains

* a proof of concept that finds collisions in the kernel regime and in the generalized regime,
* the computation of the bounds on beta and T (Tables 4, 5 and 6),
* the experiments on the distribution of Hamming weights (Figure 2),
* the experiments on the success rate of the collision attack (Figure 3).


## Dependencies

* SageMath 11.0.beta3
* Python 3.14.4 
* Python packages, installed for the Python of SageMath:
    * numpy 2.4.4
    * matplotlib 3.10.8

All scripts are run with `sage -python`. Missing packages can be installed with

```
sage -pip3 install numpy==2.4.4 matplotlib==3.10.8
```

Tested on Ubuntu 26.04.1, x86-64, 8 cores, 32 GB RAM.
Requires about 2.5 GB RAM per core to reproduce all results.

## Building

Nothing needs to be built. The output folders `results/attack`, `results/distribution`, `results/parameters`
and `results/data` are created automatically on the first run.
The scripts can be started from any folder.

## Quick check

Runs small instances of every experiment.

| Command | Time | Expected result |
|---|---|---|
| `sage -python expSimpleCollisionSearch.py` | 1 s | `Collision in kernel regime ...` and `Collision in generalized regime ...` (see below) |
| `sage -python parameters.py` | < 1 s | Tables 4, 5 and 6 printed to the console |
| `sage -python expDistribution.py -n 12 --num_exps 10` | 5 s | `Wrote file .../results/distribution/secpar-12_...pdf` |
| `sage -python expCollisionSearch.py -n 12 --num_exps 10` | 6 s | `Wrote file .../results/attack/CollisionAttack_secpar-12_...pdf` |


## Reproducing the paper's results

| Paper result | Command | Time | Output |
|---|---|---|---|
| Tables 4, 5, 6: bounds on beta and T | `sage -python parameters.py` | seconds | console and `results/parameters/` |
| Figure 2: distribution of Hamming weights | `sage -python expDistribution.py --reproduce-paper` | 7 min | `results/distribution/` |
| Figure 3: collision attack | `sage -python expCollisionSearch.py --reproduce-paper` | 3 hours | `results/attack/` |

`--reproduce-paper` runs the paper's computations with the paper's values:
secpar for each 20, 64, 128, num_exps = 100, an = 4.4, `q = 5`, T = 2^19, const = 0.9.
All other options except `--num-cpus` are ignored.

Calling `sage -python expDistribution.py --reproduce-paper` is equivalent to 
```
    sage -python expDistribution.py -n 20 --logT2 19
    sage -python expDistribution.py -n 64 --logT2 19
    sage -python expDistribution.py -n 128 --logT2 19
```
which runs in a few minutes. 

Calling `sage -python expCollisionSearch.py --reproduce-paper` is equivalent to  
```
    sage -python expCollisionSearch.py -n 20 --logT2 19
    sage -python expCollisionSearch.py -n 64 --logT2 19
    sage -python expCollisionSearch.py -n 128 --logT2 19
```
which runs in about 0:36, 1:05 and 2:10 hours.


**Memory:** SageMath's `right_kernel_matrix` leaks memory in our experiments. Every experiment therefore runs in
a new process, which frees the memory afterwards. Use `--num-cpus` to limit the number
of parallel processes if memory runs out. Larger values of `q` increase the number of kernel elements to `3**q`, 
reducing the number of calls to `right_kernel_matrix`.

**Randomness:** the experiments for Figures 2 and 3 use fresh randomness on every run,
so the numbers vary slightly between runs. The proof of concept `expSimpleCollisionSearch`
uses the fixed seed `2026433` and gives the same output on every run.


## Configuration and run modes

### `expDistribution.py` and `expCollisionSearch.py`

Without `--reproduce-paper`, the experiment runs for one security parameter given with `-n`.

| Option | Meaning | Default |
|---|---|---|
| `-n`, `--secpar` | security parameter | required, unless `--reproduce-paper` |
| `--reproduce-paper` | run the paper's computations (see above) | off |
| `--num_exps` | number of experiments | 100 |
| `--an` | n = ⌈an * secpar⌉ | 4.4 |
| `--q` | d = n − l + q | 5 |
| `--log2T` | T = 2^log2T | secpar - 1 |
| `--const` | target success probability in the success condition (Eq. 3) | 0.9 |
| `--num-cpus` | number of CPUs | all but two |

beta is not an option: `expDistribution.py` uses the largest beta that satisfies the success
condition (Eq. 3), and `expCollisionSearch.py` uses values of beta around it.

### `expSimpleCollisionSearch.py`

| Option | Meaning | Default |
|---|---|---|
| `--an` | n = ⌈an · secpar⌉, both regimes | 4 |
| `--kernel_secpar` | security parameter, kernel regime (d = n + 1) | 5 |
| `--secpar` | security parameter, generalized regime | 10 |
| `--al` | beta, l = floor(beta * secpar), generalized regime | 1.2 |
| `--q` | d = n − l + q, generalized regime | 5 |
| `--log2T` | T = 2^log2T, generalized regime | secpar − 1 |
| `--seed` | random seed, both regimes | 2026433 |

All scripts list their options with `--help`.


## Interpreting the output

### Proof of concept (`expSimpleCollisionSearch.py`)

For each regime, the script prints a collision pair x, x2 and the values
z = G * (x‖1) mod 3 (`z1_col_F3`, `z2_col_F3`) and mod 2 (`z1_col_F2`, `z2_col_F2`).

The code uses a hard-coded seed `2026433` and should produce the following output: 

```
------------------------------------------------------------------------------------------
Collision in kernel regime for n=20, d=21
-----
x        : 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 1 0 0 0 1 0
x2       : 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0
-----
z1_col_F3: 1 2 0 0 0 2 2 2 2 1 1 1 0 1 1 1 1 0 1 2
z2_col_F3: 1 2 0 0 0 2 2 2 2 1 1 1 0 1 1 1 1 0 1 2
-----
z1_col_F2: 1 1 0 0 0 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1
z2_col_F2: 1 1 0 0 0 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
Collision in generalized regime for n=40, d=33
x marks positions potentially != mod 3, but = mod 2
-----
x        : 0 0 1 1 0 0 0 0 0 1 1 0 0 0 1 0 1 1 1 0 0 1 1 1 1 1 1 1 1 1 0 1 1
x2       : 1 1 0 0 0 0 0 1 1 0 0 0 1 0 1 1 0 1 0 1 1 0 0 1 0 1 0 0 0 0 0 1 1
-----
           x         x               x x   x   x       x                     x x   x x   x
z1_col_F3: 1 1 0 0 2 2 1 2 1 2 2 2 0 2 0 0 0 0 0 0 1 1 2 0 1 2 1 0 1 2 2 1 2 1 0 1 1 2 1 1
z2_col_F3: 2 1 0 0 2 1 1 2 1 2 2 2 0 2 0 0 0 0 0 0 1 1 1 0 1 2 1 0 1 2 2 1 2 1 0 1 1 2 1 2
-----
           x         x               x x   x   x       x                     x x   x x   x
z1_col_F2: 1 1 0 0 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 1 1 0 1 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1
z2_col_F2: 1 1 0 0 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 1 1 0 1 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1
------------------------------------------------------------------------------------------

```

* Kernel regime: z1 and z2 are equal already mod 3.
* Generalized regime: the positions marked `x` may differ mod 3, but z1 and z2 are equal mod 2.

In both cases the script has checked that x and x2 give the same PRF output.

### Tables 4, 5, 6 (`parameters.py`)

Printed to the console and written as CSV to `results/parameters/` (to read into TeX/PGF):

* Table 4 (`al-T_bounds_s-128_an-4.4_q-0.tab`): smallest log_2 T_2 for each beta, secpar = 128, q = 0.
* Table 5 (`T-secpar_bounds_an-4.4_q-0.tab`): largest beta and l = floor(beta ** secpar) for T_1 = T_2 = 2^(secpar−1), q = 0.
* Table 6 (`T-19_bounds_an-4.4_q-5.tab`): largest beta and l for T_1 = T_2 = 2^19, q = 5.

### Figure 2 (`expDistribution.py`)

For each security parameter, `results/distribution/` contains a plot (`.pdf`), the same plot
with logarithmic y-axis (`_log.pdf`), and the data:

* gray bars: Hamming weights of uniform vectors in F_3^l,
* blue bars: Hamming weights hw(z_1), ..., hw(z_N) of the kernel elements (average over the experiments),
* purple bars: the R lowest of these Hamming weights (each experiment keeps its own R lowest; average over the experiments),
* red dashed line: the expected Hamming weight `E[hw(z_R)]` of the R-th lowest element.

The file `info_<...>.data` contains the parameters, the expected value `E[X_r]` and the
measured value `Max(X_r)` (average over the experiments) of the R-th lowest Hamming weight.

### Figure 3 (`expCollisionSearch.py`)

For each security parameter, `results/attack/CollisionAttack_<...>.pdf` shows the success
rate of CollisionSearch (fraction of experiments that found a collision within T_2 = T
evaluations) for values of beta around the predicted largest beta (red dashed line).
`CollisionAttack_<...>.data` contains the numbers.

### Regenerating plots

The plots can be regenerated from the saved data without rerunning the experiments,
with `distributionHW.plot_from_data(secpar, an, al, q, num_exps)` and
`collisionAttack.plot_from_file(secpar, an, q, num_exps)`.


## How the source code is organized

```
.
├── expSimpleCollisionSearch.py   # proof of concept, kernel and generalized regime
├── expDistribution.py            # Figure 2
├── expCollisionSearch.py         # Figure 3
├── parameters.py                 # Tables 4, 5, 6; parameters and success condition (Eq. 3)
├── paths.py                      # output folders
├── src/
│   ├── strongPRF.py              # the (2,3) strong PRF
│   ├── collisionAttack.py        # CollisionSearch (Fig. 1): steps (1) and (2), experiments for Figure 3
│   ├── distributionHW.py         # experiments for Figure 2
│   ├── plots.py                  # plots for Figures 2 and 3
│   └── pcommon.py                # shared helpers (randomness, kernel Hamming weights, file I/O)
└── results/                        # output, created on the first run
```


## License

MIT License, see [LICENSE](LICENSE).


## Citation

```bibtex
@article{<key>,
  author  = {<Authors>},
  title   = {<Paper title>},
  journal = {IACR Transactions on Symmetric Cryptology},
  volume  = {<Year>},
  number  = {<N>},
  pages   = {<pages>},
  doi     = {<doi>}
}
```


## Contact

Marcel Tiepelt — research (dot) mtiepelt (at) mailbox (dot) org 
https://mtiepelt.github.io/
