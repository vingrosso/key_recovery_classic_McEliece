## Description

This repository comes with the article entitled "Single-Trace Key-Recovery Template Attack on Classic McEliece Decapsulation"

  - [Installation](#installation)
  - [Usage](#usage)
  - [Reproducibility](#reproducibility)
  
## Installation

The following Python packages are required:

  - coloredlogs
  - galois
  - msgpack
  - numpy
  - tqdm

They can be installed with:

```
$ python3 -m pip install -r requirements.txt
```

### Unzipping precomputed hash tables

Precomputed hash tables can be unzipped with:

```
$ gunzip -k F2\^9_n_pow_22.msgpack.gz
$ gunzip -k F2\^12_n_pow_22.msgpack.gz
$ gunzip -k F2\^13_n_pow_22.msgpack.gz
```

## Usage
### Usage on ChipWhisperer traces

Use scripts in the `side-channel` folder

### Usage on simulated side-channel traces

```
python3 attach.py n t m --accuracy=a
```

#### Perfect accuracy

```
$ python3 attack.py 512 20 9                
[XX:XX:02] Keygen starts
[XX:XX:04] Keygen done
[XX:XX:04] Done simulating side-channel trace of shape (n, 2t) = (512, 40)
[XX:XX:04] Loading precomputed hash table for F2^9 and n_powers=22
[XX:XX:04] Loading done
[XX:XX:04] Actual attack starts !
[XX:XX:04]   HW lists matching completed after 190 iterations (190 minimum), still 322 left
[XX:XX:04]   [Alg.1 l.2] Recovered mt+δ pairs (ɑ, g(ɑ))
[XX:XX:05]   [Alg.1 l.3] Polynomial g recovered with interpolation
[XX:XX:05]     Searching for an invertible submatrix in Hpub
[XX:XX:05]     Invertible submatrix found in Hpub
[XX:XX:05]   [Alg.1 l.4] Vandermonde matrix V constructed with g and the mt+δ (ɑ, g⁻²(ɑ)) pairs
[XX:XX:05]   [Alg.1 l.5] Computed the change-of-basis S using V and H_pub
[XX:XX:05]   [Alg.1 l.6] H_priv = S⁻¹H_pub recovered 
0.000746, 0.405126, 0.052468, 0.050393, 0.010012, 0.000143, 
[XX:XX:05]   [Alg.1 l.7] Full permuted support L recovered 
[XX:XX:05] Attack successful ! Done in 0:00:00.520167
```

#### Custom accuracy

```
$ python3 attack.py 512 20 9 --accuracy=0.97
[XX:XX:12] Keygen starts
[XX:XX:15] Keygen done
[XX:XX:15] Done simulating side-channel trace of shape (n, 2t) = (512, 40)
[XX:XX:15] Loading precomputed hash table for F2^9 and n_powers=22
[XX:XX:15] Loading done
[XX:XX:15] Actual attack starts !
[XX:XX:15]   HW lists matching completed after 349 iterations (190 minimum), still 163 left
[XX:XX:15]   [Alg.1 l.2] Recovered mt+δ pairs (ɑ, g(ɑ))
[XX:XX:15]   [Alg.1 l.3] Polynomial g recovered with interpolation
[XX:XX:15]     Searching for an invertible submatrix in Hpub
[XX:XX:15]     Invertible submatrix found in Hpub
[XX:XX:15]   [Alg.1 l.4] Vandermonde matrix V constructed with g and the mt+δ (ɑ, g⁻²(ɑ)) pairs
[XX:XX:15]   [Alg.1 l.5] Computed the change-of-basis S using V and H_pub
[XX:XX:15]   [Alg.1 l.6] H_priv = S⁻¹H_pub recovered 
0.004836, 0.4443, 0.028606, 0.054372, 0.010048, 0.000141, 
[XX:XX:15]   [Alg.1 l.7] Full permuted support L recovered 
[XX:XX:15] Attack successful ! Done in 0:00:00.543559
```

```
$ python3 attack.py 512 20 9 --accuracy=0.9 
[XX:XX:17] Keygen starts
[XX:XX:20] Keygen done
[XX:XX:20] Done simulating side-channel trace of shape (n, 2t) = (512, 40)
[XX:XX:20] Loading precomputed hash table for F2^9 and n_powers=22
[XX:XX:20] Loading done
[XX:XX:20] Actual attack starts !
[XX:XX:20]   Could not match enough HW lists, only 56/190. Possible fixes are:
[XX:XX:20]     - decrease delta (currently 10),
[XX:XX:20]     - get a more accurate side-channel distinguisher.
```

## Reproducibility

### Execution time

Execution times for the 6 lines of Algorithm 1 (Table 7) are measured with:

```
$ ./measure_timings.sh
```

### Success rate

The attack success rate (Figure 3) is computed with:

```
$ python3.8 compute_success_rate.py
```