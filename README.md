# Fast-Fourier Transform

#### CSE305 Project - Duong DAC

### Build
`make`

### Executables
`test`

### Structure

All of the implemented FFT transform the sequence in-place.
Main functionalities

- `MFFT::Transform` computes the 1D FFT via the N-dimension DFT,
- `FFT::FFTParallel` computes the iterative radix-2 algorithm in parallel.

Extras

- `FFT::DFT`: 
- `FFT::FFTRec`
- `FFT::FFTSeq`

## Description

The DFT of a sequence $x$ can be defined as
$$X(k) = \sum_{n = 0}^{N-1} x(n) \omega_N^{kn}$$  
with $\omega_N = e^{-i2\pi/N}$ the root of unity. The inversion TODO.

A very common FFT is the Cooley–Tukey radix-2 algorithm. The idea is to divide the sequence $x$ into two ($x_{\text{even}}$, $x_{\text{odd}}$), performs the DFT recursively on each, and construct the transformed sequence $X$ from $X_{\text{even}}$ and $X_{\text{odd}}$. The algorithm can be implemented iteratively for better performance.

`FFT::FFTParallel` implements the iterative radix-2 algorithm with parallelized elements. The sequence has to be padded (with zeroes in our case) such that the size $N$ is a power of 2. In particular,

- `RevBitParallel` rearrange the array based on the binary representation of each index: $x_{i_1, ..., i_n} = x^{\text{old}}_{i_n,...,i_1}$.
- For each of the $\log_2(N)$ steps, `TransformParallel` compute the current FFT based on the result of the previous steps. Each partition $[x_{i}, x_{i+2^{S}}]$ can be processed individually on different threads.

### MFFT (multidimensional)
`MFFT::Transform` implements the parallelized FFT described in [(1)](<https://doi.org/10.1016/0167-8191(90)90031-4>) and [(2)](https://doi.org/10.1109/SUPERC.1994.344263). The 1-D sequence $x$ is mapped onto P dimensions $N = N_1 \times ... \times N_p$. To compute, for instance, the DFT of a 2D $N \times M$ array, we can compute $N$ DFTs along one dimension, then $M$ DFTs of the previous result along the other. The array does not have to be padded to have a size of $2^k$, but each $N_i$ should (and in our case must) be divisible by the number of threads. You can indeed 

This function is the main focus of the project. While the final implementation is straightforward, there are several elements that was hard to get right. To start with, we did not refer to the source code (if there's any) of both papers. In fact, since the paper described the algorithm adapted to the architecture of the time, our version does not closely follow the description. Some differences includes:

- The memory is shared by all threads.
- For each major step (twiddle, transform, transpose), 

The 

### Parallelization
FFT involves lots of array operations (assign, numerical operations) that fortunately are not in conflict with each other. In practice, the data and auxilliary variables 

**WIP** of notes

- computed the indicies btw. the most tricky part of mfft.
- parallelized computation of $x^n$ but might actually be redundant. IT IS REDUNDANT LMAO.
- the use of `<execution>`.

## Testing

Testing environment: Laptop, AMD 4750u, 32GB of RAM.



## Application

## Extras

Assume that $N = N_1 \times ... \times N_p$.

```
N1->NP:  x[n1][...][np] = x[np + ... + (n1)*N2*...*Np]

NP->N1: x0[np][...][n1] = x[n1][...][np]
                        = x[np + ... + (n1)*N2*...*Np]


NP->N1: xp[kp][...][k1] = X[kp][...][k1]
                        = X[k1 + N1*k2 + ... + N1*...*Np-1*(kp)]
```

$$x_p[k_p,...,k_1] = \sum_{n_1, ..., n_p} x_0[n_p,...,n_1] \omega_N^{f_p g_p}$$
where
$$f_j = [n_1,...,n_j] = n_j + ... + n_1 * N_2 * ... * N_j$$
and
$$g_j = [k_j,...,k_1] = k_1 + ... + N_1 * ... * N_{j-1} * k_j$$

Recursively, $f_j = n_j * f_{j-1}$ and $g_j = g_{j-1} + L_{j-1}k_j$.

$$
X_j[k_j][...][k_1][n_p][...][n_{j+1}] = \sum_{n_j = 0}^{N_j-1} X_{j-1}[k_{j-1}][...][k_1][n_p][...][n_{j}] \omega_{L_j}^{n_j g_{j-1}} \omega_{N_j}^{n_j k_j }
$$

```
Benchmark FFT::FFTRec:
 1 thread(s):          5,876        12,169        53,761       110,177       473,985       987,936    19,185,559    79,985,081

Benchmark FFT::FFTSeq:
 1 thread(s):          3,669         8,358        38,714        79,187       365,946       777,700    15,298,724    66,768,968

Benchmark FFT::FFTParallel:
 1 thread(s):          7,090        12,727        51,596       102,063       431,429       875,527    16,399,721    66,400,844
 2 thread(s):          6,868        12,319        49,454        97,917       410,845       875,725    10,592,994    43,867,619
 4 thread(s):          6,931        10,962        37,279        66,916       270,927       619,316     7,587,008    31,359,079
 8 thread(s):          7,497        10,972        31,914        59,430       224,654       456,690     6,224,218    27,583,350
16 thread(s):         14,503        18,095        34,172        60,801       214,588       418,912     5,794,894    24,762,653
32 thread(s):         20,216        25,326        43,795        72,407       202,037       360,616     5,856,252    25,580,457

Benchmark MFFT::Transform:
 1 thread(s):         34,472        73,671       315,031       667,066     2,969,815     6,050,454   107,953,848   465,013,701
 2 thread(s):         24,631        58,069       235,275       632,515     1,836,301     3,566,514    61,783,226   267,046,006
 4 thread(s):         15,079        22,989       110,238       234,732       798,842     1,423,182    22,274,169    92,561,386
 8 thread(s):         12,122        17,223        68,021       139,285       514,793       937,963    13,050,468    55,216,706
16 thread(s):         20,051        24,610        63,431       135,024       369,941       790,540    12,608,772    54,972,195
32 thread(s):         23,984        25,767        70,300       121,659       438,211       932,407    17,810,395    75,619,461
```
