# Fast-Fourier Transform

#### CSE305 Project - Duong DAC

### Build

`make`

### Structure

All of the implemented FFT transform the sequence in-place. No locks are used.

Main functionalities

- `MFFT::Transform` computes the 1D FFT via the N-dimension DFT,
- `FFT::FFTParallel` computes the non-recursive radix-2 algorithm in parallel.

Extras

- `FFT::DFT`
- `FFT::FFTRec`
- `FFT::FFTSeq`

## Description

The DFT of a sequence $x$ can be defined as
$$X(k) = \sum_{n = 0}^{N-1} x(n) \omega_N^{kn}$$  
with $\omega_N = e^{-i2\pi/N}$ the root of unity.

A very common FFT is the Cooley–Tukey radix-2 algorithm. The idea is to divide the sequence $x$ into two ($x_{\text{even}}$, $x_{\text{odd}}$), performs the DFT recursively on each, and construct the transformed sequence $X$ from $X_{\text{even}}$ and $X_{\text{odd}}$. The algorithm can be implemented iteratively for better performance.

`FFT::FFTParallel` implements the iterative radix-2 algorithm with parallelized elements. The sequence has to be padded (with zeroes in our case) such that the size $N$ is a power of 2. In particular,

- `RevBitParallel` rearrange the array based on the binary representation of each index: $x_{i_1, ..., i_n} = x^{\text{old}}_{i_n,...,i_1}$.
- For each of the $\log_2(N)$ steps, `TransformParallel` compute the current FFT based on the result of the previous steps. Each partition $[x_{i}, x_{i+2^{S}}]$ can be processed individually on different threads.

`MFFT::Transform` implements the parallelized FFT described in [(1)](<https://doi.org/10.1016/0167-8191(90)90031-4>) and [(2)](https://doi.org/10.1109/SUPERC.1994.344263). The 1-D sequence $x$ is mapped onto P dimensions $N = N_1 \times ... \times N_p$. To compute, for instance, the DFT of a 2D $N \times M$ array, we can compute $N$ DFTs along one dimension, then $M$ DFTs of the previous result along the other.

**WIP** of notes

- computed the indicies btw. the most tricky part of mfft.
- parallelized computation of $x^n$ but might actually be redundant.
- the use of `<execution>`.

## Testing

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
