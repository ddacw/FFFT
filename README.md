# Fast-Fourier Transform

#### CSE305 Project - Duong DAC

`Note: the generated PDF does not display the test result properly. Please check README.md on GitHub!`

### Build

`make`

### Executables

`test` `main`

### Structure

All of the implemented FFT transform the sequence in place.
Main functionalities:

- `MFFT::Transform` computes the 1D FFT via the $N^{th}$-dimension DFT,
- `FFT::FFTParallel` computes the iterative radix-2 algorithm in parallel.

Extras:

- `FFT::DFT` performs the DFT by defintion in $O(N^2)$,
- `FFT::FFTRec`: recursive radix-2,
- `FFT::FFTSeq`: iterative radix-2.

## Description

The DFT of a sequence $x$ can be defined as
$$X(k) = \sum_{n = 0}^{N-1} x(n) \omega_N^{kn}$$  
with $\omega_N = e^{-i2\pi/N}$ the root of unity. The inversion involves reapplying the DFT with $-\omega_N$ and possibly multiplication with a scalar ($\frac{1}{N}$).

A very common FFT is the Cooley–Tukey radix-2 algorithm. The idea is to divide the sequence $x$ into two ($x_{\text{even}}$, $x_{\text{odd}}$), performs the DFT recursively on each, and construct the transformed sequence $X$ from $X_{\text{even}}$ and $X_{\text{odd}}$. The algorithm can be implemented iteratively for better performance.

`FFT::FFTParallel` implements the iterative radix-2 algorithm with parallelized elements. The sequence has to be padded (with zeroes in our case) such that the size $N$ is a power of 2. In particular,

- `RevBitParallel` rearrange the array based on the binary representation of each index: $x_{i_1, ..., i_n} = x^{\text{old}}_{i_n,...,i_1}$.
- For each of the $\log_2(N)$ steps, `ComputeThread` compute the current FFT based on the result of the previous steps. Each partition $[x_{i}, x_{i+2^{S}}]$ can be processed individually on different threads.

### MFFT (multi-dimensional)

`MFFT::Transform` implements the parallelized FFT described in [(1)](<https://doi.org/10.1016/0167-8191(90)90031-4>) and [(2)](https://doi.org/10.1109/SUPERC.1994.344263). The 1-D sequence $x$ is mapped onto P dimensions $N = N_1 \times ... \times N_p$. To compute, for instance, the DFT of a 2D $N \times M$ array, we can compute $N$ DFTs along one dimension, then $M$ DFTs of the previous result along the other. The array does not have to be padded to have a size of $2^k$, but each $N_i$ should (and in our case must) be divisible by the number of threads. It is possible in our code to, for instance, set `num_threads = 5` to process an array of size $5^k$ , but for simplicity, the default is to use a power of 2 with padding.

To compute the $i^{th}$-dimensional DFT, first the array $x$ has to be transposed so that the $i^{th}$ dimension is distributed along `x[]`. This is done by transposing $x$ after each DFT computation of slices of `x[]`. Twiddle factors are needed so that the result of the $i^{th}$-D DFT matches that of the ordinary DFT. For this, I implemented the $O(N^2)$ DFT to process segments of size $N_i$. Twiddle factors computation method was gradually improved during the project, but for now, I am not sure that this DFT can be replaced with the iterative 1D-FFT.

This function is the main focus of the project. While the final implementation is straightforward, several elements were hard to get right. To start with, I did not refer to the source code (if there is any) of both papers. Since the paper described the algorithm adapted to the architecture of the time, our version does not closely follow the description and is streamlined with `std::thread`. Nevertheless, mathematical operations (twiddle factors, transposing, etc.) are computed as described in both papers. I personally enjoyed writing the formula for calculating the indices of the transposition.

### Parallelization

FFT involves lots of array operations (assign, numerical operations) that fortunately are not in conflict with each other. In practice, the data and auxilliary variables can be structured to facilitate parallelization. I employed the very first technique learned in the course: splitting the array into blocks that are processed individually by threads.

The `Array` class store the processing array as well as supporting variables. This was realized mid-development of MFFT and was gradually adapted to the previously implemented algorithm to streamline the interface. Doing so slightly impacted the performance but greatly improves readability and practicality.

`Parallelize(std::function<> f, Array& arr, int inc)` splits the array into blocks while also considering the increment of the `for` loop in `f` which process `arr`.

`Parallelize2()` uses the same technique, but for a nested `for` loop of 2. The interface choose the loop with the largest "length" to split.

For simpler operations such as copying arrays or dividing each elements by N, I experimented with `<execution>`, particularly with the execution polities `par` and `par_unseq`. However, I encountered some memory issues with `std::copy` and have removed its usage.

## Testing

Testing environment: Personal laptop, AMD 4750u (8 cores, 16 threads), 32GB of RAM.

### Correctness

- The transformed data can vary between algorithms, either because of numerical error or of different padded sizes and thus are not automatically compared with each other, although we did manual testing to calibrate the algorithms. Instead, the FFT is followed with the Inverse FFT, and the recovered result is compared with the original data within a very small error margin.

- Tests are with varying sizes that are not necessarily powers of 2, and varying numbers of threads. Arrays are of moderate sizes so that the $O(N^2)$ DFT can also be tested. The varying sizes help catches bugs such as out-of-bounds array access or mis-initialization of `Array`.

### Speed

FFT algorithms in $O(N \log N)$ are fast and will process an array of size $10^6$ in less than a second (weather data won't be enough). Given the processing power, the maximum size tested was $5\times 10^7$.

Evidently, the iterative FFT is faster than the recursive one (and is easier to parallelize). `FFT::FFTParallel` is also faster with more threads (until 16). For instance, the runtime is reduced by half using 4 threads on large tests compared to the single-threaded one.

`MFFT::Transform` is slower than `FFT::FFTParallel` most of the time, partly because my implementation is sensitive to `num_threads` (each DFT runs at approximately $O(num thread)$). It is also more complex implementation-wise and is not exactly optimized to compete with `FFT::FFTParallel` (extremely slow with 1 thread). Still, the improved runtime when more threads (especially 8) are used is promising.

`FFT::FFTRec`:
|Size | 5,000|10,000|50,000|100,000|500,000|1e6|1e7|5e7|
|:---:|---|---|---|---|---|---|---|---|
||5,720|12,219|52,428|108,835|475,401|977,845|18,182,418|79,521,071|

`FFT::FFTSeq`:
|Size | 5,000|10,000|50,000|100,000|500,000|1e6|1e7|5e7|
|:---:|---|---|---|---|---|---|---|---|
|Size | 5,000|10,000|50,000|100,000|500,000|1e6|1e7|5e7|
||3,705|8,106|40,306|80,198|365,443|769,962|14,916,680|65,636,838|

`FFT::FFTParallel`:
|num_threads / N| 5,000|10,000|50,000|100,000|500,000|1e6|1e7|5e7|
|:---:|---|---|---|---|---|---|---|---|
|1 |7,194|12,507|50,097|101,384|430,135|898,810|15,862,235|64,536,779|
|2 |7,570|14,484|49,993|105,174|430,362|794,205|10,276,082|43,075,659|
|4 |7,196|10,899|36,706|71,918|301,769|618,313|7,802,009|32,140,696|
|8 |9,969|13,656|33,594|58,673|235,739|476,824|6,263,392|26,871,738|
|16|15,021|19,166|37,262|60,317|213,435|426,420|6,004,054|24,228,015|
|32|28,341|33,373|49,111|73,230|193,107|365,076|6,285,282|24,801,308|

`MFFT::Transform`:
|num_threads / N| 5,000|10,000|50,000|100,000|500,000|1e6|1e7|5e7|
|:---:|---|---|---|---|---|---|---|---|
|1 |36,532|77,379|329,839|700,964|3,036,914|6,251,118|108,702,607|466,899,596|
|2 |43,342|65,665|249,764|544,445|1,791,600|3,704,867|62,543,744|268,254,384|
|4 |14,239|23,474|112,399|236,746|834,821|1,469,157|22,148,206|92,937,463|
|8 |13,079|19,355|69,051|127,059|469,806|818,195|12,581,343|56,490,857|
|16|18,771|24,636|64,398|128,553|334,743|683,803|13,517,104|54,627,316|
|32|28,891|31,611|71,933|121,123|424,933|914,447|16,689,929|76,115,464|

## Application

I did not focus much on the applications of FFT. You can find in `main.cpp` a simple compression/extraction algorithm and the polynomimal multiplication algorithm that uses the Fourier transform.

## Misc

I created a functional makefile and managed to separate the headers and the source code!
