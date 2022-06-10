# FFFT

## Build
`make`

## Rewritten notation

$$X(k) = \sum_{n = 0}^{N-1} x(n) \omega_N^{kn}$$  
with $\omega_N = e^{-i2\pi/N}$.

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

