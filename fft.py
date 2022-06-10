import numpy as np

with open('fft_py.txt', 'w') as out:
    x = np.loadtxt('data/temp.txt')
    X = np.fft.fft(x)
    # X = np.fft.ifft(X)
    for c in X:
        out.write('{} {}\n'.format(c.real, c.imag))
