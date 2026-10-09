#!/usr/bin/env python3
"""Least-squares half-band tap design for eon::Oversampler.

The polyphase kernel stores the even-indexed taps hE; the odd taps are zero
except the centre 0.5. The 2x-rate interpolator response is LINEAR in hE:

    H_up(W) = exp(-jW(p+1)) + 2 * sum_i hE[i] * exp(-jW(2i+1))

so a weighted LS on a dense frequency grid gives the kernel directly.
Target: pass 1 below wp (2x-rate), stop 0 above ws, don't-care between.
"""
import numpy as np

# Current embedded Kaiser taps (from Dsp/Oversampling.h), kept for comparison.
KAISER = {
 71: (17, [-1.2477987897e-06,1.1227480081e-05,-4.2213195345e-05,1.1705190499e-04,
 -2.7257107958e-04,5.6356400147e-04,-1.0666100923e-03,1.8837546556e-03,
 -3.1467476030e-03,5.0238835479e-03,-7.7341787388e-03,1.1579539945e-02,
 -1.7020470787e-02,2.4864705852e-02,-3.6794020275e-02,5.7171141772e-02,
 -1.0208649494e-01,3.1694998446e-01,3.1694998446e-01,-1.0208649494e-01,
 5.7171141772e-02,-3.6794020275e-02,2.4864705852e-02,-1.7020470787e-02,
 1.1579539945e-02,-7.7341787388e-03,5.0238835479e-03,-3.1467476030e-03,
 1.8837546556e-03,-1.0666100923e-03,5.6356400147e-04,-2.7257107958e-04,
 1.1705190499e-04,-4.2213195345e-05,1.1227480081e-05,-1.2477987897e-06]),
 47: (11, [-4.9151111759e-06,6.5122378319e-05,-2.8715247150e-04,8.7813478063e-04,
 -2.1786964615e-03,4.6960948225e-03,-9.1540676449e-03,1.6629580781e-02,
 -2.9012452245e-02,5.0748369748e-02,-9.7845787349e-02,3.1546742160e-01,
 3.1546742160e-01,-9.7845787349e-02,5.0748369748e-02,-2.9012452245e-02,
 1.6629580781e-02,-9.1540676449e-03,4.6960948225e-03,-2.1786964615e-03,
 8.7813478063e-04,-2.8715247150e-04,6.5122378319e-05,-4.9151111759e-06]),
 35: (8, [1.7121717578e-05,-2.6899266591e-04,1.2540227767e-03,-3.9265122295e-03,
 9.8371917699e-03,-2.1471419812e-02,4.3766261984e-02,-9.2881796171e-02,
 3.1366794053e-01,3.1366794053e-01,-9.2881796171e-02,4.3766261984e-02,
 -2.1471419812e-02,9.8371917699e-03,-3.9265122295e-03,1.2540227767e-03,
 -2.6899266591e-04,1.7121717578e-05]),
}

def kernel(hE, N):
    """hE -> full length-N kernel: even indices = hE, centre = 0.5."""
    h = np.zeros(N)
    c = (N - 1) // 2
    for i, v in enumerate(hE):
        h[2 * i] = v
    h[c] = 0.5
    return h

def resp(h, th):
    return np.exp(-1j * np.outer(th, np.arange(len(h)))) @ h

def evaluate(hE, N, wp=0.43 * np.pi, ws=np.pi - 0.43 * np.pi, npts=8192):
    th = np.linspace(0, np.pi, npts)
    K = np.abs(resp(kernel(hE, N), th))
    pass_err_db = 20 * np.log10(np.abs(K[th <= wp] - 1.0).max() + 1e-300)
    stop_db = 20 * np.log10(K[th >= ws].max() + 1e-300)
    return pass_err_db, stop_db

def design_ls(N, wp=0.43 * np.pi, ws=np.pi - 0.43 * np.pi,
              weight_stop=1.0):
    """Weighted LS half-band design via scipy's firls (analytic Toeplitz solve).

    firls returns the full kernel; symmetric bands about pi/2 make the
    odd-offset taps ~0 by construction. We snap them to zero, renormalise
    the centre tap to 0.5, and return the even-indexed taps hE.
    """
    from scipy.signal import firls
    h = firls(N, [0, wp / np.pi, ws / np.pi, 1.0],
              [1, 1, 0, 0], weight=[1.0, weight_stop], fs=2.0)
    # For N = 4m+3 the centre index (2m+1) is ODD, so it lives in h[1::2]
    # together with the even-offset taps that a half-band design leaves ~0.
    h[1::2] = 0.0                    # snap the near-zero even-offset taps
    h[(N - 1) // 2] = 0.5            # restore the centre tap exactly
    return h[0::2]

if __name__ == "__main__":
    for N, (p, kaiser) in KAISER.items():
        L = len(kaiser)
        kp, ks = evaluate(kaiser, N)
        best = None
        for w in (1.0, 4.0, 16.0, 64.0, 256.0):
            g = design_ls(N, weight_stop=w)
            lp, ls_ = evaluate(g, N)
            best = (w, g, lp, ls_)
            if lp > kp - 0.5: break
        w, g, lp, ls_ = best
        print(f"N={N} p={p}  Kaiser: pass {kp:7.1f} dB  stop {ks:7.1f} dB")
        print(f"N={N} LS(w={w}): pass {lp:7.1f} dB  stop {ls_:7.1f} dB")
        print("   LS hE =", ", ".join(f"{v:.10e}" for v in g))
