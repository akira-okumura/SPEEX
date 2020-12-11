# SPEEX2
Single Photo-Electron EXtractor 2

SPEEX2 is a small utility projec to extract the single photoelectron (p.e.) spectrum from a multi-p.e. charge distribution taken by PMT or SiPM. The basic idea is based on the iteration method proposed by Takahashi et al. 2018 \[1\], but we also plan to add an improved method that takes into account the pedestal distribution more accurately.

As of Dec 2020, the program is written in C++ and ROOT, and only the original iteration method, which is not applicable when the pedestal is very dirty, is implemented.

\[1\] M. Takahashi, et al. (2018) NIMA **894** 1–7, "A technique for estimating the absolute gain of a photomultiplier tube" https://dx.doi.org/10.1016/j.nima.2018.03.034

# Usage

```
$ root
root [0] .L speex2.C+
root [1] auto anas = test()
```
Please read the `test` function to see how to use the `SinglePEAnalyzer` class.
