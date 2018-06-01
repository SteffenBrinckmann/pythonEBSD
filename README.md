Read and plot electron backscatter diffraction (EBSD) data

I look for more contributors, see below

# Features:
  - File formats accepted .ang | .osc | .crc | .txt
  - can write .ang for FCC. Others could be added
  - fast plotting interaction using virtual mask (only used for plotting);
    - even though some interpolation and reading takes time
    - can be removed just before final plotting
  - verified with the OIM software and mTex
  - heavily tested for cubic
  - some educational plotting
  - examples and lots of docummentation via Doxygen
  - doctest verifies functionality


EBSD-Inverse Pole Figure (IPF) of polycrystalline Copper
![EBSD of polycrystalline Copper](HTMLInput/ebsd_py_ND.png)

corresponding pole-figure
![Pole figure](HTMLInput/ebsd_py_PF100.png)


# What features I do not envision:
  - include all crystal symmetries (materials science can mostly live with few)
  - other Euler angle definitions than Bunge; materials science does not use those


# Future features
  - introduce Kernel average misorientation
  - improve cleaning, grain identification methods

# Help wanted
 - sample files other than copper OIM files
 - feedback on tutorials
 - any feedback on functionality
 - help with cleaning and grain identification


# Build docummentation
```bash
python verifyAll.py
```
