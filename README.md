Read and plot electron backscatter diffraction (EBSD) data

Features:
  - File formats accepted .ang | .osc | .crc | .txt
  - fast plotting interaction using virtual mask (only used for plotting);
    - even though some interpolation and reading takes time
    - can be removed just before final plotting
  - verified with the OIM software and mTex
  - heavily tested for cubic

EBSD-Inverse Pole Figure (IPF) of polycrystalline Copper
![EBSD of polycrystalline Copper](HTMLInput/ebsd_py_ND.png)

corresponding pole-figure
![Pole figure](HTMLInput/ebsd_py_PF100.png)


Future features
  - improve cleaning, grain identification methods
