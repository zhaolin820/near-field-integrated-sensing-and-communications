# Near-Field Integrated Sensing and Communications

The code for the papers 

**Z. Wang, X. Mu, and Y. Liu, “Near-Field Integrated Sensing and Communications,” *IEEE Commun. Lett.*, vol. 27, no. 8, pp. 2048-2052, Aug. 2023** [[IEEE](https://ieeexplore.ieee.org/abstract/document/10135096)] [[Arxiv](https://arxiv.org/abs/2302.01153)]

Note: Although the trend observed in Fig. 4(b) holds for sensing-only systems, it may not be a general trend in near-field ISAC systems due to the randomness of beamforming caused by the communication functions. As a result, it is expected that different trends may appear in simulations. Thanks to Panpan Qian (qianpanpan@nuaa.edu.cn) from NUAA for pointing out this issue.

## Running the simulations

### Prerequisites

- [MATLAB](https://uk.mathworks.com/products/matlab.html)
- [CVX](http://cvxr.com/cvx/)

### Launch

Run `main.m`

### Expected Results

#### Spectrum of MUSIC
<img decoding="async" src="./results/MUSIC_spectrum.jpg" width="50%">

## Citing
If you in any way use this code for research, please cite our original articles listed above. The corresponding BiBTeX citation is given below:
```
@article{wang2023near,
  title={Near-field integrated sensing and communications},
  author={Wang, Zhaolin and Mu, Xidong and Liu, Yuanwei},
  journal={IEEE Wireless Commun. Lett.},
  year={2023},
  month=aug,
  volume={27},
  number={8},
  pages={2048-2052}
}
```
