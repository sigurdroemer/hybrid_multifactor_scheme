# Hybrid multifactor scheme for stochastic Volterra equations with completely monotone kernels
We implement the hybrid multifactor scheme of (Rømer, 2022) for equations of the form

![equation](https://latex.codecogs.com/svg.image?X_t%26space%3B%3D%26space%3Bg_0%28t%29%26space%3B%26plus%3B%26space%3B%5Cint_0%5Et%26space%3BK%28t-s%29b%28s%2CX_s%29ds%26space%3B%26plus%3B%26space%3B%5Cint_0%5Et%26space%3BK%28t-s%29%5Csigma%28s%2CX_s%29dW_s%2C%26space%3B%5Cphantom%7Bxxx%7D%26space%3Bt%26space%3B%5Cgeq%26space%3B0%2C)

where ![equation](https://latex.codecogs.com/svg.image?g_0%3A%5Cmathbb%7BR%7D_%26plus%3B%26space%3B%5Crightarrow%26space%3B%5Cmathbb%7BR%7D)

To get started, see the example scripts in the folder '.../get_started'.

If you decide to use this code in your research please cite the repository and the paper https://www.ssrn.com/abstract=3706253.

# Reference
- Rømer, S.E., Hybrid multifactor scheme for stochastic Volterra equations with completely monotone kernels (2022). Working paper, available at SSRN:3706253.

# External packages
The following external packages are used in the project:
- Adi Navve (2020). Pack & Unpack variables to & from structures with enhanced functionality (https://www.mathworks.com/matlabcentral/fileexchange/31532-pack-unpack-variables-to-from-structures-with-enhanced-functionality), MATLAB Central File Exchange. Retrieved March 16, 2020.
- Sheung Hun Cheng and Nicholas J. Higham (2015). Modified-cholesky downloaded from https://github.com/higham/modified-cholesky on May 17, 2020.
