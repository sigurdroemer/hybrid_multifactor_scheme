# Hybrid multifactor scheme for stochastic Volterra equations with completely monotone kernels
We implement the hybrid multifactor scheme for the simulation of stochastic Volterra equations with completely monotone kernels.

![equation](https://latex.codecogs.com/svg.image?X_t&space;=&space;g_0(t)&space;&plus;&space;\int_0^t&space;K(t-s)b(s,X_s)ds&space;&plus;&space;\int_0^t&space;K(t-s)\sigma(s,X_s)dW_s,&space;\phantom{xxx}&space;t&space;\geq&space;0,)



To get started, see the example scripts in the folder '.../get_started'.

If you decide to use this code in your research please cite the repository and the paper https://www.ssrn.com/abstract=3706253.

# Reference
- Rømer, S.E., Hybrid multifactor scheme for stochastic Volterra equations with completely monotone kernels (2022). Working paper, available at SSRN:3706253.

# External packages
The following external packages are used in the project:
- Adi Navve (2020). Pack & Unpack variables to & from structures with enhanced functionality (https://www.mathworks.com/matlabcentral/fileexchange/31532-pack-unpack-variables-to-from-structures-with-enhanced-functionality), MATLAB Central File Exchange. Retrieved March 16, 2020.
- Sheung Hun Cheng and Nicholas J. Higham (2015). Modified-cholesky downloaded from https://github.com/higham/modified-cholesky on May 17, 2020.
