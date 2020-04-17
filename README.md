DSurfTomo is a surface wave inversion program which can irectly invert surface wave dispersion data directly to 3D shear wavespeed models, without the intermediate step of constructing the phase or group velocity maps.
The fast marching method (FMM) (Rawlinson et al., 2004) is used to compute, at each period, surface wave travel times and ray paths between sources and receivers. This avoids the assumption of great-circle propagation that is used in most surface wave tomographic studies, but which is not appropriate in complex media. 
Please refer to Fang et al. (2015 , GJI) for the detail description of the method.

Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. Geophysical Journal International, 201(3), 1251-1263. 

Rawlinson, N. & Sambridge, M., 2004. Wave front evolution in strongly heterogeneous layered media using the fast marching method, Geophys. J. Int., 156(3), 631–647 

For V2.0, I incorporated the random projections based inversion using Poisson Voronoi cells. However, it has not been fully tested, so any feedbacks will be helpful. For more information about the inversion method, please refer to:

Fang, H., van der Hilst, R. D., de Hoop, M. V., Kothari, K., Gupta, S., & Dokmanić, I. (2020). Parsimonious Seismic Tomography with Poisson Voronoi Projections: Methodology and Validation. Seismological Research Letters, 91(1), 343-355.
