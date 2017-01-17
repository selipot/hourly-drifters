# hourly-drifters
Matlab scripts to interpolate to hourly intervals the GDP dataset as described in  Elipot S., R. Lumpkin, R. C. Perez, J. M. Lilly, J. J. Early, A. M. Sykulski (2016), A global surface drifter dataset at hourly resolution, J. Geophys. Res.-Oceans, 121, doi:10.1002/2016JC011716.

These scripts were written in Matlab R2016a and require some routines from the jLab Matlab toolbox available at this repository:
https://github.com/jonathanlilly/jLab

global_interpolation_GPS_1.01.m : code to apply the LOWESS plus linear interpolation for GPS-tracked drifters. This script calls argos_reduce.m and LatLonLikelihoodHessian.m.
  
global_interpolation_wmle_errors_1.01.m : cote to apply the weighted maximum likelihood estimation for Argos-tracked drifters. This script calls gps_reduce.m, LatLonLocalWess.m, and piecelinvar.m
