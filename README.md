# hourly-drifters
Matlab code/scripts to linearly interpolate the GDP dataset as described in Elipot et al, http://dx.doi.org/10.1002/2016JC011716

These scripts were last written in Matlab R2016a and may require routines from the jLab Matlab toolbox available at this repository:
https://github.com/selipot/jLab

global_interpolation_GPS_1.01.m : code to apply the LOWESS plus linear interpolation for GPS-tracked drifters
  dependencies: argos_reduce.m, LatLonLikelihoodHessian.m
  
global_interpolation_wmle_errors_1.01.m : cote to apply the weighted maximum likelihood estimation for Argos-tracked drifters
  dependencies: gps_reduce.m, LatLonLocalWess.m, piecelinvar.m
