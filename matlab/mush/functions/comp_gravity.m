function gravity = comp_gravity(LLA)
% Computation of gravitational acceleration in Earth's frame of reference based on WGS84 Gravity model
% Precision of the computation is 10^-6 up to the 20000m of height above surface
% INPUT:          LLA = [latitude (deg); longitude (deg); height (m)];
% OUTPUT:         gravitational acceleration (m/s^2)
lat = deg2rad(LLA(1)); alt = LLA(3);
% WGS84 Gravity model constants:
a1 = 9.7803267715;
a2 = 0.0052790414;
a3 = 0.0000232718;
a4 = -0.0000030876910891;
a5 = 0.0000000043977311;
a6 = 0.0000000000007211;
% Gravity computation
gn = a1*(1+a2*sin(lat)^2+a3*sin(lat)^4)+(a4+a5*sin(lat)^2)*alt+a6*alt;
gravity = gn;

 