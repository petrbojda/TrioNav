function [ geod ] = ned2geod( P_NED,initial_geod )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Geodeticke konstanty
a = 6378137.0;                  % Polomer zemske semi-major osy (m)
b = 6356752.3142;               % Polomer zemske semi-minor osy (m)
e = sqrt(1 - b^2/a^2);          % Excentricita Zeme

iniLAT =initial_geod(1);
iniLON = initial_geod(2);
iniALT = initial_geod(3);

R_N = a*(1 - e^2)/(1 - e^2*(sin(deg2rad(iniLAT))^2))^(3/2);  % Radius of curvature in the prime vertical
R_M = a/(1 - e^2*(sin(deg2rad(iniLAT))^2))^(1/2);            % Meridian radius of curvature

lat = P_NED(1)*rad2deg(atan(1/R_M))+iniLAT;
long = P_NED(2)*rad2deg(atan(1/(R_N*cos(deg2rad(iniLAT)))))+iniLON;
alt = iniALT-P_NED(3);
geod =[lat;long;alt];
end

