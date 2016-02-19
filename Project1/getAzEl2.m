function [ solAz, solEl ] = getAzEl2(LST, RA, DEC, phi)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%LST, RA, DEC, phi
rad = pi/180;
LST = LST*rad;
RA = RA*rad;
DEC = DEC*rad;
phi = phi*rad;
LHA = LST-RA;
rad = pi/180;
syms Az El
eqn1 = sin(DEC)-sin(El)*sin(phi)- cos(El)*cos(phi)*cos(Az) == 0;
eqn2 = sin(LHA)+sin(Az)*cos(El)/cos(DEC) == 0;
%eqn3 = cosd(LHA)-(cosd(phi)*sin(El)-sind(phi)*cos(Az)*cos(El))/cos(DEC) == 0;
[solAz, solEl] = solve([eqn1, eqn2], [Az El]);
solAz = double(solAz);
solEl = double(solEl);
solAz = (solAz)/rad;
solEl = (solEl)/rad;
% if solAz <0
%     solAz = 360 +solAz;
% end
% if solEl <0
%     solEl = 180 +solEl;
% end
% if solEl >180
%     solEl = solEl-180;
% end
end

