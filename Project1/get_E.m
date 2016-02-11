function [ E ] = get_E( M, ecc )
%GET_E gets the eccentric anomoly, E, from the Mean Aonomoly M and the ecc
%   @param M    Mean Anomoly, deg
%   @param ecc  eccentricity
%   @return E   Eccentric Anomoly, deg
tol = 10^-8;
ratio =1;
M = M*pi/180; %use radians for calcs
E = M-ecc;
while(abs(ratio) >tol)
    ratio = ((-E+M+ecc*sin(E))/(-1+ecc*cos(E)));
    E = E - ratio;
end
E = E*180/pi; %return degrees
end

