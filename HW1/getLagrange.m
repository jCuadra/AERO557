function [ f, g ] = getLagrange( r2, tau )
%getLagrange gets the lagrange multipliers necessary to fine r1 and r3 in
%Gauss's OD method
%   @param r2   Magnitude of the r vector from the second observation, km
%   @param tau  Time between the ith and second observation, seconds
%   @return f   f lagrange multiplier for ith observation
%   @return g   g lagrange multiplier for ith observation
mu = 398600;
u = mu/r2^3;
f = 1-u/2*tau^2;
g = tau-u/6*tau^3;


end

