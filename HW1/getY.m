function [ y ] =getY( r_0, r, A, psi, c2, c3 )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

y = norm(r_0)+norm(r)+(A*((psi*c3)-1)/sqrt(c2));
end

