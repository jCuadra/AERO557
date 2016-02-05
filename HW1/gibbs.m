function [ v2 ] = gibbs( r1, r2, r3 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
mu = 398600;
Z12 = cross(r1, r2);
Z23 = cross(r2, r3);
Z31 = cross(r3, r1);

% alpha_cop = asin(dot(Z23, r1)/norm(Z23)/norm(r1));
% c_alpha12 = dot(r1, r2)/norm(r1)/norm(r2)*180/pi;
% c_alpha23 = dot(r2, r3)/norm(r2)/norm(r3)*180/pi;
N = norm(r1)*Z23 + norm(r2)*Z31 + norm(r3)*Z12;
D = Z12 + Z23 + Z31;
S = (norm(r2)-norm(r3))*r1 + (norm(r3)-norm(r1))*r2 +(norm(r1)-norm(r2))*r3;
B = cross(D, r2);
L_g = sqrt(mu/norm(N)/norm(D));
v2 = L_g/norm(r2)*B+L_g*S;

end

