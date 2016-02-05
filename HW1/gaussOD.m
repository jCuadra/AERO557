function [ r, v, RSite, qHat, q, tau1, tau3 ] = gaussOD(alpha, delta, siteInfo)
%gaussOD Gauss method for orbit determination for an oject whose RAAN and 
%Decclination were observed at 3 different time steps within 10 degrees of 
%each other. ASSUMPTIONS: Coplanar observations (no pertibations between 
%observations)
%   @param alpha     Array of the three right ascention of ascending node 
%                    from each observation, degree (1x3)
%   @param delta     Array of the three declinations from each observation
%   @param siteInfo should containeach of the following , degree (1x3)
%                   fields (one element of the array per observation):
%       H            Altitude of the observation location, km
%       latSite      Latitude of the observation location, degrees
%       longSite     Longitude of the observation location, degrees
%       longmin      Minute unit of the longitude of the observation
%                    location, minute
%       year         Year of the observation (1x3)
%       month        Month of the year of the ovservation (1x3)
%       day          Day of the month of the observation (1x3)
%       hour         Hour of the day of the observation (1x3)
%       min          Minute of the hour of the observation (1x3)
%       sec          Second of the minute of the observation (1x3)

global NUM_OBSV
mu = 398600;
r = 0; v = 0;
% Get RSite vectors
RSite = getRSite(siteInfo);
%RSite = [4054.881 3956.224 3905.073; 2748.195 2888.232 2956.935; 4074.237 4074.364 4074.430]
% Get slant range unit vector
qHat = getQHat(alpha, delta); 
%qHat = [.9473 .5742 .3007; .0155 .5747 .7399; .3201 .5830 .6018];
% get r2
tau = getTau(siteInfo, 3, 1);
tau3 = getTau(siteInfo, 3, 2);
tau1 = getTau(siteInfo, 1, 2);

 r2 = sym('r2');
[a1, a1u] = getAs(tau, tau3);
[a3, a3u] = getAs(tau, tau1);
M = inv(qHat)*RSite;
d1 = M(2, 1)*a1 - M(2, 2) +M(2, 3)*a3;
d2 = M(2, 1)*a1u +M(2,3)*a3u;
c = dot(qHat(:, 2), RSite(:, 2));
rSite2 = norm(RSite(:, 2));
eqn = r2^8-(d1^2 +2*c*d1+rSite2^2)*r2^6-2*mu*(c*d2+d1*d2)*r2^3-mu^2*d2^2;
ans = double(solve(eqn==0, r2));
realR2 = ans(find(ans(find(imag(ans)==0))>0,1));

%Use r2 to get r1 and r3 using lagrange multipliers
[f1, g1] = getLagrange(realR2, tau1);
[f3, g3] = getLagrange(realR2, tau3);
u = mu/(realR2^3);
c1 = a1+a1u*u;
c2 = -1;
c3 = a3+a3u*u;
q = inv(qHat)*RSite*[-c1; -c2; -c3]./[c1; c2; c3];
r = RSite + qHat*diag(q);
v = 1/(f1*g3-f3*g1)*(-f3*r(:,1)+f1*r(:,3));
end

