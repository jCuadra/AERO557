function [ p, a, ecc, inc, Omega, w, theta ] = RV2COEd( r, v)
%Convert R and V vectors into the classical orbital elements in degrees.
%Created from algorithm 9 in Vallado.
%   @param r        radius vector
%   @param v        velocity vector
%   @return p       semiparameter
%   @return a       semimajor axis
%   @return ecc     eccentricity
%   @return inc     inclination
%   @return Omega   right ascention of ascending node
%   @return w       argument of perigee
%   @return theta   true anomoly
mu = 398600; %gravitational parameter of the earth
h = cross(r, v); %
n = cross([0 0 1], h);
ecc = ((norm(v)^2-mu/norm(r))*r - dot(r, v)*v)/mu;
chi = norm(v)^2/2-mu/norm(r);
if norm(ecc) ~= 1
    a = -mu/2/chi;
    p = a*(1-norm(ecc)^2);
else
    p = norm(h)^2/mu;
    a = inf;
end

inc = acosd(h(3)/norm(h));

Omega = acosd(n(1)/norm(n));
if n(2) <0
    Omega = 360-Omega;
end

w = acosd(dot(n, ecc)/norm(n)/norm(ecc));
if ecc(3) < 0
    w = 360-w;
end

theta = acosd(dot(ecc, r)/norm(ecc)/norm(r));
if dot(r, v)<0
    theta = 360-theta;
end

end

