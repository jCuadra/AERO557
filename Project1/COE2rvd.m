function [ r, v ] = COE2rvd( p, ecc, inc, RAAN, w, theta)
%COE2RV converts the classical orbital elements into radius and velocity
%vectors of the orbit
%   @param p        semi-parameter of the orbit, km
%   @param ecc      eccentricity of the orbit
%   @param inc      inclucation of the orbit, deg
%   @param RAAN     right ascention of the ascending node, deg
%   @param w        argument of perigee, deg
%   @param theta    true anomoly, deg
%   @return r       radius vector of the orbit, km
%   @return v       velocity vector of the orbit, km/s
mu = 398600;
if ecc < .0000001 && inc <.0000001 %circular equitorial
    w = 0;
    RAAN = 0;
elseif ecc <.0000001 %circular
    w = 0;
elseif inc<.0000001 %elliptical equitorial
    RAAN = 0;
end

rPQW = [p*cosd(theta)/(1+ecc*cosd(theta)); 
        p*sind(theta)/(1+ecc*cosd(theta));
        0];

vPQW = [-sqrt(mu/p)*sind(theta);
        sqrt(mu/p)*(ecc+cosd(theta));
        0];

r = R3d(-RAAN)*R1d(-inc)*R3d(-w)*rPQW;
v = R3d(-RAAN)*R1d(-inc)*R3d(-w)*vPQW;

end

