function [ RSite, LSTkeep ] = getRSite( siteInfo, begin )
%getRSite Calculate the location, in ECI (km), of the observatory location 
%   @param phi  the latitude of the observation location, degrees
%   @param H    altitude of the observation location, km
%   @param theta    Longitude, degrees
%   @return RSite   R vector of the observation location in ECI, km
Re = 6378.137; %radius of the earth
f = .003353; %oblateness factor of the earth
rad = pi/180;
phi = siteInfo.latSite;
H = siteInfo.H;
rad = pi/180;
for i = begin:begin+2
    thetaLST = LST(siteInfo.year(i), siteInfo.month(i), siteInfo.day(i),...
                    siteInfo.hour(i), siteInfo.min(i), siteInfo.sec(i),...
                    siteInfo.longSite, siteInfo.longmin);
    thetaLST = thetaLST*rad;
    A = (Re/sqrt(1-(2*f-f^2)*sin(phi*rad)^2) + H)*cos(phi*rad);
    RSite(1, i) = A*cos(thetaLST);
    RSite(2, i) = A*sin(thetaLST);
    RSite(3, i) = (Re/sqrt(1-(2*f-f^2)*sin(phi*rad)^2)*(1-f)^2 +H)*sin(phi*rad);
    LSTkeep(i) = thetaLST/rad;
end
end

