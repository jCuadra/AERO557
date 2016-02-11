function [ Az, El ] = getAzEl( rho )
%Approximation of Azimuth and Elevation from the slant range, rho, in ECI.
%Would be more accurate if slant range was in ECEF. From Vallado pg 263
%   @param rho  Slant range, km, ECI
%   @return Az  Azimuth, degs from observation site
%   @return El  Elevation, deg from observation site

El = asind(rho(2)/norm(rho));
if (El ~= 90)
    Az = asind(rho(1)/sqrt(sum(rho(1:2).^2)));
else
   disp('Somehow your elevation is exactly 90 deg and we dont wanna code that') 
end

end

