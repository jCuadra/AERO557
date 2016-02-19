function [ Az, El ] = getAzEl( rho, phi, long, lst)
%Approximation of Azimuth and Elevation from the slant range, rho, in ECI.
%Would be more accurate if slant range was in ECEF. From Vallado pg 263
%   @param rho  Slant range, km, ECI
%   @return Az  Azimuth, degs from observation site
%   @return El  Elevation, deg from observation site
rad = pi/180;
phi = phi*rad;
long = long*rad;
rot = [sin(phi)*cos(long) sin(phi)*sin(long) -cos(phi);
      -sin(long) cos(long) 0;
      cos(phi)*cos(long) cos(phi)*cos(long) sin(phi)];
  rho = rot*rho;
El = asind(rho(3)/norm(rho))-phi/rad;
if (El ~= 90)
    Azs = rho(2)/sqrt(sum(rho(1:2).^2));
    Azc = rho(1)/sqrt(sum(rho(1:2).^2));
    Az = atand(Azs/Azc);
    Az = lst-Az;
    if Az >360
        Az = Az-360;
    end
else
   disp('Somehow your elevation is exactly 90 deg and we dont wanna code that') 
end

end

