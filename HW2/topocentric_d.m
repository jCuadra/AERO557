function [ RA_DEC] = topocentric_d( r, RSite )
%Convert r and v vectors from observation site to look angles, RA and dec
%(alpha and delta). Implemented from Algorithm 26 in Vallado
%   @param r        ECI radius vector of satellite
%   @param v        ECI velocit vector of satellite
%   @param RSite    Coordinate of observation vector
%   @return RA_DEC  Column vector of RA and DEC look angles of satellite, degrees

rho = r-RSite;
rhoNorm = norm(rho);
RA_DEC(2, 1) = asind(rho(3)/rhoNorm);
rhoIJ_norm = sqrt(sum(rho(1:2).^2));
if(rhoIJ_norm~= 0)
    RA_DEC(1,1) = asind(rho(2)/rhoIJ_norm);
else
    disp('Now you have to implement the rest of this algorithm becuase its')
    disp('a hyperbolic orbit :(')
    RA_DEC = -1;
end

end

