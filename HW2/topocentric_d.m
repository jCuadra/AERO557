function [ RA_DEC] = topocentric_d( r, RSite, len )
%Convert r and v vectors from observation site to look angles, RA and dec
%(alpha and delta). Implemented from Algorithm 26 in Vallado
%   @param r        ECI radius vector of satellite
%   @param v        ECI velocit vector of satellite
%   @param RSite    Coordinate of observation vector
%   @return RA_DEC  Column vector of RA and DEC look angles of satellite, degrees
deg = 180/pi;
for i = 1:len
    rho = r(:, i)-RSite(:, i);
    rhoNorm = norm(rho);
    dec = asin(rho(3)/rhoNorm)*deg;
    if dec<0
        RA_DEC(2, i) = dec+360;
    else
        RA_DEC(2, i) = dec;
    end
    rhoIJ_norm = norm(rho(1:2));
    if(rhoIJ_norm~= 0)
        raS = rho(2)/rhoIJ_norm;
        raC = rho(1)/rhoIJ_norm;
        ra = atan(raS/raC)*deg;
        if ra <0
            RA_DEC(1,i) = 360+ra;
        else
            RA_DEC(1,i) = ra;
        end
    else
        disp('Now you have to implement the rest of this algorithm becuase its')
        disp('a hyperbolic orbit :(')
        RA_DEC = -1;
    end
end
end

