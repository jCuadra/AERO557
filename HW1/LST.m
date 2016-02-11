function [ Theta, UT ] = LST( y, m, d, hr, min, sec, long, longmin )
%LST Finds theta LST for a location at a time stamp in UT time
%   @param y    year
%   @param m    month of the year
%   @param d    day of the month
%   @param hr   hour of the day
%   @param min  minute of the hour
%   @param sec  second of the minute
%   @param lat  latitude, deg
%   @param latmin   minute unit of the latitude
%   @return Theta   The latitude in local sidereal time, degrees
J = sidereal(y, m, d);
T = (J-2451545)/36525;
gw0 = 100.4606184 + 36000.77004*T+0.000387933*T^2-2.58*10^-8 *T^3;

while gw0>360
    gw0 = gw0-360;
end

while gw0<0
    gw0 = gw0+360;
end

UT = hr + min/60 + sec/3600;

gw = gw0 + 360.98564724*(UT/24);

Theta = gw + long + longmin/60;

while Theta>360
    Theta = Theta-360;
end

while Theta<0
    Theta = Theta + 360;
end

end

