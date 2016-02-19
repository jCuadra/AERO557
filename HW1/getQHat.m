function [ qHat ] = getQHat( alpha, delta )
%getQHat get the unit vector of the slant range from RA and DEC
%   @param alpha    RA of each observation 1x3, degrees
%   @param delta    DEC of each observation 1x3, degrees
%   @return qHat    Unit vector of slant range for each of the three
%                   observations, one column vector per observation,
%                   3xNUM_OBSV
for i = 1:3
    qHat(1, i) = cosd(delta(i))*cosd(alpha(i));
    qHat(2, i) = cosd(delta(i))*sind(alpha(i));
    qHat(3, i) = sind(delta(i));
end

end

