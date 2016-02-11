function [ DCM ] = R1d( x )
%This will perform an R1 rotation through the angle x (radians)
%   Detailed explanation goes here
DCM = [1 0 0; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];

end

