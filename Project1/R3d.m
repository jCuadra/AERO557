function [ DCM ] = R3d( x)
%This will perform an R3 rotation through the angle x (radians)
%   Detailed explanation goes here
DCM = [cosd(x) sind(x) 0; -sind(x) cosd(x) 0; 0 0 1];

end

