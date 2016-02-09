function [ A ] = orb_prop( tspan, state )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mu = 398600;
v = state(4:6,1);
r = state(1:3, 1);
accel =-mu*r/norm(r)^3;
A = [v; accel];
end

