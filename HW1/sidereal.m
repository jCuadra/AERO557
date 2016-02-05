function [ J ] = sidereal( y, m, d )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

J = 367*y-floor((7*(y+floor((m+9)/12)))/4) + floor((275*m)/9)+d+1721013.5;


end

