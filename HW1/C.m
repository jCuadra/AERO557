function [c ] = C(z )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
    if z>0
        c = (1-cos(z^.5))/z;
    elseif z<0
        c = (cosh(abs(z)^.5)-1)/abs(z);
    else
        c = .5;
    end
    
end

