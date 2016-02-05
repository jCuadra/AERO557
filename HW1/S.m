function [ s ] = S( z )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
    if z>0
        s= (z^.5 - sin(z^.5))/(z^(3/2));
    elseif z< 0
            s = ((sinh(abs(z)^.5)-abs(z)^.5))/abs(z)^(3/2);
    else
        s = 1/6;
    end 

end

