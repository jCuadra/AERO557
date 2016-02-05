function [ a, au ] = getAs( tau, taui)
%getAs Calculate coefficients for obtaining scalar of splant range
%   @param tau      Total time between first and 3rd observation, sec
%   @param taui     Time between ith and second observation, sec
%   @param r2       Radius of the second observation vector, symbolic var
%   @return a       Tau coefficient
%   @return au      Tau and u coefficient
mu = 398600;
a = taui/tau;
au = taui/tau*(tau^2-taui^2)/6;
if taui<0
    a = -a;
    au = -au;
end

end

