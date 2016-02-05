function [ c2, c3 ] = lambertC( psi)
%Get c2 and c3 for use in lambert's universal variable solver. Vallado
%algorithm 1. doesn't work if psi = 0
%   @param psi  parameter in Lambert's UV X^2/a
%   @return c2
%   @return c3

if psi >1e-6
     c2 = (1-cos(sqrt(psi)))/psi;
     c3 = (sqrt(psi)-sin(sqrt(psi)))/(psi^(3/2));
else
    if psi < -1e-6
        c2 = (1-cosh(sqrt(-psi)))/psi;
        c3 = (sinh(sqrt(-psi))-sqrt(-psi))/((-psi)^(3/2));
    else
        c2 = 1/2;
        c3 = 1/6;
    end
end

end

