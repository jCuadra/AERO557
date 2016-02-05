function [ r, v ] = extendedGaussUV( RSite, v, qHat, q, tau1, tau3)
%Iterative method
%   Detailed explanation goes here

mu = 398600;
err = 10;
r1Old = [1 1 1]';
count_q = 0;
qOld = q;
while(err >1e-6)
    r = RSite + qHat*diag(q);
    r1 = r(:, 1);
    r2 = r(:, 2);
    r3 = r(:, 3);
    
    [ f1, g1] = keplerMY( r1, r2, tau1);
    [ f3, g3] = keplerMY( r2, r3, tau3);
    
    if count_q >1
        f1 = (f1+f1old)/2;
        f3 = (f3+f3old)/2;
        g1 = (g1+g1old)/2;
        g3 = (g3+g3old)/2;
    end

    v = gibbs(r1, r2, r3);
    
    C1 = g3/(f1*g3-f3*g1);
    C2 = -1;
    C3 = -g1/(f1*g3-f3*g1);
    
    q = inv(qHat)*RSite*[-C2; -C2; -C3]./[C1; C2; C3];
    
    err = abs(q(2) -qOld(2));
    qOld = q;

    f1old = f1;
    f3old = f3;
    g1old = g1;
    g3old = g3;
end

end


