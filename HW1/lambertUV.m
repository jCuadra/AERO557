function [ f, g, gdot] = lambertUV( r1, r2, tau, tm )
%Lambert's Universal variable function to find the lagrange multipliers in
%order to find velpcity vectors
%   @param r1       Initial r vector, 3x1 km
%   @param r2       Final r vector, 3x1 km
%   @param tau       Time between r vectors s
%   @param tm       1 or -1 for short way or long way around
%   @return f       Lagrange multiplier
%   @return g       Lagrange multiplier
%   @return g_dot   Lagrange multiplier
mu = 398600;
psi = 0;
c2 = 1/2;
c3 = 1/6;
psi_l = -4*pi;
psi_h = 4*pi^2;

c_dv = dot(r1, r2)/norm(r1)/norm(r2);
A = tm*sqrt(norm(r1)*norm(r2)*(1+c_dv));

count = 0;
err_dt = 10;
while(err_dt >1e-6)
    y =getY( r1, r2, A, psi, c2, c3 );
    while(A >0 && y<0)
       % disp('adjusting')
        psi_l = psi_l +.01;
        psi = (psi_l + psi_h)/2;
        y =getY( r1, r2, A, psi, c2, c3 );
    end
    
    X = sqrt(y/c2);
    dt = ((X^3)*c3 + A*sqrt(y))/sqrt(mu);
    
    if dt <= tau
        psi_l = psi;
    else
        psi_h = psi;
    end
    
    psi = (psi_l + psi_h)/2;
    
    [c2, c3] = lambertC(psi);
    err_dt = abs(dt-tau);
    count = count +1;
    if count >50
        disp('quitting')
        break;
    end
end
f = 1-(y/norm(r1));
g = A*sqrt(y/mu);
gdot = 1-(y/norm(r2));
end

