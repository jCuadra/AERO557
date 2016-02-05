function [ rf, vf ] = LambertUniversal( r, v, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u = 398600;
r0 = norm(r)
v0 = norm(v)
vr0 = dot( v, r)/r0
%[ h, i, RAAN, ecc, w, theta] = rv2COE(r, v);
%a = get_a(h, u, ecc);
%alpha = 1/a
alpha = 2/r0-(v0^2/u)
X = (sqrt(u))*t*abs(alpha)
z = alpha*X^2
tol = 1e-8;
ratio = 1;
count = 0;
while(ratio >= tol)
    if count >1000
        break;
    end
    c = C(z);
    s = S(z);

    top = (r0*vr0/sqrt(u)*(X^2)*c + (1-alpha*r0)*s*X^3 + r0*X-(sqrt(u))*t);
  
    bottom = (r0*vr0/sqrt(u)*X*(1-z*s) + (1-alpha*r0)*(X^2)*c+r0);
    ratio = top/bottom;
    X = X - ratio
    z = alpha*X^2;
   count = count +1;
end
c = C(z);
s = S(z);
f = 1 - ((X^2)/r0)*c; %1-(X^2)*C/v0;
g = t - (1/u)*s*X^3;%t-(1/(u^.5))*(X^3)*S; %t - (1/u)*S*X^3;
rf = f*r +g*v;
rfn = norm(rf);
fDot = ((u^.5)/(rfn*r0))*(X*(z*s -1));
gDot = 1-(X^2/rfn)*c;
vf = fDot*r +gDot*v;
fprintf('\n Final position vector (km):')
   %  fprintf('\n   r = (%g, %g, %g)\n', rf(1), rf(2), rf(3))
     fprintf('\n Final velocity vector (km/s):')
    % fprintf('\n   v = (%g, %g, %g)', vf(1), vf(2), vf(3))
end

