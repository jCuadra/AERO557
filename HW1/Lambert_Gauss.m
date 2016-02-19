function [ f, g, f_dot, g_dot ] = Lambert_Gauss( r_0, r, dt, tm )
%Lambert- Gauss solver to obtain the f and g lambert variables in order to
%find velocity.
%   @param r_0      Initial r vector, 3x1 km
%   @param r        Final r vector, 3x1 km
%   @param dt       Time between r vectors s
%   @param tm       1 or -1 for short way or long way around
%   @return f       Lagrange multiplier
%   @return g       Lagrange multiplier
%   @return f_dot   Lagrange multiplier
%   @return g_dot   Lagrange multiplier
mu = 398600;
r_N =  norm(r);
r_0_N = norm(r_0);

c_dv = dot(r_0, r)/(r_0_N*r_N);
dv = acos(c_dv);
c_dv_2 = cos(dv/2);
s_dv = tm*sqrt(1-c_dv^2);

l = ((r_0_N + r_N)/(4*sqrt(r_0_N*r_N)*c_dv_2)) - 1/2;
m = mu*dt^2/((2*sqrt(r_0_N*r_N)*c_dv_2)^3);
y = 1;

err = 10;
tol = 1e-6;
yOld = 1;
while err >tol
    x1 = m/y^2 - l;
    x2 = (4/3)*(1 + 6*x1/5 + 48/35*x1^2+ 480/315*x1^3);
    y = 1+x2*(l+x1);
    err = abs(y - yOld);
    yOld = y;
end
    
    c_dE_2 = 1-2*x1;
    %2*acos(c_dE_2)
   % if 2*acos(c_dE_2) > 0 % eccentric
        p = (r_0_N*r_N*(1-c_dv))/(r_0_N+r_N - 2*sqrt(r_0_N*r_N)*c_dv_2*c_dE_2);
        %end
        
        if 2*acos(c_dE_2)  < 0 % hyperbolic
            cs_dH_2 = 1-2*x1;
            p = (r_0_N*r_N*(1-c_dv))/(r_0_N+r_N-2*sqrt(r_0_N*r_N)*c_dv_2*cs_dH_2);
        end
    
f = 1-r_N/p*(1-c_dv);
g = r_N*r_0_N*s_dv/sqrt(mu*p);
f_dot = sqrt(1/p)*tan(dv/2)*((1-c_dv)/p - 1/r_N - 1/r_0_N);
g_dot = 1-r_0_N/p*(1-c_dv);
%v_0 = (r-f*r_0)/g
%v = (g_dot*r-r_0)/g

end

