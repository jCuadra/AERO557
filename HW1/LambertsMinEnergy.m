function [a_min, e_min, t_min_abs, v_0] = LambertsMinEnergy( r_0, r )
%Lambert's minimum energy problem to find the velocity a r_0. Implemented
%from algorithm 56 in Vallado.
%   @param r_0          r vector at unknown velocity
%   @param r            other known r vector within the orbit
%   @return a_min       semi-major axis of the minimum transfer orbit
%   @return e_min       eccentricty of the minimum transfer orbit
%   @return t_min_abs   time fo flight corresponding to minimum semi-major 
%                       axis
%   @return v_0         velocity vector corresponding to r_0
mu = 398600;
r_0_N = norm(r_0);
r_N =  norm(r);
c_dv = dot(r_0, r)/r_0_N/r_N;
s_dv = sin(acos(c_dv));
c = sqrt(r_0_N^2+ r_N^2-2*r_0_N*r_N*c_dv);
s = (r_0_N + r_N + c)/2;
a_min = s/2;
p_min = r_0_N*r_N/c*(1-c_dv);
e_min = sqrt(1-(2*p_min/s));
alpha_e = pi;
s_be = sqrt((s-c)/s);
beta_e = 2*asin(s_be);
t_min_a = sqrt(a_min^3/mu)*(alpha_e - (beta_e -s_be));
t_min_abs = sqrt(2/mu)*(s^(3/2)-(s-c)^(3/2))/3;
%(r-(1-(r_N/p_min)*(1-c_dv))*r_0)
v_0 = sqrt(mu*p_min)/(r_0_N*r_N*s_dv).*(r-(1-(r_N/p_min)*(1-c_dv))*r_0);

end

