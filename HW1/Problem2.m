%% Problem 2
%Morgan Yost
%AERO 557 HW1
clear all
%% Set Up
r0 = [15945.34; 0; 0];
r = [12214.83899; 10249.4731; 0];
tm = 1;
dt = 76*60;
mu = 398600;
%% Lambert's Universal
[f, g, gdot] = lambertUV(r0, r, dt, tm);
v0_universal = (r-f*r0)/g;
v_universal = (gdot*r-r0)/g;

%% Min Energy
[a_min, e_min, t_min_abs, v0_minE] = LambertsMinEnergy( r0, r );

%% Lambert's Gauss
[f, g, ~, gdot] = Lambert_Gauss( r0, r, dt, tm );
v0_LM = (r-f*r0)/g;
v_LM = (gdot*r-r0)/g;

%% Izzo Gooding
[v0_izzo, v_izzo, extremal_distances, exitflag] = IzzoGooding(...
        r0', r', dt/(60*60*24), 0, mu);
%% Result
fprintf('Universal Variable v0: %f %f %f\n', v0_universal');
fprintf('norm: %f\n', norm(v0_universal));
fprintf('Min Energy v0: %f %f %f\n', v0_minE);
fprintf('norm: %f\n', norm(v0_minE));
fprintf('Lambert Gauss v0: %f %f %f\n', v0_LM);
fprintf('norm: %f \n', norm(v0_LM));
fprintf('Izzo Gooding v0: %f %f %f\n', v0_izzo);
fprintf('norm: %f \n', norm(v0_izzo));
%% Explanation
fprintf('The variance in v vectors from observed r vectors and time are very \n')
fprintf('small. The difference comes from minor differences in implementation. \n')
fprintf('Lamberts min energy is a direct calculation scheme while the other\n')
fprintf('methods are not. It is not suprising to me that the min energy and \n')
fprintf('Lambert Gauss methods are so similar because they are both geometric \n')
fprintf('methods. \n')
