%% Problem 2
%Morgan Yost
%AERO 557 HW1
clear all
close all
%% Set Up
addpath ../Vallado
addpath ../HW1
global NUM_OBSV
NUM_OBSV = 3;
fid = fopen('singlePassInput.txt', 'r');
count = 0;
deg_sec = 306/86400;
Re = 6378;
while(count <7)
    count = count +1;
    line = fgetl(fid);
    type = strsplit(line, ' ');
    switch type{1}
        case 'LAT:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Latitude')
                return;
            end
            siteInfo.latSite = str2num(type{2});
            continue;
        case 'LONG:' 
            if(length(type) ~= 3)
                length(type)
                disp('Not enough inputs for Longitude')
                return;
            end
            siteInfo.longSite = str2num(type{2});
            siteInfo.longmin = str2num(type{3});
            continue;
        case 'ALT:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Altitude')
                return;
            end
            siteInfo.H = str2num(type{2});
            continue;
        case 'YEAR:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Year')
                return;
            end
            siteInfo.year(1:3) = str2num(type{2});
            continue;
        case 'MONTH:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Month')
                return;
            end
            siteInfo.month(1:3) = str2num(type{2});
            continue;
        case 'DAY:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Day')
                return;
            end
            siteInfo.day(1:3) = str2num(type{2});
            continue;
        case 'TIME'
            if type{2} == 'RA'
                if(length(type) ~=3)
                    disp('Not enough inputs for Day')
                    return;
                end
                for i = 1:length(type)
                    line = fgetl(fid);
                    word = strsplit(line, ' ');
                    time = strsplit(word{1}, ':');
                    siteInfo.hour(i) = str2num(time{1});
                    siteInfo.min(i) = str2num(time{2});
                    siteInfo.sec(i) = str2num(time{3});
                    alpha(i) = str2num(word{2});
                    delta(i) = str2num(word{3});
                end
            elseif type{2} == 'DEC'
                if(length(type) ~=3)
                    disp('Not enough inputs for Day')
                    return;
                end
                for i = 1:length(type)
                    line = fgetl(fid);
                    word = strsplit(line, ' ');
                    time = strsplit(word, ':');
                    siteInfo.hour(i) = str2num(time{1});
                    siteInfo.min(i) = str2num(time{2});
                    siteInfo.sec(i) = str2num(time{3});
                    alpha(i) = str2num(word{3});
                    delta(i) = str2num(word{2});
                end
            else
                disp('bad data table')
                return;
            end
        otherwise
            continue;
    end
end
fclose(fid);
if count ~= 7
    disp('Missing information for execution. Exiting.')
    return;
end

%% Gauss Angles Only
%Site Info Struct
%       H            Altitude of the observation location
%       longSite     Latitude of the observation location
%       latSite      Longitude of the observation location
%       latmin       Minute unit of the latitude of the observation
%                    location
%       year         Year of the observation
%       month        Month of the year of the ovservation
%       day          Day of the month of the observation
%       hour         Hour of the day of the observation
%       min          Minute of the hour of the observation
%       sec          Second of the minute of the observation
[r, ~, RSite, qHat, q, tau1, tau3] = gaussOD(alpha, delta, siteInfo, 1);
vAO = gibbs(r(:, 1), r(:, 2), r(:, 3)); %refine velocity
% Get COEs
[ ~, a, ecc, inc, Omega, w, theta ] = RV2COEd( r(:, 2), vAO);
fprintf('The following results are an approximation of the satelleite\n')
fprintf('at %d:%d:%f :\n', siteInfo.hour(2), siteInfo.min(2), siteInfo.sec(2));
fprintf('r: %f %f %f, \t norm(r): %f \n', r(:, 2), norm(r(:, 2)));
fprintf('v: %f %f %f, \t norm(v): %f \n', vAO, norm(vAO));
fprintf('\n')
fprintf('COES: \n a: %f km\n ecc: %f\n inc: %f deg\n RAAN: %f deg\n', a,...
    norm(ecc), inc, Omega);
fprintf(' w: %f deg\n theta: %f deg\n', w, theta);
%% TLE Decision
fprintf('\n')
fprintf('Because of the strong agreement with the RAAN and inclination vectors,\n');
fprintf('I think the object we observed was the second TLE. The RAAN and\n')
fprintf('inc angles are almost spot on, and the arguement of perigee and \n')
fprintf('true anomoly are about 180 degs off from the TLE which means the\n')
fprintf('arguement of perigee vector was just flipped in my conversion to \n')
fprintf('a COE. The only element that is of concern is the eccentricty. The\n')
fprintf('eccentricty that I calculated is off by a factor of 10. Even still,\n')
fprintf('the matching of the other COEs is enough to say with confidence that\n')
fprintf('the second TLE is the TLE for the object that was observed.\n')
%% Propogate R and V Forward
mu = 398600;
state = [r(:, 2);vAO]
numOrbs = 5;
options = odeset('RelTol', 1e-8);
T = 2*pi*sqrt(a^3/mu);
[t, y] = ode45(@orb_prop, [0 numOrbs*T], state, options);
r2 = y(end, 1:3)'
v2 = y(end, 4:6)'
[~, a2, ecc2, inc2, Omega2, w2, theta2] = RV2COEd(r2, v2);
fprintf('After propogating the R and V for %d periods of the orbit\n', numOrbs)
fprintf('COES: \n a: %f km\n ecc: %f\n inc: %f deg\n RAAN: %f deg\n', a2,...
    norm(ecc2), inc2, Omega2);
fprintf(' w: %f deg\n theta: %f deg\n', w2, theta2);

%% Propogate TLE forward
clear y state
[rTLE, vTLE, yr, time] = TLE2rv('singlePassTLE.txt');
state = [rTLE; vTLE]
norm([r(:, 2) rTLE])
[t, y] = ode45(@orb_prop, [0 numOrbs*T], state, options);
r3 = y(end, 1:3)'
v3 = y(end, 4:6)'
[~, a3, ecc3, inc3, Omega3, w3, theta3] = RV2COEd(r3, v3);
fprintf('After propogating the TLE for %d periods of the orbit\n', numOrbs)
fprintf('COES: \n a: %f km\n ecc: %f\n inc: %f deg\n RAAN: %f deg\n', a3,...
    norm(ecc3), inc3, Omega3);
fprintf('\n')
fprintf(' w: %f deg\n theta: %f deg\n', w3, theta3);
fprintf('After propogating the orbits for 5 periods, my COEs have remained\n')
fprintf('unchanged because I propogated them without pertibations. My\n')
fprintf('COEs are still off by the same amount that they were off by before\n')
fprintf('so I think that I am still correct.\n')
fprintf('\n')
%% Lambert UV
dts = [-tau1 tau3];
count = 1;
for i = 1:2
    [ f, g, g_dot ] = lambertUV( r(:, i), r(:, i+1), dts(i), 1 );
    v(:, count) = (r(:, i+1)-f*r(:, i))/g;
    count = count +1;
    v(:, count) = (g_dot*r(:, i+1)-r(:, i))/g;
    count = count+1;
end
vUV1 = v(:, 1);
vUV2 = mean([v(:, 2) v(:, 3)], 2)
vUV3 = v(:, 4);

%% Lambert-Gauss
clear v v1 v2 v3
count = 1;
for i = 1:2
    [ f, g, f_dot, g_dot ] = Lambert_Gauss( r(:, i), r(:, i+1), dts(i), 1 );
    v(:, count) = (r(:, i+1)-f*r(:, i))/g;
    count = count +1;
    v(:, count) = (g_dot*r(:, i+1)-r(:, i))/g;
    count = count +1;
end
vG1 = v(:, 1);
vG2 = mean([v(:, 2) v(:, 3)], 2)
vG3 = v(:, 4);

disp('UV error to angles only')
norm(vUV1-vAO)/norm(vAO)

disp('Gauss error to angles only')
norm(vG1-vAO)/norm(vAO)

disp('error between UV and Gauss first set')
norm(vUV1-vG1)/norm(vG1)


%% Plots for Sanity Check
% figure(1)
% hold on
% plot3([0 RSite(1, 1)], [0 RSite(2, 1)], [0 RSite(3, 1)])
% plot3([0 RSite(1, 2)], [0 RSite(2, 2)], [0 RSite(3, 2)], 'r')
% plot3([0 RSite(1, 3)], [0 RSite(2, 3)], [0 RSite(3, 3)], 'g')
% %%
% plot3([RSite(1,1) qHat(1, 1)+RSite(1,1)], [RSite(2, 1) qHat(2, 1)+RSite(2, 1)]...
%     , [RSite(3,1) qHat(3, 1)+RSite(3,1)], 'c')
% plot3([RSite(1, 2) qHat(1, 2)+RSite(1, 2)], [RSite(2, 2) qHat(2, 2)+ RSite(2, 2)],...
%     [RSite(3, 2) qHat(3, 2)+RSite(3, 2)], 'm')
% plot3([RSite(1, 3) qHat(1, 3)+RSite(1, 3)], [RSite(2, 3) qHat(2, 3)+RSite(2, 3)],...
%     [RSite(3, 3) qHat(3, 3)+RSite(3, 3)], 'y')
% plot3([0 r(1, 1)], [0 r(2, 1)], [0 r(3, 1)], 'k')
% plot3([0 r(1, 2)], [0 r(2, 2)], [0 r(3, 2)], 'k')
% plot3([0 r(1, 3)], [0 r(2, 3)], [0 r(3, 3)], 'k')
% view(-10, 0)