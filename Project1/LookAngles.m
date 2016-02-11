%Change num obsv to 3 so we can find rSite and qhat the same way and get 3
%ra and dec angles. need to change RSite...
clearvars
close all
addpath ../HW1
addpath ../HW2
global NUM_OBSV;
NUM_OBSV = 3;
%%Cal Poly Location
siteInfo.latSite = 35.30; % north 
siteInfo.longSite = -120.66; % west 
siteInfo.longmin = 0;
siteInfo.H = 105.8/1000; % km altitude
siteInfo.year = 2016;
siteInfo.month = 2;
siteInfo.day = 10;
siteInfo.hour = 12+6.5;
obs_hr = 12+9; %local
siteInfo.min = 0;
siteInfo.sec = 0;
dLocalTime = 8;
%%TLE to RV
[r, v, yr, time] = TLE2rv('TLE.txt');
[m, d, hr, min, sec] = dayOfTheYear(time); %TLE time in UTC
UT = hr + min/60 + sec/3600;
%%propogate RV to when pass happens
[lstf, UTf] = LST( siteInfo.year, siteInfo.month, siteInfo.day, siteInfo.hour,...
                    siteInfo.min, siteInfo.sec,siteInfo.longSite, siteInfo.longmin); %end of obsevation window, local
dt = (UTf+dLocalTime-UT) + (siteInfo.day-d-dLocalTime/24)*24; %add dLocalTime for universal

state = [r;v];
options = odeset('RelTol', 1e-8);
[t, y] = ode45(@orb_prop, [0 dt*3600], state, options);
t = UT*3600+t; %put time in local time
ind = find(t >= siteInfo.hour, 1, 'first');
t(ind)
count = 1;
Az = zeros(1, length(1:100:length(t)));
El = Az;
for i = ind:100:length(t)
    r = y(i, 1:3)';
    %%get Rsite
    siteInfo.sec = t(i);
    RSite = getRSite(siteInfo, 1);
    %%get RA and DEC
    RA_DEC = topocentric_d(r, RSite);
    %%get slant range
    qHat = getQHat(RA_DEC(1), RA_DEC(2));
    %rho = (r-Rsite)\qHat;
    %%get az and el
    [Az(count), El(count)] = getAzEl(qHat);
    obs_time(count) = t(i);
    count = count+1;
end
figure(1)
hold on
%plot3([0 RSite(1)], [0 RSite(2)], [0 RSite(3)])
plot3([RSite(1) qHat(1)+RSite(1)], [RSite(2) qHat(2)+RSite(2)]...
    , [RSite(3) qHat(3)+RSite(3)], 'c')
%plot3([0 r(1)], [0 r(2)], [0 r(3)], 'k')
% plot3([0 r(1, 2)], [0 r(2, 2)], [0 r(3, 2)], 'k')
% plot3([0 r(1, 3)], [0 r(2, 3)], [0 r(3, 3)], 'k')
%plot3(El, Az, obs_time)
%plot3(y(ind:end, 1), y(ind:end, 2), y(ind:end, 3))