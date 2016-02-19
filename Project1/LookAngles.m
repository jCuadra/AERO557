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
siteInfo.year(1:3) = 2016;
siteInfo.month(1:3) = 2;
siteInfo.day(1:3) = 16;
siteInfo.hour = 12+6;
obs_hr = 12+9; %local
siteInfo.min(1:3) = [20 25 30];
siteInfo.sec(1:3) = 0;
dLocalTime = 8;
%%TLE to RV
[r, v, yr, time] = TLE2rv('TLE.txt');
[m, d, hr, min, sec] = dayOfTheYear(time) %TLE time in UTC
UT = hr-dLocalTime + min/60 + sec/3600;
%%propogate RV to when pass happens
[lstf, UTf] = LST( siteInfo.year(1), siteInfo.month(1), siteInfo.day(1), siteInfo.hour(1),...
                    siteInfo.min(1), siteInfo.sec(1),siteInfo.longSite, siteInfo.longmin); %end of obsevation window, local
dt = siteInfo.hour + siteInfo.min(3)/60-UT
                %dt = (UTf-UT) + (siteInfo.day(1)-d+dLocalTime/24)*24; %add dLocalTime for universal
%%
state = [r;v];
options = odeset('RelTol', 1e-8);
[t, y] = ode45(@orb_prop, [UT*3600 (UT+dt)*3600], state, options);
%t = siteInfo.hour*3600+siteInfo.min(1)*60+t; %put time in local time
ind = find(t >= siteInfo.hour*3600, 1, 'first');
count = 1;
count1 = 1;
lookAt = [ind:3:length(t)];

figure(1)
hold on
HAaz = [328 359 5 10 15 25 30 37 44 51 56];
HAel = [0 5 8 10 10.5 11 11 10.5 10 9 8];
for j = 1:length(lookAt)
    %%get az and el
    i= lookAt(j);
    [~, ~, h, min, sec] = dayOfTheYear(t(i)/24/3600);
   if(h == 18 && min >= 18 && min <=25 && count <=length(HAaz) )
       count
        r = y(i:i+2, 1:3)';
        %%get Rsite
        [~, ~, siteInfo.hour(1:3), siteInfo.min(1:3), siteInfo.sec(1:3)] = ...
            dayOfTheYear(t(i:i+2)/24/3600);
        h
        min
        dcm = dcmeci2ecef('IAU-2000/2006', [2000+yr, m, d, h, min, sec]);
        [RSite, LSTcurrent] = getRSite(siteInfo, 1);
        %r = dcm*r;
%         RSite2 = dcm*RSite(:, 2);
%         plane = cross(r(:, 3), r(:, 1));
%         ElAngle = asind(dot(plane, RSite2)/norm(plane)/norm(RSite2));
%         ElAngleSite = acosd(dot(RSite2, [0 1 0])/norm(RSite2))
%         ElAngleSat = acosd(dot([0 1 0], r(:, 2))/norm(r(:, 2)))-siteInfo.latSite
%         El2(count) = -siteInfo.latSite+ ElAngleSite-ElAngleSat;
%         AzAngleSite = acosd(dot(RSite2, [1 0 0])/norm(RSite2))
%         AzAngleSat = acosd(dot(r(:, 2), [1 0 0])/norm(r(:, 2)))
        %%get RA and DEC
      %  RSite = -1*RSite;
        RSite = abs(RSite);
        r = abs(dcm*r);
        [ra(count), dec(count)] = azl2radc(HAaz(count)*pi/180, HAel(count)*pi/180, siteInfo.latSite*pi/180, LSTcurrent(2)*pi/180);
        ra(count) = ra(count)*180/pi;
        if ra(count)>360
            ra(count) = ra(count)-360;
        end
        dec(count) = dec(count)*180/pi;
        RA_DEC = topocentric_d(r, RSite, 3);
%         Az(count) = AzAngleSat-siteInfo.longSite;
%         if Az(count)>360
%             Az(count) = Az(count)-360;
%         end
        out(:, count) = RA_DEC(:, 2);
        %%get slant range
        qHat = getQHat(RA_DEC(1, :), RA_DEC(2, :));
        rho = (r-RSite);
        [ Az(count), El(count)] = getAzEl( rho(:, 2), siteInfo.latSite,...
            siteInfo.longSite, LSTcurrent(2));
       % [Az(count), El(count)] = getAzEl2(LSTcurrent(2), RA_DEC(1), RA_DEC(2), ...
          %  siteInfo.latSite);
          
        obs_time(count) = datenum(2000+yr, m, d, h, min, sec);
       % if(El(count)>0
            fprintf('The satellite is passing over at %d:%d:%d\n', h,  min, sec)
            fprintf('At an Azimuth of %f\n Heavens above Az: %f\n',Az(count), HAaz(count))
            fprintf('And a Elevation of %f\n Heavens Above El: %f\n\n', El(count), HAel(count))
            fprintf('Heavens above predicts RA: %f \n ', ra(count))
            fprintf('I predict the RA will be: %f\n', out(1, count));
            fprintf('Heavens above predicts Dec: %f \n ', dec(count))
            fprintf('I predict the Dec will be: %f\n', out(2, count));
        %out
        count = count+1;
        
        plot3([0 r(1)], [0 r(2)], [0 r(3)])
        plot3([0 RSite(1)], [0 RSite(2)], [0 RSite(3)], 'r')
   elseif (h == 18 && min >=14 && min <18)
        r = y(i:i+2, 1:3)';
        %%get Rsite
        [~, ~, siteInfo.hour(1:3), siteInfo.min(1:3), siteInfo.sec(1:3)] = ...
            dayOfTheYear(t(i:i+2)/24/3600);
        dcm = dcmeci2ecef('IAU-2000/2006', [2000+yr, m, d, h, min, sec]);
        [RSite, LSTcurrent] = getRSite(siteInfo, 1);
        RSite = abs(RSite);
        r = abs(dcm*r);
        RA_DEC = topocentric_d(r, RSite, 3);
        out1(:, count1) = RA_DEC(:, 2);
        qHat = getQHat(RA_DEC(1, :), RA_DEC(2, :));
        rho = (r-RSite);
        [ Az1(count1), El1(count1)] = getAzEl( rho(:, 2), siteInfo.latSite,...
            siteInfo.longSite, LSTcurrent(2));
        obs_time1(count1) = datenum(2000+yr, m, d, h, min, sec);
        count1 = count1+1;
    end
    if(h >=18 && min >26 && sec >=30)
        break;
    end
    
end
fSize = 14;
figure(2)
plot(obs_time, El, 'Linewidth', 2)
hold on
plot(obs_time, HAel, 'Linewidth', 2)
datetick('x', 'HH:MM:SS')
title('My Predicted Elevation versus Heavens Above (Interpolated)','FontSize', fSize)
xlabel('Time','FontSize', fSize)
ylabel('Elevation (deg)','FontSize', fSize)
legend('My El', 'Heavens Above')

figure(3)
plot(obs_time, Az, 'Linewidth', 2)
hold on
plot(obs_time, HAaz, 'Linewidth', 2)
datetick('x', 'HH:MM:SS')
title('My Predicted Azimuth versus Heavens Above (Interpolated)','FontSize', fSize)
xlabel('Time','FontSize', fSize)
ylabel('Azimuth (deg)','FontSize', fSize)
legend('My Az', 'Heavens Above')

figure(4)
plot([obs_time1 obs_time], [Az1 Az], 'Linewidth', 2)
datetick('x', 'HH:MM:SS')
title('My Azimuth Prediction','FontSize', fSize)
xlabel('Time','FontSize', fSize)
ylabel('Azimuth (deg)','FontSize', fSize)

figure(5)
plot([obs_time1 obs_time], [El1 El], 'Linewidth', 2)
datetick('x', 'HH:MM:SS')
title('My Elevation Prediction', 'FontSize', fSize)
xlabel('Time', 'FontSize', 14)
ylabel('Elevation (deg)', 'FontSize', fSize)

myObs = [340 345 356 13 27 32 47];
myObsmin = [22 22 23 24 24 24 25];
myObssec = [0 30 30 0 20 40 0];
for i = 1:length(myObs)
obs_time2(i) = datenum(2000+yr, m, d, h, myObsmin(i), myObssec(i) );
end
figure(6)
plot(obs_time(1:length(myObs)), Az(1:length(myObs)), 'Linewidth', 2)
hold on
plot(obs_time2, myObs, 'Linewidth', 2)
datetick('x', 'HH:MM:SS')
title('My Predicted Azimuth versus My Actual Observations','FontSize', fSize)
xlabel('Time','FontSize', fSize)
ylabel('Azimuth (deg)','FontSize', fSize)
legend('My Az', 'Heavens Above')


% figure(1)
% hold on
% plot3([0 RSite(1)], [0 RSite(2)], [0 RSite(3)])
% plot3([RSite(1) qHat(1)+RSite(1)], [RSite(2) qHat(2)+RSite(2)]...
%     , [RSite(3) qHat(3)+RSite(3)], 'c')
% plot3([0 r(1)], [0 r(2)], [0 r(3)], 'k')
% plot3([0 r(1, 2)], [0 r(2, 2)], [0 r(3, 2)], 'k')
% plot3([0 r(1, 3)], [0 r(2, 3)], [0 r(3, 3)], 'k')
% plot3(El, Az, obs_time)
% plot3(y(ind:end, 1), y(ind:end, 2), y(ind:end, 3))