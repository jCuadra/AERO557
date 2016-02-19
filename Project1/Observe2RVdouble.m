clearvars
addpath ../HW1
global NUM_OBSV NUM_VAR
NUM_OBSV = 8;
NUM_VAR = 2;
fid = fopen('doublePassInput.txt', 'r');
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
            siteInfo.latSite = str2double(type{2});
            continue;
        case 'LONG:' 
            if(length(type) ~= 3)
                length(type)
                disp('Not enough inputs for Longitude')
                return;
            end
            siteInfo.longSite = str2double(type{2});
            siteInfo.longmin = str2double(type{3});
            continue;
        case 'ALT:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Altitude')
                return;
            end
            siteInfo.H = str2double(type{2});
            continue;
        case 'YEAR:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Year')
                return;
            end
            siteInfo.year(1:NUM_OBSV) = str2double(type{2});
            continue;
        case 'MONTH:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Month')
                return;
            end
            siteInfo.month(1:NUM_OBSV) = str2double(type{2});
            continue;
        case 'DAY:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Day')
                return;
            end
            siteInfo.day(1:NUM_OBSV) = str2double(type{2});
            continue;
        case 'TIME'
            if strcmp(type{2},'RA')
                if(length(type) ~=3)
                    disp('Not enough inputs for RA')
                    return;
                end
                line = fgetl(fid);
                i = 1;
                while(line ~= -1)
                    word = strsplit(line, ' ');
                    time = strsplit(word{1}, ':');
                    siteInfo.hour(i) = str2double(time{1});
                    siteInfo.min(i) = str2double(time{2});
                    siteInfo.sec(i) = str2double(time{3});
                    alpha(i) = str2double(word{2});
                    delta(i) = str2double(word{3});
                    line = fgetl(fid);
                    i = i+1;
                end
            elseif strcmp(type{2},'DEC')
                if(length(type) ~=3)
                    disp('Not enough inputs for DEC')
                    return;
                end
                line = fgetl(fid);
                i = 1;
                while(line ~= -1)
                    word = strsplit(line, ' ');
                    time = strsplit(word, ':');
                    siteInfo.hour(i) = str2double(time{1});
                    siteInfo.min(i) = str2double(time{2});
                    siteInfo.sec(i) = str2double(time{3});
                    alpha(i) = str2double(word{3});
                    delta(i) = str2double(word{2});
                    line = fgetl(fid);
                    i = i+1;
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
len = 9;
keep = 1;
options = odeset('RelTol', 1e-8);
t = 3600*siteInfo.hour(:) + 60*siteInfo.min(:) + siteInfo.sec(:);
for i = [1 2 5 6]
    [r, v, RSite, qHat, q, tau1(keep),tau3(keep)] = gaussOD(alpha, delta, siteInfo, i);
    state = [r(:, 2); v];
    [ ~, y] = ode45(@orb_prop, [t(i+1) t(1)], state, options); %propogate middle vector
    rprop(:, keep) = y(end, 1:3)';
    vprop(:, keep) = y(end, 4:6)';
    if i == 1 || i == 5
    rkeep(:, i:i+2) = r;
    end
    clear y
    keep = keep +1;
end

%% Average vectors
rnom = mean(rprop, 2);
vnom = mean(vprop, 2);

%% Get COEs

[ ~, a, ecc, inc, Omega, w, theta ] = RV2COEd( rnom, vnom);
fprintf('The orbital elements are:\n');
fprintf('COES: \n a: %f km\n ecc: %f\n inc: %f deg\n RAAN: %f deg\n', a,...
    norm(ecc), inc, Omega);
fprintf(' w: %f deg\n theta: %f deg\n', w, theta);
fprintf('\n')
fprintf('The R and V vector at %d:%d:%f was calculated by using the gauss\n',...
    siteInfo.hour(1), siteInfo.min(1), siteInfo.sec(1))
fprintf('angles only method on observations 1 through 3, 2 through 4, 5 through\n')
fprintf('7 and 6 through 8 then propogating all the r and v vectors to the \n')
fprintf('time of the first observation and averaging the 4 vectors. The r\n')
fprintf('vector was found to be: \n r = %f i %f j %f k km\n', rnom)
fprintf('and the v vector was found to be: \n v = %f i %f j %f k km/s\n', vnom)

%% TLE Decision
fprintf('\n')
fprintf('Because of the strong agreement with the RAAN and inclination vectors,\n');
fprintf('I think the object we observed was the first TLE. The RAAN and\n')
fprintf('inc angles are very close. Once again, the ecc is of concern  The\n')
fprintf('eccentricty that I calculated is off by quite a lot, but none of\n')
fprintf('the other TLEs have a more promising value for eccentricity.\n')

%% Propogate R and V Forward
mu = 398600;
state = [rnom;vnom];
numOrbs = 5;
options = odeset('RelTol', 1e-8);
T = 2*pi*sqrt(a^3/mu);
[t, y] = ode45(@orb_prop, [0 numOrbs*T], state, options);
r2 = y(end, 1:3)';
v2 = y(end, 4:6)';
[~, a2, ecc2, inc2, Omega2, w2, theta2] = RV2COEd(r2, v2);
fprintf('After propogating the R and V for %d periods of the orbit\n', numOrbs)
fprintf('COES: \n a: %f km\n ecc: %f\n inc: %f deg\n RAAN: %f deg\n', a2,...
    norm(ecc2), inc2, Omega2);
fprintf(' w: %f deg\n theta: %f deg\n', w2, theta2);

%% Propogate TLE forward
clear y state
[rTLE, vTLE, yr, time] = TLE2rv('doublePassTLE.txt');
state = [rTLE; vTLE];
[t, y] = ode45(@orb_prop, [0 numOrbs*T], state, options);
r3 = y(end, 1:3)';
v3 = y(end, 4:6)';
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
dts = [-tau1(1) tau3(1)];
count = 1;
for i = 1:2
    [ f, g, g_dot ] = lambertUV( rkeep(:, i), rkeep(:, i+1), dts(i), 1 );
    v(:, count) = (rkeep(:, i+1)-f*rkeep(:, i))/g;
    count = count +1;
    v(:, count) = (g_dot*rkeep(:, i+1)-rkeep(:, i))/g;
    count = count+1;
end
vUV1a = v(:, 1);
vUV2a = mean([v(:, 2) v(:, 3)], 2)
vUV3a = v(:, 4);

dts = [-tau1(3) tau3(3)];
count = 1;
for i = 5:6
    [ f, g, g_dot ] = lambertUV( rkeep(:, i), rkeep(:, i+1), dts(i-4), 1 );
    v(:, count) = (rkeep(:, i+1)-f*rkeep(:, i))/g;
    count = count +1;
    v(:, count) = (g_dot*rkeep(:, i+1)-rkeep(:, i))/g;
    count = count+1;
end
vUV1 = v(:, 1);
vUV2 = mean([v(:, 2) v(:, 3)], 2)
vUV3 = v(:, 4);
%% Lambert-Gauss
clear v v1 v2 v3
count = 1;
dts = [-tau1(1) tau3(1)];
for i = 1:2
    [ f, g, f_dot, g_dot ] = Lambert_Gauss( rkeep(:, i), rkeep(:, i+1), dts(i), 1 );
    v(:, count) = (rkeep(:, i+1)-f*rkeep(:, i))/g;
    count = count +1;
    v(:, count) = (g_dot*rkeep(:, i+1)-rkeep(:, i))/g;
    count = count+1;
end
vG1a = v(:, 1);
vG2a = mean([v(:, 2) v(:, 3)], 2)
vG3a = v(:, 4);

dts = [-tau1(3) tau3(3)];
for i = 5:6
    [ f, g, f_dot, g_dot ] = Lambert_Gauss( rkeep(:, i), rkeep(:, i+1), dts(i-4), 1 );
    v(:, count) = (rkeep(:, i+1)-f*rkeep(:, i))/g;
    count = count +1;
    v(:, count) = (g_dot*rkeep(:, i+1)-rkeep(:, i))/g;
    count = count+1;
end
vG1 = v(:, 1);
vG2 = mean([v(:, 2) v(:, 3)], 2)
vG3 = v(:, 4);

disp('UV error to angles only')
norm(vUV1a-vnom)/norm(vnom)

disp('Gauss error to angles only')
norm(vG1a-vnom)/norm(vnom)

disp('error between UV and Gauss first set')
norm(vUV1a-vG1a)/norm(vG1a)


disp('error between UV and Gauss second set')
norm(vUV1-vG1)/norm(vG1)




