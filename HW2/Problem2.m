clear all
addpath ../HW1
global NUM_OBSV NUM_VAR
NUM_OBSV = 9;
NUM_VAR = 2;
fid = fopen('input.txt', 'r');
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
            siteInfo.year(1:NUM_OBSV) = str2num(type{2});
            continue;
        case 'MONTH:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Month')
                return;
            end
            siteInfo.month(1:NUM_OBSV) = str2num(type{2});
            continue;
        case 'DAY:'
            if(length(type) ~= 2)
                disp('Not enough inputs for Day')
                return;
            end
            siteInfo.day(1:NUM_OBSV) = str2num(type{2});
            continue;
        case 'TIME'
            if type{2} == 'RA'
                if(length(type) ~=3)
                    disp('Not enough inputs for RA')
                    return;
                end
                line = fgetl(fid);
                i = 1;
                while(line ~= -1)
                    word = strsplit(line, ' ');
                    time = strsplit(word{1}, ':');
                    siteInfo.hour(i) = str2num(time{1});
                    siteInfo.min(i) = str2num(time{2});
                    siteInfo.sec(i) = str2num(time{3});
                    alpha(i) = str2num(word{2});
                    delta(i) = str2num(word{3});
                    line = fgetl(fid);
                    i = i+1;
                end
            elseif type{2} == 'DEC'
                if(length(type) ~=3)
                    disp('Not enough inputs for DEC')
                    return;
                end
                line = fgetl(fid);
                i = 1;
                while(line ~= -1)
                    word = strsplit(line, ' ');
                    time = strsplit(word, ':');
                    siteInfo.hour(i) = str2num(time{1});
                    siteInfo.min(i) = str2num(time{2});
                    siteInfo.sec(i) = str2num(time{3});
                    alpha(i) = str2num(word{3});
                    delta(i) = str2num(word{2});
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
raBias = .0017; %deg
decBias = .001; %deg
len = 9;
[r, v, RSite, qHat, q, tau1, tau3] = gaussOD(alpha, delta, siteInfo);
%% WLS Method
W = [1/raBias^2 0; 0 1/decBias^2];
W = W/norm(W);

% prop to time 1
options = odeset('RelTol', 1e-8);
t = 3600*siteInfo.hour(:) + 60*siteInfo.min(:) + siteInfo.sec(:);

state1 = [r(:, 2); v(:, 1)];
[ ~, y1] = ode45(@orb_prop, [t(2) t(1)], state1, options);
r1 = y1(end, 1:3);
v1 = y1(end, 4:6);

state2 = [r(:, 5); v(:, 4)];
[ ~, y2] = ode45(@orb_prop, [t(5) t(1)], state2, options);
r2 = y2(end, 1:3);
v2 = y2(end, 4:6);

state3 = [r(:, 8); v(:, 7)];
[ ~, y3] = ode45(@orb_prop, [t(8) t(1)], state3, options);
r3 = y3(end, 1:3);
v3 = y3(end, 4:6);


rnom = mean([r1; r2; r3]);
vnom = mean([v1; v2; v3]);
Xnom = [rnom'; vnom';];
%%
obsOrig = [alpha; delta];
current_err = 10;
count = 0;
while( current_err >1e-8 && count <50)
    %Prop Xnorm to t9
    state = Xnom(:, 1);
    for i = 1:NUM_OBSV
        if i >1
            [ ~, y] = ode45(@orb_prop, [t(1) t(i)], state, options);
            Xnom(:, i) = y(end, 1:6)';
            Xnom(:, i) = Xnom(:, 1);
           % state = y(end, 1:6)';
        end
        %Finite Difference
        xdelta(:, i) = Xnom(:, i)*.001;
        XMod(:, i) = Xnom(:, i)+xdelta(:, i);
        
        %convert Xnom back to RA and DEC
        obsNom(:, i) = topocentric_d(Xnom(1:3, i), RSite(:, i));
        %Convert Xmod back to RA and DEC
        obsMod(:, i) = topocentric_d(XMod(1:3, i), RSite(:, i));

    end

    %Calculate H
    diff = (obsMod-obsNom);
    for i = 1:NUM_OBSV
        for j = 1:length(state)
            H(:, j, i) = diff(:, i)/xdelta(j,i);
        end
    end
    Havg = mean(H, 3);

    %Calculate y tilda
    y_tild = obsOrig-obsMod;
    y_tildAvg = mean(y_tild,2);
    
    %Calculate P
    [u, s, v] = svd(Havg'*W*Havg);
    P = v*pinv(s)*u';
    
    %Get RMS
    new_rms = sqrt(y_tildAvg'*W*y_tildAvg/NUM_VAR/NUM_OBSV);
    if count >0
        current_err = abs(old_rms/new_rms - new_rms);
    else
        current_err = new_rms;
    end
   % current_err = RMS_err(count+1);
    err(count+1) = current_err; %new_rms;
    deltaX = P*Havg'*W*y_tild;
    Xnom = Xnom + deltaX;
    count = count +1;
    old_rms = new_rms;
    
end


[ ~, a, ecc, inc, Omega, ~, ~ ] = RV2COEd( Xnom(1:3, 1), Xnom(4:6, 1));
fprintf('The orbital elements are:\n');
fprintf('COES: \n a: %f km\n ecc: %f\n inc: %f deg\n RAAN: %f deg\n', a,...
    norm(ecc), inc, Omega);
confidence = diag(P);
ijk = 'ijk';
fprintf('r= \n')
for i = 1:3
    fprintf(' %f+/-%f %c km\n', Xnom(i, 1), confidence(i), ijk(i));
end
fprintf('v = \n')
for i = 1:3
    fprintf(' %f+/-%f %c km/s \n',Xnom(i+3, 1), confidence(i+3), ijk(i));
end

%% Error Ellipsoid
close all

[r2_err] = sqrt(abs(eig(Xnom(1:3, 1:3))));
[r5_err] = sqrt(abs(eig(Xnom(1:3, 4:6))));
[r8_err] = sqrt(abs(eig(Xnom(1:3, 7:9))));


fprintf('The error ellipsoids for the 3 r vectors are:\n')

figure(1)
[x, y, z] = ellipsoid(0, 0, 0, r2_err(1), r2_err(2), r2_err(3));
surf(x, y, z)
title('Error ellipsoid for second observation r vector')
xlabel('km')
ylabel('km')
zlabel('km')
figure(2)
[x, y, z] = ellipsoid(0, 0, 0, r5_err(1), r5_err(2), r5_err(3));
surf(x, y, z)
title('Error ellipsoid for fifth observation r vector')
xlabel('km')
ylabel('km')
zlabel('km')
figure(3)
[x, y, z] = ellipsoid(0, 0, 0, r8_err(1), r8_err(2), r8_err(3));
surf(x, y, z)
title('Error ellipsoid for eigth observation r vector')
xlabel('km')
ylabel('km')
zlabel('km')

%% Analysis
worst_err = max(max([r2_err r5_err r8_err]));
fprintf('The worst case error is +/- %f km\n', worst_err)
fprintf('Because this is a GEO bird, this is only and error of approx %f deg\n',...
    worst_err/norm(Xnom(1:3, 1))*180/pi)
fprintf('in the observation angles. In order to improve this error, I would\n')
fprintf('probably apply a Kalman filter to my data so that the state and \n')
fprintf('covariance matrix are propogated in the predictor correct scheme \n')
fprintf('to hopefully result in more accurate mixing of data from always \n')
fprintf('updating the epoch.\n')