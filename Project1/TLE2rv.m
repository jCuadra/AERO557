function [ r, v, year, time] = TLE2rv( TLEpath)
%TLE2RV takes in a spacecraft TLE in a text file and parse the file to pass
%the TLE elements to the COE2rvd function
%   @param TLEpath  Path to TLE file, string
%   @return r       radius vector of the satellite, km
%   @return v       velocity vector of the satellite, km/s
%   @retrun t       time of the TLE
mu = 398600;
fid = fopen(TLEpath, 'r');
count = 0;
while(count <2)
    count = count +1;
    line = fgetl(fid);
    type = strsplit(line, ' ');
    if count ==1
        date = type{4};
        year = str2double(date(1:2));
        time = str2double(date(3:end));
    end
end
fclose(fid);
% extract things from second line
inc = str2double(type{3});
RAAN = str2double(type{4});
ecc = str2double(strcat(['.', type{5}]));
w = str2double(type{6});
M = str2double(type{7});
E = get_E(M, ecc);
theta = acosd((cos(E)-ecc)/(1-ecc*cos(E)));
mm = str2double(type{8})*2*pi/3600/24; % convert from rev/day to rad/s
a = (mu/mm^2)^(1/3);
if ecc ~=1
    p = a*(1-ecc^2);
else
    disp('You did not write TLE to rv to handle hyperbolic orbits')
end

[r, v ] = COE2rvd(p, ecc, inc, RAAN, w, theta);

end

