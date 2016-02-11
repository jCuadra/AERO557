function [ month, day, hr, min, sec ] = dayOfTheYear( dayIn )
%DAYOFTHEYEAR converts the date in day of the year form to month, day, min
%and sec
%   @param  dayIn    number between 1 an 365 to represent the day of
%                           of the year
%   @return month           
%   @return day
%   @return hr
%   @return min
%   @return sec
m = [0 31 59 90 120 151 181 212 243 273 304 224 365];
for i = 2:length(m)
    if(dayIn <= m(i))
        month = i-1;
        day = floor(dayIn-m(i-1));
        time = (dayIn-day-m(i-1))*24*3600; %s
        break;
    end
end
hr = floor(time/3600);
min = floor((time-hr*3600)/60);
sec = floor(time -hr*3600-min*60);

end

