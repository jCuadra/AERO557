function [ tau ] = getTau( siteInfo, a, b )
%getTau Calculated the time between observation a and b (a - b)
%   @param siteInfo    Previouslly used siteInfo structure
%   @param a           First time index
%   @param b           Second time index, to be subtracted from a
%   @return tau        Time between a and b, seconds
tau = 3600*siteInfo.hour(a) + 60*siteInfo.min(a) + siteInfo.sec(a)-...
       (3600*siteInfo.hour(b) + 60*siteInfo.min(b) + siteInfo.sec(b));


end

