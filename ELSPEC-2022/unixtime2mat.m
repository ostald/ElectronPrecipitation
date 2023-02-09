
function matlab_time = unixtime2mat(unix_time);
% unixtime2mat  Converts unix time stamps (seconds since Jan 1, 1970) to
%               Matlab serial date number (decimal days since Jan 1 0000).
%               
%               USAGE:
%                      unixtime2mat(unix_time)
%
%
%               The function may not handle leap years or leap seconds
%               appropriately.  
%               
%               Originally written by:
%               Val Schmidt 
%               Center for Coastal and Ocean Mapping
%               2007
%               
%               Original Link: 
%               http://www.mathworks.com/matlabcentral/fileexchange/24024-convert-unix-time-seconds-since-jan-1-1970-to-matlab-serial-time/content/unixtime2mat.m
%               
%               Tweaked by:
%               João Nogueira
%               Portugal Telecom Inovação, SA
%               Instituto de Telecomunicações de Aveiro
%               Universidade de Aveiroo
%               2012
%
matlab_time = unix_time./86400 + datenummx(1970,1,1,0,0,0);
