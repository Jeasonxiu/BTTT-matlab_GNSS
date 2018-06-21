function [jd] = date2jd(yy, mo, dd, h, m, s)
% 
%function [jd] = date2jd(yy,mo,dd,h,m,s)
%
% DO: Computes JD and MJD from calendar day and time
%
%       - JD(Julian Day), MJD(Modified Julian Day)
%       - This function is valid from 3/1/1900 through 2/28/2100
%
%   Fundamentals of Astrodynamics and Applications (David A. Vallado)
%   pp. 67-68 (Algorithm 2: Julian Date)
%
% Copyright: Kwan-Dong Park, August 2001 @CfA
%

jd = 367*yy - floor(7*(yy + floor((mo + 9)/12))/4) + floor(275*mo/9) + dd + 1721013.5 + ((s/60 + m)/60 + h)/24;
