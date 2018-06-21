function [dis] = disbw2deg(Lat1, Lat2, Lon1, Lon2)
%
% function [dis] = disbw2deg(Lat1, Lat2, Lon1, Lon2)
%
%   Finding the distance in meters between 2 Decimal degrees Coordinates
%   
%   input Lat1, Lat2, Lon1, Lon2 : Decimal degree
%
%   Example : [dis] = disbw2deg(36.7, 36.71,126.1, 126.2)
%
%   coded by Joonseong Gim, Feb 19, 2016


PolarRadius = 6356750; EquatorialRadius = 6378200;
La_dis = abs(Lat1 - Lat2)*(pi/180)*PolarRadius;
Lo_dis = abs(Lon1 - Lon2)*(pi/180)*(EquatorialRadius*cos(((Lat1 + Lat2)/2)*pi/180));
dis =sqrt(La_dis^2 + Lo_dis^2);

            