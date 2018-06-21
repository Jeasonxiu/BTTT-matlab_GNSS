function [real, real_dg] = choiceP(point)
%
%function [real, real_dg] = choiceP(point)
%
%   Choice reference Point A or B
%   
%   input point : reference Point A or B
%
%   Example : [real, real_dg] = choiceP(point)
%
%   coded by Joonseong Gim, Jan 12, 2016
%

point = char(point);
switch point
    case 'A'
        real_x = -3041235.578;  % 대성디폴리스 옥상 A 지점 x
        real_y = 4053941.677;   % 대성디폴리스 옥상 A 지점 y
        real_z = 3859881.013;   % 대성디폴리스 옥상 A 지점 z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % 대성디폴리스 옥상 A 지점 gd
    case 'a'
        real_x = -3041235.578;  % 대성디폴리스 옥상 A 지점 x
        real_y = 4053941.677;   % 대성디폴리스 옥상 A 지점 y
        real_z = 3859881.013;   % 대성디폴리스 옥상 A 지점 z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % 대성디폴리스 옥상 A 지점 gd
    case 'B'
        real_x = -3041241.741;  % 대성디폴리스 옥상 B 지점 x
        real_y = 4053944.143;   % 대성디폴리스 옥상 B 지점 y
        real_z = 3859873.640;   % 대성디폴리스 옥상 B 지점 z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % 대성디폴리스 옥상 B 지점 gd
    case 'b'
        real_x = -3041241.741;  % 대성디폴리스 옥상 B 지점 x
        real_y = 4053944.143;   % 대성디폴리스 옥상 B 지점 y
        real_z = 3859873.640;   % 대성디폴리스 옥상 B 지점 z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % 대성디폴리스 옥상 B 지점 gd
end