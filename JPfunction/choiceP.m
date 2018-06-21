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
        real_x = -3041235.578;  % �뼺�������� ���� A ���� x
        real_y = 4053941.677;   % �뼺�������� ���� A ���� y
        real_z = 3859881.013;   % �뼺�������� ���� A ���� z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % �뼺�������� ���� A ���� gd
    case 'a'
        real_x = -3041235.578;  % �뼺�������� ���� A ���� x
        real_y = 4053941.677;   % �뼺�������� ���� A ���� y
        real_z = 3859881.013;   % �뼺�������� ���� A ���� z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % �뼺�������� ���� A ���� gd
    case 'B'
        real_x = -3041241.741;  % �뼺�������� ���� B ���� x
        real_y = 4053944.143;   % �뼺�������� ���� B ���� y
        real_z = 3859873.640;   % �뼺�������� ���� B ���� z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % �뼺�������� ���� B ���� gd
    case 'b'
        real_x = -3041241.741;  % �뼺�������� ���� B ���� x
        real_y = 4053944.143;   % �뼺�������� ���� B ���� y
        real_z = 3859873.640;   % �뼺�������� ���� B ���� z
        real = [real_x real_y real_z];  
        real_dg = xyz2gd(real); % �뼺�������� ���� B ���� gd
end