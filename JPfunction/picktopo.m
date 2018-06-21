function [pick,  value, direction] = picktopo(topo)
%
%function [] = picktopo(topo)
%
%   Read the topo matrix and find a north or east pick
%   
%   input topo : topo matrix
%
%   output pick : absolute pick value
%   output value : [pick value, north or east]
%
%   Example : [pick, direction] = picktopo(topo)
%
%   coded by Joonseong Gim, Jan 14, 2016
%

max_north = max(topo(:,1)); % north error�� �ִ밪
max_east = max(topo(:,2));  % east error�� �ִ밪
min_north = min(topo(:,1)); % north error�� �ִ밪
min_east = min(topo(:,2));  % east error�� �ִ밪

abs_pick = [abs(max_north) abs(max_east) abs(min_north) abs(min_east)];
pick = max(abs_pick);

if pick == abs(max_north)
    value = max_north; direction = 'max_north';
elseif pick == abs(max_east)
    value = max_east; direction = 'max_east';
elseif pick == abs(min_north)
    value = min_north; direction = 'min_north';
else
    value = min_east; direction = 'min_east';
end
