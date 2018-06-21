function [h24] = gs2h24(gs)
%function [h24] = gs2h24(gs)
%
% DO: Converts GPS Week-Seconds to 24H 
%
% <input> gs: gps week seconds
%
% <output> h24: hour scale
%
% Copyright: Kwan-Dong Park, October 2013 @LDEO
%
h24 = mod(gs,86400)/3600;
