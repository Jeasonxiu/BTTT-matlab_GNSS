function [out]= KPVeloTrans(in, to)
%
% DO: Transform coordinates determined relatively e.g. VRS 
%       - VRS surveys are currently based on the epoch time 2002.0
%       - KP velocity is taken from the SUWN site
%           V = [-0.0271 -0.0102 -0.0095]
%
% <in>  in: XYZ coordinates to translate from [meters; 1x3 row vector]
%       to: Epoch time to which coordinates should be translated [yrscale]
%               eg. 2005.5 for June 30, 2005
%
% <out> out: XYZ coordinates tranlated [meters]
%
% Copyright: Kwan-Dong Park, June 29, 2015
% Modified by : Jihye Won, Aug 16, 2014 @ Jipyong Space
%               Reference epoch was changed from 2000 to 2002

%% 함수 작성 이전 입출력 테스트
% in = [-3026443.3260  4067351.6671  3857212.8707];
% to = 2001.0;
% to = 2015.5;
%% 
KPVelo = [-0.0271 -0.0102 -0.0095];
dt = to - 2002.0;
out = in + KPVelo*dt;
