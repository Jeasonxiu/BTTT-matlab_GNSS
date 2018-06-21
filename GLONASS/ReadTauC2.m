function [tau_c] = ReadTauC2(GLON_file)

% function [tau_c] = ReadTauC(GLON_file)
% GLONASS 방송궤도력 파일에서 TauC값을 읽어들임
% 2014.10.12 김미소

% GLON_file = 'brdm2200.15p';
fid = fopen(GLON_file, 'r');

ready = 0;

while ~ready
    s = fgetl(fid);
    if length(s) <= 80
        if s(01 : 04) == 'GLUT'
            tau_c = str2num(s(06 : 22));
            break;
        end
    end
end