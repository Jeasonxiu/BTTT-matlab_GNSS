function [s] = Get2ENDsp3(fid_sp3)
%
%function [] = Get2ENDsp3(sp3_file)
% 
% DO: Get to the last part of SP3 Header
%
% <input> fid_sp3: SP3 File Handler
%
% Copyright: Kwan-Dong Park@LDEO October 25, 2013
%

%% SP3 파일의 헤더를 #, +, %, /로 판단해서 제거
ready = 0;
while ready == 0
    s = fgets(fid_sp3);
    switch s(1,1)
        case '#'
        case '+'
        case '%'
        case '/'
        otherwise
            ready = 1;  
    end
end
%% 한 줄 되돌리는 방법이 없어서 출력으로 s를 뱉아 내야함 ㅠ
