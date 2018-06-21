function DATA = getRTCM(filename, mode)
%
%function DATA = getRTCM(filename, mode)
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global 초기화
global now_word;    now_word = uint32(1);
global sync;        sync = 0;
global start;       start = 0;
global glo_hour; glo_hour=0;
global glo_min; glo_min=0;
global glo_sec; glo_sec=0;
%% 초기 설정
fid = fopen(filename);
f41 = fopen('message41.txt', 'w');
DATA = cell(86400,1);
iter = 0;
%%
while ~feof(fid)
    iter = iter + 1;
    switch(mode)
        case 'jprt'
            DATA{iter,1} = rtcm_mode_jprt(fid, f41);xss
        case 'ntrip'
            DATA{iter,1} = rtcm_mode_ntrip(fid);
        otherwise
            error('error');
    end    
    fprintf('%d %d %d %d \n', iter,glo_hour,glo_min,glo_sec);
    
end
%% 마무리
DATA = DATA(1:iter);
fclose(fid);
fclose(f41);