function [bin] = rtcm_get_bin(numofword)

%% Global
global now_word;
%%
bin = [];
for m=1:numofword
    bin = [bin; dec2bin(now_word,32)];
    rtcm_next_word();
end
1;