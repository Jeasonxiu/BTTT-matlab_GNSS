function data = rtcm_get_type41(numofword, Zcount, time, f41)
%
%function [] = rtcm_get_type1()
%
%   Read type1(Differential GPS Corrections)
%
%   <input>
%           numofword : 해당 프레임의 총 워드수( header(6) )
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global
global now_word;
%% Function
numofbits = @(n_) 2^n_ - 1;     % eg) numofbits(4) = 15(10) = 1111(2)
%%
bin = []; S = [0 200 0 0 0 100]+100;
for m=1:numofword
    bin = [bin; dec2bin(now_word,32)];
    rtcm_next_word();
end

bin = reshape(bin(:,3:end-6)', 1, []);
sfactor= bin2dec(bin(1));
system = S( bin2dec(bin(1:4)) );
frontf = floor( length(bin(9:end)) / 36 ) * 36 + 8;
bin = reshape(bin(9:frontf), 36, [])';
[r, c] = size(bin);

[dum gs] = date2gwgs(time(1), time(2), time(3), time(4), time(5), time(6));

warning off
for i=1:r
    prn = system + bin2dec( bin(i,2:6) );
    prc = bin2dec( bin(i,21:end) );
    prc = bitcmp2( prc, 16 ) * ( 0.02  + sfactor*0.3 );
    fprintf(f41, '%8.1f %7d %5d %10.2f\n', gs, round(Zcount*0.6), prn, prc);
end
warning on