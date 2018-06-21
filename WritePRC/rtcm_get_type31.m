function data = rtcm_get_type31(numofword)
%
%function [] = rtcm_get_type31()
%
%   Read type31(Differential GLONASS Corrections)
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
k = floor( double(numofword) / 5);
j = mod(numofword, 5);
count = k*3+j/2;
data = zeros(count, 7);
i=1;

for m=1:k
    data(i,1) = bitand( bitshift( now_word, -29 ), numofbits(1) );
    data(i,2) = bitand( bitshift( now_word, -27 ), numofbits(2) );
    data(i,3) = bitand( bitshift( now_word, -22 ), numofbits(5) );
    data(i,4) = bitand( bitshift( now_word, -6  ), numofbits(16));
    rtcm_next_word();
    data(i,5) = bitand( bitshift( now_word, -22 ), numofbits(8) );
    data(i,6) = bitand( bitshift( now_word, -21 ), numofbits(1) );
    data(i,7) = bitand( bitshift( now_word, -14 ), numofbits(7) );
    i=i+1;
    data(i,1) = bitand( bitshift( now_word, -13 ), numofbits(1) );
    data(i,2) = bitand( bitshift( now_word, -11 ), numofbits(2) );
    data(i,3) = bitand( bitshift( now_word, -6  ), numofbits(5) );
    rtcm_next_word();
    data(i,4) = bitand( bitshift( now_word, -14 ), numofbits(16));
    data(i,5) = bitand( bitshift( now_word, -6  ), numofbits(8) );
    rtcm_next_word();
    data(i,6) = bitand( bitshift( now_word, -29 ), numofbits(1) );
    data(i,7) = bitand( bitshift( now_word, -22 ), numofbits(7) );
    i=i+1;
    data(i,1) = bitand( bitshift( now_word, -21 ), numofbits(1) );
    data(i,2) = bitand( bitshift( now_word, -19 ), numofbits(2) );
    data(i,3) = bitand( bitshift( now_word, -14 ), numofbits(5) );
    data(i,4) = bitand( bitshift( now_word, -6  ), numofbits(8) );
    rtcm_next_word();
    data(i,4) = bitshift(data(i,4), 8) + ...
                bitand( bitshift( now_word, -22 ), numofbits(8) );
    data(i,5) = bitand( bitshift( now_word, -14 ), numofbits(8) );
    data(i,6) = bitand( bitshift( now_word, -13 ), numofbits(1) );
    data(i,7) = bitand( bitshift( now_word, -6  ), numofbits(7) );
    rtcm_next_word();
    i=i+1;
end

if j==2
    data(i,1) = bitand( bitshift( now_word, -29 ), numofbits(1) );
    data(i,2) = bitand( bitshift( now_word, -27 ), numofbits(2) );
    data(i,3) = bitand( bitshift( now_word, -22 ), numofbits(5) );
    data(i,4) = bitand( bitshift( now_word, -6  ), numofbits(16));
    rtcm_next_word();
    data(i,5) = bitand( bitshift( now_word, -22 ), numofbits(8) );
    data(i,6) = bitand( bitshift( now_word, -21 ), numofbits(1) );
    data(i,7) = bitand( bitshift( now_word, -14 ), numofbits(7) );
    rtcm_next_word();
elseif j==4
    data(i,1) = bitand( bitshift( now_word, -29 ), numofbits(1) );
    data(i,2) = bitand( bitshift( now_word, -27 ), numofbits(2) );
    data(i,3) = bitand( bitshift( now_word, -22 ), numofbits(5) );
    data(i,4) = bitand( bitshift( now_word, -6  ), numofbits(16));
    rtcm_next_word();
    data(i,5) = bitand( bitshift( now_word, -22 ), numofbits(8) );
    data(i,6) = bitand( bitshift( now_word, -21 ), numofbits(1) );
    data(i,7) = bitand( bitshift( now_word, -14 ), numofbits(7) );
    i=i+1;
    data(i,1) = bitand( bitshift( now_word, -13 ), numofbits(1) );
    data(i,2) = bitand( bitshift( now_word, -11 ), numofbits(2) );
    data(i,3) = bitand( bitshift( now_word, -6  ), numofbits(5) );
    rtcm_next_word();
    data(i,4) = bitand( bitshift( now_word, -14 ), numofbits(16));
    data(i,5) = bitand( bitshift( now_word, -6  ), numofbits(8) );
    rtcm_next_word();
    data(i,6) = bitand( bitshift( now_word, -29 ), numofbits(1) );
    data(i,7) = bitand( bitshift( now_word, -22 ), numofbits(7) );
    rtcm_next_word();
end