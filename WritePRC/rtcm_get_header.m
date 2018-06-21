function header = rtcm_get_header()
%
%function [] = rtcm_get_header()
%
%   Read header
%
%   Copyright: taeil Kim, February 20, 2015@INHA University

%% Global
global now_word;
global sync;
%% Function
numofbits = @(n_) 2^n_ - 1;     % eg) numofbits(4) = 15(10) = 1111(2)
%%
%--- preamble check -------------------------------------------------------
if ~rtcm_preamble(now_word)
    sync = 0;
    return;
end
%--- first word of each message -------------------------------------------
preamble = bitand( bitshift( now_word, -22 ), numofbits(8) );
msgtype  = bitand( bitshift( now_word, -16 ), numofbits(6) );
stationID= bitand( bitshift( now_word, -6 ), numofbits(10) );

rtcm_next_word();
%--- second word of each message ------------------------------------------
M_Z_count= bitand( bitshift( now_word, -17 ),numofbits(13) );
seqno    = bitand( bitshift( now_word, -14 ), numofbits(3) );
numofword= bitand( bitshift( now_word, -9  ), numofbits(5) );
sv_health= bitand( bitshift( now_word, -6  ), numofbits(3) );

rtcm_next_word();
%--- ¸¶¹«¸® ---------------------------------------------------------------
header = [preamble msgtype stationID M_Z_count seqno numofword sv_health];