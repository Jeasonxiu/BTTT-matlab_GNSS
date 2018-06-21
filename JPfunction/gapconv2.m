function [target,gsUBLOX] = gapconv2(VRS_text, gap)
%
%
%   input VRS_text : VRS logged text file's name
%                   ex) 'SDT1_VRS_16237.txt'
%   input gap : distance between reference receiver and target receiver
%                   ex) 0.43 (meter)
%
%   output target : N X 4 matric = [gs or sec, x, y, z]

                    
YY = str2num(VRS_text(10:11)); 
DOY = str2num(VRS_text(12:14)); 

VRS = load(VRS_text);
[dd, mm, yy] = ydoy2ymd(YY, DOY) ;     % 입력 연도와 DOY로 년 월 일 계산

for i =2:length(VRS(:,1))
    utc = num2str(VRS(i,1));            % VRS UTC
    h = str2num(utc(1)); m = str2num(utc(2:3)); s = str2num(utc(4:5));
    [gw, gs(i,1)] = date2gwgs(yy, mm, dd, h, m, s);     % UTC 2 Gps second
    
    A = VRS(i,5:7);
    B = VRS(i-1,5:7);
    dis(i,:) = norm(A-B); Dis(i,1) = VRS(i,1); Dis(i,2) = dis(i,1);
    target(i,1) = VRS(i,1);
    target(i,2:4) = [0,0,0]; C = [0, 0, 0];
    gsUBLOX(i,1) = round(gs(i,1));

    if dis(i) > 0.1                 % VRS 움직임 시작
        m = (dis(i) - gap) / dis(i); M(i,1) = m;
        n = gap / dis(i); M(i,2) = n;
        x = (m * A(1) + n * B(1)) ;
        y = (m * A(2) + n * B(2)) ;
        z = (m * A(3) + n * B(3)) ;
        target(i,1) = VRS(i,1);
        target(i,2:4) = [x,y,z]; C = [x, y, z];

    elseif dis(i) < gap && target(i-1,2) ~= 0       % target의 이전좌표가 존재하나 VRS 이동거리가 gap보다 작을때 VRS 차이값 만큼 target 좌표 생성
        AA(i-1,1) = 1;
        C = target(i-1,2) + A - B;
        target(i,2:4) = target(i-1,2:4) + A - B;

    end

end

i = 1;
udVRS = flipud(VRS);
udtarget = flipud(target);
while target(i,2) == 0
    i = i + 1;
    start = target(i,1);
end
for j = find(udVRS(:,1) == start) : length(VRS)-1
    invT = udVRS(j,1);
    A = udVRS(j,5:7);
    B = udVRS(j+1,5:7);
    C = A - B;
    udtarget(j+1,2:4) = udtarget(j,2:4) - C;
    udtarget(j+1,1) = invT;
end
target = flipud(udtarget);
gsUBLOX(1,1) = gsUBLOX(2,1) -1;
gsUBLOX(:,2:4) = target(:,2:4);
