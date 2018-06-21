function Skymap(cut_ele)

% -- skyplot 을 그리는 함수
%
%     Usage:       skyplot(Az, El, cut_ele) 
%
%          eg:       skyplot(120, 30, 10)
% 
% -- 방위각, 고도각, 임계고도각을 degree로 입력한다(개별적 수치 또는 배열형태).
% -- Hye-In, Kim

hold on;

circlea = 0:pi/30:2*pi;
cx  = cos(circlea);
cy  = sin(circlea);

% 30도, 60도, 90도를 나타내는 원을 실선으로 그린다(고도각을 나타내는 부분).
for i= [30 60 90]
    plot(cx*i,cy*i,'-','color','k','linewidth',1);
end

% 15도, 45도, 75도를 나타내는 원을 점선으로 그린다(고도각을 나타내는 부분).
for i=[15 45 75]
    plot(cx*i,cy*i,':','color','k','linewidth',1);
end

% 임계고도각을 파란색 실선으로 나타낸다. 
plot(cx*(cut_ele-90),cy*(cut_ele-90),'-','color','b','linewidth',1);

lenmax = 90;

% 원을 12등분으로 나눈다(방위각을 나타내는 부분).
circleTick = (1:6)*pi/6;
cosct = cos(circleTick); 
sinct = sin(circleTick);
cax = [-cosct; cosct];
say = [-sinct; sinct];
plot(lenmax*cax,lenmax*say,'-','color','k','linewidth',1);

%filename = strcat('[ "', filename,'" Invisibility Sky Plot ]');
% title(site, 'position',[-60 95]);

% 고도각을 나타내는 텍스트를 삽입한다(30, 60, 90, 120도 ....).
rlen = 1.06*lenmax;
for i = 1:length(circleTick) 
    ticm1 = int2str(i*30);
    ticm2 = int2str(180+i*30);
    
    % 방위각 360도를 나타내는 부분은 북쪽을 나타내는 'N'으로 표기한다.
    if ticm2 == '360'
        ticm2 ='N';
    end
    
    text( rlen*sinct(i), rlen*cosct(i),ticm1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','horizontalalignment','center');     
    text(-rlen*sinct(i),-rlen*cosct(i),ticm2,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','horizontalalignment','center');
end

axis('equal');
axis('off');
