function [] = lineCross(x,y,color,width)
% usage : [] = lineCross(x,y,color,width)
% ex) lineCross(0,0,'k',3);
% 축 그래프 그리는 함수
% 작성자 : 원지혜 (Mar 10, 2016)

if nargin == 3, width = 2; end
if nargin == 2, color='k'; width = 2; end


X_lim=get(gca,'xlim');
Y_lim=get(gca,'ylim');

hold on;
plot([x x],Y_lim,'-','color',color,'lineWidth',width);
plot(X_lim,[y y],'-','color',color,'lineWidth',width);

