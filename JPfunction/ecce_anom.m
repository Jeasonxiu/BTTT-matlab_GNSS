function [Ek]=ecce_anom(Mk,e,d)
Ek=Mk; % Ek 의 초기값 선정.
for i=1:d
    fE=Mk-Ek+e*sin(Ek);
    fpE=-1+e*cos(Ek);
    Ek_next=Ek-fE/fpE; % 뉴튼-랩슨법
    Ek=Ek_next; %구한 값을 다시 초기값으로 해줌.
end