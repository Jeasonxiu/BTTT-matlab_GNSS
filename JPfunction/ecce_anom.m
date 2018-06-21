function [Ek]=ecce_anom(Mk,e,d)
Ek=Mk; % Ek �� �ʱⰪ ����.
for i=1:d
    fE=Mk-Ek+e*sin(Ek);
    fpE=-1+e*cos(Ek);
    Ek_next=Ek-fE/fpE; % ��ư-������
    Ek=Ek_next; %���� ���� �ٽ� �ʱⰪ���� ����.
end