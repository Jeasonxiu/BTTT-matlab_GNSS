function [dNE dV d3] = MEAN(dNEV)
%
%function [dNE dV d3] = MEAN(dNEV)
%
%   MEAN estimation form dNEV
%   
%   input dNEV : delta N, E, V
%   output dNE : delta North East
%          dV  : delta V
%          d3  : delta #d
%   
%   Example : [dNE dV d3] = MEAN(dNEV)
%
%   coded by Joonseong Gim, Jan 18, 2016
%

dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
% dNE = mean(sqrt(dN.^2 + dE.^2));
% dV = mean(sqrt(dV.^2));
% d3 = mean(sqrt(dN.^2 + dE.^2 + dV.^2)) ;
dNE = mean(sqrt(dNEV(:,1).^2 + dNEV(:,2).^2));
dV = mean(dNEV(:,3));
d3 = mean(sqrt(dNEV(:,1).^2 + dNEV(:,2).^2 + dNEV(:,3).^2)) ;

% for i = 1:96
%     d3(i,1) = (sqrt(dN(i)^2 + dE(i)^2 + dV(i)^2)) ;
%     d3(i,2) = (sqrt(dNEV(i,1)^2 + dNEV(i,2)^2 + dNEV(i,3)^2)) ;
% end
