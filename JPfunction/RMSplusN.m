function [dNE dV d3] = RMS(dNEV)
%
%function [dNE dV d3] = RMS(dNEV)
%
%   RMS estimation form dNEV
%   
%   input dNEV : delta N, E, V
%   output dNE : delta North East
%          dV  : delta V
%          d3  : delta #d
%   
%   Example : [dNE dV d3] = RMS(dNEV)
%
%   coded by Joonseong Gim, Jan 18, 2016
%

dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dNE = rms(sqrt(dN.^2 + dE.^2));
dV = rms(sqrt(dV.^2));
d3 = rms(sqrt(dN.^2 + dE.^2 + dV.^2)) ;

