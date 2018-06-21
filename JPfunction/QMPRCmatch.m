% function [QM PRC] = QMPRCmatch(QM, PRC)
%
%function [] = picktopo(topo)
%
%   Read the topo matrix and find a north or east pick
%   
%   input topo : topo matrix
%
%   output pick : absolute pick value
%   output value : [pick value, north or east]
%
%   Example : [pick, direction] = picktopo(topo)
%
%   coded by Joonseong Gim, Jan 14, 2016
%
PRCgs = unique(PRC(:,1));
QMgs = unique(QM(:,1));

if length(QMgs) >= length(PRCgs)
    gs = PRCgs;
    for i = 1:length(gs)
        gsindex(i,1) = find(QM(:,1) == PRC(i,1));
    end
else
    gs = QMgs;
end

    