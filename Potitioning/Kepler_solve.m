function [E] = Kepler_solve(M, e, iter);
%
%
% Equations for Newton's Method
% 
% 
E =pi/6;
for i = 1:iter;
    E_new = E - (M - E + e*sin(E))/(e*cos(E) - 1);
    E = E_new;
end

