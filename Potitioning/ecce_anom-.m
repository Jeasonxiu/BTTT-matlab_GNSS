function E = ecce_anom(M,e,n_iter)

E = M + e.*sin(M);
for i=1:n_iter
  E = M + e.*sin(E);
end


