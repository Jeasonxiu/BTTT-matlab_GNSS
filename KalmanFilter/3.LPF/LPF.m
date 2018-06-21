function xlpf = LPF(x)
%
%
persistent prevX
persistent firstRun


if isempty(firstRun)
  prevX = x;
  
  firstRun = 1;  
end


alpha = 0.5;
xlpf = alpha*prevX + (1 - alpha)*x;

prevX = xlpf;