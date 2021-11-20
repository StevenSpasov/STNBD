function STNBDFull(B,C)
% displays the full factorization from STN
n = length(B)-1;
for i = -n:n
  disp (STNBDFactor(B,C,i));
end
