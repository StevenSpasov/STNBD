function A = STNBDCheck(B,C)
n = length(B)-1;
A = STNBDFactor(B,C,-n);
for i = -(n-1):n
    A = A*STNBDFactor(B,C,i);
end
end
