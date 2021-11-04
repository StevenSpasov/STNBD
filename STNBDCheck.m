function A = STNBDCheck(B,C)
n = length(B)-1;
A=zeros(n+1);
if isa(B,'sym') || isa(C,'sym')
    A=sym(A);
end
for i = -n:n
    A = A*STNBDFactor(B,C,i);
end
end