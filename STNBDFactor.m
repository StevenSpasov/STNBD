function A = STNBDFactor(B,C,i)

% if A=L1 L2 ... D U_{n-1} ... U_1 this code extracts D for i=0
% extracts L_{n-1} for i=-1, L_{n-2} for i=-2, D_{n-1} for i=1, etc.

if i==0
    A=diag(diag(B));
else
    n=size(B,1);
    D=diag([ones(abs(i)-1,1); diag(C,i);1]);
    A=D+diag([zeros(abs(i)-1,1); diag(B,i)],sign(i));
end