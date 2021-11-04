function L=LM(q,t)
% computes the Lupas matrix
n=length(t)-1;
L=zeros(n+1);
if isa(q,'sym') || isa(t,'sym')
    L=sym(L);
end
for i=1:n+1
    for j=1:n+1
        L(i,j)=LP(q,n,j-1,t(i));
    end
end

function l = LP(q,n,i,t)
%computes the Lupas polynomial
l = qbc(q,n,i);
l=l*q^((i^2-i)/2)*t^i*(1-t)^(n-i);
l=l/prod(1-t+q.^(1:n-1)*t);
end

end