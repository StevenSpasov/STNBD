function V=QBVM(q,x)
% computes the q-Bernstein-Vandermonde matrix
n=length(x)-1;
V=zeros(n+1);
if isa(q,'sym') || isa(x,'sym')
    V=sym(V);
end
for i=1:n+1
    for j=1:n+1
        V(i,j)=qBVP(q,n,j-1,x(i));
    end
end

function s=qBVP(q,n,i,x)
% computes the q-BernsteinPolynomial b_{i,j}^n(x)
s=qbc(q,n,i)*x^i*prod(1-q.^(0:n-i-1)*x);
end

end