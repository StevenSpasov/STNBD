function V = VM(x)
%computes the Vandermonde matrix
n=length(x)-1;
V=zeros(n+1);
if isa(x,'sym')
   V=sym(V);
end
for i = 1:n+1
    for j = 1:n+1
        V(i,j)=x(i)^(j-1);
    end
end