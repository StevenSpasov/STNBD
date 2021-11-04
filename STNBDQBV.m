% function [M,C] = STNBDQBV(q, x)

% Computes the bidiagonal decomposition of a potentially singular q-Bernstein-Vandermonde matrix

function [M,C] = STNBDQBV(q,x)
n=length(x)-1;
M=zeros(n+1);
C=ones(n+1);
if class(x) == 'sym'
    M=sym(M);
    C=sym(C);
end
% first loop: m_ij
for i=0:n
    for j=0:i-1
        M(i+1,j+1)=(1-q^(n-j)*x(i-j)) ...
            *prod(1-q.^(0:n-j-1)*x(i+1)) ...
            /prod(1-q.^(0:n-j)*x(i));
    end
end
for i=2:n+1
    for j=1:i-2
        C(i,j+1)=x(i-1)-x(i-j-1);
    end
end
% second loop: m~_ij
for i=0:n
    for j=0:i-1
       M(j+1,i+1)=qi(q,n-i+1)/qi(q,i)*x(j+1)/(1-q^(n-i)*x(j+1));
       M(j+1,i+1)=M(j+1,i+1)*prod(1-q^(n-i+1)*x(1:j)) ...
            /prod(1-q^(n-i)*x(1:j));
    end
end
% third loop: the diagonal entries p_ii
for i=0:n
    M(i+1,i+1)=qbc(q,n,i);
    M(i+1,i+1)=M(i+1,i+1)*prod(1-q.^(0:n-i-1)*x(i+1))/prod(1-q^(n-i)*x(1:i));
end
for i=2:n+1
    M(n+1,i)=M(n+1,i)*prod(x(n+1)-x(n-i+2:n));
end
end