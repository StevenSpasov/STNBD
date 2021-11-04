%function [M,C] = STNBDL(q,x)
%computes the bidiagonal decomposition of a potentially singular Lupas matrix

function [M,C] = STNBDL(q,x)

[m1,n1]=size(x);

n=n1-1;
M=zeros(n1);
A=zeros(n,n1);
C=ones(n+1);

if class(x) == 'sym'
    M=sym(M);
    A=sym(A);
    C = sym(C);
end

for i=1:n
    for j=1:n1
	A(i,j)=(1-x(j))+q^(i-1)*x(j);
    end
end

v=prod(A);

% first loop: m_ij
for i=2:n+1
    auxM=((1-x(i))^n*v(i-1))/((1-x(i-1))^(n+1)*v(i));
    M(i,1)=(1-x(i-1))*auxM;
    for j=1:i-2
        auxM=auxM*((1-x(i-1)))/(1-x(i));
        M(i,j+1)=(1-x(i-j-1))*auxM;
    end
end

 for i=2:n+1
    for j=1:i-2
      C(i,j+1)=x(i-1)-x(i-j-1);
   end
end   

% second loop: m~_ij
for j=1:n
    cj=x(j)/(1-x(j)); 
    for i=j+1:n+1
        M(j,i)=cj*q^(i-2)*qi(q,n-i+2)/qi(q,i-1);
    end
end

% third loop: the diagonal entries p_ii
r=1;
M(1,1)=(1-x(1))^n/v(1);
for i=1:n
    r=(qi(q,n-i+1)/qi(q,i))*(1/(1-x(i)))*r;
    M(i+1,i+1)=r*q^(i*(i-1)/2)*(1-x(i+1))^(n-i)/v(i+1);
end
for i=2:n+1
    M(n+1,i)=M(n+1,i)*prod(x(n+1)-x(n-i+2:n));
end