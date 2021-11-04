function s=qbc(q,n,r)

% function s=qbc(q,n,r)
% computes the q-binomial coefficient

s=1;
for i=1:r
    s=s*qi(q,n-i+1)/qi(q,i);
end
end