function s=qi(q,r)

% function s=qi(q,r)
% computes q-integer

if q==1 
    s=r;
else
    s=(1-q^r)/(1-q);
end
end