function [res] = nck(n,k)
%self defined nchoosek
if n<0
    error('n<0');
end
if k<0
    error('k<0');
end
if k>n
    error('k>n')
end
res = 1;
if n>0 && k==0
    return;
end

for i=n-k+1:n
    res = res*i;
end
for i=1:k
    res = res / i;
end

end

