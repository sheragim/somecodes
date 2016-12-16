function ind = transient_remove(x)
% INPUT:
% x: is the signal which after some transient part the signal becomes stationary
% OUTPUT:
% ind: is the index of signal that after that the signal could be
% considered as a stationry signal.
% Author: sheragim parsa

N = length(x);
for kk = 1:floor(length(x)/2)
    y(kk) = var(x(kk+1:end))./(N-kk);
end
[~,ind] = min(y);
end