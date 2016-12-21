function [ Cost, gradCost ] = DiagonalScaleCost( a, xi, yi )
%The Loss (cost) function that  how smooth the interpolating funciton is.
inter_par= interpolateparametarization(xi.*repmat(a,1,size(xi,2)),yi,1);
w = inter_par{2};
% Cost =sum(w.^2)*norm(a)^3;
Cost =sum(w.^2);
% keyboard
if nargout>1
%The gradient of Loss (cost) function that indicates how smooth the interpolating funciton is.
inter_par= interpolateparametarization(xi.*repmat(a,1,size(xi,2)),yi,1);
w = inter_par{2};
Dw = Jacobi_a(  a,  xi, yi );
if length(a) ==1
    a = [1;a];
end
% gradCost =2*norm(a)^3 * Dw*w + 3 * norm(a) * norm(w)^2*a;
gradCost =2*Dw*w;
end
end


function [ Dw, Dv ] = Jacobi_a(  a, xi, yi )
% if length(a) ==1
%     a = [1;a];
% end
% gradCost =2*norm(a)^3 * Dw*w + 3 * norm(a) * norm(w)^2*a;
if nargin < 3
    a  = ones(size(xi,1),1);
end
inter_par= interpolateparametarization(xi.*repmat(a,1,size(xi,2)),yi,1);
 N = size(xi,2);n = size(xi,1); A = zeros(N,N); dA = zeros(N,N,n);
for ii = 1 : 1 : N
    for jj = 1 : 1 : N
%         A(ii,jj) = ((a.*(xi(:,ii) - xi(:,jj)))' * a.*(xi(:,ii) - xi(:,jj))) ^ (3 / 2);
            A(ii,jj) = (norm(a.*(xi(:,ii) - xi(:,jj))))^3;
%         for kk = 1 : 1 : n
            dA(ii,jj,:) = 3* a .* (xi(:,ii) - xi(:,jj)).^2 *norm(a.*(xi(:,ii) - xi(:,jj)));  
%         end
    end
end
V = [ones(1,N); xi];
A = [A V'; V zeros(n+1,n+1)];
Dw = []; Dv=[];
wv = [inter_par{2} ; inter_par{3}];
for kk=1:n
b{kk} = -[dA(:,:,kk) zeros(size(V')); zeros(size(V)) zeros(n+1,n+1)]*wv;
Dwv = A \ b{kk}; % solve the associated linear system
Dw = [Dw; Dwv(1:N)'];  
Dv = [Dv; Dwv(N+1:end)'];
end
% inter_par{1}=1;
% inter_par{2} = wv(1:N); inter_par{3} = wv(N+1:N+n+1); 
% inter_par{4}= xi;
end
