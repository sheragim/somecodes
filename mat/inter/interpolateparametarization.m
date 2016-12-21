function inter_par= interpolateparametarization(xi1,yi1,inter_method,interpolate_index)
% keyboard
global xi yi y0 w 
xi=xi1; yi=yi1;
n=size(xi,1);
% keyboard
% polyharmonic spline interpolation
if inter_method==1
    N = size(xi,2); A = zeros(N,N);
for ii = 1 : 1 : N
    for jj = 1 : 1 : N
        A(ii,jj) = ((xi(:,ii) - xi(:,jj))' * (xi(:,ii) - xi(:,jj))) ^ (3 / 2);
    end
end
%keyboard
V = [ones(1,N); xi];
A = [A V'; V zeros(n+1,n+1)];
wv = A \ [yi.'; zeros(n+1,1)]; % solve the associated linear system
inter_par{1}=1;
inter_par{2} = wv(1:N); inter_par{3} = wv(N+1:N+n+1); 
inter_par{4}= xi;
end

% new 
if inter_method == 7 || inter_method == 8
 keyboard
a = ones(size(xi,1),1);
lambda = 1;
% keyboard
for i=1:20
% a0=[1;1];  [inter_par, a]  = Scale_interpar( xi,yi,a0, lambda); % method1
[inter_par,a]  = Scale_interpar( xi,yi,a, lambda); %method2
inter_par{1}=inter_method;
lambda = lambda/2;

for jj=1:numel(yi)
ygps(jj)= interpolate_val(xi(:,jj),inter_par);
% equation 19 in MAPS
deltaPx  = abs(ygps(jj)-yi(jj));
DeltaFun = abs(yi(jj)-y0);
% keyboard
if deltaPx/DeltaFun > 0.1
break;
elseif jj==numel(yi)
    return;
end

end


end
if inter_method==8
%      keyboard
epsFun = yi-y0;
inter_par{8}=epsFun;
end


end





end

%% Scaled Polyharmonic Spline
function [ inter_par,a ] = Scale_interpar( xi,yi,a0, lambda0)
% This function is for spinterpolation and finds the scaling factor for
% polyharmonic spline interpolation
global lambda
lambda = lambda0;
n=size(xi,1);
if nargin <3
a0=ones(n,1);
lambda =1;
end
% keyboard

% options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed' );
options = optimoptions(@fmincon,'Algorithm','sqp');
options = optimoptions(options,'GradObj','on');
lb = zeros(n,1); ub = ones(n,1)*n;   % No upper or lower bounds
% for ii=1:20
fung = @(a)DiagonalScaleCost(a,xi,yi);
% keyboard
[a,fval] = fmincon(fung,a0,[],[],ones(1,n),n,lb,ub,[],options);
% end
[ff,gf,inter_par] = DiagonalScaleCost(a,xi,yi);

end
%
function [ Cost, gradCost, inter_par ] = DiagonalScaleCost( a, xi, yi)
%The Loss (cost) function that  how smooth the interpolating funciton is.
global lambda
% keyboard
inter_par= interpolateparametarization_scaled(xi,yi,a,1, lambda);
w = inter_par{2};
Cost =sum(w.^2);
% keyboard
if nargout>1
%The gradient of Loss (cost) function that indicates how smooth the interpolating funciton is.
Dw = inter_par{5};
gradCost =2*Dw*w;
end
inter_par{7}=a;
inter_par{1}=7;
%%%%%%%%%%%%%
end
%
function inter_par= tabehSP_SCALED(xi1,yi1,a, inter_method,lambda,interpolate_index)
global xi yi y0 w 
H= diag(a);
if nargin < 4
lambda = 0;
%lambda = 1e-3;
end
xi= xi1;
yi=yi1;
n=size(xi,1);
% keyboard
% polyharmonic spline interpolation
if inter_method==1
    N = size(xi,2); A = zeros(N,N);
for ii = 1 : 1 : N
    for jj = 1 : 1 : N
        A(ii,jj) = ((xi(:,ii) - xi(:,jj))' *H* (xi(:,ii) - xi(:,jj)))^(3 / 2);
        dA(ii,jj,:) =3/2.* (xi(:,ii) - xi(:,jj)).^2 *  ((xi(:,ii) - xi(:,jj))' *H* (xi(:,ii) - xi(:,jj)))^(1/2);
    end
end
% keyboard
V = [ones(1,N); xi1];
A = A + eye(N)*lambda;
A = [A V'; V zeros(n+1,n+1)];
%%%wv = pinv(A)* [yi.'; zeros(n+1,1)]; % solve the associated linear system
% keyboard
wv = A\[yi.'; zeros(n+1,1)];
%
% bb=[yi.'; zeros(n+1,1)], AA= A*A',  WV = AA\bb,   XX = A'*WV
% err = A*XX-[yi.'; zeros(n+1,1)]
inter_par{1}=1;
inter_par{2} = wv(1:N); inter_par{3} = wv(N+1:N+n+1); 
inter_par{4}= xi1;
% calculating the gradient  
Dw = []; Dv=[];
for kk=1:n
b{kk} = -[dA(:,:,kk) zeros(size(V')); zeros(size(V)) zeros(n+1,n+1)]*wv;
% Dwv = A \ b{kk}; % solve the associated linear system
% keyboard
Dwv = pinv(A)* b{kk}; % solve the associated linear system
Dw = [Dw; Dwv(1:N)'];  
Dv = [Dv; Dwv(N+1:end)'];
end
inter_par{5} = Dw;
inter_par{6} = Dv; 
inter_par{7} = a; 
end

end
