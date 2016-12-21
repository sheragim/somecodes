clear all; close all;
load sample_data
a = [1;1]; n = size(xi,1);
xi = [0.5000    0.8000    0.5000    0.2000    0.5000; ...
    0.5000    0.5000    0.8000    0.5000    0.2000]

fun=  @(x)(x(1,:)-0.45).^2+0.001*(x(2,:)-0.45).^2;
yi = fun(xi);
%  keyboard
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed' );
options = optimoptions(options,'GradObj','on');
lb = [0;0]; ub = [n;n];   % No upper or lower bounds
fung = @(x)DiagonalScaleCost(x,xi,yi);
[a,fval] = fmincon(fung,x0,[],[],ones(1,n),n,lb,ub,[],options)

%%
% a=[1.8;0.2]
a=sqrt(a);
% a1=a0;  
% a=a1;
%  a=[1;1]
inter_par= interpolateparametarization(xi.*repmat(a,1,size(xi,2)),yi,1);
xv=0:0.05:1;
for ii=1:length(xv)
   %  for jj=1:length(xv)
  yr(ii)=fun([xv(ii);0.5]);
  yp(ii)=interpolate_val([xv(ii);0.5].*a,inter_par);
  yr1(ii)=fun([0.5;xv(ii)]);
  yp1(ii)=interpolate_val([0.5;xv(ii)].*a,inter_par);
     for jj=1:length(xv)
  Ur(ii,jj)=fun([xv(ii);xv(jj)]);
  Up(ii,jj)=interpolate_val([xv(ii);xv(jj)],inter_par);
     end
end

% figure(1)
% plot(xv,yr,'r-',xv,yp,'b--')
% 
% figure(2)
% plot(xv,yr1,'r-',xv,yp1,'b--')
figure(3)
hold on
subplot(2,1,1)
plot(xv,yr,'r-',xv,yp,'g--','linewidth',3)

% figure(4)
subplot(2,1,2)
hold on
plot(xv,yr1,'r-',xv,yp1,'g--','linewidth',3)
% %%
% % a=[2/5;1-2/5]
% for ii=.01:0.1:5
%     a=[2/ii;2-2/ii];
%      a=[1;1]
% %     a=a0
%   
% % a0=a;
% % a1=a0;  
% % a=a1;
% %  a=[1;1]
% inter_par= interpolateparametarization(xi.*repmat(a,1,size(xi,2)),yi,1);
% xv=0:0.05:1;
% for ii=1:length(xv)
%    %  for jj=1:length(xv)
%   yr(ii)=fun([xv(ii);0.5]);
%   yp(ii)=interpolate_val([xv(ii);0.5].*a,inter_par);
%   yr1(ii)=fun([0.5;xv(ii)]);
%   yp1(ii)=interpolate_val([0.5;xv(ii)].*a,inter_par);
%      for jj=1:length(xv)
%   Ur(ii,jj)=fun([xv(ii);xv(jj)]);
%   Up(ii,jj)=interpolate_val([xv(ii);xv(jj)],inter_par);
%      end
% end
% a
% figure(4)
% hold on
% subplot(2,1,1)
% plot(xv,yr,'r-',xv,yp,'b--')
% 
% % figure(4)
% subplot(2,1,2)
% hold on
% plot(xv,yr1,'r-',xv,yp1,'b--')
% end
% 
