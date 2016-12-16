xx= load('data1.txt');
x1=xx(1:10^4);
% x1=xx(1:,log(1+length(xx))*length(xx)/1000)
plot(x1)

ind1 = transient_remove(x1)

figure(5);clf
plot(x1); hold on
plot([ind1,ind1],[min(x1)/2, max(x1)], '--r')
