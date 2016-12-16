function [y,x1,index] = mindistance(x,xi)
% calculates the minimum distance from all the existing points
% xi all the previous points
% x the new point
%Author: sheragim parsa

y=inf;
for i=1:size(xi,2)
   y1=norm(x-xi(:,i));
    if y1<y
        y=y1;
        x1=xi(:,i);
        index=i;
    end
end  
    


end

