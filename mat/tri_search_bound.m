function [xm cse indm ym] = tri_search_bnd(inter_par,xiT)

global n bnd1 bnd2 y0
% Update global interpolation

tri=delaunayn(xiT.');
ym=inf; cse=2; 
    for ind=1:size(tri,1)
      if min(tri(ind,:))>n+1
      [xc,R2]=circhyp(xiT(:,tri(ind,:)), n);
       x=xiT(:,tri(ind,:))*ones(n+1,1)/(n+1);
              %%% skipping these functions
          x=A_K_Search(x,inter_par,xc,R2);
          y=(interpolate_val(x,inter_par)-y0)/(R2-norm(x-xc)^2);
              %%%
      if (y<ym)
          ym=y; xm=x; indm=ind;
      end
      %end
      end
    end
    
end




