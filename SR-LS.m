function [X] = GTRS(D, dim, pars)
PP=pars.PP;
n=size(D,2);
m=pars.m;
allnodes.num=m;
allnodes.dim=dim;
allnodes.loc=pars.PP(:,1:pars.m)';
dis=sqrt(D);

X(1:m,:)=pars.PP(:,1:pars.m)';
for i=m+1:n
    X(i,:)=SRLS(dis(i,1:m),allnodes);
end
X=X';

RMSD = sum(sum((X(:,m+1:n)-PP(:,m+1:n)).^2))/(n-m);
RMSD = sqrt(RMSD);

markersize=5;
if (pars.plotyes)
plotEMBED(X, dim, markersize, m, PP, RMSD);
end


% function plotEMBED(X, dim, markersize, m, PP, RMSD)
% 
% if nargin < 3
%     markersize = 5;
%     m = 0;
%     np = 0;
%     PP = [];
% end
% if nargin < 4
%     m = 0;
%     np = 0;
% end
% if nargin >= 5
%    np = size(PP, 2);
% end
% 
% if nargin < 6
%     RMSD = 0;
% end
% 
% if (dim == 2)
%  % plot the recovered anchors and sensors
%   figure
%  plot(X(1, :), X(2, :), '*r', 'markersize',markersize);
%  hold on
%  
%  % plot the anchors
%   if (m > 0) 
%     h = plot(PP(1,1:m),PP(2,1:m),'d','markersize',markersize); 
%     set(h,'linewidth',3,'color','b');
%   end
%  % plot the sensors
%  if (np > m) 
%      plot(PP(1, (m+1):end), PP(2, (m+1):end), 'ok', 'markersize',markersize); % plot the position of AS in R^2
%  % link the known anchor/sensors and the recoverd points
%      xy = [X(:, 1:np)'; PP'];
%      I  = (1:np)'; J = (np+1:2*np)';
%      a  = ones(np,1);
%      E  = sparse(I,J,a, 2*np, 2*np);
%      gplot(E, xy);
%  end
%  %xie add%    axis([-0.6 0.6 -0.6 0.6]);
%    xlabel(['EMBED: RMSD = ', sprintf('%4.2e', RMSD)]);
% end
% 
% % plot in R^3
% % 
% if (dim == 3)
%     figure
%     plot3(X(1, :), X(2, :), X(3, :), '*r', 'markersize',markersize);
%     hold on; grid on;
%    if (m > 0)
%        h = plot3(PP(1,1:m),PP(2,1:m),PP(3,1:m),'d','markersize',markersize); 
%        set(h,'linewidth',3,'color','b');     
%    end  
% 
%     if (np > m)
%        plot3(PP(1,(m+1):end),PP(2,(m+1):end),PP(3,(m+1):end),'ok','markersize',markersize);
%        hold on; grid on;
%        plot3([X(1,(m+1):end); PP(1,(m+1):end)],[X(2,(m+1):end); PP(2,(m+1):end)],[X(3,(m+1):end); PP(3,(m+1):end)],'b');
%        legend('GTRS',3);
%     end
%      
%       %axis('square'); 
%       %axis(0.6*BoxScale*[-1,1,-1,1,-1,1])
%       %xlabel(['GTRS: RMSD = ', sprintf('%4.2e', RMSD)]);
%       xlabel('x轴');
%       ylabel('y轴');
%       zlabel('z轴');
%      % pause(0.1);
%      legend('GTRS','参考节点','真实值',3);
%       hold off
% end
% 



return