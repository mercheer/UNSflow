
figure
set(gcf,'Units','Inches','Position',[1 1 3.25 2.6]);
set(gcf,'DefaultAxesFontName','Helvetica');
set(gcf,'DefaultTextFontName','Helvetica');
set(gcf,'DefaultAxesFontSize',8);
set(gcf,'DefaultTextFontSize',8);

set(gcf,'PaperUnits',get(gcf,'Units'));
pos = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 pos(3) pos(4)]);


theo=load('force_25pr.dat');
exp=load('exp_25pr.dat');
comp=load('comp_25pr.dat');

[ax, h1, h2] = plotyy(exp(:,1),exp(:,3)*2,...
		      theo(:,1),theo(:,2));

set(ax(1),'YColor','k');
set(ax(2),'YColor','k');

set(h1,'Color','k');

set(h2,'Color',[0.5 0.5 0.5]);
set(h2,'LineWidth',1.2);

axes(ax(1));

hold on
plot(comp(:,1)+0.1,comp(:,5),'-b','linewidth',0.1)
hcomp=plot(comp(1:18:end,1)+0.1,comp(1:18:end,5),'.b','linewidth',1.2)
plot(theo(:,1),smooth(theo(:,10)),'--','Color',[0 0.4 0],'linewidth',1.5) 


plot(2.2882,interp1(theo(:,1),theo(:,10),2.2882,'pchip'),'^k','markersize',4,'markerfacecolor','w','linewidth',1)
plot(4.0679,interp1(theo(:,1),theo(:,10),4.0679,'pchip'),'^k','markersize',4,'markerfacecolor','k','linewidth',1)


%hleg = legend('Experiment','CFD','LOM','Location','NorthWest');

%set(hleg,'Fontsize',7,'Box','off')
xlabel('t^*');
%ylabel('C_l')
  text(-0.5,1.5,'C_l');
set(gca,'XLim',[0 7],'YLim',[-1 3]);
set(gca,'XTick',[0:1:7],'YTick',[-1:1:3]);
% -grid on;

axes(ax(2));
set(gca,'YLim',[-10 30],'XLim',[0 7],'YTick',[-10:10:30]);
%h = text(4.4,41,'\alpha (right axis)');
%set(h,'Color',[0.5 0.5 0.5]);
%ylabel('\alpha, deg');
h = text(7.05,15,'\alpha (deg)');
set(h,'Color','k');
  

hold on;
plot(linspace(3.5,4.2,4),linspace(27.5,27.5,4),'.b','linewidth',1.5)
line([3.5 4.2],[27.5 27.5],'linewidth',0.1,'color','b')
     line([0.2 0.9],[27.5 27.5],'Color','k')
text(1.0,27.5,'Exp','Fontsize',7)
    
text(4.3,27.5,'CFD','Fontsize',7)
 line([5.2 5.9],[27.5 27.5],'Linestyle','--','Color',[0 0.4 0],'linewidth',1.5) 
text(6.0,27.5,'LDVM','Fontsize',7)

%plot(4.3,22,'ws','markersize',8,'markerfacecolor','w','markeredgecolor','w')

text(4.3,22,'\alpha','Fontsize',7,'background','w')
%set(gcf,'color','white')
%set(gcf,'inverthardcopy','off')
%line([5.7 6.4],[27.5 27.5],'Color',[0.5 0.5 0.5],'Linewidth',1.2)
%text(6.5,25,'\alpha','Fontsize',7)


set(ax(1),'Units','Inches');
set(ax(1),'Position',[0.35 0.45 2.5 2.0]);
set(ax(2),'Units','Inches');
set(ax(2),'Position',[0.35 0.45 2.5 2.0]);

print -depsc -loose ../figs/cl_25pr.eps

