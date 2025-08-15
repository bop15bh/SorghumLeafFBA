%% net CO2 uptake vs SWC taken from ogbaga et al., 2014 physiologia plantarum
clear
close all
co2=[ 13.75 11.75 7.5 2.5]
%swc=[ 100 40 25 15 10]
swc=[10 15 25 40 100]

full=16.5;

dec=abs((co2-full)/full)

doo=1-dec
doodoo=vertcat((dec')*100,100)
moon=doo*294.6979;
lolo=vertcat(moon',294.6979)
%popo=[53 159.1 249.2 291.7 350]
popo=[44.6512 133.9536 209.8606 245.5816 294.6979]
  
%drip=[40 270]
%soo=[10 100]
%po=10:20:270
 figure(1)            
            plot(swc,(doodoo),'color','r','LineWidth',10);
           % plot(soo,drip,'color','k','LineWidth',4);
          %  legend('data','v2','v3' ,'Location','Best');
            xlabel('Soil water content (%)','FontSize',60)
            ylabel('Net assimilation rate (%)','FontSize',60)
              set(gca,'LineWidth',8,'FontSize',60)
                           set(gcf, 'PaperUnits', 'inches'); 

                             x_width=16 ;y_width=16;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
  print('co2perc_vs_SWC','-djpeg','-loose');
  figure(2)            
            plot(swc,sort(lolo),'color','r','LineWidth',10);
            %hold on, drawnow 
           % plot(swc,popo,'color','k','LineWidth',10);
             hold on, drawnow 
           % plot(soo,drip,'color','k','LineWidth',4);
           % legend('data','Location','Best');
            xlabel('Soil water content (%)','FontSize',60)
            ylabel('Net assimilation rate {\mu}mol/g/hr','FontSize',60)
              set(gca,'LineWidth',8,'FontSize',60)
                           set(gcf, 'PaperUnits', 'inches'); 

                             x_width=16 ;y_width=16;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
  print('co2_vs_SWC','-djpeg','-loose');
