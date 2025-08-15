clear
close all
% FVA pipeline
%addpath /mnt/home/holla293/Documents/cobratoolbox
%    initCobraToolbox;

%changeCobraSolver('glpk')
%load('SC_constrained_unblocked.mat');
%load('constrained_leaf_noredundancies.mat')
%load('300deadsMar21.mat')
load('Leaf_balanced_FINAL0625.mat')
model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 10, 'l');
model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 40, 'l');
%model = changeRxnBounds(model,'D-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN_2[M]', 0, 'l');
%model = changeRxnBounds(model,'ATR_L-ASPARTATE_[cb]_[cm]', 100, 'u');
model=removeRxns(model,{'4.1.1.32-RXN[M]','4.1.1.32-RXN[B]'})
form=' ADP[cm] + PHOSPHO-ENOL-PYRUVATE[cm] + PROTON[cm] -> ATP[cm] + PYRUVATE[cm] ';
model=addReaction(model,'PEPDEPHOS-RXN_1[M]',form,[],0,0,1000);
form=' ADP[cb] + PHOSPHO-ENOL-PYRUVATE[cb] + PROTON[cb] -> ATP[cb] + PYRUVATE[cb] ';
model=addReaction(model,'PEPDEPHOS-RXN_1[B]',form,[],0,0,1000);
% load('increase_withCO219nov24.mat')
% load('decrease_withCO219nov24.mat')
% %% 
% load('mixed_withCO219nov24.mat')
%load('lis_movedate.mat')
%load('lis_March25.mat')
%model=model5;


changeCobraSolver('glpk');

%load('SC_constrained_unblocked.mat');
%load('SC_constrainedmodelFINAL.mat');


sol_og=[]; pepc_og=[];co2_og=[];rub_og=[];h2o_og=[];
for n=0:10:290
 model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -n, 'l'); 
 model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -n, 'u'); 
og=optimizeCbModel(model);
sol_og=[sol_og,og.f*24/1000];
pepc_og1=find(contains(model.rxns,'PEPCARBOX-RXN[M]'));
co2_og1=find(contains(model.rxns,'EX_CARBON-DIOXIDE_EXTRACELLULAR'));
rub_og1=find(contains(model.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
h2o_og1=find(contains(model.rxns,'EX_WATER_EXTRACELLULAR'));
h2o_og=[h2o_og,og.v(h2o_og1)];

pepc_og=[pepc_og,og.v(pepc_og1)];
co2_og=[co2_og,og.v(co2_og1)];
rub_og=[rub_og,og.v(rub_og1)];
end

  figure(1)            
            plot(-co2_og,pepc_og,'color','r','LineWidth',4);
            hold on, drawnow
            plot(-co2_og,rub_og,'color','b','LineWidth',4);
            hold on, drawnow
                        plot(-co2_og,-h2o_og,'color','g','LineWidth',4);


               %         legend('100% N','0 % N' ,'Soybean N fixing','Location','Best');

            legend('PEPC','RuBisCO','Water uptake' ,'Location','Best');
            xlabel('CO2 uptake mmolg^{-1}hr^{-1}','FontSize',40)
            ylabel('Flux {\mu}molg^{-1}hr^{-1}','FontSize',40)
              set(gca,'LineWidth',2,'FontSize',40)
                             x_width=30 ;y_width=30;

  print('fluxes_co2','-djpeg','-loose');

 figure(3)            
            plot(-co2_og,sol_og,'color','k','LineWidth',4);
            hold on, drawnow
       %     plot(co2_sc,sol_sc,'color','k','LineWidth',4);

               %         legend('100% N','0 % N' ,'Soybean N fixing','Location','Best');

          %  legend('unconstrained','constrained' ,'Location','Best');
            xlabel('CO2 uptake mmolg^{-1}hr^{-1}','FontSize',40)
            ylabel('Relative growth rate (gg^{-1}d^{-1})','FontSize',40)
              set(gca,'LineWidth',2,'FontSize',40)  
                             x_width=30 ;y_width=30;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
              print('rgr_co2','-djpeg','-loose');
