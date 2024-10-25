%% comparing FVA results at different drought strengths multiple data points 
clear
close all
% FVA pipeline
%addpath /mnt/home/holla293/Documents/cobratoolbox
%    initCobraToolbox;

%changeCobraSolver('glpk')
%load('SC_constrained_unblocked.mat');
load('constrained_leaf_noredundancies.mat')
model=model5;

load('importantlist_oneCAleft.mat')
spread=un;
%spread={'RXN-12300[B]'}
%spread={'TRIOSEPISOMERIZATION-RXN_1[M]','RXN-1347[B]','RXN-12300[B]','NADH-DEHYDROG-A-RXN_1[B]','ADENYL-KIN-RXN_2[M]'}
%spread={'CYTOCHROME-C-OXIDASE-RXN_1[B]'}
%spread={'HYDROXYPYRUVATE-REDUCTASE-RXN_2[B]'}
%spread={'RXN-14182[B]'}
%spread=  {'5.1.1.18-RXN[M]','RXN-13470_1[M]'}
%load('fluct.mat')
%spread=lis;
%spread={'biomass','1TRANSKETO-RXN_2[B]','1TRANSKETO-RXN_2[M]','2TRANSKETO-RXN[B]','2TRANSKETO-RXN[M]','RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]' }
%spread={'XYLISOM-RXN_1[M]'}
%spread={'TRIOSEPISOMERIZATION-RXN_2[B]','TRIOSEPISOMERIZATION-RXN_1[M]','TRIOSEPISOMERIZATION-RXN_1[B]','TRIOSEPISOMERIZATION-RXN_2[M]','TRANS-RXN-TPT-Plastid[B]','TRANS-RXN-TPT-Plastid[M]'}
%spread={'GLYCEROL-DEHYDROGENASE-NADP+-RXN[B]','GLYCEROL-DEHYDROGENASE-NADP+-RXN[M]'}
%spread={'TRIOKINASE-RXN[B]','TRIOKINASE-RXN[M]'}
%spread={'TRIOSEPISOMERIZATION-RXN_1[B]','TRIOSEPISOMERIZATION-RXN_1[M]','TRIOSEPISOMERIZATION-RXN_2[B]','TRIOSEPISOMERIZATION-RXN_2[M]'}
%spread={'SERINE--PYRUVATE-AMINOTRANSFERASE-RXN[B]','SERINE--PYRUVATE-AMINOTRANSFERASE-RXN[M]'}
%spread={'HYDROXYPYRUVATE-REDUCTASE-RXN_1[B]','HYDROXYPYRUVATE-REDUCTASE-RXN_2[B]'}
%spread={'SHIKIMATE-KINASE-RXN[M]','RXN-1106[B]','PHENYLALANINE-AMMONIA-LYASE-RXN[M]','MALIC-NADP-RXN_1[B]','CELLULOSE-SYNTHASE-UDP-FORMING-RXN_1[M]','CARBOXYCYCLOHEXADIENYL-DEHYDRATASE-RXN[M]','RXN-9531_1[M]','RXN-9659_1[M]','RXN-9528_1[M]','RXN-2621[B]','RXN-2602[B]'}
%spread={'RXN-8473[B]','RXN-8473[M]'};%,'ADENYL-KIN-RXN_2[M]','ADENYL-KIN-RXN_2[B]','SHIKIMATE-5-DEHYDROGENASE-RXN[M]','SHIKIMATE-5-DEHYDROGENASE-RXN[B]'}
%spread={'ADENYL-KIN-RXN_1[B]','ADENYL-KIN-RXN_1[M]','ADENYL-KIN-RXN_2[B]','ADENYL-KIN-RXN_2[M]'}
%spread={'SHIKIMATE-5-DEHYDROGENASE-RXN[M]','SHIKIMATE-5-DEHYDROGENASE-RXN[B]'}
%spread={'dihydroorotate dehydrogenase(NAD+)[B]','dihydroorotate dehydrogenase(NAD+)[M]'};model=model5;
%pos_spread=find(contains(model.rxns,spread));
%rxns=model.rxns(pos_spread);
% spread={'TRIOSEPISOMERIZATION-RXN_2[B]','TRIOSEPISOMERIZATION-RXN_1[M]','TRIOSEPISOMERIZATION-RXN_2[M]','TRIOSEPISOMERIZATION-RXN_1[B]','TRIOKINASE-RXN[B]','TRIOKINASE-RXN[M]',...
% 'TRANS-RXN-TPT-Plastid[B]','SHIKIMATE-5-DEHYDROGENASE-RXN[M]','SHIKIMATE-5-DEHYDROGENASE-RXN[B]',...
% 'GLYCEROL-DEHYDROGENASE-NADP+-RXN[B]','GLYCEROL-DEHYDROGENASE-NADP+-RXN[M]','SERINE--PYRUVATE-AMINOTRANSFERASE-RXN[B]','SERINE--PYRUVATE-AMINOTRANSFERASE-RXN[M]','RXN-8473[M]', 'RXN-8473[B]',...
% 'HYDROXYPYRUVATE-REDUCTASE-RXN_1[B]','HYDROXYPYRUVATE-REDUCTASE-RXN_2[B]','GLYCEROL-KIN-RXN[B]','GLYCEROL-KIN-RXN[M]',...
% 'ADENYL-KIN-RXN_2[M]', 'ADENYL-KIN-RXN_2[B]', 'ADENYL-KIN-RXN_1[B]', 'ADENYL-KIN-RXN_1[M]'}


model = changeRxnBounds(model,'EX_Light_Compound_EXTRACELLULAR', 3000, 'u');
model=removeRxns(model,{'ATR_CARBON-DIOXIDE_CYTOSOL_PLASTID-STR[B]'})
 model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u'); 
   model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -270.3777, 'l'); 

model = changeRxnBounds(model,'RXN-15479_1[M]', 2000, 'u');
model = changeRxnBounds(model,'ATPSYN-RXN-Plastid[M]', 2000, 'u');
light=find(contains(model.rxns,'ATR_Light'))
model.lb(light)=-5000;
model.ub(light)=5000;
model=removeRxns(model,{'RXN0-5224[B]','RXN0-5224[M]'})
pos_spread=[];
for k=1:length(spread)
pos_spread1=find(contains(model.rxns,spread(k)));
    if length(pos_spread1)>1
        for lom=1:length(pos_spread1)
        poos=strmatch(spread{k},model.rxns{pos_spread1(lom)},'exact')
        if ~isempty(poos)
        pos_spread=[pos_spread,pos_spread1(lom)]
        else
        end
        end

    else
         poos=strmatch(spread{k},model.rxns(pos_spread1))
       if ~isempty(poos)
        pos_spread=[pos_spread,pos_spread1]
        else
        end
    end
end

rxns=model.rxns(pos_spread);
og=optimizeCbModel(model);
og_flux=og.v;
og_flux1=og_flux(pos_spread);

position_ca=find(contains(model.rxns,'CARBODEHYDRAT-RXN[M]'));
lis={};
%po=10:20:270
popo=[40.9091 122.7 192.3 225 265.3777];
%popo=[40 97.5 155 212.5 270]
po=popo;
%popo=po;
step_size = 2.5;  % Step change size
num_steps = 5;    % Number of steps
range = [-2, -1, 0, 1, 2];  % Relative steps from the original point


for n=1:5
    if n<5   
    sub_po = po(n) + step_size * range  % Create 5 points around original point
    flux_results = zeros(length(spread), num_steps);  % Store flux for each point
    else
        num_steps=1 
        sub_po = po(n) + step_size * 2;  % Create 5 points around original point
    flux_results = zeros(length(spread), num_steps);  % Store flux for each point

    end
    for s=1:num_steps

   model1=model;
   model1 = changeRxnBounds(model1,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -sub_po(s), 'l'); 
pep=optimizeCbModel(model1);


%jeff=cell2table(horzcat(model.rxns,num2cell(og.v),num2cell(minFlux),num2cell(maxFlux),num2cell(pep.v),num2cell(minFlux_pep),num2cell(maxFlux_pep)));
% rxns=table2cell(jeff(:,1))
% rxns(1)='';

dro=pep.v(pos_spread);
flux_results(:, s) = dro;  % Store each sub-point's fluxes
%rxns=model.rxns;
%rxns=spread;
% og_flux_min=table2cell(jeff(:,3));
% og_flux_max=table2cell(jeff(:,4));
% og_flux_min(1)='';
% og_flux_max(1)='';
% og_flux_min=cell2mat(og_flux_min);
% og_flux_max=cell2mat(og_flux_max);
% dro_min=table2cell(jeff(:,6));
% dro_max=table2cell(jeff(:,7));
% dro_min(1)='';
% % dro_max(1)='';
% dro_min=cell2mat(dro_min);
% dro_max=cell2mat(dro_max);

%dro1=dro(pos_spread);
%
%gas(:,n)=sub_po(s)
%floopy(:,n)=flux_results

end

  % Calculate median flux across the 5 sub-points for each reaction
    median_flux = median(flux_results, 2);

FC=abs(median_flux./og_flux1)

for i=1:length(spread)
whole_list(i,n)=rxns(i);
og_flurx(i,n)=og_flux1(i);
dro_flurx(i,n)=median_flux(i);
end
end

%jeff=cell2table(horzcat(model.rxns,num2cell(og.v),num2cell(minFlux),num2cell(maxFlux),num2cell(pep.v),num2cell(minFlux_pep),num2cell(maxFlux_pep)));
% rxns=table2cell(jeff(:,1))
% rxns(1)='';
% og_flux=table2cell(jeff(:,2));
% og_flux(1)='';
% dro=table2cell(jeff(:,5));
% % dro(1)='';

% for i=1:length(spread)
%     pos=strmatch(spread(i),model1.rxns,'exact');
%     whole_list(i,n)=model1.rxns(pos);
%     og_flurx(i,n)=og.v(pos);
%     dro_flurx(i,n)=pep.v(pos);
% end

rxns=whole_list(:,1);
control=og_flurx(:,1);
%co2=10:20:270;
co2=popo;
tempFC = {}; % Temporary cell array to collect valid rows
log={};tempoo={};
for p = 1:length(rxns)
    if sum(dro_flurx(p, :) > 0) < 5 && sum(dro_flurx(p,:)<0) < 5
        % Do nothing if condition is not met
   
   else
        log{end+1,1}=rxns(p);
        tempFC{end+1, 2} = abs(dro_flurx(p, :) ./ control(p)); % Collect valid rows
    end
end   
FC = cell2mat(tempFC);


co2=co2';

for n=1:length(log)
     if contains(log{n},'ATR_') || contains(log{n},'EX_')

     else

    fold=log2(FC(n,:));
   if  max(fold)>0.1315 || contains(log{n},'biomass')%|| min(fold)<-1 
lis=[lis,log{n}]
fold=fold'
X = [ones(length(co2),1) co2];
b = X\fold;
yCalc2 = X*b;
Rsq1 = 1 - sum((fold - yCalc2).^2)/sum((fold - mean(fold)).^2)
noo = length(co2);  % Number of data points
r = sum((co2 - mean(co2)) .* (fold - mean(fold))) / sqrt(sum((co2 - mean(co2)).^2) * sum((fold - mean(fold)).^2));

% Calculate t-statistic
t_stat = r * sqrt((noo-2) / (1 - r^2));

% Calculate p-value
p_value = 2 * (1 - tcdf(abs(t_stat), noo-2));
id=log{n};
   % if b(2)<0 || contains(log{n},'XYLISOM-RXN_1[M]')%|| min(fold)<-1 


%if Rsq1>0.
            figure(n)
            plot(co2,fold,'color','r','LineWidth',4);
         % plot(co2,abs(dro_flurx(n, :)),'color','r','LineWidth',4);
           hold on, drawnow
        %    plot(co2,yCalc2,'--')


            legend(['R^2=' num2str(round(Rsq1,3)) ',p=' num2str(p_value) ],'Location','Best');
            %  legend(['R^2=' num2str(Rsq1)],'Location','Best');

            xlabel('CO2 uptake umol/g/hr','FontSize',40)
            ylabel('Log 2 fold change','FontSize',40)
            title(log{n})
              set(gca,'LineWidth',2,'FontSize',40)
               x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
               print(['og5points' id{:} ''],'-djpeg','-loose');

 else
 end
   % else
   % end
     end
end

