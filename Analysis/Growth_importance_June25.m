clear
close all
% drought at 53 ~10% SWC
% FVA pipeline
%addpath /mnt/home/holla293/Documents/cobratoolbox
%    initCobraToolbox;
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

changeCobraSolver('glpk');
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -44, 'l');
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');
og=optimizeCbModel(model);
lis=model.rxns;
pos=find(contains(model.rxns,lis));
ogfl=og.v(pos)
lis=model.rxns(pos);
% model2=model;
%    model2 = changeRxnBounds(model2,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -40.9091 , 'l');
% contrl rgr is  0.1275
%og=optimizeCbModel(model);
model_limits=model;
for n=1:length(model.rxns)
if  og.v(n)<0
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n), 'l');
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n)-og.v(n)*0.1, 'u');
elseif og.v(n)>0
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n), 'u');
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n)-og.v(n)*0.1, 'l');
elseif  og.v(n)==0
model_limits = changeRxnBounds(model_limits,model.rxns(n), 0, 'u');
model_limits = changeRxnBounds(model_limits,model.rxns(n), 0, 'l');
end
end
% limits look in line with max fluxes at co2 limitation
dro=optimizeCbModel(model_limits)
og_flux=og.v
lb=model_limits.lb
ub=model_limits.ub
pos=find(contains(model.rxns,lis))
model_limits = changeRxnBounds(model_limits,'biomass', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass', 1000000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 100000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 100000, 'u');
% model_limits = changeRxnBounds(model_limits,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');
% RGR is same in this version to drought condition
sol=[];
for n=1:length(model_limits.rxns)
model2=model_limits;
if  og.v(n)<0
model2 = changeRxnBounds(model2,model_limits.rxns(n),og.v(n)-og.v(n)*0.1 , 'l');
%model2 = changeRxnBounds(model2,lis(n), og.v(pos(n))*0.95  , 'u');
else
model2 = changeRxnBounds(model2,model_limits.rxns(n), og.v(n)-og.v(n)*0.1 , 'u');
%    model2 = changeRxnBounds(model2,lis(n), og.v(pos(n))*0.95  , 'l');
end
new_dro=optimizeCbModel(model2);
sol=[sol,new_dro.f*24/1000];
end
cond=sol<(og.f*24/1000)*0.99;
sens_dro=model_limits.rxns(cond)
clearvars -except sens_dro

%%

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
changeCobraSolver('glpk');
%load('SC_constrained_unblocked.mat');
%load('SC_constrainedmodelFINAL.mat');
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -294.6979, 'l');
og=optimizeCbModel(model);

lis=model.rxns;
pos=find(contains(model.rxns,lis))
ogfl=og.v(pos)
lis=model.rxns(pos);
% model2=model;
%    model2 = changeRxnBounds(model2,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -40.9091 , 'l');
% contrl rgr is  0.1275
%og=optimizeCbModel(model);
model_limits=model;
for n=1:length(model.rxns)
if  og.v(n)<0
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n), 'l');
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n)-og.v(n)*0.1, 'u');
elseif og.v(n)>0
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n), 'u');
model_limits = changeRxnBounds(model_limits,model.rxns(n), og.v(n)-og.v(n)*0.1, 'l');
elseif  og.v(n)==0
model_limits = changeRxnBounds(model_limits,model.rxns(n), 0, 'u');
model_limits = changeRxnBounds(model_limits,model.rxns(n), 0, 'l');
end
end
% limits look in line with max fluxes at co2 limitation
%dro=optimizeCbModel(model_limits)
og_flux=og.v
lb=model_limits.lb
ub=model_limits.ub
pos=find(contains(model.rxns,lis))
model_limits = changeRxnBounds(model_limits,'biomass', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass', 1000000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 100000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 100000, 'u');
% model_limits = changeRxnBounds(model_limits,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');
% RGR is same in this version to drought condition
sol=[];
for n=1:length(model_limits.rxns)
model2=model_limits;
if  og.v(n)<0
model2 = changeRxnBounds(model2,model_limits.rxns(n),og.v(n)-og.v(n)*0.1 , 'l');
%model2 = changeRxnBounds(model2,lis(n), og.v(pos(n))*0.95  , 'u');
else
model2 = changeRxnBounds(model2,model_limits.rxns(n), og.v(n)-og.v(n)*0.1 , 'u');
%    model2 = changeRxnBounds(model2,lis(n), og.v(pos(n))*0.95  , 'l');
end
new_dro=optimizeCbModel(model2);
sol=[sol,new_dro.f*24/1000];
end
cond=sol<(og.f*24/1000)*0.99;
sens_control=model_limits.rxns(cond)

dro_only=setdiff(sens_dro,sens_control)
con_only=setdiff(sens_control,sens_dro)
save('sens_controlJune25.mat','sens_control')
save('sens_droJune25.mat','sens_dro')
