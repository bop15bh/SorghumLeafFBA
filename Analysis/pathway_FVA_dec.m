clear
close all
% FVA pipeline
%addpath /mnt/home/holla293/Documents/cobratoolbox
%    initCobraToolbox; 

load('Leaf_model_biomass.mat')
 trans=readcell('old_transport.xlsx');
 trans=trans(:,1);% nbefore 19 was close 26 closest, 28,29,30,35, 38, 41,42
 % 49 best 
 %  nitrate only: 17 close, 19 bit weird, 22 weirdish, 23 but min mdh, got up to 35
% 3 remo: 19 weirdish, 20 is very good. 
% 3rd try 31 maybe 33 looks reasonable 
%4th try no dups 22 weird, 25 also weird 40 is v close, 42 only1, 49
%mental, 54 close (2 inc), 55 malo 0, 57 wild
 model=removeRxns(model,{'ATR_NITRITE_[cb]_[cm]','RXN-14932_1[M]','RXN-14932_2[M]'})
%model=removeRxns(model,{'ATR_NITRITE_[cb]_[cm]'});
model=removeRxns(model,{'RXN-13697[B]','RXN-13697[M]','TSA-REDUCT-RXN[B]','TSA-REDUCT-RXN[M]'});

% 
%  model=removeRxns(model,trans(40));
% deads = model.mets(detectDeadEnds(model));
% [rxns,rxnforms]=findRxnsFromMets(model,deads);
%  model=removeRxns(model,rxns);
%  deads = model.mets(detectDeadEnds(model))
% deads
% 19, 20, 26,30 broken 39 1 rxn, 40 good 1 rxn, 46 sort of broken
% 49 sort of broken, 54 malo broken 1, 55 close 
%model=removeRxns(model,{'TSA-REDUCT-RXN[B]','TSA-REDUCT-RXN[M]','5-METHYLTHIORIBOSE-KINASE-RXN[B]','5-METHYLTHIORIBOSE-KINASE-RXN[M]'})
%40 after dead removal, only 1 inc, 25 close 
dups=readcell('duplicates1.xlsx');
dups2=readcell('dupdup.xlsx');
model=removeRxns(model,dups(:,1))
model=removeRxns(model,dups2(:,1))
emets=model.mets(find(contains(model.mets,'E-[')))
[rxns,rxnforms]=findRxnsFromMets(model,emets)
model=removeRxns(model,rxns)
%%
% setting bounds to all rxns to be 1000 mmol/g/hr
model.ub(model.ub>1000)=1000;
model.lb(model.lb<-1000)=-1000;
%% Setting limit to Rubisco
model = changeRxnBounds(model,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]', 6.5, 'u');
%%
%NGAM 
model = changeRxnBounds(model,'ATPASE-RXN[B]', 0.0229, 'l');
model = changeRxnBounds(model,'ATPASE-RXN[M]', 0.0229, 'l');



%% ensuring c4 cycle 
model=removeRxns(model,{'4.1.1.32-RXN[M]','4.1.1.32-RXN[B]'});
%model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.2, 'l');
%model=removeRxns(model,{'ATR_G3P_[cb]_[cm]','ATR_CPD-1777_[cb]_[cm]'})
%model=removeRxns(model,{'ATR_GLYCERATE_[cb]_[cm]'})

%model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.5, 'l');
%model = changeRxnBounds(model,'MALATE-DEH-RXN_2[M]', 0.2, 'u');

%model=removeRxns(model,{'MALATE-DEH-RXN_2[M]','MALATE-DEH-RXN_2[B]'})

model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.3, 'l');
%model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.2, 'l');
model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 0.6, 'l');
%model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 10, 'u');
%model = changeRxnBounds(model,'ATR_MAL_[cb]_[cm]', 0.6, 'l');
%model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 0.9, 'l');
%%model = changeRxnBounds(model,'D-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN_2[M]', 0, 'l');
%model = changeRxnBounds(model,'ATR_L-ASPARTATE_[cb]_[cm]', 0.2, 'u');
%model = changeRxnBounds(model,'ATR_G3P_[cb]_[cm]', 0.3, 'u');
%model = changeRxnBounds(model,'ATR_CPD-1777_[cb]_[cm]', 0.3, 'u');


%model=removeRxns(model,{'4.1.1.32-RXN[M]','4.1.1.32-RXN[B]'})%,'RXN-13697[M]','RXN-13697[B]'})
%model=removeRxns(model,{'RXN-14932_2[M]','MALATE-DEH-RXN_1[M]','ATR_PLASTOQUINONE-9_[cb]_[cm]','ATR_G3P_[cb]_[cm]'})
%model=removeRxns(model,{'ATR_CPD-12829_[cb]_[cm]','ATR_PLASTOQUINONE-9_[cb]_[cm]'})
%model=removeRxns(model,{'RXN-14932_2[M]','MALATE-DEH-RXN_1[M]','ATR_G3P_[cb]_[cm]'})
%model=removeRxns(model,{'MALATE-DEH-RXN_1[M]','ATR_PLASTOQUINONE-9_[cb]_[cm]','ATR_G3P_[cb]_[cm]'})

 form=' ADP[cm] + PHOSPHO-ENOL-PYRUVATE[cm] + PROTON[cm] -> ATP[cm] + PYRUVATE[cm] ';
 model=addReaction(model,'PEPDEPHOS-RXN_1[M]',form,[],0,0,1000);
 form=' ADP[cb] + PHOSPHO-ENOL-PYRUVATE[cb] + PROTON[cb] -> ATP[cb] + PYRUVATE[cb] ';
 model=addReaction(model,'PEPDEPHOS-RXN_1[B]',form,[],0,0,1000);
form=' ADP[sm] + PHOSPHO-ENOL-PYRUVATE[sm] + PROTON[sm] -> ATP[sm] + PYRUVATE[sm] ';
 model=addReaction(model,'PEPDEPHOS-RXN_2[M]',form,[],0,0,1000);
 form=' ADP[sb] + PHOSPHO-ENOL-PYRUVATE[sb] + PROTON[sb] -> ATP[sb] + PYRUVATE[sb] ';
 model=addReaction(model,'PEPDEPHOS-RXN_2[B]',form,[],0,0,1000);

changeCobraSolver('glpk');

% load('increase_withCO219nov24.mat')
% load('decrease_withCO219nov24.mat')
% %% 
% load('mixed_withCO219nov24.mat')
%load('lis_movedate.mat')
%model=model5;
%load('sens_droMay25.mat')
%load('sens_controlMay25.mat')
%load('sens_dro350.mat')

%load('sens_control350.mat')
%load("sens_dro2PG.mat")
%load('sens_con2PG.mat')
%load('sens_droNN.mat')
%load('sens_controlNN.mat')
%load('sens_droJune25.mat');
%load('sens_controlJune25.mat');
load('sens_droDec25.mat');
load('sens_controlDec25.mat');
sens_dro_nocell=erase(sens_dro,'[M]')
sens_dro_nocell=erase(sens_dro_nocell,'[B]');
sens_dro_nocell=unique(sens_dro_nocell)
sens_con_nocell=erase(sens_control,'[M]');
sens_con_nocell=erase(sens_con_nocell,'[B]');
sens_con_nocell=unique(sens_con_nocell)
dro_nocell=setdiff(sens_dro_nocell,sens_con_nocell)
con_nocell=setdiff(sens_con_nocell,sens_dro_nocell)
dro_only=setdiff(sens_dro,sens_control);
con_only=setdiff(sens_control,sens_dro);
same=intersect(sens_dro,sens_control)
same_nocell=intersect(sens_con_nocell,sens_dro_nocell);
sames=same(find(contains(same,same_nocell)))
drodro=dro_only(find(contains(dro_only,dro_nocell)))
control=con_only(find(contains(con_only,con_nocell)))


%% control only 
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -5.587, 'l'); 
 model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u'); 
og=optimizeCbModel(model);
dro_version = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -0.847, 'l'); 
newfl=optimizeCbModel(dro_version);
pos_con=[];
for n=1:length(control)
pos2=find(strcmp(model.rxns,control{n}));
pos_con=[pos_con,pos2]
end
ogfl_con=og.v(pos_con)
drought_flux_con=newfl.v(pos_con);
lis_con=model.rxns(pos_con);
con_subs=model.subSystems(pos_con);

[minFlux_con, maxFlux_con] = fluxVariability(model, 99, 'max', lis_con);
conTAB = cell2table(horzcat(lis_con,num2cell(ogfl_con),num2cell(drought_flux_con), num2cell(minFlux_con),num2cell(maxFlux_con),con_subs));
writetable(conTAB, 'Control_only_rxnsP_analysisDEC25.txt', 'Delimiter', '\t') % tab-delimited

%% dro only 
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -0.847, 'l'); 
 model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u'); 
og=optimizeCbModel(model);
con_version = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -5.587, 'l'); 
newfl=optimizeCbModel(con_version);
pos_dro=[];
for n=1:length(drodro)
pos2=find(strcmp(model.rxns,drodro{n}));
pos_dro=[pos_dro,pos2]
end
ogfl_dro=og.v(pos_dro)
lis_dro=model.rxns(pos_dro);
dro_subs=model.subSystems(pos_dro);
con_flux_dro=newfl.v(pos_dro);

[minFlux_dro, maxFlux_dro] = fluxVariability(model, 99, 'max', lis_dro);
droTAB = cell2table(horzcat(lis_dro,num2cell(con_flux_dro),num2cell(ogfl_dro), num2cell(minFlux_dro),num2cell(maxFlux_dro),dro_subs));
writetable(droTAB, 'Drought_only_rxnsP_analysisDEC25.txt', 'Delimiter', '\t') % tab-delimited

%% same only 
pos_same=[];
for n=1:length(sames)
pos2=find(strcmp(model.rxns,sames{n}));
pos_same=[pos_same,pos2]
end
ogfl_same=og.v(pos_same)
same_flux_dro=newfl.v(pos_same);

lis_same=model.rxns(pos_same);
same_subs=model.subSystems(pos_same);

[minFlux_same, maxFlux_same] = fluxVariability(model, 99, 'max', lis_same);
sameTAB = cell2table(horzcat(lis_same,num2cell(same_flux_dro),num2cell(ogfl_same), num2cell(minFlux_same),num2cell(maxFlux_same),same_subs));
writetable(sameTAB, 'Same_only_rxnsP_analysisDEC25.txt', 'Delimiter', '\t') % tab-delimited
