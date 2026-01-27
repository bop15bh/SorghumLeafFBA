clear
close all
% drought at 53 ~10% SWC
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

model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -0.847, 'l');
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');

changeCobraSolver('glpk');

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
model_limits = changeRxnBounds(model_limits,'biomass', 1000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 1000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 1000, 'u');
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
%%
clearvars -except sens_dro

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

model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -5.587, 'l');
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');
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
model_limits = changeRxnBounds(model_limits,'biomass', 1000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[M]', 1000, 'u');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 0, 'l');
model_limits = changeRxnBounds(model_limits,'biomass[B]', 1000, 'u');
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
save('sens_controlDec25.mat','sens_control')
save('sens_droDec25.mat','sens_dro')
