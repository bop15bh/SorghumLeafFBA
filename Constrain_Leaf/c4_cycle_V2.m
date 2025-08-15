%% fixing c4 cycle V2
clear
close all
changeCobraSolver('glpk');
load('constrained_leaf_noredundancies.mat')
%load('Constrained_model_Mar25.mat')
load('Constrained_model_MAY25.mat')
model=removeRxns(model,{'RXN-4144[M]','S-NORLAUDANOSOLINE-SYNTHASE-RXN[B]'})

pos=find(contains(model.rxns,'RXN-19748'))
model=removeRxns(model,model.rxns(pos))

%model=model5;
%load('SC_constrained_unblocked.mat');
%load('SC_constrainedmodelFINAL.mat');
model = changeRxnBounds(model,'EX_Light_Compound_EXTRACELLULAR', 3000, 'u');
model=removeRxns(model,{'ATR_CARBON-DIOXIDE_CYTOSOL_PLASTID-STR[B]','ATR_CARBON-DIOXIDE_[cb]_[cm]'})
model = changeRxnBounds(model,'RXN-15479[M]', 2000, 'u');
model = changeRxnBounds(model,'ATPSYN-RXN-Plastid[M]', 2000, 'u');
model=removeRxns(model,{'ATR_Light_Compound_[s]_[t][B]' ,'ATR_Light_Compound_[s]_[t][M]'})
light=find(contains(model.rxns,'ATR_Light'))
model.lb(light)=-5000;
model.ub(light)=5000;
% removing extra CA's stops growth
model1=model;
to_rem={'3.6.3.25-RXN_1[B]','3.6.3.25-RXN_1[M]','3.6.3.26-RXN[B]','3.6.3.26-RXN[M]','3.6.3.30-RXN[B]','3.6.3.30-RXN[M]',...
'3.6.3.4-RXN[B]','3.6.3.4-RXN[M]','3.6.3.6-RXN[B]','3.6.3.6-RXN[M]','ABC-27-RXN_1[B]','ABC-27-RXN_1[M]','ATR_SUCROSE_CYTOSOL_EXTRACELLULAR[B]',...
'ATR_SUCROSE_CYTOSOL_EXTRACELLULAR[M]','RXN-9615[B]','RXN-9615[M]','TRANS-RXN-206[B]','TRANS-RXN-206[M]','TRANS-RXN-207[B]','TRANS-RXN-207[M]','TRANS-RXN-213[B]','TRANS-RXN-213[M]'};
rem=to_rem;
%rem=setdiff(to_rem,'3.6.3.25-RXN_1[M]')
soly=[];model2=model1;
for n=1:length(rem)
model2=removeRxns(model2,rem(n));
sol=optimizeCbModel(model2);
soly=[soly,sol.f]
end
form='SULFATE[e] -> SULFATE[cb]';
form2='SULFATE[e] -> SULFATE[cm]';
model2 = addReaction(model2,'ATR_SULFATE_CYTOSOL_EXTRACELLULAR[M]',form2,[],0,0,1000);
model2 = addReaction(model2,'ATR_SULFATE_CYTOSOL_EXTRACELLULAR[B]',form,[],0,0,1000);
form2='NITRATE[e] -> NITRATE[cb]';
form3='NITRATE[e] -> NITRATE[cm]';
model2 = addReaction(model2,'ATR_NITRATE_CYTOSOL_EXTRACELLULAR[B]',form2,[],0,0,1000);
model2 = addReaction(model2,'ATR_NITRATE_CYTOSOL_EXTRACELLULAR[M]',form3,[],0,0,1000);
form4=' Pi[e] -> Pi[cb] ';
form5=' Pi[e] -> Pi[cm] ';
model2 = addReaction(model2,'ATR_Pi_CYTOSOL_EXTRACELLULAR[B]',form4,[],0,0,1000);
model2 = addReaction(model2,'ATR_Pi_CYTOSOL_EXTRACELLULAR[M]',form5,[],0,0,1000);
sol=optimizeCbModel(model2)

rub=find(contains(model2.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
pepc=find(contains(model2.rxns,'PEPCARBOX'))
CA1=find(contains(model2.rxns,'CARBODEHYDRAT-RXN'))
%   960.1139 for CA1
CAM=find(contains(model2.rxns,'RXN0-5224[M]'))
% -1000 CAM
CAB=find(contains(model2.rxns,'RXN0-5224[B]'))
% -13 CAB
[rxns,form]=findRxnsFromMets(model2,'HCO3[cm]','ConsumersOnly',1)
pos=find(contains(model2.rxns,rxns))
flux=sol.v(pos)
names=model2.rxns(pos);
forms=printRxnFormula(model2,names)

[rxns5,form]=findRxnsFromMets(model5,'HCO3[cm]','ConsumersOnly',1)
pos=find(contains(model5.rxns,rxns5))
deads=model2.mets(detectDeadEnds(model2))
model3=removeRxns(model2,{'RXN0-5224[M]','RXN-22735_1[M]','RXN-19659[M]','CARBAMATE-KINASE-RXN[M]'})


% rub=find(contains(model3.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
% pepc=find(contains(model3.rxns,'PEPCARBOX'))
% CA1=find(contains(model3.rxns,'CARBODEHYDRAT-RXN'))
% 
% [rxns,form]=findRxnsFromMets(model3,'CARBON-DIOXIDE[cm]','ConsumersOnly',1)
% pos=find(contains(model3.rxns,rxns))
% flux=sol.v(pos)
% names=model3.rxns(pos);
% forms=printRxnFormula(model3,names)
% make following two reactions irreversible to force PEPC!
model4=changeRxnBounds(model3,'1.2.1.2-RXN[M]',0,'l');
model4=changeRxnBounds(model4,'4.1.1.32-RXN[M]',0,'l')
rub=find(contains(model4.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
pepc=find(contains(model4.rxns,'PEPCARBOX'))
CA1=find(contains(model4.rxns,'CARBODEHYDRAT-RXN'))
mal=find(contains(model4.rxns,'ATR_MAL_[cb]_[cm]'))


%% removing final CA in BS and getting malate to be used properly
% malate is being sent from M to BS
model6=model4;
form1='MAL[sb] + NADP[sb]   ->   CARBON-DIOXIDE[sb] + NADPH[sb] + PYRUVATE[sb]';
model6 = addReaction(model6,'MALIC-NADP-RXN[B]',form1,[],0,0,1000);
%model4=removeRxns(model4,{'1.1.1.39-RXN[B]','RXN0-5224[B]'})

model6=changeRxnBounds(model6,'1.2.1.2-RXN[B]',0,'l');
model6=changeRxnBounds(model6,'4.1.1.32-RXN[B]',0,'l')
% now rubisco is being used but i still have to remove extra CA and now MAL
% is being sent from BS to M
% so lets get rid of some of these MDH's 

%% removing extra MDHs
model7=removeRxns(model6,{'1.1.1.39-RXN[B]','1.1.1.39-RXN[M]'})
ro=optimizeCbModel(model7)
rub=find(contains(model7.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
pepc=find(contains(model7.rxns,'PEPCARBOX'))
CA1=find(contains(model7.rxns,'CARBODEHYDRAT-RXN'))
mal=find(contains(model7.rxns,'ATR_MAL_[cb]_[cm]'))
ro.v(rub)
ro.v(mal)
ro.v(pepc)

%% comparing M<->B transport
forms=printRxnFormula(model7,model7.rxns);
cb=forms(find(contains(forms,'[cb]')));
rxnscb=model7.rxns(find(contains(forms,'[cb]')));
cbcm=cb(find(contains(cb,'[cm]')))
rxnscbcm=rxnscb(find(contains(cb,'[cm]')))
pos=find(contains(model7.rxns,rxnscbcm));
names=model7.rxns(pos);
forms=printRxnFormula(model7,names);
flux=ro.v(pos)
%%
model8=removeRxns(model7,'RXN0-5224[B]');
sol=optimizeCbModel(model8)
diff=setdiff(model5.rxns,model7.rxns)

for n=1:length(diff)
pos=strmatch(diff(n),model5.rxns,'exact');
formy = printRxnFormula(model5, model5.rxns(pos), false);
lb=model5.lb(pos);
ub=model5.ub(pos);
model8=addReaction(model8,diff{n},formy{:},[],0,lb,ub);
end
ro=optimizeCbModel(model8);
test=diff;
test=setdiff(diff,{'RIB5PISOM-RXN_3[M]','RXN-10763[M]','RXN0-1441_2[M]','RXN0-5398_1[M]','RXN-15346[M]','RXN0-5199[M]'});%,'RXN0-1441_2[M]','RXN0-5199[M]'})
model9=removeRxns(model8,test); 
%model9=model8;
%   soly=[];     
% for n=401:1:length(test)
% model9=removeRxns(model9,test(n));
% rop=optimizeCbModel(model9);
% soly=[soly,rop.f]
% end
form='L-ALPHA-ALANINE[cb] <=> L-ALPHA-ALANINE[cm]';
model9=addReaction(model9,'ATR_L-ALPHA-ALANINE[cb]_[cm]',form,[],0,-1000,1000);

% all are going right way up to here except alanine and aspartate
%% comparing M<->B transport
forms=printRxnFormula(model9,model9.rxns);
cb=forms(find(contains(forms,'[cb]')));
rxnscb=model9.rxns(find(contains(forms,'[cb]')));
cbcm=cb(find(contains(cb,'[cm]')))
rxnscbcm=rxnscb(find(contains(cb,'[cm]')))
pos=find(contains(model9.rxns,rxnscbcm));
names=model9.rxns(pos);
forms=printRxnFormula(model9,names);
flux=ro.v(pos)

%% testing if removing serine transport helps
model10=changeRxnBounds(model9,'ATR_MAL_[cb]_[cm]',0,'l')

form='PYRUVATE[cb] -> PYRUVATE[cm] '
model10=addReaction(model10,'ATR_PYRUVATE_[cb]_[cm]',form,[],0,0,1000);
form='G3P[cb] -> G3P[cm] '
model10=addReaction(model10,'ATR_G3P_[cb]_[cm]',form,[],0,0,1000);

model10=changeRxnBounds(model10,'ATR_L-ASPARTATE_[cb]_[cm]',0,'l')

model10=changeRxnBounds(model10,'ATR_L-ALPHA-ALANINE[cb]_[cm]',0,'l')

model10=removeRxns(model10,{'ATR_SER_[cb]_[cm]','ATR_GLY_[cb]_[cm]','ATR_PHOSPHO-ENOL-PYRUVATE_[cb]_[cm]','ATR_GLT_[cb]_[cm]','ATR_2-KETOGLUTARATE_[cm]_[cb]'} )
model=model10;
 save('Constrained_unblocked_MAY25.mat','model')
%% Unblocking remaining mets
clearvars -except model
deads=model.mets(detectDeadEnds(model))
[rxns5,form]=findRxnsFromMets(model,'1-183-2-183-SN-GLYCEROL-PHOSPHOCHOLINE[cm]')
model1=removeRxns(model,rxns5);
deadsnew=model1.mets(detectDeadEnds(model1));
extra=setdiff(deadsnew,deads);

[rxns5,form]=findRxnsFromMets(model1,'CPD-8090[cm]')
[rxns,form]=findRxnsFromMets(model1,'CPD-8093[cm]')
to_rem=vertcat(rxns5,rxns)
model2=removeRxns(model1,to_rem);
deads2=model2.mets(detectDeadEnds(model2));
extra=setdiff(deads2,deads);
[rxns5,form]=findRxnsFromMets(model2,extra)
model3=removeRxns(model2,rxns5);
deads3=model3.mets(detectDeadEnds(model3));
extra=setdiff(deads3,deads);

[rxns5,form]=findRxnsFromMets(model3,extra);
model4=removeRxns(model3,rxns5)
deads4=model4.mets(detectDeadEnds(model4));
extra=setdiff(deads4,deads);
m=setdiff(model.rxns,model4.rxns);
m=strrep(m,'[M]','[B]');

model5=removeRxns(model4,m);
deads5=model5.mets(detectDeadEnds(model5));
extra=setdiff(deads5,deads)
model6=model5;
load('constrained_leaf_noredundancies.mat')

[rxns5,form]=findRxnsFromMets(model6,'12-OXO-TRANS-10-DODECENOATE[cb]');
model7=removeRxns(model6,rxns5);
deads6=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads6,deads)
[rxns5,form]=findRxnsFromMets(model7,extra);
model8=removeRxns(model7,rxns5);
deads7=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads7,deads)
[rxns5,form]=findRxnsFromMets(model8,extra);
model9=removeRxns(model8,rxns5);
deads8=model9.mets(detectDeadEnds(model9));
extra=setdiff(deads8,deads);

b=setdiff(model6.rxns,model9.rxns);
b=strrep(b,'[B]','[M]');
model9=removeRxns(model9,b);
deads9=model9.mets(detectDeadEnds(model9))
extra=setdiff(deads9,deads)

[rxns5,form]=findRxnsFromMets(model9,'16-EPIVELLOSIMINE[cb]');

model10=removeRxns(model9,rxns5)
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads);

[rxns5,form]=findRxnsFromMets(model10,'2-ALPHA-HYDROXYETHYL-THPP[cb]');
model11=removeRxns(model10,rxns5);
deads11=model11.mets(detectDeadEnds(model11));
extra=setdiff(deads11,deads);
[rxns5,form]=findRxnsFromMets(model11,extra);
model12=removeRxns(model11,rxns5);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads);
b=setdiff(model9.rxns,model12.rxns);
b=strrep(b,'[B]','[M]');
model12=removeRxns(model12,b);
deads12=model12.mets(detectDeadEnds(model12));

model12=removeRxns(model12,{'ISOFLAVONE-2-HYDROXYLASE-RXN[B]','ISOFLAVONE-2-HYDROXYLASE-RXN[M]'})
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads11)
model12=removeMetabolites(model12,{'2-HYDROXYFORMONONETIN[cm]','2-HYDROXYFORMONONETIN[cb]'})
deads12=model12.mets(detectDeadEnds(model12));

model13=removeRxns(model12, {'RXN-8172[B]','RXN-8172[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
%% this may need to be removed based on RNA seq 
form='2-PHOSPHO-4-CYTIDINE-5-DIPHOSPHO-2-C-MET[cb] <=>   2C-METH-D-ERYTHRITOL-CYCLODIPHOSPHATE[cb] + CMP[cb] ';
model14=addReaction(model13,'RXN0-302[B]',form,[],0,-1000,1000);
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
form1='2-PHOSPHO-4-CYTIDINE-5-DIPHOSPHO-2-C-MET[cm] <=>   2C-METH-D-ERYTHRITOL-CYCLODIPHOSPHATE[cm] + CMP[cm] ';
model14=addReaction(model14,'RXN0-302[M]',form1,[],0,-1000,1000);
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
%%
model15=removeRxns(model14,{'RXN-2209[B]','RXN-2209[M]'})
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)

[rxns5,form]=findRxnsFromMets(model15,'4-hydroxybenzoate[cb]');

model16=removeRxns(model15,rxns5)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads);

[rxns5,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns5)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads);

[rxns5,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns5)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns5,form]=findRxnsFromMets(model16,extra)
model16=removeRxns(model16,rxns5)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns5,form]=findRxnsFromMets(model16,extra)

model16=removeRxns(model16,rxns5)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

B=setdiff(model15.rxns,model16.rxns);
B=strrep(B,'[B]','[M]');
model16=removeRxns(model16,B);

form='NADPH[cb] + OXYGEN-MOLECULE[cb] + PALMITATE[cb] + PROTON[cb]  <=>   16-HYDROXYPALMITATE[cb] + NADP[cb] + WATER[cb]';
model17=addReaction(model16,'RXN-16398[B]',form,[],0,-1000,1000);
form1='NADPH[cm] + OXYGEN-MOLECULE[cm] + PALMITATE[cm] + PROTON[cm]  <=>   16-HYDROXYPALMITATE[cm] + NADP[cm] + WATER[cm]';
model17=addReaction(model17,'RXN-16398[M]',form1,[],0,-1000,1000);

deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

form2='NADPH[cb] + OLEATE-CPD[cb] + OXYGEN-MOLECULE[cb] + PROTON[cb]  <=>   18-HYDROXYOLEATE[cb] + NADP[cb] + WATER[cb] '
model18=addReaction(model17,'RXN-16400[B]',form2,[],0,-1000,1000);
form3='NADPH[cm] + OLEATE-CPD[cm] + OXYGEN-MOLECULE[cm] + PROTON[cm]  <=>   18-HYDROXYOLEATE[cm] + NADP[cm] + WATER[cm] '
model18=addReaction(model18,'RXN-16400[M]',form3,[],0,-1000,1000);

deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)

form4='3-HYDROXY-3-METHYL-GLUTARYL-COA[cb]   <=>   3-KETOBUTYRATE[cb] + ACETYL-COA[cb] ';
model19=addReaction(model18,'HYDROXYMETHYLGLUTARYL-COA-LYASE-RXN[B]',form4,[],0,-1000,1000);
form5='3-HYDROXY-3-METHYL-GLUTARYL-COA[cm]   <=>   3-KETOBUTYRATE[cm] + ACETYL-COA[cm] ';
model19=addReaction(model19,'HYDROXYMETHYLGLUTARYL-COA-LYASE-RXN[M]',form5,[],0,-1000,1000);

deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)

[rxns5,form]=findRxnsFromMets(model19,'44-DIMETHYL-CHOLESTA-814-24-TRIENOL[cb]')
model20=removeRxns(model19,rxns5);
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads);
[rxns5,form]=findRxnsFromMets(model20,extra);
model20=removeRxns(model20,rxns5);
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)
[rxns5,form]=findRxnsFromMets(model20,extra);
model20=removeRxns(model20,rxns5);
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)

B=setdiff(model19.rxns,model20.rxns);
B=strrep(B,'[B]','[M]');
model20=removeRxns(model20,B)
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)

model20=removeRxns(model20,{'DOPACHROME-DELTA-ISOMERASE-RXN[B]','DOPACHROME-DELTA-ISOMERASE-RXN[M]'})
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads);
model20=removeRxns(model20,{'RXN-7653[B]','RXN-7653[M]'})
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)

model20=removeRxns(model20, {'PROSTAGLANDIN-E-SYNTHASE-RXN[B]', 'PROSTAGLANDIN-E-SYNTHASE-RXN[M]'})
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)
model20=removeRxns(model20, {'RXN-14704[B]','RXN-14704[M]','7KAPSYN-RXN[B]','7KAPSYN-RXN[M]'})
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)




% 
 [rxns5,form]=findRxnsFromMets(model20,'91018-TRIHYDROXYSTEARATE[cb]')
 model21=removeRxns(model20,rxns5);
 deads21=model21.mets(detectDeadEnds(model21));
 extra=setdiff(deads21,deads)
  [rxns5,form]=findRxnsFromMets(model21,extra);
  model21=removeRxns(model21,rxns5);
 deads21=model21.mets(detectDeadEnds(model21));
 extra=setdiff(deads21,deads)
  [rxns5,form]=findRxnsFromMets(model21,extra);
model21=removeRxns(model21,rxns5);
 deads21=model21.mets(detectDeadEnds(model21));
 extra=setdiff(deads21,deads)

 B=setdiff(model20.rxns,model21.rxns);
B=strrep(B,'[B]','[M]');
model21=removeRxns(model21,B)
deads21=model21.mets(detectDeadEnds(model21));
extra=setdiff(deads21,deads)

model21=removeRxns(model21,{'ALLYL-ALCOHOL-DEHYDROGENASE-RXN[B]','ALLYL-ALCOHOL-DEHYDROGENASE-RXN[M]'})
deads21=model21.mets(detectDeadEnds(model21));
extra=setdiff(deads21,deads)

model21 = changeRxnBounds(model21,'ALLANTOICASE-RXN[M]', -1000, 'l');
model21 = changeRxnBounds(model21,'ALLANTOICASE-RXN[B]', -1000, 'l');

deads21=model21.mets(detectDeadEnds(model21));

 [rxns5,form]=findRxnsFromMets(model21,'ALPHA-TOCOPHEROL[cb]')
model22=removeRxns(model21,rxns5);
 deads22=model22.mets(detectDeadEnds(model22));
 extra=setdiff(deads22,deads)
 [rxns5,form]=findRxnsFromMets(model22,extra)
 model22=removeRxns(model22,rxns5);
 deads22=model22.mets(detectDeadEnds(model22));
 extra=setdiff(deads22,deads)
 [rxns5,form]=findRxnsFromMets(model22,extra)
 model22=removeRxns(model22,rxns5);
 deads22=model22.mets(detectDeadEnds(model22));
 extra=setdiff(deads22,deads)

 B=setdiff(model21.rxns,model22.rxns);
 B=strrep(B,'[B]','[M]')
 model22=removeRxns(model22,B);
 deads22=model22.mets(detectDeadEnds(model22));
 extra=setdiff(deads22,deads)

% BENZALDEHYDE[cm] is blocked based on RNA seq?
[rxns5,form]=findRxnsFromMets(model22,'BETA-TOCOPHEROL[cb]')
model23=removeRxns(model22,rxns5);
 deads23=model23.mets(detectDeadEnds(model23));
 extra=setdiff(deads23,deads);
 [rxns5,form]=findRxnsFromMets(model23,extra)
model23=removeRxns(model23,rxns5);
 deads23=model23.mets(detectDeadEnds(model23));
 extra=setdiff(deads23,deads);
 [rxns5,form]=findRxnsFromMets(model23,extra)
model23=removeRxns(model23,rxns5);
 deads23=model23.mets(detectDeadEnds(model23));
 extra=setdiff(deads23,deads)

 B=setdiff(model22.rxns,model23.rxns)
  B=strrep(B,'[B]','[M]')
  model23=removeRxns(model23,B)
   deads23=model23.mets(detectDeadEnds(model23));
 extra=setdiff(deads23,deads)

 [rxns5,form]=findRxnsFromMets(model23,'CARBON-MONOXIDE[cb]')
model24=removeRxns(model23,rxns5);
 deads24=model24.mets(detectDeadEnds(model24));
 extra=setdiff(deads24,deads)
[rxns5,form]=findRxnsFromMets(model24,extra)
model24=removeRxns(model24,rxns5);
 deads24=model24.mets(detectDeadEnds(model24));
 extra=setdiff(deads24,deads)
 [rxns5,form]=findRxnsFromMets(model24,extra)
model24=removeRxns(model24,rxns5);
 deads24=model24.mets(detectDeadEnds(model24));
 extra=setdiff(deads24,deads)
 [rxns5,form]=findRxnsFromMets(model24,extra)
model24=removeRxns(model24,rxns5);
 deads24=model24.mets(detectDeadEnds(model24));
 extra=setdiff(deads24,deads)
 [rxns5,form]=findRxnsFromMets(model24,extra)
model24=removeRxns(model24,rxns5);
 deads24=model24.mets(detectDeadEnds(model24));
 extra=setdiff(deads24,deads);

 B=setdiff(model23.rxns,model24.rxns)
  B=strrep(B,'[B]','[M]')
  model24=removeRxns(model24,B)
   deads24=model24.mets(detectDeadEnds(model24));
 extra=setdiff(deads24,deads)

 [rxns5,form]=findRxnsFromMets(model24,'CHOLESTEROL[cb]')
model25=removeRxns(model24,rxns5);
 deads25=model25.mets(detectDeadEnds(model25));
 extra=setdiff(deads25,deads)
  [rxns5,form]=findRxnsFromMets(model25,extra)
model25=removeRxns(model25,rxns5);
 deads25=model25.mets(detectDeadEnds(model25));
 extra=setdiff(deads25,deads)
   [rxns5,form]=findRxnsFromMets(model25,extra)
model25=removeRxns(model25,rxns5);
 deads25=model25.mets(detectDeadEnds(model25));
 extra=setdiff(deads25,deads)
   [rxns5,form]=findRxnsFromMets(model25,extra)
model25=removeRxns(model25,rxns5);
 deads25=model25.mets(detectDeadEnds(model25));
 extra=setdiff(deads25,deads)
   [rxns5,form]=findRxnsFromMets(model25,extra)
model25=removeRxns(model25,rxns5);
 deads25=model25.mets(detectDeadEnds(model25));
 extra=setdiff(deads25,deads)
   [rxns5,form]=findRxnsFromMets(model25,extra)
model25=removeRxns(model25,rxns5);
 deads25=model25.mets(detectDeadEnds(model25));
 extra=setdiff(deads25,deads)

 B=setdiff(model24.rxns,model25.rxns)
  B=strrep(B,'[B]','[M]')
  model25=removeRxns(model25,B)
   deads25=model25.mets(detectDeadEnds(model25));
 extra=setdiff(deads25,deads)

 model25 = changeRxnBounds(model25,'RXN-17473[B]', -1000, 'l');
model25 = changeRxnBounds(model25,'RXN-17473[M]', -1000, 'l');


 [rxns5,form]=findRxnsFromMets(model25,'Behenoyl-ACPs_Compound[cb]')
model26=removeRxns(model25,rxns5);
 deads26=model26.mets(detectDeadEnds(model26));
 extra=setdiff(deads26,deads)
 [rxns5,form]=findRxnsFromMets(model26,extra)
model26=removeRxns(model26,rxns5);
 deads26=model26.mets(detectDeadEnds(model26));
 extra=setdiff(deads26,deads)
 [rxns5,form]=findRxnsFromMets(model26,extra)
model26=removeRxns(model26,rxns5);
 deads26=model26.mets(detectDeadEnds(model26));
 extra=setdiff(deads26,deads)
 [rxns5,form]=findRxnsFromMets(model26,extra)
model26=removeRxns(model26,rxns5);
 deads26=model26.mets(detectDeadEnds(model26));
 extra=setdiff(deads26,deads)
 [rxns5,form]=findRxnsFromMets(model26,extra)
model26=removeRxns(model26,rxns5);
 deads26=model26.mets(detectDeadEnds(model26));
 extra=setdiff(deads26,deads)
 [rxns5,form]=findRxnsFromMets(model26,extra)
model26=removeRxns(model26,rxns5);
 deads26=model26.mets(detectDeadEnds(model26));
 extra=setdiff(deads26,deads)
 
 B=setdiff(model25.rxns,model26.rxns)
  B=strrep(B,'[B]','[M]')
  model26=removeRxns(model26,B)
   deads26=model26.mets(detectDeadEnds(model26));
 extra=setdiff(deads26,deads)

 model26=removeRxns(model26,{'ATR_CU+2_CYTOSOL_EXTRACELLULAR[M]','ATR_CU+2_CYTOSOL_EXTRACELLULAR[B]','ATR_TYR_CYTOSOL_PLASTID-STR[M]'})
   deads26=model26.mets(detectDeadEnds(model26));
%% fixing glucopyranose
   glu=find(contains(model26.mets,'Glucopyranose'))
    [rxns,form]=findRxnsFromMets(model26,'Glucopyranose[cb]');
    for n=1:length(rxns)
        form1=strrep(form{n},'Glucopyranose[cb]','ALPHA-GLUCOSE[cb]')
        if contains(form1,'<=>')==1
        model26=addReaction(model26,rxns{n},form1,[],0,-1000,1000);
        else
        model26=addReaction(model26,rxns{n},form1,[],0,0,1000);
        end
    end
    [rxns,form]=findRxnsFromMets(model26,'Glucopyranose[cm]');
    for n=1:length(rxns)
        form1=strrep(form{n},'Glucopyranose[cm]','ALPHA-GLUCOSE[cm]')
        if contains(form1,'<=>')==1
        model26=addReaction(model26,rxns{n},form1,[],0,-1000,1000);
        else
        model26=addReaction(model26,rxns{n},form1,[],0,0,1000);
        end
    end
   deads26=model26.mets(detectDeadEnds(model26));


   model26=removeRxns(model26,{'ATR_NITRITE_PLASTID-STR_CYTOSOL[B]','ATR_LYS_CYTOSOL_PLASTID[M]', 'RXN-18031[B]','ATR_GLYCOLLATE_CYTOSOL_PLASTID-STR[M]', ...
       'RXN-18232[B]','RXN-18232[M]', 'D-RIBULOKIN-RXN[B]', 'DCYSDESULF-RXN[M]', 'DCYSDESULF-RXN[B]'})
   deads26=model26.mets(detectDeadEnds(model26));
model26=removeMetabolites(model26,{'Glucopyranose[cb]','Glucopyranose[cm]'})


model27=removeRxns(model26,{'TRANSALDOL-RXN[B]','TRANSALDOL-RXN[M]'})
   deads27=model27.mets(detectDeadEnds(model27));
extra=setdiff(deads27,deads)
model27=removeRxns(model27,{'ATR_MAL_[s]_[c][M]'})
   deads27=model27.mets(detectDeadEnds(model27));

       [rxns,form]=findRxnsFromMets(model27,'MANNITOL[cb]');
model28=removeRxns(model27,rxns)
   deads28=model28.mets(detectDeadEnds(model28));
extra=setdiff(deads28,deads)
       [rxns,form]=findRxnsFromMets(model28,extra);
      model28=removeRxns(model28,rxns)
    deads28=model28.mets(detectDeadEnds(model28));
extra=setdiff(deads28,deads)
B=setdiff(model27.rxns,model28.rxns)
B=strrep(B,'[B]','[M]')
model28=removeRxns(model28,B)
    deads28=model28.mets(detectDeadEnds(model28));
extra=setdiff(deads28,deads)

       [rxns,form]=findRxnsFromMets(model28,'MYRICETIN[cb]');
model29=removeRxns(model28,rxns)
deads29=model29.mets(detectDeadEnds(model29));
extra=setdiff(deads29,deads)
       [rxns,form]=findRxnsFromMets(model29,extra);
       model29=removeRxns(model29,rxns)
deads29=model29.mets(detectDeadEnds(model29));
extra=setdiff(deads29,deads)
B=setdiff(model28.rxns,model29.rxns)
B=strrep(B,'[B]','[M]')
model29=removeRxns(model29,B)
deads29=model29.mets(detectDeadEnds(model29));
extra=setdiff(deads29,deads)
model29=removeRxns(model29,{'RXN-13076[B]','RXN-13076[M]','RXN-8037[B]','RXN-8037[M]'})
deads29=model29.mets(detectDeadEnds(model29));
extra=setdiff(deads29,deads)       
[rxns,form]=findRxnsFromMets(model29,extra);
model30=removeRxns(model29,rxns)
deads30=model30.mets(detectDeadEnds(model30));
extra=setdiff(deads30,deads) 
[rxns,form]=findRxnsFromMets(model30,extra);
model30=removeRxns(model30,rxns)
deads30=model30.mets(detectDeadEnds(model30));
extra=setdiff(deads30,deads) 
[rxns,form]=findRxnsFromMets(model30,extra);
model30=removeRxns(model30,rxns)
deads30=model30.mets(detectDeadEnds(model30));
extra=setdiff(deads30,deads) 
[rxns,form]=findRxnsFromMets(model30,extra);
model30=removeRxns(model30,rxns)
deads30=model30.mets(detectDeadEnds(model30));
extra=setdiff(deads30,deads) 

form='ISOCHORISMATE[cb] <=> CHORISMATE[cb]';
model31=addReaction(model30,'ISOCHORSYN-RXN[B]',form,[],0,-1000,1000);
deads31=model31.mets(detectDeadEnds(model31));
extra=setdiff(deads31,deads) 

[rxns,form]=findRxnsFromMets(model31,'BETAINE[cb]');

model32=removeRxns(model31,rxns)
deads32=model32.mets(detectDeadEnds(model32));
extra=setdiff(deads32,deads) 

B=setdiff(model31.rxns,model32.rxns);
B=strrep(B,'[B]','[M]')
model32=removeRxns(model32,B)
deads32=model32.mets(detectDeadEnds(model32));
extra=setdiff(deads32,deads)

model32=removeRxns(model32,{ 'RXN-4702[B]','RXN-4702[M]'})
deads32=model32.mets(detectDeadEnds(model32));
extra=setdiff(deads32,deads)

[rxns,form]=findRxnsFromMets(model32,'CINNAMALDEHYDE[cb]');

model33=removeRxns(model32,rxns)
deads33=model33.mets(detectDeadEnds(model33));
extra=setdiff(deads33,deads)
model33=removeRxns(model33,{'CINNAMOYL-COA-REDUCTASE-RXN[M]'})
deads33=model33.mets(detectDeadEnds(model33));
extra=setdiff(deads33,deads)

model33=removeRxns(model33,{'RXN-9619[B]','RXN-9619[M]'})
deads33=model33.mets(detectDeadEnds(model33));
extra=setdiff(deads33,deads)
[rxns,form]=findRxnsFromMets(model33,extra);

model34=removeRxns(model33,rxns)
deads34=model34.mets(detectDeadEnds(model34));
extra=setdiff(deads34,deads)
[rxns,form]=findRxnsFromMets(model34,extra);

model34=removeRxns(model34,rxns)
deads34=model34.mets(detectDeadEnds(model34));
extra=setdiff(deads34,deads)
[rxns,form]=findRxnsFromMets(model34,extra);

model34=removeRxns(model34,rxns)
deads34=model34.mets(detectDeadEnds(model34));
extra=setdiff(deads34,deads)
[rxns,form]=findRxnsFromMets(model34,extra);

model34=removeRxns(model34,rxns)
deads34=model34.mets(detectDeadEnds(model34));
extra=setdiff(deads34,deads)

 model34 = changeRxnBounds(model34,'RXN-13939[B]', -1000, 'l');
model34 = changeRxnBounds(model34,'RXN-13939[M]', -1000, 'l');
deads34=model34.mets(detectDeadEnds(model34));
extra=setdiff(deads34,deads)

model35=removeRxns(model34,{'RXN-7645[M]','RXN-7645[B]'})
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)
[rxns,form]=findRxnsFromMets(model35,extra);
model35=removeRxns(model35,rxns)
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)

[rxns,form]=findRxnsFromMets(model35,extra);
model35=removeRxns(model35,rxns)
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)

model35=removeRxns(model35,{ 'RXN-8502[B]', 'RXN-8502[M]'})
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)
[rxns,form]=findRxnsFromMets(model35,extra);
model35=removeRxns(model35,rxns)
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)

model35=removeRxns(model35, {'ATR_ZN+2_CYTOSOL_EXTRACELLULAR[M]','SIROHEME-FERROCHELAT-RXN[B]','SIROHEME-FERROCHELAT-RXN[M]'})
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)
[rxns,form]=findRxnsFromMets(model35,extra);
model35=removeRxns(model35,rxns)
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)
[rxns,form]=findRxnsFromMets(model35,extra);
model35=removeRxns(model35,rxns)
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)
[rxns,form]=findRxnsFromMets(model35,extra);
model35=removeRxns(model35,rxns)
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)

[rxns,form]=findRxnsFromMets(model35,{'SECOLOGANIN-CPD[cm]','SECOLOGANIN-CPD[cb]'});
model35=removeRxns(model35,rxns)
deads35=model35.mets(detectDeadEnds(model35));
extra=setdiff(deads35,deads)

[rxns,form]=findRxnsFromMets(model35,{'SARCOSINE[cm]','SARCOSINE[cb]'});
model36=removeRxns(model35,rxns)
deads36=model36.mets(detectDeadEnds(model36));
extra=setdiff(deads36,deads)


 model36 = changeRxnBounds(model36,'RXN-2881[B]', -1000, 'l');
model36 = changeRxnBounds(model36,'RXN-2881[M]', -1000, 'l');
deads36=model36.mets(detectDeadEnds(model36));
extra=setdiff(deads36,deads)

model36=removeRxns(model36,{'ATR_SACCHAROPINE_CYTOSOL_PLASTID[M]'})
deads36=model36.mets(detectDeadEnds(model36));
extra=setdiff(deads36,deads)

model37=removeRxns(model36,{'S-NORLAUDANOSOLINE-SYNTHASE-RXN[B]','S-NORLAUDANOSOLINE-SYNTHASE-RXN[M]'})
deads37=model37.mets(detectDeadEnds(model37));
extra=setdiff(deads37,deads)
[rxns,form]=findRxnsFromMets(model37,extra);
model37=removeRxns(model37,rxns)
deads37=model37.mets(detectDeadEnds(model37));
extra=setdiff(deads37,deads)

[rxns,form]=findRxnsFromMets(model37,'RIBITOL[cm]');
model38=removeRxns(model37,rxns)
deads38=model38.mets(detectDeadEnds(model38));
extra=setdiff(deads38,deads)
[rxns,form]=findRxnsFromMets(model38,extra);
model38=removeRxns(model38,rxns)
deads38=model38.mets(detectDeadEnds(model38));
extra=setdiff(deads38,deads)

model38=removeRxns(model38,{'RXN-17809[B]','RXN-17809[M]'})
deads38=model38.mets(detectDeadEnds(model38));
extra=setdiff(deads38,deads)
[rxns,form]=findRxnsFromMets(model38,extra);
model38=removeRxns(model38,rxns)
deads38=model38.mets(detectDeadEnds(model38));
extra=setdiff(deads38,deads)
[rxns,form]=findRxnsFromMets(model38,extra);
model38=removeRxns(model38,rxns)
deads38=model38.mets(detectDeadEnds(model38));
extra=setdiff(deads38,deads)

[rxns,form]=findRxnsFromMets(model38,{'PHTYOSPHINGOSINE-1-P[cb]','PHTYOSPHINGOSINE-1-P[cm]'});
model39=removeRxns(model38,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra);
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)

[rxns,form]=findRxnsFromMets(model39,{'PALMITALDEHYDE[cm]','PALMITALDEHYDE[cb]'})
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra);
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra);
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra);
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)

[rxns,form]=findRxnsFromMets(model39,{'OH[cm]','OH[cb]'})
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra)
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra)
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)

[rxns,form]=findRxnsFromMets(model39,{'NITRIC-OXIDE[cb]','NITRIC-OXIDE[cm]'})
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra)
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)

[rxns,form]=findRxnsFromMets(model39,{'MEK[cb]','MEK[cm]'})
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra)
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)
[rxns,form]=findRxnsFromMets(model39,extra)
model39=removeRxns(model39,rxns)
deads39=model39.mets(detectDeadEnds(model39));
extra=setdiff(deads39,deads)

[rxns,form]=findRxnsFromMets(model39,{'DTDP-RHAMNOSE[cb]','DTDP-RHAMNOSE[cm]'})
model40=removeRxns(model39,rxns)
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)

[rxns,form]=findRxnsFromMets(model40,extra)
model40=removeRxns(model40,rxns)
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)
[rxns,form]=findRxnsFromMets(model40,extra)
model40=removeRxns(model40,rxns)
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)
[rxns,form]=findRxnsFromMets(model40,extra)
model40=removeRxns(model40,rxns)
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)
[rxns,form]=findRxnsFromMets(model40,extra)
model40=removeRxns(model40,rxns)
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)
[rxns,form]=findRxnsFromMets(model40,extra)
model40=removeRxns(model40,rxns)
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)

[rxns,form]=findRxnsFromMets(model40,extra)
model40=removeRxns(model40,rxns)
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)

form='7-8-DIHYDROPTEROATE[cb] + ATP[cb] + GLT[cb]   ->   ADP[cb]   + DIHYDROFOLATE-GLU-N[cb] + PROTON[cb] + Pi[cb] ';
model40=addReaction(model40,'DIHYDROFOLATESYNTH-RXN[B]',form,[],0,0,1000);
form='7-8-DIHYDROPTEROATE[cm] + ATP[cm] + GLT[cm]   ->   ADP[cm]   + DIHYDROFOLATE-GLU-N[cm] + PROTON[cm] + Pi[cm] ';
model40=addReaction(model40,'DIHYDROFOLATESYNTH-RXN[M]',form,[],0,0,1000);
model40=removeMetabolites(model40,{'DIHYDROFOLATE[cb]','DIHYDROFOLATE[cm]'})
deads40=model40.mets(detectDeadEnds(model40));
extra=setdiff(deads40,deads)

[rxns,form]=findRxnsFromMets(model40,{'ETHYLENE-CMPD[cb]','ETHYLENE-CMPD[cm]'})
model41=removeRxns(model40,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)
[rxns,form]=findRxnsFromMets(model41,extra)
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)

[rxns,form]=findRxnsFromMets(model41,{'DEOXYCYTIDINE[cb]','DEOXYCYTIDINE[cm]'})
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)
[rxns,form]=findRxnsFromMets(model41,extra)
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)

[rxns,form]=findRxnsFromMets(model41,{'CH33ADO[cb]','CH33ADO[cm]'})
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)
[rxns,form]=findRxnsFromMets(model41,extra)
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)

[rxns,form]=findRxnsFromMets(model41,{'CPD-10279[cb]','CPD-10279[cm]'});
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)

[rxns,form]=findRxnsFromMets(model41,extra)
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)
[rxns,form]=findRxnsFromMets(model41,extra)
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)
[rxns,form]=findRxnsFromMets(model41,extra)
model41=removeRxns(model41,rxns)
deads41=model41.mets(detectDeadEnds(model41));
extra=setdiff(deads41,deads)

[rxns,form]=findRxnsFromMets(model41,{'CPD-10511[cm]','CPD-10511[cb]'});
model42=removeRxns(model41,rxns)
deads42=model42.mets(detectDeadEnds(model42));
extra=setdiff(deads42,deads)
[rxns,form]=findRxnsFromMets(model42,extra)
model42=removeRxns(model42,rxns)
deads42=model42.mets(detectDeadEnds(model42));
extra=setdiff(deads42,deads)

model42=removeRxns(model42,{'RXN-9868[B]','RXN-9868[M]'})
deads42=model42.mets(detectDeadEnds(model42));
extra=setdiff(deads42,deads)

[rxns,form]=findRxnsFromMets(model41,{'CPD-10687[cb]','CPD-10687[cm]'});
model42=removeRxns(model42,rxns)
deads42=model42.mets(detectDeadEnds(model42));
extra=setdiff(deads42,deads)

[rxns,form]=findRxnsFromMets(model41,extra)
model42=removeRxns(model42,rxns)
deads42=model42.mets(detectDeadEnds(model42));
extra=setdiff(deads42,deads)

model=model42;
clearvars -except model
 save('300deadsMar21.mat','model')
%% continuing dead end removal
model1=model;
deads=model1.mets(detectDeadEnds(model1));
form1='1 CPD-734[cm] + 1 ILE[cm] + 1 ATP[cm] -> 1 PROTON[cm] + 1 CPD-11232[cm] + 1 AMP[cm] + 1 PPI[cm]'
model1=addReaction(model1,'RXN-10435[M]',form1,[],0,0,1000);
deads1=model1.mets(detectDeadEnds(model1));
extra=setdiff(deads1,deads)
fixed=setdiff(deads,deads1)
form2='CPD-11232[cm] -> CPD-11259[cm]';
model1=addReaction(model1,'RXN-10454[M]',form2,[],0,0,1000);
deads1=model1.mets(detectDeadEnds(model1));
extra=setdiff(deads1,deads)
fixed=setdiff(deads,deads1)
form3='1 CPD-734[cb] + 1 ILE[cb] + 1 ATP[cb] -> 1 PROTON[cb] + 1 CPD-11232[cb] + 1 AMP[cb] + 1 PPI[cb]'
model1=addReaction(model1,'RXN-10435[B]',form3,[],0,0,1000);
form4='CPD-11232[cb] -> CPD-11259[cb]';
model1=addReaction(model1,'RXN-10454[B]',form4,[],0,0,1000);
deads1=model1.mets(detectDeadEnds(model1));
extra=setdiff(deads1,deads)
fixed=setdiff(deads,deads1)
model2=removeRxns(model1,{'RXN-12879[B]','RXN-12879[M]'})
deads2=model2.mets(detectDeadEnds(model2));
extra=setdiff(deads2,deads)
fixed=setdiff(deads,deads2)

[rxns,form]=findRxnsFromMets(model2,{'CPD-11712[cb]','CPD-11712[cm]'});
model3=removeRxns(model2,rxns)
deads3=model3.mets(detectDeadEnds(model3));
extra=setdiff(deads3,deads)
[rxns,form]=findRxnsFromMets(model3,extra);
model3=removeRxns(model3,rxns)
deads3=model3.mets(detectDeadEnds(model3));
extra=setdiff(deads3,deads)
[rxns,form]=findRxnsFromMets(model3,extra);
model3=removeRxns(model3,rxns)
deads3=model3.mets(detectDeadEnds(model3));
extra=setdiff(deads3,deads)

model4=changeRxnBounds(model3,'RXN-10974[B]', -1000, 'l');
model4=removeRxns(model4, {'RXN-10974[M]'})
deads4=model4.mets(detectDeadEnds(model4));
extra=setdiff(deads4,deads)

model4=removeRxns(model4,{'RXN-11524[B]','RXN-11524[M]'})
deads4=model4.mets(detectDeadEnds(model4));
extra=setdiff(deads4,deads)

model4=removeRxns(model4,{'RXN-11582[B]','RXN-11582[M]'})
deads4=model4.mets(detectDeadEnds(model4));
extra=setdiff(deads4,deads)

model4=removeRxns(model4,{'RXN-12325[B]','RXN-12325[M]'})
deads4=model4.mets(detectDeadEnds(model4));
extra=setdiff(deads4,deads)
[rxns,form]=findRxnsFromMets(model4,extra)
model5=removeRxns(model4,rxns)
deads5=model5.mets(detectDeadEnds(model5));
extra=setdiff(deads5,deads)

model6=model5;
model6=removeRxns(model6,{'MANDELONITRILE-LYASE-RXN[B]','MANDELONITRILE-LYASE-RXN[M]'})
deads6=model6.mets(detectDeadEnds(model6));
extra=setdiff(deads6,deads)
[rxns,form]=findRxnsFromMets(model6,extra);
model6=removeRxns(model6,rxns)
deads6=model6.mets(detectDeadEnds(model6));
extra=setdiff(deads6,deads)
[rxns,form]=findRxnsFromMets(model6,extra);
model6=removeRxns(model6,rxns)
deads6=model6.mets(detectDeadEnds(model6));
extra=setdiff(deads6,deads)

model6=removeRxns(model6,{'RXN-14237[B]','RXN-14237[M]'})
deads6=model6.mets(detectDeadEnds(model6));
extra=setdiff(deads6,deads)

model7=removeRxns(model6,{'RXN-7569[B]','RXN-7569[M]'})
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
[rxns,form]=findRxnsFromMets(model7,extra);
model7=removeRxns(model7,rxns)
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
[rxns,form]=findRxnsFromMets(model7,extra);
model7=removeRxns(model7,rxns)
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)

model7=removeRxns(model7,{'RXN-11895[B]','RXN-11895[M]'})
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
model7=removeRxns(model7,{'RXN-11907[B]','RXN-11907[M]'})
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
model7=removeRxns(model7,{'RXN-12142[B]','RXN-12142[M]'})

deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
[rxns,form]=findRxnsFromMets(model7,extra);
model7=removeRxns(model7,rxns)
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)

model7=removeRxns(model7,{'RXN-12140[B]','RXN-12140[M]'})
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
[rxns,form]=findRxnsFromMets(model7,extra);
model7=removeRxns(model7,rxns)
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
model7=removeRxns(model7,{'RXN-12141[B]','RXN-12141[M]'})
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)
[rxns,form]=findRxnsFromMets(model7,extra);
model7=removeRxns(model7,rxns)
deads7=model7.mets(detectDeadEnds(model7));
extra=setdiff(deads7,deads)

model8=removeRxns(model7,{'1.1.1.271-RXN[B]','1.1.1.271-RXN[M]'})
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
[rxns,form]=findRxnsFromMets(model8,extra);
model8=removeRxns(model8,rxns)
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
model8=removeRxns(model8,{'RXN-13492[B]','RXN-13492[M]'})
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
[rxns,form]=findRxnsFromMets(model8,extra);
model8=removeRxns(model8,rxns)
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
[rxns,form]=findRxnsFromMets(model8,extra);
model8=removeRxns(model8,rxns)
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
model8=removeRxns(model8,{'RXN-12872[B]','RXN-12872[M]'})
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
[rxns,form]=findRxnsFromMets(model8,extra);
model8=removeRxns(model8,rxns)
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)

model8=removeRxns(model8,{'RXN-13150[B]','RXN-13150[M]'})
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
model8=removeRxns(model8,{'RXN-13459[B]','RXN-13459[M]'})
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)

model8=removeRxns(model8, {'RXN-13642[B]','RXN-13642[M]'});
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)
[rxns,form]=findRxnsFromMets(model8,extra);
model8=removeRxns(model8,rxns)
deads8=model8.mets(detectDeadEnds(model8));
extra=setdiff(deads8,deads)

model9=removeRxns(model8,{'RXN-12882[B]','RXN-12882[M]'})
deads9=model9.mets(detectDeadEnds(model9));
extra=setdiff(deads9,deads)
[rxns,form]=findRxnsFromMets(model9,extra);
model9=removeRxns(model9,rxns)
deads9=model9.mets(detectDeadEnds(model9));
extra=setdiff(deads9,deads)
[rxns,form]=findRxnsFromMets(model9,extra);
model9=removeRxns(model9,rxns)
deads9=model9.mets(detectDeadEnds(model9));
extra=setdiff(deads9,deads)

[rxns,form]=findRxnsFromMets(model9,{'CPD-14649[cb]','CPD-14649[cm]'});
model10=removeRxns(model9,rxns)
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)
[rxns,form]=findRxnsFromMets(model10,extra);
model10=removeRxns(model10,rxns)
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)


model10=removeRxns(model10,{'RXN-14424[B]','RXN-14424[M]'})
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)
[rxns,form]=findRxnsFromMets(model10,extra);
model10=removeRxns(model10,rxns)
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)
[rxns,form]=findRxnsFromMets(model10,extra);
model10=removeRxns(model10,rxns)
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)
model10=removeRxns(model10,{'RXN-15372[B]','RXN-15372[M]'})
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)
model10=removeRxns(model10,{'RXN-14427[B]','RXN-14427[M]'})
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)
[rxns,form]=findRxnsFromMets(model10,extra);
model10=removeRxns(model10,rxns)
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)
[rxns,form]=findRxnsFromMets(model10,extra);
model10=removeRxns(model10,rxns)
deads10=model10.mets(detectDeadEnds(model10));
extra=setdiff(deads10,deads)

form='ATP[cb] + CPD-15317[cb]   ->   AMP[cb] + PROTON[cb] + PRPP[cb]';
model11 = addReaction(model10,'PRPPSYN-RXN_2[B]',form,[],0,0,1000);
form1='ATP[cm] + CPD-15317[cm]   ->   AMP[cm] + PROTON[cm] + PRPP[cm]';
model11 = addReaction(model11,'PRPPSYN-RXN_2[M]',form1,[],0,0,1000);
deads11=model11.mets(detectDeadEnds(model11));
extra=setdiff(deads11,deads)

model11=removeRxns(model11,{'RXN-14661[B]','RXN-14661[M]'})
deads11=model11.mets(detectDeadEnds(model11));
extra=setdiff(deads11,deads)
[rxns,form]=findRxnsFromMets(model11,extra);
model11=removeRxns(model11,rxns)
deads11=model11.mets(detectDeadEnds(model11));
extra=setdiff(deads11,deads)

model11=removeRxns(model11,{'RXN-14816[B]','RXN-14816[M]'})
deads11=model11.mets(detectDeadEnds(model11));
extra=setdiff(deads11,deads)
[rxns,form]=findRxnsFromMets(model11,extra);
model11=removeRxns(model11,rxns)
deads11=model11.mets(detectDeadEnds(model11));
extra=setdiff(deads11,deads)

model11=removeRxns(model11,{'RXN-8205[B]','RXN-8205[M]'})
deads11=model11.mets(detectDeadEnds(model11));
extra=setdiff(deads11,deads)

[rxns,form]=findRxnsFromMets(model11,{'CPD-16676[cb]','CPD-16676[cm]'});
model12=removeRxns(model11,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

[rxns,form]=findRxnsFromMets(model12,{'CPD-16681[cb]','CPD-16681[cm]'});
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

model12=removeRxns(model12,{'RXN-15532[B]','RXN-15532[M]'});
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)


model12=removeRxns(model12,{'RXN-16117[B]','RXN-16117[M]'});
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

model12=removeRxns(model12,{'RXN-16200[B]','RXN-16200[M]'});
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

model12=removeRxns(model12,{'RXN-16418[B]','RXN-16418[M]'});
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

model12=removeRxns(model12,{'RXN-16403[B]','RXN-16403[M]'});
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

model12=removeRxns(model12,{'RXN-16389[B]','RXN-16389[M]'});
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

model12=removeRxns(model12,{'RXN-16405[B]','RXN-16405[M]'});
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)
[rxns,form]=findRxnsFromMets(model12,extra);
model12=removeRxns(model12,rxns);
deads12=model12.mets(detectDeadEnds(model12));
extra=setdiff(deads12,deads)

model13=removeRxns(model12,{'RXN-16406[B]','RXN-16406[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)

model13=removeRxns(model13,{'RXN-18392[B]','RXN-18392[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)



model13=removeRxns(model13,{'RXN-18480[B]','RXN-18480[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
model13=removeRxns(model13,{'RXN-18482[B]','RXN-18482[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)

model13=removeRxns(model13,{'RXN-115[B]','RXN-115[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)

model13=removeRxns(model13,{'RXN-18668[B]','RXN-18668[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
model13=removeRxns(model13,{'RXN-18670[B]','RXN-18670[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)

[rxns,form]=findRxnsFromMets(model13,{'CPD-2183[cb]','CPD-2183[cm]'});
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)


[rxns,form]=findRxnsFromMets(model13,{'CPD-2187[cb]','CPD-2187[cm]'});
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)

model13=removeRxns(model13,{'RXN-3221[B]','RXN-3221[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
model13=removeRxns(model13,{'RXN-4725[B]','RXN-4725[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)

model13=removeRxns(model13,{'RXN-12861[B]','RXN-12861[M]'})
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
[rxns,form]=findRxnsFromMets(model13,extra);
model13=removeRxns(model13,rxns)
deads13=model13.mets(detectDeadEnds(model13));
extra=setdiff(deads13,deads)
form5='ATP[cb] + CO-A[cb] + CPD-3617[cb]   ->   AMP[cb] + CPD-10267[cb]   + PPI[cb]';
form6='ATP[cm] + CO-A[cm] + CPD-3617[cm]   ->   AMP[cm] + CPD-10267[cm]   + PPI[cm]';
model14 = addReaction(model13,'ACYLCOASYN-RXN_2[B]',form5,[],0,0,1000);
model14 = addReaction(model14,'ACYLCOASYN-RXN_2[M]',form6,[],0,0,1000);
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
model14=removeRxns(model14,{'RXN-8352[B]','RXN-8352[M]'})
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
model14=removeRxns(model14,{'RXN-4243[B]','RXN-4243[M]'})
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)

model14=removeRxns(model14,{'TYRAMINE-N-FERULOYLTRANSFERASE-RXN[B]','TYRAMINE-N-FERULOYLTRANSFERASE-RXN[M]'})
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)

model14=removeRxns(model14,{'RXN-4736[B]','RXN-4736[M]'})
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)

model14=removeRxns(model14,{'PSEUDOURIDINE-KINASE-RXN[B]','PSEUDOURIDINE-KINASE-RXN[M]'})
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)

model14=removeRxns(model14,{'RXN-6425[B]','RXN-6425[M]'})
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)

model14=removeRxns(model14,{'RXN-6687[B]','RXN-6687[M]'})
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)
[rxns,form]=findRxnsFromMets(model14,extra);
model14=removeRxns(model14,rxns)
deads14=model14.mets(detectDeadEnds(model14));
extra=setdiff(deads14,deads)

model15=removeRxns(model14,{'1.2.3.14-RXN[B]','1.2.3.14-RXN[M]'})
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)
[rxns,form]=findRxnsFromMets(model15,extra);
model15=removeRxns(model15,rxns)
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)
[rxns,form]=findRxnsFromMets(model15,extra);
model15=removeRxns(model15,rxns)
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)
[rxns,form]=findRxnsFromMets(model15,extra);
model15=removeRxns(model15,rxns)
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)
[rxns,form]=findRxnsFromMets(model15,extra);
model15=removeRxns(model15,rxns)
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)

model15=removeRxns(model15,{'RXN-7633[B]','RXN-7633[M]'})
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)
[rxns,form]=findRxnsFromMets(model15,extra);
model15=removeRxns(model15,rxns)
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)
[rxns,form]=findRxnsFromMets(model15,extra);
model15=removeRxns(model15,rxns)
deads15=model15.mets(detectDeadEnds(model15));
extra=setdiff(deads15,deads)

form7=' CPD-7025[cb] + ATP[cb] ->  PHYTYL-PYROPHOSPHATE[cb] + ADP[cb] + PROTON[cb]';
model16 = addReaction(model15,'RXN-7763[B]',form7,[],0,0,1000);

form8=' CPD-7025[cm] + ATP[cm] ->  PHYTYL-PYROPHOSPHATE[cm] + ADP[cm] + PROTON[cm]';
model16 = addReaction(model16,'RXN-7763[M]',form8,[],0,0,1000);
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model16=removeRxns(model16,{'RXN-7741[B]','RXN-7741[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model16=removeRxns(model16,{'RXN-8228[B]','RXN-8228[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

form9=' P-COUMAROYL-COA[cb] +  CPD-7138[cb] ->  CPD-7714[cb] +  CO-A[cb]';
model16 = addReaction(model16,'RXN-8170[B]',form9,[],0,0,1000);
form10=' P-COUMAROYL-COA[cm] +  CPD-7138[cm] ->  CPD-7714[cm] +  CO-A[cm]';
model16 = addReaction(model16,'RXN-8170[M]',form10,[],0,0,1000);
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model16=removeRxns(model16,{'RXN-715[B]','RXN-715[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
model16=removeRxns(model16,{'RXN-8139[B]','RXN-8139[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model16=removeRxns(model16,{'RXN-7982[B]','RXN-7982[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model16=removeRxns(model16,{'RXN-12292[B]','RXN-12292[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model16=removeRxns(model16,{'METHIONINE-GAMMA-LYASE-RXN[B]','METHIONINE-GAMMA-LYASE-RXN[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)


model16=removeRxns(model16,{'RXN-8171[B]','RXN-8171[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
model16=removeRxns(model16,{'RXN-8271[B]','RXN-8271[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)
[rxns,form]=findRxnsFromMets(model16,extra);
model16=removeRxns(model16,rxns)
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model16=removeRxns(model16,{'RXN-8313[B]','RXN-8313[M]','RXN-8310[B]','RXN-8310[M]'})
deads16=model16.mets(detectDeadEnds(model16));
extra=setdiff(deads16,deads)

model17=removeRxns(model16,{'RXN-8360[B]','RXN-8360[M]'})
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

model17=removeRxns(model17,{'RXN-8363[B]','RXN-8363[M]'})
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
[rxns,form]=findRxnsFromMets(model17,extra);
model17=removeRxns(model17,rxns)
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

model17=removeRxns(model17,{'RXN-8366[B]','RXN-8366[M]','RXN-8367[B]','RXN-8367[M]'})
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

model17=removeRxns(model17,{'RXN-15043[B]','RXN-15043[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

model17=removeRxns(model17,{'RXN-8471[B]','RXN-8471[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
model17=removeRxns(model17,{'RXN-8476[B]','RXN-8476[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
model17=removeRxns(model17,{'RXN-8491[B]','RXN-8491[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
[rxns,form]=findRxnsFromMets(model17,extra);
model17=removeRxns(model17,rxns)
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
[rxns,form]=findRxnsFromMets(model17,extra);
model17=removeRxns(model17,rxns)
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
[rxns,form]=findRxnsFromMets(model17,extra);
model17=removeRxns(model17,rxns)
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

model17=removeRxns(model17,{'RXN-8497[B]','RXN-8497[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
model17=removeRxns(model17,{'RXN-8508[B]','RXN-8508[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
model17=removeRxns(model17,{'RXN-11575[B]','RXN-11575[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
model17=removeRxns(model17,{'RXN-8684[B]','RXN-8684[M]'});
    deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

[rxns,form]=findRxnsFromMets(model17,extra);
model17=removeRxns(model17,rxns)
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)
[rxns,form]=findRxnsFromMets(model17,extra);
model17=removeRxns(model17,rxns)
deads17=model17.mets(detectDeadEnds(model17));
extra=setdiff(deads17,deads)

model18=removeRxns(model17,{'RXN-13433[B]','RXN-13433[M]', 'RXN-8681[B]','RXN-8681[M]'});
    deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)

model18=removeRxns(model18,{'RXN-8696[B]','RXN-8696[M]'});
    deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)

model18=removeRxns(model18,{'RXN-8909[B]','RXN-8909[M]'});
    deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
model18=removeRxns(model18,{'RXN-9021[B]','RXN-9021[M]'});
    deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
model18=removeRxns(model18,{'RXN-9193[B]','RXN-9193[M]'});
    deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)

model18=removeRxns(model18,{'RXN-5962[B]','RXN-5962[M]'});
    deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,{'CPD-5661[cb]','CPD-5661[cm]'});
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,{'CPD1F-118[cb]','CPD1F-118[cm]'});
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,{'CPD1F-115[cb]','CPD1F-115[cm]'});
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
deads=vertcat(deads,extra')
model18=removeRxns(model18,{'RXN1F-165[B]','RXN1F-165[M]'});
    deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)
[rxns,form]=findRxnsFromMets(model18,extra);
model18=removeRxns(model18,rxns)
deads18=model18.mets(detectDeadEnds(model18));
extra=setdiff(deads18,deads)

model19=removeRxns(model18,{'RXN1F-461[B]','RXN1F-461[M]'});
    deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)

[rxns,form]=findRxnsFromMets(model19,{'CPD1F-98[cb]','CPD1F-98[cm]'});
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)

model19=removeRxns(model19,{'RXNQT-4329[B]','RXNQT-4329[M]'});
    deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
model19=removeRxns(model19,{'RXNQT-4330[B]','RXNQT-4330[M]'});
    deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
model19=removeRxns(model19,{'RXNQT-4332[B]','RXNQT-4332[M]'});
    deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
model19=removeRxns(model19,{'RXNQT-4333[B]','RXNQT-4333[M]'});
    deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
model19=removeRxns(model19,{'RXNQT-4331[B]','RXNQT-4331[M]'});
    deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
model19=removeRxns(model19,{'RXNQT-4397[B]','RXNQT-4397[M]'});
    deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)
[rxns,form]=findRxnsFromMets(model19,extra);
model19=removeRxns(model19,rxns)
deads19=model19.mets(detectDeadEnds(model19));
extra=setdiff(deads19,deads)

form11=' CPD-11259[cb] +  WATER[cb] ->  CPD-731[cb] + ILE[cb]';
model20 = addReaction(model19,'RXN-18484[B]',form11,[],0,0,1000);
form12=' CPD-11259[cm] +  WATER[cm] ->  CPD-731[cm] + ILE[cm]';
model20 = addReaction(model20,'RXN-18484[M]',form12,[],0,0,1000);
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)
form13=' CPD-506[cb] +  WATER[cb] ->  INOSITOL-1-3-4-TRIPHOSPHATE[cb] +  Pi[cb]';
form14=' CPD-506[cm] +  WATER[cm] ->  INOSITOL-1-3-4-TRIPHOSPHATE[cm] +  Pi[cm]';
model20 = addReaction(model20,'RXN-8730[B]',form13,[],0,0,1000);
model20 = addReaction(model20,'RXN-8730[M]',form14,[],0,0,1000);
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)

model20=removeRxns(model20,{'2.7.1.151-RXN[B]','2.7.1.151-RXN[M]','3.1.3.56-RXN[B]','3.1.3.56-RXN[M]'});
    deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)
[rxns,form]=findRxnsFromMets(model20,extra);
model20=removeRxns(model20,rxns)
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)
[rxns,form]=findRxnsFromMets(model20,extra);
model20=removeRxns(model20,rxns)
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)
model20=removeRxns(model20,{'ATR_CPD-7866_[cb]_[cm]','ATR_CPD-6955_[cb]_[cm]','ATR_CPD-591_[cb]_[cm]','ATR_CPD-15924_[cb]_[cm]'});
deads20=model20.mets(detectDeadEnds(model20));
extra=setdiff(deads20,deads)
%%
to_rem=setdiff(model.rxns,model20.rxns);
added=setdiff(model20.rxns,model.rxns);
added=added(find(contains(added,'[B]')))
added_forms=printRxnFormula(model20,added)
completed=[added added_forms]
%% save('deads_to_remove.mat','to_rem')


%% save('to_add_completeness.mat','completed')
model=model20
clearvars -except model
 save('Constrained_unblocked_MAY25_Final.mat','model')
%%
% clear
% load('SC_constrained_Mar25_BIO.mat');
% model=model10;
% deads=model.mets(detectDeadEnds(model));
% model10=removeRxns(model9,'RXN-10763[M]')
% deadsnew=model10.mets(detectDeadEnds(model10))
% extra=setdiff(deadsnew,deads)
% [rxns9,form]=findRxnsFromMets(model9,'CPD-16551[cm]')
% pos=find(contains(model9.rxns,rxns9))
% names=model9.rxns(pos)
% flux=ro.v(pos)
% forms=printRxnFormula(model9,names);
% jeffyNEW=cell2table(horzcat(names,forms,num2cell(flux)))
% 
% 
% model10=removeRxns(model9,'RXN0-5398_1[M]')
% 
% 
