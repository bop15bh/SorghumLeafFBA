% pipeline to add old model reactions to get the model to grow 
SC_multigene_constrain
diff=setdiff(model.rxns,final_model.rxns)
final=final_model;
for n=1:length(diff)
pos=strmatch(diff(n),model.rxns,'exact');
formy = printRxnFormula(model, model.rxns(pos), false);
lb=model.lb(pos);
ub=model.ub(pos);
final=addReaction(final,diff{n},formy{:},[],0,lb,ub);
end
ro=optimizeCbModel(final)
transporters_add
%diff=setdiff(diff,broken)


% test=setdiff(diff,{'ATPSYN-RXN-Plastid[M]','DIHYDROOROT-RXN[M]','DUTP-PYROP-RXN[M]','OROPRIBTRANS-RXN[M]' ...
% 'OROTPDECARB-RXN[M]','RXN-1143_2[M]','RXN-14014_2[M]','RXN-6883_1[B]','RXN0-5224[M]' ...
% 'RXN0-723[M]','RXN0-882[M]','RXN0-884[M]','SEDOBISALDOL-RXN[M]','SEDOHEPTULOSE-BISPHOSPHATASE-RXN_1[M]' ...
% 'SUCROSE-SYNTHASE-RXN[M]','SULFITE-REDUCTASE-FERREDOXIN-RXN_1[M]','THREONINE-ALDOLASE-RXN[M]' ...
% 'UDPREDUCT-RXN_1[M]','dihydroorotate dehydrogenase(NAD+)[M]'});
% test=setdiff(diff,{'2-ISOPROPYLMALATESYN-RXN[B]','2TRANSKETO-RXN[M]','ADPREDUCT-RXN_1[B]','DIHYDROOROT-RXN[M]','DUTP-PYROP-RXN[M]' ...
%     'GUANYL-KIN-RXN[B]','HOMOSERKIN-RXN[M]','MALONYL-COA-ACP-TRANSACYL-RXN_1[M]','OROPRIBTRANS-RXN[M]','OROTPDECARB-RXN[M]' ...
%     'RXN-1101[M]','RXN-1106[M]','RXN-13398[B]','RXN-14014_2[M]','RXN-5682[M]','SEDOHEPTULOSE-BISPHOSPHATASE-RXN_1[M]' ...
%     'SULFITE-REDUCTASE-FERREDOXIN-RXN_1[M]','TRANS-CINNAMATE-4-MONOOXYGENASE-RXN[B]','UDP-GLUCURONATE-4-EPIMERASE-RXN[B]' ...
%     'UGD-RXN[B]','TRIOSEPISOMERIZATION-RXN_2[B]','RXN-14182[M]'})
% soly=[];    model1=final;    
% for n=1:1:length(test)
% model1=removeRxns(model1,test(n));
% rop=optimizeCbModel(model1);
% soly=[soly,rop.f]
% end



%% group for 2sd stringency
test=setdiff(diff,{'HOMOSERKIN-RXN[B]','HOMOSERKIN-RXN[M]','RXN-1106[M]','UGD-RXN[B]'})
soly=[];   model1=final;
for n=1:1:length(test)
     
model1=removeRxns(model1,test(n));
rop=optimizeCbModel(model1);
soly=[soly,rop.f]
end


%% comparing list of known rxns with constrained 1.3sd
%list=readcell('list.xlsx');% 75 rxns 

%missing=setdiff(list,model1.rxns); % 8 rxns missing from the list

%% 
%% comparing list of known rxns with constrained 1.3sd
% list=readcell('list.xlsx');% 75 rxns 
% 
% missing=setdiff(list,model1.rxns); % 8 rxns missing from the list

%% only {'GLY3KIN-RXN[B]'} missing when using 2sd, 2.1 also still 1 
%% 4 rxns missing at 1.6 
%% 2 rxns missing at 1.8 

not_list=readcell('not_list.xlsx');% 
not_list{5}='RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[M]';
missing=intersect(not_list,model1.rxns); % 8 rxns missing from the list
%% 14 of the ones which shouldn't be present are present at 2sd
  %  maybekeep={'CATAL-RXN[M]','CARBODEHYDRAT-RXN[B]' };
  maybekeep={'CATAL-RXN[M]'};
missing=setdiff(missing,maybekeep)

soly=[];        
 model2=model1;
for n=1:1:length(missing)
   
model2=removeRxns(model2,missing(n));
rop=optimizeCbModel(model2);
soly=[soly,rop.f]
end
model2=removeRxns(model2,{'PSII-RXN_1[B]','RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[M]'})
save('SC_constrainedmodelFINAL.mat','model2')
%% removing genes from rxns which have 1 which work
% pos=find(contains(model.rxns,'GLYOHMETRANS-RXN_1'));
% model.rules(pos)={'x(117)'};
% model.grRules(pos)={'Sobic.001G097100'};
% pos1=find(contains(model.rxns,'GPH-RXN'));
% model.rules(pos1)={'x(2395)'};
% model.grRules(pos1)={'Sobic.006G130300'};
% 
% pos2=find(contains(model.rxns,'PEPCARBOX-RXN'));
% model.rules(pos2)={'x(3466)'};
% model.grRules(pos2)={'Sobic.010G160700'};
% pos3=find(contains(model.rxns,'PHOSPHORIBULOKINASE-RXN'));
% model.rules(pos3)={'x(1894)'};
% model.grRules(pos3)={'Sobic.004G272100'};
% pos4=find(contains(model.rxns,'RXN-961['));
% model.rules(pos4)={'x(2048)'};
% model.grRules(pos4)={'Sobic.005G042000'};
