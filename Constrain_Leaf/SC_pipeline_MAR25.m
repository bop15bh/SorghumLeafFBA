SC_multigene_constrain_FEB25
load('constrained_leaf_noredundancies.mat')
final=final_model;
transporters_add
%diff=setdiff(model5.rxns,final.rxns)
% 
% for n=1:length(diff)
% pos=strmatch(diff(n),model5.rxns,'exact');
% formy = printRxnFormula(model5, model5.rxns(pos), false);
% lb=model5.lb(pos);
% ub=model5.ub(pos);
% final=addReaction(final,diff{n},formy{:},[],0,lb,ub);
% end
ro=optimizeCbModel(final)

%% 
not_list=readcell('not_list.xlsx');% 
list=readcell('list_Mar25.xlsx')
intersect(list,final.rxns)
not_list{5}='RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[M]';
missing=intersect(not_list,final.rxns); % 10 present which shouldnt be
%% 14 of the ones which shouldn't be present are present at 2sd
  %  maybekeep={'CATAL-RXN[M]','CARBODEHYDRAT-RXN[B]' };
  maybekeep={'CATAL-RXN[M]'};
missing=setdiff(missing,maybekeep)

soly=[];      model2=final;
      
for n=1:1:length(missing)

model2=removeRxns(model2,missing(n));
rop=optimizeCbModel(model2);
soly=[soly,rop.f]
end
model2=removeRxns(model2,{'PSII-RXN[B]','RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[M]'})

%%
% test=diff;
% test=setdiff(diff,{'DUTP-PYROP-RXN[B]'})%,'RXN-1142[B]'})
% % test=setdiff(diff,{'1TRANSKETO-RXN_2[M]','RXN-1106[M]','UGD-RXN[B]','DUTP-PYROP-RXN[B]','HOMOSERKIN-RXN[B]','HOMOSERKIN-RXN[M]','RXN-1101[B]','RXN-1142[B]','TRANSALDOL-RXN_2[M]'});
% %test=setdiff(diff,{'1TRANSKETO-RXN_2[M]','RXN-1106[M]','UGD-RXN[B]','DUTP-PYROP-RXN[B]','HOMOSERKIN-RXN[B]','HOMOSERKIN-RXN[M]','RXN-1101[B]','RXN-1142[B]','RXN0-5224[B]','TRANSALDOL-RXN_2[M]'})
% %test=setdiff(diff,{'1TRANSKETO-RXN_2[M]','RXN-1106[M]','RXN-7577[B]','RXN-7577[M]','RXN-7578[B]','RXN-7578[M]', ...
% %    'RXN0-5224[B]','TRANSALDOL-RXN_2[M]','UGD-RXN[B]','DUTP-PYROP-RXN[B]','HOMOSERKIN-RXN[B]','HOMOSERKIN-RXN[M]','RXN-1101[B]','RXN-1142[B]'})
% %test=setdiff(diff,{'RXN-1106[M]','RXN-7577[B]','RXN-7577[M]','RXN-7578[B]','RXN-7578[M]','UGD-RXN[B]'})
% soly=[];    model1=model2;  
% model1=removeRxns(model1,test)
% for n=1:1:length(test)
% model1=removeRxns(model1,test(n));
% rop=optimizeCbModel(model1);
% soly=[soly,rop.f]
% end


model=model2;
%save('Constrained_model_Mar25.mat','model')
save('Constrained_model_MAY25.mat','model')
% %%
% load('Constrained_model_Mar25.mat')
% load('constrained_leaf_noredundancies.mat')
% e=model5.mets(find(contains(model5.mets,'[e]')));
% [rxnso,rxnformo]=findRxnsFromMets(model5,e);
% e_rem=model.mets(find(contains(model.mets,'[e]')));
% [rxns_rem,form_rem]=findRxnsFromMets(model,e_rem);
% model=removeRxns(model,rxns_rem);
% for n=1:length(rxnso)
% if contains(rxnformo{n},'<=>')
%     model = addReaction(model,rxnso{n},rxnformo{n},[],0,-1000,1000);
% else
%         model = addReaction(model,rxnso{n},rxnformo{n},[],0,0,1000);
% end
% end
% 
% 
% %%
% model=removeRxns(model, {'RXN-19748_1[B]'})
% form='MAL[sb] + NADP[sb]   ->   CARBON-DIOXIDE[sb] + NADPH[sb] + PYRUVATE[sb] ';
%         model = addReaction(model,'MALIC-NADP-RXN[B]',form,[],0,0,1000);
% 
% %%
% allform=printRxnFormula(model,model.rxns);
% 
% cb=allform(find(contains(allform,'[cb]')));
% cbpos=model.rxns(find(contains(allform,'[cb]')));
% 
% cbcm=cb(find(contains(cb,'[cm]')));
% cbcmpos=cbpos(find(contains(cb,'[cm]')));
% 
% pos=find(contains(model.rxns,cbcmpos));
% names=model.rxns(pos)
% flux=sol.v(pos)
% 
