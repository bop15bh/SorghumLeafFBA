%% Varied co2 analysis with additional CA's removed to get  new list of important rxns

% FVA pipeline
%addpath /mnt/home/holla293/Documents/cobratoolbox
    %initCobraToolbox;

%changeCobraSolver('glpk')
%load('SC_constrained_unblocked.mat');
load('constrained_leaf_noredundancies.mat')
model=model5;
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
model=removeRxns(model, {'RXN-11334_1[B]','RXN-11334_1[M]','RXN-11334_3[B]','RXN-11334_3[M]','TRANS-RXN-194[B]','TRANS-RXN-194[M]','TRANS-RXN-206[B]','TRANS-RXN-206[M]','TRANS-RXN-213[B]','TRANS-RXN-213[M]','OXALOACETASE-RXN[B]','OXALOACETASE-RXN[M]'})
cpd1=find(contains(model.mets,'CPD-16551[cb]'));
cpd2=find(contains(model.mets,'CPD-16551[cm]'));
model.mets(cpd1)={'CPD-15317[cb]'}
model.mets(cpd2)={'CPD-15317[cm]'}
og=optimizeCbModel(model);
position_ca=find(contains(model.rxns,'CARBODEHYDRAT-RXN[M]'));
popo=[40.9091 122.7 192.3 225 270.3777];

for n=1:8
   model1=model;
   model1 = changeRxnBounds(model1,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -270.3777*(n/8), 'l'); 
pep=optimizeCbModel(model1);


%jeff=cell2table(horzcat(model.rxns,num2cell(og.v),num2cell(minFlux),num2cell(maxFlux),num2cell(pep.v),num2cell(minFlux_pep),num2cell(maxFlux_pep)));
% rxns=table2cell(jeff(:,1))
% rxns(1)='';

dro=pep.v;

og_flux=og.v;
rxns=model.rxns;
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
og_flux_pos=og_flux>1
rxns1=rxns(og_flux_pos);
og_flux1=og_flux(og_flux_pos);
dro1=dro(og_flux_pos);
% og_flux_min1=og_flux_min(og_flux_pos);
% og_flux_max1=og_flux_max(og_flux_pos);
% dro_min1=dro_min(og_flux_pos);
% dro_max1=dro_max(og_flux_pos);
dro_pos=dro1>1;
rxns2=rxns1(dro_pos);
og_flux2=og_flux1(dro_pos);
dro2=dro1(dro_pos);
% og_flux_min2=og_flux_min1(dro_pos);
% og_flux_max2=og_flux_max1(dro_pos);
% dro_min2=dro_min1(dro_pos);
% dro_max2=dro_max1(dro_pos);

neg_pos=og_flux<-1;

rxnsneg=rxns(neg_pos);
og_fluxneg=og_flux(neg_pos);
droneg=dro(neg_pos);
% og_flux_minneg=og_flux_min(neg_pos);
% og_flux_maxneg=og_flux_max(neg_pos);
% dro_minneg=dro_min(neg_pos);
% dro_maxneg=dro_max(neg_pos);
dro_posneg=droneg<-1;
rxnsneg2=rxnsneg(dro_posneg);
droneg2=droneg(dro_posneg);
og_fluxneg2=og_fluxneg(dro_posneg);
% og_flux_minneg2=og_flux_minneg(dro_posneg);
% og_flux_maxneg2=og_flux_maxneg(dro_posneg);
% dro_minneg2=dro_minneg(dro_posneg);
% dro_maxneg2=dro_maxneg(dro_posneg);

rxns=vertcat(rxns2,rxnsneg2);
og_flux=vertcat(og_flux2,og_fluxneg2);
% og_flux_min=vertcat(og_flux_min2,og_flux_minneg2);
% og_flux_max=vertcat(og_flux_max2,og_flux_maxneg2);
dro=vertcat(dro2,droneg2);
% dro_min=vertcat(dro_min2,dro_minneg2);
% dro_max=vertcat(dro_max2,dro_maxneg2);
FC=abs(dro./og_flux)

%jeff=cell2table(horzcat(rxns,num2cell(og_flux),num2cell(og_flux_min),num2cell(og_flux_max),num2cell(dro),num2cell(dro_min),num2cell(dro_max),num2cell(FC)));

fcinc=FC>1.1;
og_1=og_flux(fcinc);
dro_1=dro(fcinc);
list_1=rxns(fcinc);
fcdec=FC<0.9;
list_2=rxns(fcdec);
og_2=og_flux(fcdec);
dro_2=dro(fcdec);
%whole_list(:,n)=vertcat(list_2,list_1);
boop=vertcat(list_1,list_2);
boop2=vertcat(og_1,og_2);
boop3=vertcat(dro_1,dro_2);
for i=1:length(boop)
whole_list(i,n)=boop(i);
og_flurx(i,n)=boop2(i);
dro_flurx(i,n)=boop3(i);
end
end
%%
list1={};list2={};list3={};list4={};list5={};list6={};list7={};list8={};
for n=1:length(whole_list(:,1))
if ~isempty(whole_list{n,1})
    list1=[list1,whole_list(n,1)];
else 
end
end
for n=1:length(whole_list(:,2))
    if ~isempty(whole_list{n,2})
    list2=[list2,whole_list(n,2)];
    else
    end
end
 for n=1:length(whole_list(:,3))
if   ~isempty(whole_list{n,3})
    list3=[list3,whole_list(n,3)];
else
end
 end
 for n=1:length(whole_list(:,4))
if ~isempty(whole_list{n,4})
    list4=[list4,whole_list(n,4)];
else
end
 end
 for n=1:length(whole_list(:,5))
    if ~isempty(whole_list{n,5})
    list5=[list5,whole_list(n,5)];
    else
    end
 end
 for n=1:length(whole_list(:,6))
    if ~isempty(whole_list{n,6})
    list6=[list6,whole_list(n,6)];
    else
    end
 end
 for n=1:length(whole_list(:,7))
    if ~isempty(whole_list{n,7})
     
    list7=[list7,whole_list(n,7)];
    else
    end
 end
%  for n=1:length(whole_list(:,8))
%     if ~isempty(whole_list{n,8})
%     list8=[list8,whole_list(n,8)];
% else
% end
%  end
 list1=transpose(list1);
  list2=transpose(list2)
 list3=transpose(list3)
 list4=transpose(list4)
 list5=transpose(list5)
 list6=transpose(list6)
 list7=transpose(list7)
 %list8=transpose(list8)
all=vertcat(list1,list2,list3,list4,list5,list6,list7);%,list8);

un=unique(all)
%save('importantlist_oneCAleft.mat','un')
save('important_list_Nov24.mat','un')
%save('workspace_FVA_varied.mat')
%clear
%% exporting the rxns and FC for each iteration
%jeff=cell2table(horzcat(whole_list(:,1),num2cell(og_flurx(:,1)),num2cell(dro_flurx(:,1)),whole_list(:,2),num2cell(og_flurx(:,2)),num2cell(dro_flurx(:,2)),whole_list(:,3),num2cell(og_flurx(:,3)),num2cell(dro_flurx(:,3)),whole_list(:,4),num2cell(og_flurx(:,4)),num2cell(dro_flurx(:,4)),whole_list(:,5),num2cell(og_flurx(:,5)),num2cell(dro_flurx(:,5)),whole_list(:,6),num2cell(og_flurx(:,6)),num2cell(dro_flurx(:,6)),whole_list(:,7),num2cell(og_flurx(:,7)),num2cell(dro_flurx(:,7)),whole_list(:,8),num2cell(og_flurx(:,8)),num2cell(dro_flurx(:,8))));
%writetable(jeff,'varied_FVA_withflux.txt','Delimiter','\t')