%% Generating M/B ratio for every reaction in the model which has gene expression data from the spreadsheet
% Feb 2025 with updated base sorghum model 
 addpath /Users/holla293/Documents/MSU_MacbookPro/Tn-Core-master/Most-Recent-Version
%model=readCbModel('model_0601.mat');
%model=readCbModel('model_balanced_databases_subs.mat');
%model=readCbModel('SorghumBicolorPMNV8_genes.mat')
%deads=model.mets(detectDeadEnds(model))
%model=readCbModel('SorghumBicolor_lessdeads.mat');
load('Sorghum_base_Mar25.mat')
 load('deads_to_remove.mat')
load('to_add_completeness.mat')
to_rem=erase(to_rem,'[B]');
to_rem=erase(to_rem,'[M]');
to_rem=unique(to_rem)
model=removeRxns(model,to_rem);
names=erase(completed(:,1),'[B]')
forms=strrep(completed(:,2),'[cb]','[c]');
for n=1:length(names)
  model = addReaction(model,names{n},forms{n},[],0,0,1000);
end
  model = tncore_fix(model)
model.rev=[];
for i=1:length(model.rxns)
    if model.lb(i) < 0
        model.rev(i,1)=1;
    else
        model.rev(i,1)=0;
    end
end

pos_rubpc=find(contains(model.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN'));
pos_ox=find(contains(model.rxns,'RXN-961'));
model.rules(pos_ox)=model.rules(pos_rubpc);
model.grRules(pos_ox)=model.grRules(pos_rubpc);

%
pos_trans=find(contains(model.rxns,'2TRANSKETO-RXN'));
pos_1trans=find(contains(model.rxns,'1TRANSKETO'));
model.rules(pos_trans)=model.rules(pos_1trans(1));
model.grRules(pos_trans)=model.grRules(pos_1trans(1));
% adding genes to NADP ME
genies=readcell('NADPMEgenes.xlsx');
rule={};
for n=1:length(genies)
    pos=find(contains(model.genes,genies(n)));
    if n==1
        rule=['(' genies{n}]
    elseif n==9
        rule=[rule ' or ' genies{n} ')' ]
    else
            rule=[rule ' or ' genies{n} ]

    end
end
pos=find(contains(model.rxns,'MALIC-NADP-RXN'));
model.grRules{pos(1)}=rule;
for n=1:length(genies)
    possy=find(contains(model.genes,genies(n)));
    if ~isempty(possy)
        rule=strrep(rule,genies{n},['x(' num2str(possy) ')'])
    else
        le=length(model.genes);
        model.genes{le+1}=genies{n}
        rule=strrep(rule,genies{n},['x(' num2str(le+1) ')'])

    end
end

rule=strrep(rule,'or','|');
model.rules{pos(1)}=rule;
% adding genes to PYRUVATEORTHOPHOSPHATE-DIKINASE-RXN
genies={'Sobic.001G326900','Sobic.009G132900'}
rule1={};
for n=1:length(genies)
    pos=find(contains(model.genes,genies(n)));
    if n==1
        rule1=['(' genies{n}]
    else
        rule1=[rule1 ' or ' genies{n} ')']
    end

   
end
pos=find(contains(model.rxns,'PYRUVATEORTHOPHOSPHATE-DIKINASE-RXN'));
model.grRules{pos(1)}=rule1;
for n=1:length(genies)
    possy=find(contains(model.genes,genies(n)));
    if ~isempty(possy)
        rule1=strrep(rule1,genies{n},['x(' num2str(possy) ')'])
    else
        le=length(model.genes);
        model.genes{le+1}=genies{n}
        rule1=strrep(rule1,genies{n},['x(' num2str(le+1) ')'])

    end
end

rule1=strrep(rule1,'or','|');
model.rules{pos(1)}=rule1;

% RXN-1106
genies={'Sobic.010G195500','Sobic.007G141200','Sobic.001G500100','Sobic.002G250100'}
rule3={};
for n=1:length(genies)
    pos=find(contains(model.genes,genies(n)));
    if n==1
        rule3=['(' genies{n}]
    elseif n==4
        rule3=[rule3 ' or ' genies{n} ')' ]
    else
            rule3=[rule3 ' or ' genies{n} ]

    end
end

pos=find(contains(model.rxns,'RXN-1106'));
model.grRules{pos(1)}=rule3;
for n=1:length(genies)
    possy=find(contains(model.genes,genies(n)));
    if ~isempty(possy)
        rule3=strrep(rule3,genies{n},['x(' num2str(possy) ')'])
    else
        le=length(model.genes);
        model.genes{le+1}=genies{n}
        rule3=strrep(rule3,genies{n},['x(' num2str(le+1) ')'])

    end
end

rule3=strrep(rule3,'or','|');
model.rules{pos(1)}=rule3;

%%

% % PRPPSYN-RXN
% genies={'Sobic.010G193600','Sobic.003G281200','Sobic.004G025800','Sobic.004G261500','Sobic.001G544600'};
% rule3={};
% for n=1:length(genies)
%     pos=find(contains(model.genes,genies(n)));
%     if n==1
%         rule3=['(' genies{n}]
%     elseif n==4
%         rule3=[rule3 ' or ' genies{n} ')' ]
%     else
%             rule3=[rule3 ' or ' genies{n} ]
% 
%     end
% end
% 
% pos=find(contains(model.rxns,'PRPPSYN-RXN_2'));
% model.grRules{pos(1)}=rule3;
% for n=1:length(genies)
%     possy=find(contains(model.genes,genies(n)));
%     if ~isempty(possy)
%         rule3=strrep(rule3,genies{n},['x(' num2str(possy) ')'])
%     else
%         le=length(model.genes);
%         model.genes{le+1}=genies{n}
%         rule3=strrep(rule3,genies{n},['x(' num2str(le+1) ')'])
% 
%     end
% end
% 
% rule3=strrep(rule3,'or','|');
% model.rules{pos(1)}=rule3;
% RXN-10435
% genies={'Sobic.010G241200','Sobic.002G093500','Sobic.005G108300','Sobic.009G249900','Sobic.003G015100','Sobic.003G015000','Sobic.003G224100'};
% rule3={};
% for n=1:length(genies)
%     pos=find(contains(model.genes,genies(n)));
%     if n==1
%         rule3=['(' genies{n}]
%     elseif n==4
%         rule3=[rule3 ' or ' genies{n} ')' ]
%     else
%             rule3=[rule3 ' or ' genies{n} ]
% 
%     end
% end
% 
% pos=find(contains(model.rxns,'RXN-10435'));
% model.grRules{pos(1)}=rule3;
% for n=1:length(genies)
%     possy=find(contains(model.genes,genies(n)));
%     if ~isempty(possy)
%         rule3=strrep(rule3,genies{n},['x(' num2str(possy) ')'])
%     else
%         le=length(model.genes);
%         model.genes{le+1}=genies{n}
%         rule3=strrep(rule3,genies{n},['x(' num2str(le+1) ')'])
% 
%     end
% end
% 
% rule3=strrep(rule3,'or','|');
% model.rules{pos(1)}=rule3;
% 
% % RXN-18484
% genies={'Sobic.002G093500','Sobic.002G093300','Sobic.002G093400','Sobic.006G150600','Sobic.001G021900','Sobic.001G022000','Sobic.008G016600','Sobic.003G190500','Sobic.003G273000','Sobic.010G241200'};
% rule3={};
% for n=1:length(genies)
%     pos=find(contains(model.genes,genies(n)));
%     if n==1
%         rule3=['(' genies{n}]
%     elseif n==4
%         rule3=[rule3 ' or ' genies{n} ')' ]
%     else
%             rule3=[rule3 ' or ' genies{n} ]
% 
%     end
% end
% 
% pos=find(contains(model.rxns,'RXN-18484'));
% model.grRules{pos(1)}=rule3;
% for n=1:length(genies)
%     possy=find(contains(model.genes,genies(n)));
%     if ~isempty(possy)
%         rule3=strrep(rule3,genies{n},['x(' num2str(possy) ')'])
%     else
%         le=length(model.genes);
%         model.genes{le+1}=genies{n}
%         rule3=strrep(rule3,genies{n},['x(' num2str(le+1) ')'])
% 
%     end
% end
% 
% rule3=strrep(rule3,'or','|');
% model.rules{pos(1)}=rule3;


% RXN-8170
genies={'Sobic.010G034800'};
rule3={};
for n=1:length(genies)
    pos=find(contains(model.genes,genies(n)));
    if n==1
        rule3=['(' genies{n}]
    elseif n==4
        rule3=[rule3 ' or ' genies{n} ')' ]
    else
            rule3=[rule3 ' or ' genies{n} ]

    end
end

pos=find(contains(model.rxns,'RXN-8170'));
model.grRules{pos(1)}=rule3;
for n=1:length(genies)
    possy=find(contains(model.genes,genies(n)));
    if ~isempty(possy)
        rule3=strrep(rule3,genies{n},['x(' num2str(possy) ')'])
    else
        le=length(model.genes);
        model.genes{le+1}=genies{n}
        rule3=strrep(rule3,genies{n},['x(' num2str(le+1) ')'])

    end
end

rule3=strrep(rule3,'or','|');
model.rules{pos(1)}=rule3;
% 
% 
% % RXN-8730
% genies={'Sobic.009G211400','Sobic.009G211600','Sobic.007G198300','Sobic.002G254700'};
% rule3={};
% for n=1:length(genies)
%     pos=find(contains(model.genes,genies(n)));
%     if n==1
%         rule3=['(' genies{n}]
%     elseif n==4
%         rule3=[rule3 ' or ' genies{n} ')' ]
%     else
%             rule3=[rule3 ' or ' genies{n} ]
% 
%     end
% end
% 
% pos=find(contains(model.rxns,'RXN-8730'));
% model.grRules{pos(1)}=rule3;
% for n=1:length(genies)
%     possy=find(contains(model.genes,genies(n)));
%     if ~isempty(possy)
%         rule3=strrep(rule3,genies{n},['x(' num2str(possy) ')'])
%     else
%         le=length(model.genes);
%         model.genes{le+1}=genies{n}
%         rule3=strrep(rule3,genies{n},['x(' num2str(le+1) ')'])
% 
%     end
% end
% 
% rule3=strrep(rule3,'or','|');
% model.rules{pos(1)}=rule3;

model = buildRxnGeneMat(model)

%% identifying rxns which break model if removed 
% load('constrained_simple.mat');
% const=final_model;
% diff=setdiff(model.rxns,const.rxns)
%  soly=[];
% for n=1:length(diff)
%     model1=model;
%     model1=removeRxns(model1,diff(n));
%     rop=optimizeCbModel(model1);
% soly=[soly,rop.f];
% end
% broke=soly<0.001;
% broken=diff(broke);
% 
% clear const final_model model1


%%
mes=readcell('mesophyllEXP.xlsx');
bundle=readcell('bundlesheathEXP.xlsx');

bundle_names=bundle(:,1);
mes_names=mes(:,1);
[A,B] = ismember(string(mes_names),model.genes);

genes_in_model=mes_names(A);
mes_exp=mes(A,2);
bundle_exp=bundle(A,2);
%[Mes_results,MesReacs] = findRxnsFromGenes(model,genes_in_model,1,1);

% isolating unique rxns which arent M or B
tr=find(contains(model.rxns,'ATR_'));
ex=find(contains(model.rxns,'EX_'))
bio=find(contains(model.rxns,'biomass'))
more=find(contains(model.rxns,'TRANS-RXN-'))
extr=vertcat(tr,ex,bio,more)
keep=model.rxns(extr);
%keep=vertcat(keep,{'PHOSPHORIBULOKINASE-RXN[B]'},broken)
mode={'PHOSPHORIBULOKINASE-RXN','CARBODEHYDRAT-RXN','HOMOSERKIN-RXN','RXN-1101','DUTP-PYROP-RXN'};%,'HOMOSERKIN-RXN','RXN-1101','RXN-1142','DUTP-PYROP-RXN'};
keep=vertcat(keep,transpose(mode))
extr=find(contains(model.rxns,keep))

left=setdiff(1:length(model.rxns),extr);
rxn=model.rxns(left);
% rxn=erase(rxn,'[M]');
% rxn=erase(rxn,'[B]');
% total number of unique rxns in the model 
rxns=unique(rxn);

mesmean=sum(cell2mat(mes_exp))/length(mes_exp);
bumean=sum(cell2mat(bundle_exp))/length(bundle_exp)
rxns_with_expression={};MB_ratio=[];zer={};rat=[];rxnno=[];no_genes={};
toolong=[];
for n=1:length(rxns)
    ma=find(strcmp(rxns(n), model.rxns));
    rxn_check=model.rxns(ma)
   % rxn_check1=erase(rxn_check,'[M]');
    %rxn_check1=erase(rxn_check1,'[B]');
    real_rxn=strmatch(rxns(n),rxn_check,'exact')
    geneList = findGenesFromRxns(model,rxn_check(real_rxn));
    mop=geneList{:}
    if isempty(mop) 
        %no_genes=[no_genes,rxns(n)]
    else
    end
mes_no=find(contains(genes_in_model,geneList{1,1}))
mesophyll_sum=sum(cell2mat(mes_exp(mes_no)));
bundle_sum=sum(cell2mat(bundle_exp(mes_no)));
mb=mesophyll_sum/(mesophyll_sum+bundle_sum);
if ~isempty(mes_no)    
    rxnno=[rxnno,ma(1)];
    %boop=geneList{1,1}
    %geneNameList=[geneNameList,boop{:}]
    rxns_with_expression=[rxns_with_expression,rxns(n)];
    if bundle_sum == 0 && mesophyll_sum == 0
    zer=[zer,rxns(n)]
    %MB_ratio=[MB_ratio,0]
    rat=[rat,0];
    elseif mesophyll_sum == 0 && bundle_sum > 0
    MB_ratio=[MB_ratio,0.0001]
        rat=[rat,0.0001];
    elseif mesophyll_sum > 0 && bundle_sum == 0
    MB_ratio=[MB_ratio,1]
    rat=[rat,1];
    else
    MB_ratio=[MB_ratio,mb]   
    rat=[rat,mb];
    end
else
     no_genes=[no_genes,rxns(n)]
 %  rat=[rat,0] 
end
    
end
no_genes=transpose(no_genes)
bubs=find(contains(model.rxns,no_genes));
nog=model.rxns(bubs)
nogene_rxns={};
for n=1:length(nog)
nogene_rxns=[nogene_rxns,[nog{n} '[M]']];
nogene_rxns=[nogene_rxns,[nog{n} '[B]']];
end

keepo=setdiff(keep,{'EX_Bio'})
keepoB=strcat(keepo,'[B]');
keepoM=strcat(keepo,'[M]');
keepy=vertcat(keepoM,keepoB);

keep=vertcat(keepy,nogene_rxns')
rxns_with_expression=transpose(rxns_with_expression)
MB_ratio=transpose(MB_ratio);
MB_ratio=sort(MB_ratio)
mu=mean(MB_ratio)
sd=std(MB_ratio)
y=pdf('Normal',MB_ratio,mu,sd)
% figure(1)
% plot(MB_ratio,y,'LineWidth',6);
% histfit(MB_ratio)
% set(gca,'LineWidth',4,'FontSize',60)
%  x_width=30 ;y_width=24;
%  set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%  print('histogram','-depsc','-loose');
%qqplot(MB_ratio,y)

%UB=mu+1.9*sd
%LB=mu-1.9*sd
% UB=mu+2*sd
% LB=mu-2*sd
UB=mu+2*sd
LB=mu-2*sd
(sum(MB_ratio< LB)+sum(MB_ratio>  UB))/length(MB_ratio)
constrained_rxns={};
for n=1:length(rxns_with_expression)
    if rat(n)>UB
    constrained_rxns=[constrained_rxns,[rxns_with_expression{n} '[M]']]
    elseif rat(n)<UB && rat(n)>LB
     constrained_rxns=[constrained_rxns,[rxns_with_expression{n} '[M]']]
     constrained_rxns=[constrained_rxns,[rxns_with_expression{n} '[B]']]
    elseif rat(n)<LB && rat(n)>0
    constrained_rxns=[constrained_rxns,[rxns_with_expression{n} '[B]']]
    else
    end
end

% testing how many are M only, B only
moo=constrained_rxns(find(contains(constrained_rxns,'[M]')));
boo=constrained_rxns(find(contains(constrained_rxns,'[B]')));
moo1=erase(moo,'[M]');
boo1=erase(boo,'[B]');

monly=setdiff(moo1,boo1);
bonly=setdiff(boo1,moo1);
% 28 M only, 37 B only. 
    

constrained_rxns=transpose(constrained_rxns)
all_rxns=vertcat(constrained_rxns,keep)
all_rxns=unique(all_rxns)
all_rxns=setdiff(all_rxns,'biomass');
%% creating cell type reactions


no_letter={};rxns_in_model=[];noo=[];rxns_unique={};form_1={}; rxnNameList_1={};rev_1=[];lowerBoundList_1=[];upperBoundList_1=[];subSystemList_1={};
grRuleList_1={};
for n=1:length(all_rxns)
       B=find(contains(all_rxns(n),'[B]'));
    M=find(contains(all_rxns(n),'[M]'));
   % rxns_unique=[rxns_unique,all_rxns(n)]
   if B==1 
    rxnname=erase(all_rxns(n),'[B]');
    pos1=find(strcmp(rxnname, model.rxns));
    mesFormulas = printRxnFormula(model,model.rxns(pos1), false);
    % adding cell location to mets
    mesFormulas=strrep(mesFormulas,'[c]','[cb]'); 
    mesFormulas=strrep(mesFormulas,'[m]','[mb]'); 
    mesFormulas=strrep(mesFormulas,'[s]','[sb]'); 
    mesFormulas=strrep(mesFormulas,'[t]','[tb]'); 
    mesFormulas=strrep(mesFormulas,'[i]','[ib]'); 

    form_1=[form_1,mesFormulas];
    rxnNameList_1 = [rxnNameList_1,all_rxns{n}];
    rev_1=[rev_1,model.rev(pos1)];
    lowerBoundList_1 = [lowerBoundList_1,model.lb(pos1)];
    upperBoundList_1= [upperBoundList_1,model.ub(pos1)];
    subSystemList_1= [subSystemList_1,model.subSystems(pos1)];
    grRuleList_1 = [grRuleList_1,model.grRules(pos1)];
   elseif M==1
       rxnname=erase(all_rxns(n),'[M]');
    pos1=find(strcmp(rxnname, model.rxns));
    mesFormulas = printRxnFormula(model,model.rxns(pos1), false);
    mesFormulas=strrep(mesFormulas,'[c]','[cm]'); 
    mesFormulas=strrep(mesFormulas,'[m]','[mm]'); 
    mesFormulas=strrep(mesFormulas,'[s]','[sm]'); 
    mesFormulas=strrep(mesFormulas,'[t]','[tm]'); 
    mesFormulas=strrep(mesFormulas,'[i]','[im]'); 
    form_1=[form_1,mesFormulas];
    rxnNameList_1 = [rxnNameList_1,all_rxns(n)];
    rev_1=[rev_1,model.rev(pos1)];
    lowerBoundList_1 = [lowerBoundList_1,model.lb(pos1)];
    upperBoundList_1= [upperBoundList_1,model.ub(pos1)];
    subSystemList_1= [subSystemList_1,model.subSystems(pos1)];
    grRuleList_1 = [grRuleList_1,model.grRules(pos1)];
   else 
       no_letter=[no_letter,all_rxns(n)];
   end

end


rxns_in_model=transpose(rxns_in_model);
rxnNameList_1=transpose(rxnNameList_1);
rxns_unique=transpose(rxns_unique);
form_1=transpose(form_1);
rev_1=transpose(rev_1);
lowerBoundList_1=transpose(lowerBoundList_1);
upperBoundList_1=transpose(upperBoundList_1);
subSystemList_1=transpose(subSystemList_1);
grRuleList_1=transpose(grRuleList_1);
% keep is names of kept rxns
% extr is rxn nos of keep
%% Create the model
% Get the relevant information for bundle sheath
% ormulas_keep = printRxnFormula(model, model.rxns(extr), false);
% rxnNameList_keep = model.rxnNames(extr);
% rxnNameList=vertcat(rxnNameList_keep,rxnNameList_1);
% rxnList = vertcat(ormulas_keep,form_1);
% revFlagList_keep = model.rev(extr);
% revFlagList = vertcat(revFlagList_keep,rev_1);
% lowerBoundList_keep = model.lb(extr);
% lowerBoundList = vertcat(lowerBoundList_keep,lowerBoundList_1);
% upperBoundList_keep= model.ub(extr);
% upperBoundList = vertcat(upperBoundList_keep,upperBoundList_1);
% subSystemList_keep= model.subSystems(extr);
% subSystemList = vertcat(subSystemList_keep,subSystemList_1)
% grRuleList_keep = model.grRules(extr);
% grRuleList = vertcat(grRuleList_keep,grRuleList_1);
% %geneNameList = model.genes;
% geneNameList=model.genes;

%% fine up to 2000
% rxnAbrList2=rxnAbrList(1:2490);
% rxnNameList2=rxnNameList(1:2490);
% rxnList2=rxnList(1:2490);
% revFlagList2=revFlagList(1:2490);
% lowerBoundList2=lowerBoundList(1:2490);
% upperBoundList2=upperBoundList(1:2490);
% subSystemList2=subSystemList(1:2490);
% grRuleList2=grRuleList(1:2490);
% 


final_model = createModel(rxnNameList_1,rxnNameList_1,form_1,rev_1, ...
lowerBoundList_1, upperBoundList_1, subSystemList_1, grRuleList_1);
final_model=tncore_remove(final_model);
 %final_model = changeObjective(final_model, 'biomass')
%%
% save('fmfeb18.mat','final_model')
%load('fmfeb18.mat')

 deads=final_model.mets(detectDeadEnds(final_model))
% 1623 deads
e=final_model.mets(find(contains(final_model.mets,'[e]')))
exchange=findRxnsFromMets(final_model,e,'ProducersOnly',1)
formex=printRxnFormula(final_model,exchange);

load('constrained_leaf_noredundancies.mat')
e5=model5.mets(find(contains(model5.mets,'[e]')))
e_torem=setdiff(e,e5);
exchange=findRxnsFromMets(final_model,e_torem)
formex=printRxnFormula(final_model,exchange);

final_model=removeRxns(final_model,exchange)
exchange5=findRxnsFromMets(model5,e5)
formex5=printRxnFormula(model5,exchange5)
pos5=find(contains(model5.rxns,exchange5));
exer=model5.rxns(pos5);
lb=model5.lb(pos5);
ub=model5.ub(pos5);
formex5=printRxnFormula(model5,exer)

for n=1:length(formex5)
  final_model = addReaction(final_model,exer{n},formex5{n},[],0,lb(n),ub(n));
end

 deads=final_model.mets(detectDeadEnds(final_model))
% 1594 now 
pos1=find(contains(model5.rxns,'_[cb]_[cm]'));
pos2=find(contains(model5.rxns,'_[cm]_[cb]'));
names=model5.rxns(pos1);
lb=model5.lb(pos1);
ub=model5.ub(pos1);
formex5=printRxnFormula(model5,names)
for n=1:length(formex5)
  final_model = addReaction(final_model,names{n},formex5{n},[],0,lb(n),ub(n));
end
%1584 deads now- this may be adding new deads
 deads=final_model.mets(detectDeadEnds(final_model))

% all the transporters from previous models have been added
%final_model=removeMetabolites(final_model)

e=final_model.mets(find(contains(final_model.mets,'[e]')))
intersect(deads,e)
% no exchange mets are dead

bio=findMetsFromRxns(final_model,{'biomass[M]','biomass[B]'})
debi=intersect(deads,bio)
pos=find(contains(final_model.mets,'biomass[cm]'))
final_model.mets{pos}='biomass[m]';
pos=find(contains(final_model.mets,'biomass[cb]'))
final_model.mets{pos}='biomass[b]';
form='0.5 biomass[m] + 0.5 biomass[b] ->';
  final_model = addReaction(final_model,'biomass',form,[],0,0,1000);
 final_model = changeObjective(final_model, 'biomass')

 deads=final_model.mets(detectDeadEnds(final_model))


%%
% final_model = addReaction(final_model,'ATR_GLYCOLLATE_CYTOSOL_PLASTID-STR[M]',{'GLYCOLLATE[cm]','GLYCOLLATE[sm]'},[-1 1],0,-1000,1000);
% rxnformula='ATP[sm] + RIBULOSE-5P[sm]   ->   ADP[sm] + PROTON[sm] + D-RIBULOSE-15-P2[sm]';
% final_model = addReaction(final_model,'PHOSPHORIBULOKINASE-RXN[M]',rxnformula,[],0,0,1000);
% rxnformula='GAP[sm] + CPD-19339[sm]   <=>   CPD-15318[sm] + XYLULOSE-5-PHOSPHATE[sm]';
%final_model = addReaction(final_model,'1TRANSKETO-RXN_2[M]',rxnformula,[],0,-1000,1000);

%needed for biomass rxn
%rxnFormula='2-KETOGLUTARATE[cm] + VAL[cm]   <=>   GLT[cm] + 2-KETO-ISOVALERATE[cm]';
%final_model = addReaction(final_model,'BRANCHED-CHAINAMINOTRANSFERVAL-RXN[M]',rxnFormula,[],0,-1000,1000);
%rxnFormula1='DATP[cm] + 2 WATER[cm]   ->   DAMP[cm] + 2 PROTON[cm] + 2 Pi[cm]';
%final_model = addReaction(final_model,'RXN-14195[M]',rxnFormula1,[],0,0,1000);
%rxnFormula2='WATER[cm] + CPD-18[cm]   ->   CO-A[cm] + LINOLEIC_ACID[cm]  + PROTON[cm]';
%final_model = addReaction(final_model,'LINOLEOYL-RXN[M]',rxnFormula2,[],0,0,1000);
%form='CPD-18[cm] + Glycerolipids_Compound[cm]   <=>   CO-A[cm] + Linoleoyl-groups_Compound[cm]';
%final_model = addReaction(final_model,'RXN-16045_1[M]',form,[],0,-1000,1000);
%form2='WATER[cm] + Alpha-linolenoyl-groups_Compound[cm] ->   LINOLENIC_ACID[cm] + PROTON[cm] + Glycerolipids_Compound[cm]'
%final_model = addReaction(final_model,'Lipase_Alpha-linolenoyl-groups_1[M]',form2,[],0,0,1000);
% form3='WATER[cm] + CPD-1777[cm]   ->   CONIFERYL-ALCOHOL[cm]  + ALPHA-GLUCOSE[cm]';
% final_model = addReaction(final_model,'CONIFERIN-BETA-GLUCOSIDASE-RXN_1[M]',form3,[],0,0,1000);
% form4='WATER[cm] + CPD-1777[cm]   ->   CONIFERYL-ALCOHOL[cm] + GLC[cm]';
% final_model = addReaction(final_model,'CONIFERIN-BETA-GLUCOSIDASE-RXN_2[M]',form4,[],0,0,1000);
% form5='CONIFERYL-ALCOHOL[cm] + CPD-12575[cm]   <=>   CPD-1777[cm]   + PROTON[cm] + UDP[cm]'
% final_model = addReaction(final_model,'2.4.1.111-RXN[M]',form5,[],0,-1000,1000);
% form6='WATER[cm] + MANNOSE-6P[cm]   ->   CPD-12601[cm] + Pi[cm]';
% final_model = addReaction(final_model,'SUGAR-PHOSPHATASE-RXN_2[M]',form6,[],0,0,1000);
%form7='GDP[cm] + Red-Thioredoxin_Compound[cm]   ->   DGDP[cm] + WATER[cm]   + Ox-Thioredoxin_Compound[cm]';
%final_model = addReaction(final_model,'GDPREDUCT-RXN_1[M]',form7,[],0,0,1000);
%form8='DGDP[cm] + ATP[cm]   <=>   DGTP[cm] + ADP[cm]';
%final_model = addReaction(final_model,'DGDPKIN-RXN[M]',form8,[],0,-1000,1000);
%form9='DGTP[cm] + 2 WATER[cm]   ->   DGMP[cm] + 2 PROTON[cm] + 2 Pi[cm]';
%final_model = addReaction(final_model,'RXN-14208[M]',form9,[],0,0,1000);
%form10='2-KETOGLUTARATE[cm] + ILE[cm]   <=>   GLT[cm]   + 2-KETO-3-METHYL-VALERATE[cm] ';
%final_model = addReaction(final_model,'BRANCHED-CHAINAMINOTRANSFERILEU-RXN[M]',form10,[],0,-1000,1000);
% form11='2-KETO-3-METHYL-VALERATE[cb] <=> 2-KETO-3-METHYL-VALERATE[cm]';
% final_model = addReaction(final_model,'ATR_2-KETO-3-METHYL-VALERATE[cm]',form11,[],0,-1000,1000);

% photosynthesis
% form12='WATER[sb]   <=>   WATER[tb]';
% final_model = addReaction(final_model,'ATR_WATER_PLASTID-STR_THY-LUM[B]',form12,[],0,-1000,1000);
% zn uptake 
% form13='ATP[cm] + WATER[cm] + ZN+2[cm]   ->   ADP[cm] + PROTON[cm] + Pi[cm]   + ZN+2[e]';
% final_model = addReaction(final_model,'RXN0-5205[M]',form12,[],0,0,1000);
