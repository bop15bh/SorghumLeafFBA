%% Generating M/B ratio for every reaction in the model which has gene expression data from the spreadsheet
% Oct 2023 
% addpath /Users/holla293/Documents/MSU_MacbookPro/Tn-Core-master/Most-Recent-Version
%model=readCbModel('model_0601.mat');
model=readCbModel('model_balanced_databases_subs.mat');
%%
pos_rubpc=find(contains(model.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
pos_ox=find(contains(model.rxns,'RXN-961[B]'));
model.rules(pos_ox)=model.rules(pos_rubpc);
model.grRules(pos_ox)=model.grRules(pos_rubpc);

%
pos_trans=find(contains(model.rxns,'2TRANSKETO-RXN'));
pos_1trans=find(contains(model.rxns,'1TRANSKETO'));
model.rules(pos_trans)=model.rules(pos_1trans);
model.grRules(pos_trans)=model.grRules(pos_1trans);
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
pos=find(contains(model.rxns,'MALIC-NADP-RXN_1'));
model.grRules{pos(1)}=rule;
model.grRules{pos(2)}=rule;
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
model.rules{pos(2)}=rule;
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
model.grRules{pos(2)}=rule1;

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
model.rules{pos(2)}=rule1;

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
model.grRules{pos(2)}=rule3;
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
model.rules{pos(2)}=rule3;
model = buildRxnGeneMat(model)

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

%%
% fixing Rubisco genes

% %load('essential.mat')
% ro=optimizeCbModel(model);
% nonzer_flux=ro.v~=0;
% flux_rxns=model.rxns(nonzer_flux);
model = tncore_fix(model)
model.rev=[];
for i=1:length(model.rxns)
    if model.lb(i) < 0
        model.rev(i,1)=1;
    else
        model.rev(i,1)=0;
    end
end

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
mode={'PHOSPHORIBULOKINASE-RXN[B]','CARBODEHYDRAT-RXN[B]'};
keep=vertcat(keep,transpose(mode))

left=setdiff(1:length(model.rxns),extr);
rxn=model.rxns(left);
rxn=erase(rxn,'[M]');
rxn=erase(rxn,'[B]');
% total number of unique rxns in the model 
rxns=unique(rxn);


rxns_with_expression={};MB_ratio=[];zer={};rat=[];rxnno=[];no_genes={};
toolong=[];
for n=1:length(rxns)
    ma=find(contains(model.rxns,rxns(n)));
    rxn_check=model.rxns(ma)
    rxn_check1=erase(rxn_check,'[M]');
    rxn_check1=erase(rxn_check1,'[B]');
    real_rxn=strmatch(rxns(n),rxn_check1,'exact')
    geneList = findGenesFromRxns(model,rxn_check(real_rxn));
    mop=geneList{:}
    if isempty(mop) 
        no_genes=[no_genes,rxns(n)]
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
    if bundle_sum==0 && mesophyll_sum==0
    zer=[zer,rxns(n)]
    %MB_ratio=[MB_ratio,0]
    rat=[rat,0];
    elseif mesophyll_sum==0 && bundle_sum>0
    MB_ratio=[MB_ratio,0.0001]
        rat=[rat,0.0001];
    elseif mesophyll_sum>0 && bundle_sum==0
    MB_ratio=[MB_ratio,1]
    rat=[rat,1];
    else
    MB_ratio=[MB_ratio,mb]   
    rat=[rat,mb];
    end
else
 %  rat=[rat,0] 
end
    
end
no_genes=transpose(no_genes)
bubs=find(contains(model.rxns,no_genes));
nog=model.rxns(bubs)
keep=vertcat(keep,nog)
rxns_with_expression=transpose(rxns_with_expression)
MB_ratio=transpose(MB_ratio);
MB_ratio=sort(MB_ratio)
mu=mean(MB_ratio)
sd=std(MB_ratio)
y=pdf('Normal',MB_ratio,mu,sd)
figure(1)
plot(MB_ratio,y,'LineWidth',6);
histfit(MB_ratio)
set(gca,'LineWidth',4,'FontSize',60)
 x_width=30 ;y_width=24;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 print('histogram','-depsc','-loose');
%qqplot(MB_ratio,y)

%UB=mu+1.3*sd
%LB=mu-1.3*sd
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
    elseif rat(n)<LB
    constrained_rxns=[constrained_rxns,[rxns_with_expression{n} '[B]']]
    else
    end
end

constrained_rxns=transpose(constrained_rxns)
all_rxns=vertcat(constrained_rxns,keep)
all_rxns=unique(all_rxns)
%% isolating rxns which match original model
% if they dont match, write new formula in new compartment
rxns_in_model=[];noo=[];rxns_unique={};form_1={}; rxnNameList_1={};rev_1=[];lowerBoundList_1=[];upperBoundList_1=[];subSystemList_1={};
grRuleList_1={};
for n=1:length(all_rxns)
    pos=strmatch(all_rxns(n),model.rxns,'exact')
    %pos=find(contains(model.rxns,all_rxns(n)));
if length(pos)==1
   noo=[noo,n];
   rxns_in_model=[rxns_in_model,pos];
else  
    B=find(contains(all_rxns(n),'[B]'));
    M=find(contains(all_rxns(n),'[M]'));
    rxns_unique=[rxns_unique,all_rxns(n)]
   if B==1 
    rxnname=erase(all_rxns(n),'[B]');
    pos1=find(contains(model.rxns,rxnname));
    mesFormulas = printRxnFormula(model,model.rxns(pos1), false);
    mesFormulas=strrep(mesFormulas,'m]','b]');   
    form_1=[form_1,mesFormulas];
    rxnNameList_1 = [rxnNameList_1,model.rxnNames(pos1)];
    rev_1=[rev_1,model.rev(pos1)]
    lowerBoundList_1 = [lowerBoundList_1,model.lb(pos1)];
    upperBoundList_1= [upperBoundList_1,model.ub(pos1)];
    subSystemList_1= [subSystemList_1,model.subSystems(pos1)];
    grRuleList_1 = [grRuleList_1,model.grRules(pos1)];
   elseif M==1
       rxnname=erase(all_rxns(n),'[M]');
    pos1=find(contains(model.rxns,[rxnname{:},'[']));
    mesFormulas = printRxnFormula(model,model.rxns(pos1), false);
    mesFormulas=strrep(mesFormulas,'b]','m]'); 
    form_1=[form_1,mesFormulas];
    rxnNameList_1 = [rxnNameList_1,model.rxnNames(pos1)];
    rev_1=[rev_1,model.rev(pos1)]
    lowerBoundList_1 = [lowerBoundList_1,model.lb(pos1)];
    upperBoundList_1= [upperBoundList_1,model.ub(pos1)];
    subSystemList_1= [subSystemList_1,model.subSystems(pos1)];
    grRuleList_1 = [grRuleList_1,model.grRules(pos1)];
   end
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
ormulas = printRxnFormula(model, model.rxns(rxns_in_model), false);
ormulas_keep = printRxnFormula(model, model.rxns(extr), false);
rxnNameList = model.rxnNames(rxns_in_model);
rxnNameList_keep = model.rxnNames(extr);
rxnNameList=vertcat(rxnNameList,rxnNameList_1);
rxnAbrList = model.rxns(rxns_in_model);
rxnAbrList = vertcat(rxnAbrList,rxns_unique);
rxnList = ormulas;
rxnList = vertcat(ormulas,form_1);
revFlagList = model.rev(rxns_in_model);
revFlagList_keep = model.rev(extr);
revFlagList = vertcat(revFlagList,rev_1);
lowerBoundList = model.lb(rxns_in_model);
lowerBoundList_keep = model.lb(extr);
lowerBoundList = vertcat(lowerBoundList,lowerBoundList_1);
upperBoundList= model.ub(rxns_in_model);
upperBoundList_keep= model.ub(extr);
upperBoundList = vertcat(upperBoundList,upperBoundList_1);
subSystemList= model.subSystems(rxns_in_model);
subSystemList_keep= model.subSystems(extr);
subSystemList = vertcat(subSystemList,subSystemList_1)
grRuleList = model.grRules(rxns_in_model);
grRuleList_keep = model.grRules(extr);
grRuleList = vertcat(grRuleList,grRuleList_1);
%geneNameList = model.genes;
geneNameList=model.genes;

final_model = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList, ...
lowerBoundList, upperBoundList, subSystemList, grRuleList);
final_model=tncore_remove(final_model);
 final_model = changeObjective(final_model, 'biomass')

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
