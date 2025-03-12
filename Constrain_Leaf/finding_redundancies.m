%% comparing FVA results at different drought strengths

% FVA pipeline
%addpath /mnt/home/holla293/Documents/cobratoolbox
%    initCobraToolbox;

%changeCobraSolver('glpk')
load('SC_constrained_unblocked.mat');
%spread=readcell('spreadrxns.xlsx');
%load('listofimportantrxns.mat');
%load('importantlist_CAlimit.mat')
load('importantlist_oneCAleft.mat')
spread=un;
model=model1;
model = changeRxnBounds(model,'EX_Light_Compound_EXTRACELLULAR', 3000, 'u');
model=removeRxns(model,{'ATR_CARBON-DIOXIDE_CYTOSOL_PLASTID-STR[B]'})
model = changeRxnBounds(model,'RXN-15479_1[M]', 2000, 'u');
model = changeRxnBounds(model,'ATPSYN-RXN-Plastid[M]', 2000, 'u');
light=find(contains(model.rxns,'ATR_Light'))
model.lb(light)=-5000;
model.ub(light)=5000;


bo2=find(contains(model.rxns,'_2['));
bo3=find(contains(model.rxns,'_3['));
bo4=find(contains(model.rxns,'_4['));
bo5=find(contains(model.rxns,'_5['));
bo6=find(contains(model.rxns,'_6['));
bo7=find(contains(model.rxns,'_7['));
bo8=find(contains(model.rxns,'_8['));
bo9=find(contains(model.rxns,'_9['));
bo10=find(contains(model.rxns,'_10['));
bo11=find(contains(model.rxns,'_11['));
bob=vertcat(bo2,bo3,bo4,bo5,bo6,bo7,bo8,bo9,bo10,bo11);
% isolate which reactions in bob list have _1 version and is present in M
% and B, and which is only in one cell type
repeats=strrep(model.rxns(bob),'_2[B]','[B]');
repeats=strrep(repeats,'_2[M]','[M]');
repeats=strrep(repeats,'_3[M]','[M]');
repeats=strrep(repeats,'_4[M]','[M]');
repeats=strrep(repeats,'_5[M]','[M]');
repeats=strrep(repeats,'_6[M]','[M]');
repeats=strrep(repeats,'_7[M]','[M]');
repeats=strrep(repeats,'_8[M]','[M]');
repeats=strrep(repeats,'_9[M]','[M]');
repeats=strrep(repeats,'_10[M]','[M]');
repeats=strrep(repeats,'_11[M]','[M]');
repeats=strrep(repeats,'_3[B]','[B]');
repeats=strrep(repeats,'_4[B]','[B]');
repeats=strrep(repeats,'_5[B]','[B]');
repeats=strrep(repeats,'_6[B]','[B]');
repeats=strrep(repeats,'_7[B]','[B]');
repeats=strrep(repeats,'_8[B]','[B]');
repeats=strrep(repeats,'_9[B]','[B]');
repeats=strrep(repeats,'_10[B]','[B]');
repeats=strrep(repeats,'_11[B]','[B]');
repeats=unique(repeats)


% THIOESTER-RXN has 11 versions

nocell=erase(repeats,'[B]');
nocell=erase(nocell,'[M]')
oncellrepeat=[];
for n=1:length(nocell);
    pos=find(contains(nocell,nocell(n)));
    if length(pos)<2
    oncellrepeat=[oncellrepeat,n];
    else
    end
end
onecell=nocell(oncellrepeat);
list=unique(nocell);
list=setdiff(list,onecell);
nop=find(contains(model.rxns,onecell))
% removing 'NITRATE-REDUCTASE-NADPORNOPH-RXN_2[B]' to keep only reactions
% which are repeats -- onecell is list to fix which is only one cell type
onecell=setdiff(onecell,{'NITRATE-REDUCTASE-NADPORNOPH-RXN' });

% does the list have _1 version of all in the list? 
hasone=[];
for n=1:length(list)
    pos=find(contains(model.rxns,[list{n} '_1']))
    if ~isempty(pos)
        hasone=[hasone,n]
    else
    end
end
% 122 out of 135 reactions have original _1 and a repeat
notone=setdiff(list,list(hasone));
onlyrepeats=setdiff(list,notone); % this is the list to use 
% testing if rxns in notone list aren't actually repeated
rep=[];
for n=1:length(notone)
    pos=find(contains(model.rxns,notone(n)));
    model.rxns(pos)
    if length(pos)>2
        rep=[rep,n]
    else
    end
end
% all rxns in notone aren't repeated except for in both M and B, so we can
% ignore those 

% isolating if any reactions only have a repeat because of location change
onlyrepeats=setdiff(onlyrepeats,{'RXN-11334','ACONITATEDEHYDR-RXN','ACONITATEHYDR-RXN','RXN66-526'})
big_redundant={};
for n=1:length(onlyrepeats)
    unmets={};rxnlist={};blem={};
pos=find(contains(model.rxns,onlyrepeats(n)))
% B only since all rxns are in B and M
roo=model.rxns(pos);
Bonly=roo(find(contains(roo,'[B]')));
pos=find(contains(model.rxns,Bonly));
mets=findMetsFromRxns(model,model.rxns(pos(1)));
mets1=erase(mets,'[cb]')
mets1=erase(mets1,'[mb]')
mets1=erase(mets1,'[sb]')
mets1=erase(mets1,'[sm]')
mets1=erase(mets1,'[cm]')
mets1=erase(mets1,'[mm]')



for i=2:length(pos)
   current=findMetsFromRxns(model,model.rxns(pos(i)));
current1=erase(current,'[cb]');
current1=erase(current1,'[mb]');
current1=erase(current1,'[sb]');
current1=erase(current1,'[sm]');
current1=erase(current1,'[cm]');
current1=erase(current1,'[mm]')

un1=setdiff(mets1,current1);
un2=setdiff(current1,mets1);
if i==2
    for j=1:length(un1)
unmets(j,1)=un1(j)
blem=[blem,un1(j)]
    end
else
end
if isempty(un2)
    unmets(1,i)={''};
else
for j=1:length(un2)
 unmets(j,i)=un2(j)
 blem=[blem,un2(j)];
end
end
end

if ~isempty(blem)
    
for m=1:length(pos)
    metab={};
    for p=1:length(unmets(:,m))
        if ~isempty(unmets{p,m})
        metab=[metab,unmets{p,m}]
        else
        end
    end
    %metab=unmets(:,m);
    %bloo=metab{:};
    if length(metab)>0
    metabolites=findMetsFromRxns(model,model.rxns(pos(m))) 
    mm=metabolites(find(contains(metabolites,metab)));
    rxns1=findRxnsFromMets(model,mm);
    rxn1=erase(rxns1,'[B]');
for j=1:length(rxns1)
 rxnlist(j,m)=rxns1(j)
end
    else 
        rxnlist(1,m)={''};
    end
end
% generated list of rxns for each set of mets, now need to identify if the
% list contains unique rxns and if they are linked to biomass mets 
sizes=[];
for k=1:length(pos)
    rxnsize=[];
    for mi=1:length(rxnlist(:,k))
       if ~isempty(rxnlist{mi,k})
    rxnsize=[rxnsize,mi];
       else
       end
    end
    if isempty(rxnsize)
        listlength=0;
    else
    listlength=max(rxnsize);
    end
sizes=[sizes,listlength]
end
[M,I]=max(sizes);
% I is location of largest list
mainlist=rxnlist(:,I);
listy=erase(mainlist,'[B]')
listy=erase(listy,'[M]');
listy=erase(listy,'_1');
listy=erase(listy,'_2');
listy=erase(listy,'_3');
listy=erase(listy,'_4');
listy=erase(listy,'_5');
listy=erase(listy,'_6');
listy=erase(listy,'_7');
listy=erase(listy,'_8');
listy=erase(listy,'_9');
listy=erase(listy,'_10');
listy=erase(listy,'_11');
%%
redundant=[];
for k=1:length(pos)
if k==I
else
     match=[];
    for mi=1:length(rxnlist(:,k))
       if ~isempty(rxnlist{mi,k})
    reaction=rxnlist{mi,k};
    reaction=erase(reaction,'[B]')
reaction=erase(reaction,'[M]');
reaction=erase(reaction,'_1');
reaction=erase(reaction,'_2');
reaction=erase(reaction,'_3');
reaction=erase(reaction,'_4');
reaction=erase(reaction,'_5');
reaction=erase(reaction,'_6');
reaction=erase(reaction,'_7');
reaction=erase(reaction,'_8');
reaction=erase(reaction,'_9');
reaction=erase(reaction,'_10');
reaction=erase(reaction,'_11');
    
    pos_list=strmatch(reaction,listy,'exact');
    if ~isempty(pos_list)
        match=[match,mi]
       else
    end
       else


       end

    end
  if length(match)==sizes(k) && sizes(k)~=0
        redundant=[redundant,k]
    else
    end  
end
end
red=model.rxns(pos(redundant));
metsforrxn=findMetsFromRxns(model,red)

actmet=metsforrxn(find(contains(metsforrxn,unmets(:,redundant))));
reac1=findRxnsFromMets(model,actmet);
big_redundant=[big_redundant,reac1']
else
    % for k=2:length(pos)
    %     big_redundant=[big_redundant,model.rxns(pos(k))]
    % end
end
end
big_redundant=horzcat(big_redundant,{'RXN66-526_1[B]','RXN66-526_2[B]','RXN66-526_3[B]','GCVMULTI-RXN_1[B]','GCVMULTI-RXN_3[B]','1.5.1.20-RXN_1[B]', ...
    '3-CH3-2-OXOBUTANOATE-OH-CH3-XFER-RXN_2[B]','GCVMULTI-RXN_1[B]','GLYOHMETRANS-RXN_1[B]','HOMOCYSMET-RXN[B]','RXN-2881_1[B]','GCVMULTI-RXN_3[B]','GLYOHMETRANS-RXN_3[B]'  });
big_redundant=unique(big_redundant);

% only 59 rxns which are redundant, check if they have M or B version
big=erase(big_redundant,'[B]');
hmm=find(contains(model.rxns,big))
to_remove=model.rxns(hmm)
to_remove=setdiff(to_remove,{'PHOSMANMUT-RXN_2[B]','PHOSMANMUT-RXN_2[M]','SUGAR-PHOSPHATASE-RXN_2[B]','SUGAR-PHOSPHATASE-RXN_2[M]', ...
    'ALDOSE-1-EPIMERASE-RXN[B]','ALDOSE-1-EPIMERASE-RXN[M]','GLUCOSE-6-PHOSPHATE-1-EPIMERASE-RXN[B]','GLUCOSE-6-PHOSPHATE-1-EPIMERASE-RXN[M]','RXN-14815[B]','RXN-14815[M]','TREHALA-RXN[B]','TREHALA-RXN[M]'})
% soly=[];model1=model;
% for n=1:length(to_remove)
% model1=removeRxns(model1,to_remove(n));
% sol=optimizeCbModel(model1);
% soly=[soly,sol.f]
% end

model2=removeRxns(model,to_remove)
blocked1=model2.mets(detectDeadEnds(model2));
rxns_extra_remove=findRxnsFromMets(model2,blocked1);
model3=removeRxns(model2,rxns_extra_remove)
blocked2=model3.mets(detectDeadEnds(model3));
rxns_extra_remove2=findRxnsFromMets(model3,blocked2);
model4=removeRxns(model3,rxns_extra_remove2)
blocked3=model4.mets(detectDeadEnds(model4))
rxns_extra_remove3=findRxnsFromMets(model4,blocked3);

model5=removeRxns(model4,rxns_extra_remove3)

save('constrained_leaf_noredundancies.mat','model5')

% match identifies how many rxns in list is matched with largest list
% if this is same length as size then all rxns in that list can be added to
% a remove list, still need to add a condin
% tional if the reaction is
% involved in biomass 


% 
% 
% %% old version of code below 
% for n=1:length(onecell)
% pos=find(contains(model.rxns,onecell(n)))
% mets=findMetsFromRxns(model,model.rxns(pos(1)))
% if sum(contains(mets,'[cb]'))==length(mets) 
% mets_unique={};
% for i=2:length(pos)
%     current=findMetsFromRxns(model,model.rxns(pos(i)))
%     if sum(contains(current,'[cb]'))==length(current) || sum(contains(current,'[cm]'))==length(current)
% 
%     mets1=setdiff(mets,current);
% metsy=vertcat(mets1,setdiff(findMetsFromRxns(model,model.rxns(pos(i))),mets));
% 
%     mets_unique=[mets_unique,metsy]
%     else
%     end
% end
% firstmets=intersect(mets,mets_unique);
% firstrxns=findRxnsFromMets(model,firstmets);
% restmets=setdiff(mets_unique,firstmets);
% restrxns=findRxnsFromMets(model,restmets)
%  unique_mets=()
% if sum(contains(mets,'[cb]'))==length(mets) || sum(contains(mets,'[cm]'))==length(mets)
%   mets=findMetsFromRxns(model,model.rxns(pos(1)))
% mets_unique={};
% for 2:length(pos)
%      mets1=setdiff(mets,findMetsFromRxns(model,model.rxns(pos(i))));
% metsy=vertcat(mets1,setdiff(findMetsFromRxns(model,model.rxns(pos(i))),mets));
% 
%     mets_unique=[mets_unique,metsy]
% else
%     % rxn which has component other than cytosol
%     i
% end
% end
% 
% end
% 
% %% checking if the repeated reaction is linked to unique reactions
% % lets test with onecell first since it is only in one cell type
% for n=1:length(onecell)
%     pos=find(contains(model.rxns,onecell(n)));
%          mets=findMetsFromRxns(model,model.rxns(pos(1)))
%     for i=1:length(pos)
%      mets2=findMetsFromRxns(model,model.rxns(pos(i)));
%      diff=setdiff(mets2,mets)
