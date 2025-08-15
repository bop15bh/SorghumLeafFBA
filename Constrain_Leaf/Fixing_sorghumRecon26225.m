%% 
clear
load('SorghumBicolor_lessdeads.mat') 
model1=changeRxnBounds(model1,'RXN-15468',-1000,'l')
model1=removeRxns(model1,'RXN-15815_2')
form='PROTON[c] + NADP[c] + 2 E-[c]   ->   NADPH[c]';
model1 = addReaction(model1,'RXN-19772',form,[],0,0,1000);
form1='2 PROTON[c] + FUM[c] + 2 E-[c]   ->   SUC[c]';
model1 = addReaction(model1,'RXN0-5245_2',form1,[],0,0,1000);

olddead=model1.mets(detectDeadEnds(model1))
model2=removeRxns(model1,{'ABC-34-RXN','ABC-25-RXN','ABC-25-RXN ','ABC-27-RXN'})
% form='-> DIHYDROXY-ACETONE-PHOSPHATE[c]'
% model2 = addReaction(model2,'TEST',form,[],0,0,1000);
% that fixes it, so clearly cannot be produced

% [rxnList, rxnFormulaList] = findRxnsFromMets(model1, 'DIHYDROXY-ACETONE-PHOSPHATE[c]')


vModel = load('CobraModel_sorghumbicolor_04-Jun-2019.mat');
  rxns=vModel.rxns;
 
extra=setdiff(rxns,model2.rxns)
extra_pos=find(contains(vModel.rxns,extra))
extra_to_check=vModel.rxns(extra_pos)
formo=printRxnFormula(vModel,extra_to_check)
formo=strrep(formo,'_e0','[e]')
formo=strrep(formo,'_c0','[c]')
formo=strrep(formo,'_m0','[m]')
formo=strrep(formo,'_s0','[s]')
formo=strrep(formo,'_t0','[t]')
formo=strrep(formo,'_i0','[i]')
lbe=vModel.lb(extra_pos);
ube=vModel.ub(extra_pos);

model4=model2;
for n=1:length(extra_to_check)
     model4 = addReaction(model4,extra_to_check{n},formo{n},[],0,lbe(n),ube(n));
end
   soly=[]; model5=model4;
   rem=extra_to_check;
   rem=setdiff(extra_to_check,{'RXN-1142','RXN-11354' ...
       'RXN-15468_1'})%,'PUTRESCINE-OXIDASE-RXN','RXN-11354'})
% only PUTRESCINE-OXIDASE-RXN isnt in sorghumcyc 8
   model5=removeRxns(model5,rem);

%    for n=1:length(rem)
% model5=removeRxns(model5,rem(n));
% sol=optimizeCbModel(model5);
% soly=[soly,sol.f]
% end

fdcomp=model5.mets(find(contains(model5.mets,'ferredoxins')));
comp=model5.mets(find(contains(model5.mets,'_Compound[')));
fd=intersect(fdcomp,comp);
 [rxnList, rxnFormulaList] = findRxnsFromMets(model5, fd)
model5=removeRxns(model5,rxnList)
 %model5=removeRxns(model5,{'RXN-15468_1','ATR_Oxidized-ferredoxins_Compound_PLASTID-STR','ATR_Reduced-ferredoxins_Compound_PLASTID-STR'})
 %[rxnList, rxnFormulaList] = findRxnsFromMets(model5, fd)
form1=erase(rxnFormulaList(1),'_Compound');
form2=strrep(rxnFormulaList(2),'Oxidized-ferredoxins_Compound','Oxidized-ferredoxins');
form2=strrep(form2,'Reduced-ferredoxins_Compound','Reduced-ferredoxins');
 
model5 = addReaction(model5,rxnList{1},form1{:},[],0,0,1000);
     model5 = addReaction(model5,rxnList{2},form2{:},[],0,0,1000);
form1='Oxidized-ferredoxins[s]  <=> Oxidized-ferredoxins[c] ';                                                                       
form2='Reduced-ferredoxins[s]  <=> Reduced-ferredoxins[c] ';                                                                         

model5 = addReaction(model5,'ATR_Oxidized-ferredoxins_PLASTID-STR',form1,[],0,-1000,1000);
model5 = addReaction(model5,'ATR_Reduced-ferredoxins_PLASTID-STR',form2,[],0,-1000,1000);
deads=model5.mets(detectDeadEnds(model5))


%% fixing remaining dead ends
norxns=[];
for n=1:length(deads)
    [rxnList, rxnFormulaList] = findRxnsFromMets(model5, deads(n));
    if isempty(rxnList)
    norxns=[norxns,n]
    else
    end
end
to_rem=deads(norxns);
model5=removeMetabolites(model5,to_rem);
deads=model5.mets(detectDeadEnds(model5))
% mets that need s to c exchange 
s_mets={'AMMONIA[s]','CPD1F-130[s]','HOMO-CYS[s]','HS[s]' 'CPD-15056[s]','Stearoyl-ACPs[s]','STEARIC_ACID[s]','METHYL-GLYOXAL[s]','D-LACTATE[s]','Cytochromes-C-Oxidized[s]','Cytochromes-C-Reduced[s]','SACCHAROPINE[s]',...
    'LYS[s]','ACP[s]','Stearoyl-ACPs[s]','GERANYLGERANYL-PP[s]','CPD-4211[s]'   };
model6=model5;
for n=1:length(s_mets)
    newmet=strrep(s_mets{n},'[s]','[c]');
    form=[s_mets{n} ' <=> ' newmet];
    model6 = addReaction(model6,['ATR_' erase(s_mets{n},'[s]') '_CYTOSOL_PLASTID'],form,[],0,-1000,1000);
end
% 266 deads up to here

% mets that need m to c
m_mets={'PPI[m]','2-KETOGLUTARATE[m]','AMMONIUM[m]','CPD-9957[m]','L-DELTA1-PYRROLINE_5-CARBOXYLATE[m]'};
for n=1:length(m_mets)
    newmet=strrep(m_mets{n},'[m]','[c]');
    form=[m_mets{n} ' <=> ' newmet];
    model6 = addReaction(model6,['ATR_' erase(m_mets{n},'[m]') '_CYTOSOL_MIT'],form,[],0,-1000,1000);
end
form1='FE+2[e] <=> FE+2[c]';
model6 = addReaction(model6,'ATR_FE+2_CYTOSOL_EXTRACELLULAR',form1,[],0,-1000,1000);

% rxns to remove 
ro_rem_rxns={'TRANS-RXN4LZ-42','TRANS-RXN4LZ-6860','3.1.2.21-RXN_2',...
    'TRANS-RXN-328','HEXOKINASE-RXN_4','ATR_CPD-20699_[s]_[c]','RXN-15029',...
    'NITRITE-REDUCTASE-CYTOCHROME-RXN_1','ATR_HYDROGEN-PEROXIDE_[s]_[c]','RXN-12940',...
     'RXN-7578','RXN-7577','RXN-14503','RXN-12276','RXN-12277','RXN-12204','SELENOCYSTEINE-LYASE-RXN_1',...
     'SINAPATE-1-GLUCOSYLTRANSFERASE-RXN','TROPINE-DEHYDROGENASE-RXN','RXN-8014','RXN-12299','RXN-1346',...
     'RXN-11541_1'};


% 279 deads atm 
model6=removeRxns(model6,ro_rem_rxns)
deads=model6.mets(detectDeadEnds(model6))

form11='XYLITOL[c] + NADP[c]   <=>   XYLOSE[c] + PROTON[c] + NADPH[c]';
model6 = addReaction(model6,'RXN-8773',form11,[],0,-1000,1000);

form33='WATER[c] + TRYPTAMINE[c] + OXYGEN-MOLECULE[c]  ->   HYDROGEN-PEROXIDE[c] + AMMONIUM[c] + INDOLE_ACETALDEHYDE[c]'; 
model6 = addReaction(model6,'RXN-1401',form33,[],0,0,1000);

form44='WATER[c] + INDOLE_ACETALDEHYDE[c] + OXYGEN-MOLECULE[c] ->   HYDROGEN-PEROXIDE[c] + PROTON[c] + INDOLE_ACETATE_AUXIN[c]';
model6 = addReaction(model6,'INDOLE-3-ACETALDEHYDE-OXIDASE-RXN',form44,[],0,0,1000);

form55='ATP[c] + WATER[c] + MET[c]   <=>   PPI[c] + Pi[c] + CPD0-2554[c] ';
model6 = addReaction(model6,'S-ADENMETSYN-RXN_2',form55,[],0,-1000,1000);


form66='CPD-9965[c] + NADP[c]   <=>   PROTON[c] + CPD-14280[c] + NADPH[c]';
model6 = addReaction(model6,'RXN-13306',form66,[],0,-1000,1000);
deads=model6.mets(detectDeadEnds(model6));
model=model6;
save('Sorghum_base_Mar25.mat','model')
%% extra fixes 5th March
met=model.mets(find(contains(model.mets,'D-glucopyranose-6-phosphate')));

[rxns,rxnform]=findRxnsFromMets(model,met);
rxnform=strrep(rxnform,met,'ALPHA-GLC-6-P[c]');
for n=1:length(rxns)
    if contains(rxnform{n},'<=>')
model = addReaction(model,rxns{n},rxnform{n},[],0,-1000,1000);
    else
        model = addReaction(model,rxns{n},rxnform{n},[],0,0,1000);

    end
end

form=printRxnFormula(model,'MALIC-NADP-RXN');
form=strrep(form,'[c]','[s]')
        model = addReaction(model,'MALIC-NADP-RXN',form,[],0,0,1000);
form=printRxnFormula(model,'MALATE-DEHYDROGENASE-NADP+-RXN');
form=strrep(form,'[c]','[s]')
        model = addReaction(model,'MALATE-DEHYDROGENASE-NADP+-RXN',form,[],0,-1000,1000);

form=printRxnFormula(model,'PYRUVATEORTHOPHOSPHATE-DIKINASE-RXN');
form=strrep(form,'[c]','[s]')
        model = addReaction(model,'PYRUVATEORTHOPHOSPHATE-DIKINASE-RXN',form,[],0,-1000,1000);
form='WATER[s] + CARBON-DIOXIDE[s] + D-RIBULOSE-15-P2[s] -> 2 PROTON[s] + 2 G3P[s]'
        model = addReaction(model,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN',form,[],0,0,1000);

form='CPD-19339[s] + GAP[s]   <=>   CPD-15318[s] + XYLULOSE-5-PHOSPHATE[s] '
        model = addReaction(model,'1TRANSKETO-RXN_2',form,[],0,-1000,1000);


form='D-SEDOHEPTULOSE-7-P[s] + GAP[s]   <=>   RIBOSE-5P[s]  + XYLULOSE-5-PHOSPHATE[s] '
        model = addReaction(model,'1TRANSKETO-RXN',form,[],0,-1000,1000);
pos=find(contains(model.rxns,'1TRANSKETO-RXN'))
model.rules(pos)=model.rules(pos(1))

form='CPD-19339[s] + GAP[s]   <=>   ERYTHROSE-4P[s]  + FRUCTOSE-6P[s]'; 
        model = addReaction(model,'TRANSALDOL-RXN_2',form,[],0,-1000,1000);
pos=find(contains(model.rxns,'TRANSALDOL-RXN'))
model.rules(pos)=model.rules(pos(1))



model=buildRxnGeneMat(model)
save('Sorghum_base_Mar25.mat','model')


% %% 
% model7=model6;
% no_backup={};added={};left={};extra={};
%  for n=1:length(deads)
%  [rxnListP, rxnFormulaListP] = findRxnsFromMets(model6, deads(n),'ProducersOnly',1);    
%   [rxnListC, rxnFormulaListC] = findRxnsFromMets(model6, deads(n),'ConsumersOnly',1); 
%        pos=find(contains(rxnFormulaListP,'<=>'));
%        metC=strrep(deads{n},'[c]','_c0');
% if isempty(rxnListP) & ~isempty(rxnListC)
%     % met can be consumed only
%          [rxnList_cheng, rxnFormulaList_cheng] = findRxnsFromMets(vModel, metC,'ProducersOnly',1);    
%     if isempty(rxnList_cheng)
%         no_backup=[no_backup,deads{n}];
%     else
%         % is the rxn(s) already in model?
%         for m=1:length(rxnList_cheng)
%           check=find(strcmp(model6.rxns,rxnList_cheng{m}));
% 
%         if isempty(check)
% 
%             form=strrep(rxnFormulaList_cheng{m},'_c0','[c]')
%             if contains(rxnFormulaList_cheng{m},'<=>')
%             model7=addReaction(model7,rxnList_cheng{m},form,[],0,-1000,1000);
%             deads7=model7.mets(detectDeadEnds(model7));
%             extra=setdiff(deads7,deads);
%             else
%              model7=addReaction(model7,rxnList_cheng{m},form,[],0,0,1000);
%                          deads7=model7.mets(detectDeadEnds(model7));
% 
%              extra=setdiff(deads7,deads);
% 
%             end
%             added=[added,rxnList_cheng{m}]
%         else
%             % do nothing
%         end
%         end
% 
%     end
% elseif ~isempty(rxnListP) && isempty(rxnListC)
%     % met can be produced only, look for consumers
%      [rxnList_cheng, rxnFormulaList_cheng] = findRxnsFromMets(vModel, metC,'ConsumersOnly',1);    
%     if isempty(rxnList_cheng)
%         no_backup=[no_backup,deads{n}];
%     else
%         % is the rxn(s) already in model?
%         for j=1:length(rxnList_cheng)
%           check=find(strcmp(model6.rxns,rxnList_cheng{j}));
% 
%         if isempty(check)
% 
%             form=strrep(rxnFormulaList_cheng{j},'_c0','[c]')
%             if contains(rxnFormulaList_cheng{j},'<=>')
%             model7=addReaction(model7,rxnList_cheng{j},form,[],0,-1000,1000);
%             deads7=model7.mets(detectDeadEnds(model7));
%             extra=setdiff(deads7,deads);
%             added=[added,rxnList_cheng{j}]
%             else
%              model7=addReaction(model7,rxnList_cheng{j},form,[],0,0,1000);
%                 deads7=model7.mets(detectDeadEnds(model7));
%             extra=setdiff(deads7,deads);
%             end
%             added=[added,rxnList_cheng{j}]
%         else
%             % do nothing
%         end
%         end
% 
%     end
% 
% elseif length(intersect(rxnListP,rxnListC))==1
%     % only 1 reaction thats reversible
%     left=[left,deads{n}];
% else
% end
%  end
% %% sorting the extras
% for n=1:length(added)
% 
% core4=removeRxns(model7,added{n})
% coredead= model7.mets(detectDeadEnds(model7));
% corenew= core4.mets(detectDeadEnds(core4))
% extra=setdiff(corenew,coredead)
% met=strrep(extra,'[c]','_c0');
% sol=optimizeCbModel(core4)
% added4={};
% while ~isempty(extra)
%          [rxnList, rxnFormulaList] = findRxnsFromMets(vModel, met)
%          rxnFormulaList=strrep(rxnFormulaList,'_c0','[c]');
% for n=1:length(rxnList)
%      if contains(rxnFormulaList{n},'<=>')
%      core4 = addReaction(core4,rxnList{n},rxnFormulaList{n},[],0,-1000,1000);
%      else
%               core4 = addReaction(core4,rxnList{n},rxnFormulaList{n},[],0,0,1000);
%      end
% added4=[added4,rxnList{n}]
% end
% corenew1= core4.mets(detectDeadEnds(core4))
% extra=setdiff(corenew1,corenew)
% sol=optimizeCbModel(core4)
% end