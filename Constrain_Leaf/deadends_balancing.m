clear
close ALL
load('model_balanced_databases_subs.mat')
cheng=model;
load('Constrained_unblocked_Leaf_FINAL0525.mat')

%%
% More aggressive conversion approach
if isfield(model, 'subSystems')
    % Convert everything to string first, then to cell array of char
    if ~isempty(model.subSystems)
        % Handle different input types
        if iscell(model.subSystems)
            % Convert each element to string, then to char
            model.subSystems = cellfun(@(x) char(string(x)), model.subSystems, 'UniformOutput', false);
        else
            % Convert non-cell to cell array of char
            model.subSystems = {char(string(model.subSystems))};
        end
        
        % Ensure column vector (some tools expect this)
        model.subSystems = model.subSystems(:);
    else
        % If empty, make it an empty cell array
        model.subSystems = {};
    end
end
changeCobraSolver('glpk');
form='NADPH[sm] + PROTON[sm] + OXALACETIC_ACID[sm] -> NADP[sm] + MAL[sm] ';
model=addReaction(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]',form,[],0,0,1000);
model=removeRxns(model,{'FUMHYDR-RXN_1[M]', 'ATR_2-PG_[cb]_[cm]','MALIC-NAD-RXN[B]','MALIC-NAD-RXN[M]'})
deads=model.mets(detectDeadEnds(model))

[rxns,rxnforms]=findRxnsFromMets(model,deads)

model=removeRxns(model,rxns)
deads=model.mets(detectDeadEnds(model))
[rxns,rxnforms]=findRxnsFromMets(model,deads)

model=removeRxns(model,rxns)
%model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.1, 'l');

%% testing C4 cycle
% co2=find(contains(model.rxns,'EX_CARBON-DIOXIDE_EXTRACELLULAR'));
% CA1=find(contains(model.rxns,'CARBODEHYDRAT-RXN'))
% pepc=find(contains(model.rxns,'PEPCARBOX'))
% mal=find(contains(model.rxns,'ATR_MAL_[cb]_[cm]'))
% malo=find(contains(model.rxns,'MALATE-DEHYDROGENASE-NADP+-RXN[M]'));
% me=find(contains(model.rxns,'MALIC-NADP-RXN[B]'))
% rub=find(contains(model.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
% pyr=find(contains(model.rxns,'ATR_PYRUVATE_[cb]_[cm]'))
% nit=find(contains(model.rxns,'EX_NITRATE_EXTRACELLULAR'));
% 
% model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');
% model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -350.4454, 'l');
% og=optimizeCbModel(model);
% 
% 
% dro = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -53, 'l');
% dro = changeRxnBounds(dro,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u');


%%

pos=find(contains(model.mets,'MAL[sm]'))
model.metFormulas(pos)={'C4H4O5'};
model.metCharges(pos)=-2
pos1=find(contains(model.mets,'OXALACETIC_ACID[sm]'))
model.metFormulas(pos1)={'C4H2O5'};
model.metCharges(pos1)=-2
%%
results = verifyModel(model,'massBalance', true)
imbalanced=model.rxns(results.massBalance.imBalancedRxnBool)
ex=model.rxns(find(contains(model.rxns,'EX_')))
imbalanced=setdiff(imbalanced,ex)
atr=imbalanced(find(contains(imbalanced,'ATR_')))
imbalanced=setdiff(imbalanced,atr)
missing=model.mets(results.massBalance.missingFormulaeBool);
[rxns,rxnforms]=findRxnsFromMets(model,missing)
imbalanced=setdiff(imbalanced,rxns)
elementsList = {'H', 'C', 'O', 'P', 'S', 'N', 'Mg', 'X', 'Fe', 'Zn', 'Co', 'Ca', 'Y', 'I', 'Na', 'Cl', 'K', 'R', 'FULLR'};
oners=[];elements={};proton={};
for n=1:length(imbalanced)
    pos=find(strcmp(model.rxns,imbalanced{n}))
    ele=results.massBalance.imBalancedMass(pos)
 numElements = strfind(ele, ',');
 if length(numElements{:})==0
     oners=[oners,n]
     elements=[elements,ele]
     if contains(ele,'H')
         proton=[proton,imbalanced{n}]
     else
     end
 else
 end
end
fixed=[];not_fixed=[];
model1=model;
for n=1:length(proton)

pos=find(strcmp(model1.rxns,proton{n}))
    ele=results.massBalance.imBalancedMass(pos)
    % water or proton
% -ve is  extra on left     
mets=findMetsFromRxns(model1,proton(n));
if sum(contains(mets,'PROTON'))==1
prot=mets(contains(mets,'PROTON'));
met_pos=find(strcmp(model1.mets,prot));
stoich=model1.S(met_pos,pos);
stoich_value=full(stoich);
ele=str2double(erase(ele,' H'));
new_stoich=stoich-ele;
model1.S(met_pos,pos)=new_stoich;
results1 = verifyModel(model1,'massBalance', true);
imbalanced1=model1.rxns(results1.massBalance.imBalancedRxnBool);
if isempty(find(strcmp(imbalanced1,proton{n})))
     fixed=[fixed,n]
else
    not_fixed=[not_fixed,n]
    
end
elseif  sum(contains(mets,'PROTON'))==0
    splitter=splitString(mets{1},'[');
    new_met=['PROTON[' splitter{2}];
    met_pos=find(strcmp(model.mets,new_met));
    ele=str2double(erase(ele,' H'));
    new_stoich=ele*-1;
    model1.S(met_pos,pos)=new_stoich;
results1 = verifyModel(model1,'massBalance', true);
imbalanced1=model1.rxns(results1.massBalance.imBalancedRxnBool);
if isempty(find(strcmp(imbalanced1,proton{n})))
        fixed=[fixed,n]
else
    not_fixed=[not_fixed,n];

end
else
        not_fixed=[not_fixed,n];

end
end
results1 = verifyModel(model1,'massBalance', true);
imbalanced1=model1.rxns(results1.massBalance.imBalancedRxnBool);
ex=model.rxns(find(contains(model.rxns,'EX_')))
imbalanced1=setdiff(imbalanced1,ex)
atr=imbalanced1(find(contains(imbalanced1,'ATR_')))
imbalanced1=setdiff(imbalanced1,atr)
missing=model.mets(results1.massBalance.missingFormulaeBool);
[rxns,rxnforms]=findRxnsFromMets(model,missing)
imbalanced1=setdiff(imbalanced1,rxns)


oners=[];elements={};proton={};
for n=1:length(imbalanced1)
   
    pos=find(strcmp(model1.rxns,imbalanced{n}))
    ele=results1.massBalance.imBalancedMass(pos)
 numElements = strfind(ele, ',');
% if length(numElements{:})==0
     oners=[oners,n]
     elements=[elements,ele]
     if contains(ele,'H')
         proton=[proton,imbalanced{n}]
     else
     end
 %else
 %end
end
model2=model1;
results2 = verifyModel(model2,'massBalance', true);
imbalanced2=model2.rxns(results2.massBalance.imBalancedRxnBool);
ex=model2.rxns(find(contains(model2.rxns,'EX_')))
imbalanced2=setdiff(imbalanced2,ex)
atr=imbalanced2(find(contains(imbalanced2,'ATR_')))
imbalanced2=setdiff(imbalanced2,atr)
missing=model2.mets(results2.massBalance.missingFormulaeBool);
[rxns,rxnforms]=findRxnsFromMets(model2,missing)
imbalanced2=setdiff(imbalanced2,rxns)
pos=find(contains(model2.rxns,imbalanced2))
ele=results2.massBalance.imBalancedMass(pos)

form='5-METHYL-THF-GLU-N[cb] + NAD[cb]   <=>   METHYLENE-THF-GLU-N[cb] + NADH[cb] + PROTON[cb]';
model2 = addReaction(model2,'1.5.1.20-RXN[B]',form,[],0,-1000,1000);
form1='5-METHYL-THF-GLU-N[cm] + NAD[cm]   <=>   METHYLENE-THF-GLU-N[cm] + NADH[cm] + PROTON[cm]';
model2 = addReaction(model2,'1.5.1.20-RXN[M]',form1,[],0,-1000,1000);
form2='NAD[cm] + GLYCERATE[cm]   <=>   CPD-26279[cm] + NADH[cm]   + PROTON[cm]';
model2 = addReaction(model2,'TSA-REDUCT-RXN[M]',form2,[],0,-1000,1000);
form3='NAD[cb] + GLYCERATE[cb]   <=>   CPD-26279[cb] + NADH[cb]   + PROTON[cb]';
model2 = addReaction(model2,'TSA-REDUCT-RXN[B]',form3,[],0,-1000,1000);


%%
% 
% pos=find(contains(model2.mets,'Dodec-2-enoyl-ACPs[cb]'))
% model2.metFormulas(pos)={'C26H45N3O9PS'};
% pos1=find(contains(model2.mets,'Dodec-2-enoyl-ACPs[cm]'))
% model2.metFormulas(pos1)={'C26H45N3O9PS'};
% pos2=find(contains(model2.mets,'Acetoacetyl-ACPs[cm]'))
% model2.metCharges(pos2)=-2;
% 
% pos3=find(contains(model2.mets,'Acetoacetyl-ACPs[cb]'))
% model2.metCharges(pos3)=-2;
% pos4=find(contains(model2.mets,'Oleoyl-ACPs[cb]'))
% model2.metFormulas(pos4)={'C18H34OS'};
% pos5=find(contains(model2.mets,'Oleoyl-ACPs[cm]'))
% model2.metFormulas(pos5)={'C18H34OS'};
% pos6=find(contains(model2.mets,'Dodecanoyl-ACPs[cb]'))
% model2.metFormulas(pos6)={'C12H24OS'};
% pos7=find(contains(model2.mets,'Dodecanoyl-ACPs[cm]'))
% model2.metFormulas(pos7)={'C12H24OS'};
% pos8=find(contains(model2.mets,'CELLULOSE_Compound[cb]'))
% model2.metFormulas(pos8)={'C6H11O5'};
% model2.metCharges(pos8)=1;
% 
% pos9=find(contains(model2.mets,'CELLULOSE_Compound[cm]'))
% model2.metFormulas(pos9)={'C6H11O5'};
% model2.metCharges(pos9)=1;
% pos10=find(contains(model2.mets,'Amylose_Compound[sb]'))
% model2.metFormulas(pos10)={'C6H10O5'};
% pos11=find(contains(model2.mets,'Amylose_Compound[sm]'))
% model2.metFormulas(pos11)={'C6H10O5'};
% 
% form4='PLASTOQUINONE-9[sb] + 4 PROTON[sb] + 2 Reduced-ferredoxins[sb]   ->   CPD-12829[sb] + 2 Oxidized-ferredoxins[sb] + 2 PROTON[tb]' 
%   model2 = addReaction(model2,'PSI-CYCLIC-RXN[B]',form4,[],0,0,1000);
% form45='PLASTOQUINONE-9[sm] + 4 PROTON[sm] + 2 Reduced-ferredoxins[sm]   ->   CPD-12829[sm] + 2 Oxidized-ferredoxins[sm] + 2 PROTON[tm]' 
%   model2 = addReaction(model2,'PSI-CYCLIC-RXN[M]',form5,[],0,0,1000);
% form6='HYDROGEN-PEROXIDE[cb] + ASCORBATE[cb]  ->    + L-DEHYDRO-ASCORBATE[cb]  + 2 WATER[cb] ';
% model2 = addReaction(model2,'RXN-12440[B]',form6,[],0,0,1000);
% form7='HYDROGEN-PEROXIDE[cm] + ASCORBATE[cm]  ->    + L-DEHYDRO-ASCORBATE[cm]  + 2 WATER[cm] ';
% model2 = addReaction(model2,'RXN-12440[M]',form7,[],0,0,1000);
% pos12=find(contains(model2.mets,'Stearoyl-ACPs[sb]'))
% model2.metFormulas(pos12)={'C18H36OS'};
% pos13=find(contains(model2.mets,'Stearoyl-ACPs[sm]'))
% model2.metFormulas(pos13)={'C18H36OS'};
% pos14=find(contains(model2.mets,'Stearoyl-ACPs[cb]'))
% model2.metFormulas(pos14)={'C18H36OS'};
% pos15=find(contains(model2.mets,'Stearoyl-ACPs[cm]'))
% model2.metFormulas(pos15)={'C18H36OS'};
% 
% pos16=find(contains(model2.mets,'Phosphorylated-phosphoglucomutase[cb]'))
% model2.metFormulas(pos16)={'C3H4NO5P'};
% pos17=find(contains(model2.mets,'Phosphorylated-phosphoglucomutase[cm]'))
% model2.metFormulas(pos17)={'C3H4NO5P'};
% pos18=find(contains(model2.mets,'CPD-8124'))
% model2.metFormulas(pos18)={'C10H10N5O7P1S3Mo1'};
% pos19=find(contains(model2.mets,'Beta-3-hydroxybutyryl-ACPs'))
% model2.metFormulas(pos19)={'C10H10N5O7P1S3Mo1'};
% pos20=find(contains(model2.mets,'R-3-hydroxyhexanoyl-ACPs'))
% model2.metFormulas(pos20)={'C17H31N2O9PS'};
% 
% met_pos=find(contains(model2.mets,'PROTON[cb]'))
% rxn_pos=find(contains(model2.rxns,'RXN-9535[B]'))
% model2.S(met_pos,rxn_pos)=0
% met_pos2=find(contains(model2.mets,'PROTON[cm]'))
% rxn_pos2=find(contains(model2.rxns,'RXN-9535[M]'))
% model2.S(met_pos2,rxn_pos2)=0


results2 = verifyModel(model2,'massBalance', true);
imbalanced2=model2.rxns(results2.massBalance.imBalancedRxnBool);
ex=model2.rxns(find(contains(model2.rxns,'EX_')))
imbalanced2=setdiff(imbalanced2,ex)
atr=imbalanced2(find(contains(imbalanced2,'ATR_')))
imbalanced2=setdiff(imbalanced2,atr)
missing=model2.mets(results2.massBalance.missingFormulaeBool);
[rxns,rxnforms]=findRxnsFromMets(model2,missing)
imbalanced2=setdiff(imbalanced2,rxns)
pos=find(contains(model2.rxns,imbalanced2))
ele=results2.massBalance.imBalancedMass(pos)

%DIHYDROFOLATE-GLU-N[cb]

% pos=find(strcmp(model1.rxns,'1.5.1.20-RXN[B]'));
% met_pos=find(strcmp(model.mets,'PROTON[cb]'));
%   model2.S(met_pos,pos)=-1;
% results2 = verifyModel(model2,'massBalance', true);
% imbalanced2=model2.rxns(results2.massBalance.imBalancedRxnBool);
% ele=results1.massBalance.imBalancedMass(results2.massBalance.imBalancedRxnBool)
% find(strcmp(imbalanced2,proton{n}))


%%

list=erase(imbalanced2,'[B]')
list=erase(list,'[M]')
results_cheng = verifyModel(cheng,'massBalance', true)
imbalanced_cheng=cheng.rxns(results_cheng.massBalance.imBalancedRxnBool)
match=cheng.rxns(find(contains(cheng.rxns,list(7))))
intersect(imbalanced_cheng,match)

%%
pos11=find(contains(model2.mets,'Oleoyl-ACPs'))
model2.metFormulas(pos11)={'C18H33O2'};

pos12=find(contains(model2.mets,'ACP['))
mal=find(contains(model2.mets,'MALONYL-ACP'));
pos12=setdiff(pos12,mal)
model2.metFormulas(pos12)={'HO'};

pos13=find(contains(model2.mets,'R-3-Hydroxypalmitoyl-ACPs'))
model2.metFormulas(pos13)={'C16H31O3'};
pos14=find(contains(model2.mets,'2-Hexadecenoyl-ACPs'))
model2.metFormulas(pos14)={'C16H29O2'};

pos9=find(contains(model2.mets,'CELLULOSE_Compound'))
model2.metFormulas(pos9)={'C6H10O5'};
met_pos=find(contains(model2.mets,'PROTON[cb]'))
met_pos2=find(contains(model2.mets,'PROTON[cm]'))

rxn_pos1=find(contains(model2.rxns,'CELLULOSE-SYNTHASE-UDP-FORMING-RXN[B]'))
rxn_pos2=find(contains(model2.rxns,'CELLULOSE-SYNTHASE-UDP-FORMING-RXN[M]'))
model2.S(met_pos,rxn_pos1)=1;
model2.S(met_pos2,rxn_pos2)=1
pos10=find(contains(model2.mets,'Amylose_Compound'))
model2.metFormulas(pos10)={'C6H10O5'};
pos15=find(contains(model2.mets,'MALONYL-ACP'))
model2.metFormulas(pos15)={'C3H2O4'};
form6='HYDROGEN-PEROXIDE[cb] + ASCORBATE[cb]  ->    + L-DEHYDRO-ASCORBATE[cb]  + 2 WATER[cb] ';
model2 = addReaction(model2,'RXN-12440[B]',form6,[],0,0,1000);
form7='HYDROGEN-PEROXIDE[cm] + ASCORBATE[cm]  ->    + L-DEHYDRO-ASCORBATE[cm]  + 2 WATER[cm] ';
model2 = addReaction(model2,'RXN-12440[M]',form7,[],0,0,1000);

pos16=find(contains(model2.mets,'Phosphorylated-phosphoglucomutase'))
model2.metFormulas(pos16)={'O4P'};
pos17=find(contains(model2.mets,'Phosphoglucomutase'))
model2.metFormulas(pos17)={'HO'};
pos18=find(contains(model2.mets,'Stearoyl-ACPs'))
model2.metFormulas(pos18)={'C18H35O2'};
pos19=find(contains(model2.mets,'CPD-8124'))
 model2.metFormulas(pos19)={'C10H10N5O7P1S3Mo1'};
pos20=find(contains(model2.mets,'3-oxo-hexanoyl-ACPs'))
 model2.metFormulas(pos20)={'C6H9O3'};

 pos21=find(contains(model2.mets,'Butanoyl-ACPs'))
 model2.metFormulas(pos21)={'C4H7O2'};
 pos22=find(contains(model2.mets,'R-3-hydroxyhexanoyl-ACPs'))
 model2.metFormulas(pos22)={'C6H11O3'};

  pos23=find(contains(model2.mets,'Hex-2-enoyl-ACPs'))
 model2.metFormulas(pos23)={'C6H9O2'};
   model2.metCharges(pos23)=0;

 pos24=find(contains(model2.mets,'Hexanoyl-ACPs'))
 model2.metFormulas(pos24)={'C6H11O2'};
 model2.metCharges(pos24)=0;

 pos25=find(contains(model2.mets,'3-oxo-octanoyl-ACPs'))
 model2.metFormulas(pos25)={'C8H13O3'};
  model2.metCharges(pos25)=0;

  

  pos26=find(contains(model2.mets,'3-Hydroxy-octanoyl-ACPs'))
 model2.metFormulas(pos26)={'C8H15O3'};  
 model2.metCharges(pos26)=0;

  pos27=find(contains(model2.mets,'3-oxo-dodecanoyl-ACPs'))
 model2.metFormulas(pos27)={'C12H21O3'}; 
   pos28=find(contains(model2.mets,'Decanoyl-ACPs'))
 model2.metFormulas(pos28)={'C10H19O2'}; 
  pos29=find(contains(model2.mets,'R-3-hydroxydodecanoyl-ACPs'))
 model2.metFormulas(pos29)={'C12H23O3'}; 
pos30=find(contains(model2.mets,'Dodec-2-enoyl-ACPs'))
 model2.metFormulas(pos30)={'C12H21O2'}; 
  pos31=find(contains(model2.mets,'Dodecanoyl-ACPs'))
 model2.metFormulas(pos31)={'C12H23O2'}; 
pos32=find(contains(model2.mets,'3-oxo-myristoyl-ACPs'))
 model2.metFormulas(pos32)={'C14H25O3'}; 
  pos33=find(contains(model2.mets,'R-3-hydroxymyristoyl-ACPs'))
 model2.metFormulas(pos33)={'C14H27O3'}; 
 pos34=find(contains(model2.mets,'Tetradec-2-enoyl-ACPs'))
 model2.metFormulas(pos34)={'C14H25O2'}; 
 pos35=find(contains(model2.mets,'3-oxo-palmitoyl-ACPs'))
 model2.metFormulas(pos35)={'C16H29O3'}; 
  pos36=find(contains(model2.mets,'Myristoyl-ACPs'))
 model2.metFormulas(pos36)={'C14H27O2'}; 

  pos37=find(contains(model2.mets,'3-oxo-stearoyl-ACPs'))
 model2.metFormulas(pos37)={'C18H33O3'}; 
  pos38=find(contains(model2.mets,'Palmitoyl-ACPs'))
 model2.metFormulas(pos38)={'C16H31O2'}; 
  pos39=find(contains(model2.mets,'R-3-hydroxystearoyl-ACPs'))
 model2.metFormulas(pos39)={'C18H35O3'}; 
  pos40=find(contains(model2.mets,'Octadec-2-enoyl-ACPs'))
 model2.metFormulas(pos40)={'C18H33O2'}; 
 pos41=find(contains(model2.mets,'Trans-D2-decenoyl-ACPs'))
 model2.metFormulas(pos41)={'C10H17O2'};  
 pos42=find(contains(model2.mets,'Crotonyl-ACPs'))
 model2.metFormulas(pos42)={'C4H5O2'};  
 pos43=find(contains(model2.mets,'Octanoyl-ACPs'))
 model2.metFormulas(pos43)={'C8H15O2'};  
   model2.metCharges(pos43)=0;

 pos44=find(contains(model2.mets,'2-Octenoyl-ACPs'))
 model2.metFormulas(pos44)={'C8H13O2'};  
  model2.metCharges(pos44)=0;

 pos45=find(contains(model2.mets,'Beta-hydroxydecanoyl-ACPs'))
 model2.metFormulas(pos45)={'C10H19O3'};  
 pos46=find(contains(model2.mets,'Acetoacetyl-ACPs'))
 model2.metFormulas(pos46)={'C4H5O3'};  
met_pos=find(contains(model2.mets,'PROTON[cb]'))
met_pos2=find(contains(model2.mets,'PROTON[cm]'))

rxn_pos1=find(contains(model2.rxns,'2.3.1.180-RXN[B]'))
rxn_pos2=find(contains(model2.rxns,'2.3.1.180-RXN[M]'))
model2.S(met_pos,rxn_pos1)=-1;
model2.S(met_pos2,rxn_pos2)=-1

pos47=find(contains(model2.mets,'Beta-3-hydroxybutyryl-ACPs'))
 model2.metFormulas(pos47)={'C4H7O3'};  
 

 model2=removeRxns(model2,{'RXN-1161[B]','RXN-1161[M]'})
 form8=' NITRITE[sm] + 6 Reduced-ferredoxins[sm] 8 PROTON[sm]  ->    AMMONIUM[sm]  + 6 Oxidized-ferredoxins[sm] + 2 WATER[sm]';
 model2 = addReaction(model2,'FERREDOXIN--NITRITE-REDUCTASE-RXN[M]',form8,[],0,0,1000);


  pos48=find(contains(model2.mets,'3-oxo-decanoyl-ACPs'))
 model2.metFormulas(pos48)={'C10H17O3'};  
   model2.metCharges(pos48)=0;

met_pos=find(contains(model2.mets,'PROTON[cb]'))
met_pos2=find(contains(model2.mets,'PROTON[cm]'))

rxn_pos1=find(contains(model2.rxns,'RXN-9549[B]'))
rxn_pos2=find(contains(model2.rxns,'RXN-9549[M]'))
model2.S(met_pos,rxn_pos1)=1;
model2.S(met_pos2,rxn_pos2)=1

model=model2;
save('Leaf_balanced_FINAL0625.mat','model')
