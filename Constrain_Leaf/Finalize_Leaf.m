%% finalize sorghum model
% add subsystems, formulas etc

clear
close all
% FVA pipeline
%addpath /mnt/home/holla293/Documents/cobratoolbox
%    initCobraToolbox;

%changeCobraSolver('glpk')
%load('SC_constrained_unblocked.mat');
%load('constrained_leaf_noredundancies.mat')
%load('300deadsMar21.mat')
load('Constrained_unblocked_MAY25_Final.mat')
model.metFormulas={};

%% adding subsystems and EC numbers
rxns=erase(model.rxns,'[B]');
rxns=erase(rxns,'[M]');
rxns=unique(rxns)
rxns=erase(rxns,'_10');
rxns=erase(rxns,'_21');
rxns=erase(rxns,'_22');
rxns=erase(rxns,'_23');
rxns=erase(rxns,'_24');
rxns=erase(rxns,'_25');
rxns=erase(rxns,'_26');
rxns=erase(rxns,'_2');
rxns=erase(rxns,'_1');
rxns=erase(rxns,'_3');
rxns=erase(rxns,'_4');
rxns=erase(rxns,'_5');
rxns=erase(rxns,'_6');
rxns=erase(rxns,'_7');
rxns=erase(rxns,'_8');
rxns=erase(rxns,'_9');

rxns=unique(rxns)


ex=rxns(find(contains(rxns,'EX_')));
rxns=setdiff(rxns,ex)
atr=rxns(find(contains(rxns,'ATR_')));
rxns=setdiff(rxns,atr)

%data=readcell('Sorghumcyc8_EC_subsystems.txt');
charge=readtable('sorghum_charge (1).txt')
comp={};charge1=[];
for n=1:height(charge(:,1))
        entry=table2cell(charge(n,2));
        name=table2cell(charge(n,1));
    charge1=[charge1,entry{:}];
    comp=[comp,name{:}];
end
nans1=comp(isnan(charge1))
pos=find(strcmp(comp,'XYLOSE'))
charge1(pos)=0;
extraforms(1,1)={'XYLOSE'};
extraforms(1,2)={'C5H10O5'};
pos1=find(strcmp(comp,'S_CPD-14916'))
charge1(pos1)=-4;
extraforms(2,1)={'S_CPD-14916'};
extraforms(2,2)={'C29H46N7O18P3S'};




data=readcell('Sorghumcyc8_subsystems_names.txt');
model.metCharges=[];
model.EC={};
model.subSystems={};
for n=1:length(model.rxns)

    rxns1=erase(model.rxns{n},'[B]');
rxns1=erase(rxns1,'[M]');
rxns1=erase(rxns1,'_10');
rxns1=erase(rxns1,'_21');
rxns1=erase(rxns1,'_22');
rxns1=erase(rxns1,'_23');
rxns1=erase(rxns1,'_24');
rxns1=erase(rxns1,'_25');
rxns1=erase(rxns1,'_26');
rxns1=erase(rxns1,'_2');
rxns1=erase(rxns1,'_1');
rxns1=erase(rxns1,'_3');
rxns1=erase(rxns1,'_4');
rxns1=erase(rxns1,'_5');
rxns1=erase(rxns1,'_6');
rxns1=erase(rxns1,'_7');
rxns1=erase(rxns1,'_8');
rxns1=erase(rxns1,'_9');
pos=find(strcmp(data(:,4),rxns1));
if ~isempty(pos)
if ~ismissing(data{pos,2})
model.subSystems{n,1}=data{pos,2};
else
        model.subSystems{n,1}=' - ';

end
if ~ismissing(data{pos,3})
model.EC{n,1}=data{pos,3};
else
    model.EC{n,1}=' - ';

end
if ~ismissing(data{pos,5})
model.rxnNames{n,1}=data{pos,5};
else
end
else
    model.subSystems{n,1}=' - ';
model.EC{n,1}=' - ';
end
end



%% adding formulas to metabolites

formulas=readcell('met_formulas.txt')
for n=1:length(model.mets)
    met=erase(model.mets{n},'[cm]');
met=erase(met,'[cb]');
met=erase(met,'[sb]');
met=erase(met,'[sm]');
met=erase(met,'[tm]');
met=erase(met,'[tb]');
met=erase(met,'[mm]');
met=erase(met,'[mb]');
met=erase(met,'[e]');
pos=find(strcmp(formulas(:,1),met))
if ~isempty(pos)
if ~ismissing(formulas{pos,2})
model.metFormulas{n,1}=formulas{pos,2};
else
        model.metFormulas{n,1}=' - ';

end 
else
            model.metFormulas{n,1}=' - ';

end
end

missing_forms=model.mets(find(strcmp(model.metFormulas,' - ')));
missing_forms=erase(missing_forms,'[cb]');
missing_forms=erase(missing_forms,'[cm]');
missing_forms=erase(missing_forms,'[mb]');
missing_forms=erase(missing_forms,'[mm]');
missing_forms=erase(missing_forms,'[sm]');
missing_forms=erase(missing_forms,'[sb]');
missing_forms=erase(missing_forms,'[tm]');
missing_forms=erase(missing_forms,'[tb]');
missing_forms=erase(missing_forms,'[e]');
missing_forms=unique(missing_forms)

new_forms=readcell('formulas.txt')


for n=1:length(model.mets)
    met=erase(model.mets{n},'[cm]');
met=erase(met,'[cb]');
met=erase(met,'[sb]');
met=erase(met,'[sm]');
met=erase(met,'[tm]');
met=erase(met,'[tb]');
met=erase(met,'[mm]');
met=erase(met,'[mb]');
met=erase(met,'[e]');
pos=find(strcmp(new_forms(:,1),met))
if ~isempty(pos)
if ~ismissing(new_forms{pos,2})
model.metFormulas{n,1}=new_forms{pos,2};
else

end 
else

end
end

%% extra forms 
extraforms(3,1)={'Red-Thioredoxin'};
extraforms(3,2)={'C10H15N4O4S2'};
extraforms(9,1)={'Ox-Thioredoxin'};
extraforms(9,2)={'C10H13N4O4S2'};
extraforms(4,1)={'Oxidized-ferredoxins'};
extraforms(4,2)={'Fe2S2'};
extraforms(5,1)={'Reduced-ferredoxins'};
extraforms(5,2)={'Fe2S2'};
extraforms(6,1)={'Plastocyanin-Reduced'};
extraforms(6,2)={'Cu'};
extraforms(7,1)={'Phosphorylated-phosphoglucomutase'};
extraforms(7,2)={'C4H7N2O6P'};
extraforms(8,1)={'Phosphoglucomutase'};
extraforms(8,2)={'C3H5NO2'};
extraforms(10,1)={'S_CPD-14275'};
extraforms(10,2)={'C41H70N7O18P3S'};
extraforms(10,3)={-4};
extraforms(11,1)={'R-3-hydroxystearoyl-ACPs'};
extraforms(11,2)={'C32H59N3O10PS'};
extraforms(12,1)={'R-3-hydroxymyristoyl-ACPs'};
extraforms(12,2)={'C28H51N3O10PS'};
extraforms(13,1)={'R-3-hydroxyhexanoyl-ACPs'};
extraforms(13,2)={'C20H35N3O10PS'};
extraforms(14,1)={'R-3-hydroxydodecanoyl-ACPs'};
extraforms(14,2)={'C26H47N3O10PS'};
extraforms(15,1)={'R-3-Hydroxypalmitoyl-ACPs'};
extraforms(15,2)={'C30H55N3O10PS'};
extraforms(16,1)={'Palmitoyl-ACPs'};
extraforms(16,2)={'C16H31OS'};
extraforms(17,1)={'Oxidized-Plastocyanins'};
extraforms(17,2)={'Cu'};
extraforms(18,1)={'2-Hexadecenoyl-ACPs'};
extraforms(18,2)={'C16H29OS'};
extraforms(19,1)={'2-Octenoyl-ACPs'};
extraforms(19,2)={'C19H33N2O8PS'};
extraforms(20,1)={'3-Hydroxy-octanoyl-ACPs'};
extraforms(20,2)={'C19H35N2O9PS'};
extraforms(21,1)={'3-oxo-decanoyl-ACPs'};
extraforms(21,2)={'C10H17O2S'};
extraforms(22,1)={'3-oxo-dodecanoyl-ACPs'};
extraforms(22,2)={'C23H41N2O9PS'};
extraforms(23,1)={'3-oxo-hexanoyl-ACPs'};
extraforms(23,2)={'C17H29N2O9PS'};

extraforms(24,1)={'3-oxo-myristoyl-ACPs'};
extraforms(24,2)={'C14H25O2S'};
extraforms(25,1)={'3-oxo-octanoyl-ACPs'};
extraforms(25,2)={'C19H33N2O9PS'};
extraforms(25,3)={-2};
extraforms(26,1)={'3-oxo-palmitoyl-ACPs'};
extraforms(26,2)={'C27H49N2O9PS'};
extraforms(27,1)={'3-oxo-stearoyl-ACPs'};
extraforms(27,2)={'C32H57N3O10PS'};
extraforms(28,1)={'5-10-METHENYL-THF-GLU-N'};
extraforms(28,2)={'C20H20N7O5'};
extraforms(29,1)={'5-METHYL-THF-GLU-N'};
extraforms(29,2)={'C20H23N7O5'};

extraforms(30,1)={'ACP'};
extraforms(30,2)={'H2S'};
extraforms(31,1)={'Acetoacetyl-ACPs'};
extraforms(31,2)={'C4H5O2S'};
extraforms(32,1)={'Alpha-linolenoyl-groups'};
extraforms(32,2)={'C18H29O'};
extraforms(33,1)={'Amylose_Compound'};
extraforms(33,2)={'C16H30O10'};
extraforms(33,3)={0};
extraforms(34,1)={'Beta-3-hydroxybutyryl-ACPs'};
extraforms(34,2)={'C18H31N3O10PS'};
extraforms(35,1)={'Beta-hydroxydecanoyl-ACPs'};
extraforms(35,2)={'C10H19O2S'};
extraforms(36,1)={'Butanoyl-ACPs'};
extraforms(36,2)={'C4H7OS'};
extraforms(37,1)={'CELLULOSE_Compound'};
extraforms(37,2)={'C24H40O19'};
extraforms(37,3)={0};

extraforms(38,1)={'Charged-GLT-tRNAs'};
extraforms(38,2)={'C33H41N12O24P3'};
extraforms(39,1)={'Crotonyl-ACPs'};
extraforms(39,2)={'C18H29N3O9PS'};
extraforms(40,1)={'Cytochromes-C-Oxidized'};
extraforms(40,2)={'C34H32FeN4O4S2'};
extraforms(41,1)={'Cytochromes-C-Reduced'};
extraforms(41,2)={'C34H32FeN4O4S2'};
extraforms(42,1)={'DIHYDROFOLATE-GLU-N'};
extraforms(42,2)={'C19H19N7O5'};
extraforms(43,1)={'Decanoyl-ACPs'};
extraforms(43,2)={'C10H19OS'};
extraforms(44,1)={'Dodec-2-enoyl-ACPs'};
extraforms(44,2)={'C12H21OS'};
extraforms(45,1)={'Dodecanoyl-ACPs'};
extraforms(45,2)={'C26H47N3O9PS'};
extraforms(46,1)={'ETF-Oxidized'};
extraforms(46,2)={'C13H10N4O2'};
extraforms(47,1)={'ETF-Reduced'};
extraforms(47,2)={'C13H13N4O2'};
extraforms(48,1)={'FERRICYTOCHROME-B5'};
extraforms(48,2)={'C34H30FeN4O4'};
extraforms(49,1)={'FERROCYTOCHROME-B5'};
extraforms(49,2)={'C34H30FeN4O4'};
extraforms(50,1)={'FORMYL-THF-GLU-N'};
extraforms(50,2)={'C20H21N7O6'};
extraforms(51,1)={'GLT-tRNAs'};
extraforms(51,2)={'C28H34N11O21P3'};
extraforms(52,1)={'Hex-2-enoyl-ACPs'};
extraforms(52,2)={'C17H29N2O8PS'};
extraforms(53,1)={'Hexanoyl-ACPs'};
extraforms(53,2)={'C6H11OS'};
extraforms(54,1)={'Linoleoyl-groups'};
extraforms(54,2)={'C18H31O'};
extraforms(55,1)={'MALONYL-ACP'};
extraforms(55,2)={'C3H3O3S'};
extraforms(56,1)={'METHYLENE-THF-GLU-N'};
extraforms(56,2)={'C20H21N7O5'};
extraforms(57,1)={'Myristoyl-ACPs'};
extraforms(57,2)={'C14H27OS'};
extraforms(58,1)={'Octadec-2-enoyl-ACPs'};
extraforms(58,2)={'C32H57N3O9PS'};
extraforms(59,1)={'Octanoyl-ACPs'};
extraforms(59,2)={'C8H15OS'};
extraforms(60,1)={'Oleoyl-ACPs'};
extraforms(60,2)={'C32H57N3O9PS'};
extraforms(61,1)={'Oleoyl-lipid'};
extraforms(61,2)={'C18H33O'};
extraforms(62,1)={'Ox-NADPH-Hemoprotein-Reductases_Compound'};
extraforms(62,2)={'C12H8N4O2'};
extraforms(62,3)={-1};
extraforms(63,1)={'Red-NADPH-Hemoprotein-Reductases_Compound'};
extraforms(63,2)={'C12H11N4O2'};
extraforms(63,3)={0};

extraforms(64,1)={'Stearoyl-ACPs'};
extraforms(64,2)={'C32H59N3O9PS'};
extraforms(65,1)={'THF-GLU-N'};
extraforms(65,2)={'C19H21N7O5'};
extraforms(66,1)={'Tetradec-2-enoyl-ACPs'};
extraforms(66,2)={'C14H25OS'};
extraforms(67,1)={'Trans-D2-decenoyl-ACPs'};
extraforms(67,2)={'C24H41N3O9PS'};
extraforms(68,1)={'biomass'};
extraforms(68,3)={0};
extraforms(69,1)={'Light_Compound'};
extraforms(69,3)={0};
comp=comp';
charge1=charge1';
for n=1:length(model.mets)
    met=erase(model.mets{n},'[cm]');
met=erase(met,'[cb]');
met=erase(met,'[sb]');
met=erase(met,'[sm]');
met=erase(met,'[tm]');
met=erase(met,'[tb]');
met=erase(met,'[mm]');
met=erase(met,'[mb]');
met=erase(met,'[b]');
met=erase(met,'[m]');

met=erase(met,'[e]');
pos=find(strcmp(extraforms(:,1),met))
pos1=find(strcmp(comp(:,1),met))
if ~isempty(pos) && ~isempty(pos1)
if ~ismissing(extraforms{pos,2}) 
model.metFormulas{n,1}=extraforms{pos,2};
model.metCharges{n,1}=extraforms{pos,3};
elseif ~ismissing(extraforms{pos,3}) 
model.metCharges{n,1}=extraforms{pos,3};

else

end 
else

end
if ~isempty(pos1)
if ~isnan(charge1(pos1))
model.metCharges{n,1}=charge1(pos1);
else
end 

elseif ~isempty(pos)
model.metCharges{n,1}=extraforms{pos,3};
else
end

end

missing_forms=model.mets(find(strcmp(model.metFormulas,' - ')));
missing_forms=erase(missing_forms,'[cb]');
missing_forms=erase(missing_forms,'[cm]');
missing_forms=erase(missing_forms,'[mb]');
missing_forms=erase(missing_forms,'[mm]');
missing_forms=erase(missing_forms,'[sm]');
missing_forms=erase(missing_forms,'[sb]');
missing_forms=erase(missing_forms,'[tm]');
missing_forms=erase(missing_forms,'[tb]');
missing_forms=erase(missing_forms,'[e]');
missing_forms=unique(missing_forms)
missing_charge=[];
for n=1:length(model.mets)
if isempty(model.metCharges{n,1})
 missing_charge=[missing_charge,n]
else
end
end
model.mets(missing_charge)
model.metFormulas(strcmp(model.metFormulas, ' - ')) = {''};
model.metCharges(missing_charge(1))={0};
model.metCharges(missing_charge(2))={0};
model.metCharges(missing_charge(3))={1};
model.metCharges(missing_charge(4))={0};

charges=cell2mat(model.metCharges)
model.metCharges = charges;
model.metCharges(missing_charge(1))
results = verifyModel(model,'massBalance', true)
rxns=model.rxns(results.massBalance.imBalancedRxnBool)
save('35deads_withsubsystems.mat','model')

