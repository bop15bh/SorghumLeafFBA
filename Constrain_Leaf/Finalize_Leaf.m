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
load('300deadsMar21.mat')

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

data=readcell('Sorghumcyc8_EC_subsystems.txt');
model.metFormulas={};
model.metCharges={};
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
