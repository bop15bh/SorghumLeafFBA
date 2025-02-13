compnames=readcell("common_name_to_metacycID.xlsx");
forms1=readcell("Sorghum Reaction Equations.xlsx");
forms=readcell("sorghumbicolorcyc-8_0_0-rxns .txt")
 % pool = parpool
forms(1265,1)={'--MENTHOL-DEHYDROGENASE-RXN'}
forms(:,2)=strrep(forms(:,2),' x ',' ')
rev=[];lb=[];ub=[];
for n=1:length(forms(:,2))
if contains(forms(n,2),'<-->') || contains(forms(n,2),'=')
    rev=[rev,1];
    lb=[lb,-1000];
    ub=[ub,1000];
else
    rev=[rev,0];
    lb=[lb,0];
    ub=[ub,1000];
end
end
for n=1:length(forms(:,2))
   if  contains(forms(n,2),'<-')
       formy=strrep(forms(n,2),'<-','<')
       mets=splitString(formy,'<');
       mets2=mets{:};
       newform=[mets2{2} ' -> ' mets2{1}] 
       forms{n,2}=newform;
   else
   end
end
Formulas=forms(:,2)
 Formulas=strrep(Formulas,' = ',' <=> ')
 positi=find(contains(Formulas,'Light'))
 Formulas(positi)=strrep( Formulas(positi),'Light','Light_Compound')
positi1=find(contains(Formulas,'Red-NADPH-Hemoprotein-Reductases'))
 Formulas(positi1)=strrep( Formulas(positi1),'Red-NADPH-Hemoprotein-Reductases','Red-NADPH-Hemoprotein-Reductases_Compound')
 positi2=find(contains(Formulas,'Ox-NADPH-Hemoprotein-Reductases'))
 Formulas(positi2)=strrep( Formulas(positi2),'Ox-NADPH-Hemoprotein-Reductases','Ox-NADPH-Hemoprotein-Reductases_Compound')

rxnAbrList = forms(:,1);
revFlagList = rev;
lowerBoundList = lb;
upperBoundList= ub;
model = createModel(rxnAbrList,rxnAbrList,Formulas,revFlagList, ...
lowerBoundList, upperBoundList);
model=removeRxns(model,{'Reaction'})


%% testing which reactions formulas are completely made up of correct metabolites
%load('Correct_metabolite_names.mat')
load('Correct_metabolite_names_sorghum25.mat')
Met_list_V2=met_list
%Metabolite_list
complete={};
for n=1:length(model.rxns)
    mets=findMetsFromRxns(model,model.rxns(n))
    mets=erase(mets,'[c]');
    noon=[];
    for j=1:length(mets)
    pos=strmatch(mets(j),Met_list_V2,'exact');
    if ~isempty(pos)
    noon=[noon,pos];
    else
    end
    end
    if length(noon)==length(mets)
        complete=[complete,model.rxns(n)]
    else
    end
end

left=setdiff(model.rxns,complete)

load('base10rxns.mat') % common rxns between cheng and sorghumcyc
blue=readcell('blumatch.xlsx')
baser=setdiff(base_rxns10,complete)
% lol is green rxns from the 137 so should be kept 
lol={'ABC-70-RXN','TRANS-RXN-178','BENZALDEHYDE-DEHYDROGENASE-NAD+-RXN','RXN-161','2.7.7.13-RXN','NITRATE-REDUCTASE-NADH-RXN','2.4.1.67-RXN','RXN66-3'};
base1=vertcat(baser,lol')
base1=vertcat(base1,blue)
left=setdiff(left,unique(base1))
% model2 is base reactions including blue confirmed and missed off base 10
% rxns
model2=removeRxns(model,left)
model=removeRxns(model,model2.rxns)

% ensuring rxns from multi aren't repeated in model2
tor=readcell('comps_to_remove.xlsx');
model2=removeRxns(model2,tor);

% adding extra rxns
cheng_add=readcell('cheng_add.xlsx');
comps=readcell('multi.xlsx')
load('tcell_essentNAMES.mat')
load('tcell_essentFORMS.mat')
% removing cell type and c<->c rxns
cc=~contains(keep,'_[cb]_[cm]')
keep1=keep(cc)
keep_form1=keep_form(cc)
keep1=erase(keep1,'[B]');
keep1=erase(keep1,'[M]')
pos=find(contains(keep1,'_[cm]_[cb]'));
keep1(pos)=''
keep_form1(pos)=''
keep_form1=strrep(keep_form1,'[cb]','[c]')
keep_form1=strrep(keep_form1,'[sb]','[s]')
keep_form1=strrep(keep_form1,'[sm]','[s]')
keep_form1=strrep(keep_form1,'[cm]','[c]')
keep_form1=strrep(keep_form1,'[mm]','[m]')
keep_form1=strrep(keep_form1,'[mb]','[m]')

keep_form1=strrep(keep_form1,'[tb]','[t]')
keep_form1=strrep(keep_form1,'[tm]','[t]')
keep_form1=strrep(keep_form1,'[ib]','[i]')
keep_form1=strrep(keep_form1,'[im]','[i]')

keep1=strrep(keep1,'_[cb]_[mb]','_[c]_[m]')
keep1=strrep(keep1,'_[cm]_[mm]','_[c]_[m]')
keep1=strrep(keep1,'_[sb]_[cb]','_[s]_[c]')
keep1=strrep(keep1,'_[sm]_[cm]','_[s]_[c]')

keep1=strrep(keep1,'_[cm]_[im]','_[c]_[i]')
keep1=strrep(keep1,'_[mm]_[cm]','_[m]_[c]')
keep1=strrep(keep1,'_[mb]_[cb]','_[m]_[c]')
keep1=strrep(keep1,'_[cb]_[ib]','_[c]_[i]')
keep1=strrep(keep1,'_[sb]_[tb]','_[s]_[t]')
keep1=strrep(keep1,'_[cb]_[sb]','_[c]_[s]')
keep1=strrep(keep1,'_[cm]_[sm]','_[c]_[s]')

keep1=strrep(keep1,'[sb]_[cb]','_[s]_[c]')
keep1=strrep(keep1,'[sm]_[cm]','_[s]_[c]')
keep1=strrep(keep1,'[mb]_[cb]','_[m]_[c]')
keep1=strrep(keep1,'[mm]_[cm]','_[m]_[c]')
cc2=~(contains(keep1,'[cb]'));
keep1=keep1(cc2);
keep_form1=keep_form1(cc2);
cc3=~(contains(keep1,'[cm]'));
keep2=keep1(cc3);
keep_form2=keep_form1(cc3);
keep_form2(200)='';
keep2(200)=''
keep_form2(200)='';
keep2(200)=''
rxns_to_add=vertcat(cheng_add(:,1),comps(:,1),keep2);
forms_to_add=vertcat(cheng_add(:,2),comps(:,2),keep_form2);

rev1=[];lb1=[];ub1=[];
for n=1:length(forms_to_add)
if contains(forms_to_add(n),'<-->') || contains(forms_to_add(n),'=')
    rev1=[rev1,1];
    lb1=[lb1,-1000];
    ub1=[ub1,1000];
else
    rev1=[rev1,0];
    lb1=[lb1,0];
    ub1=[ub1,1000];
end
end
pos=find(contains(forms_to_add, 'Light[t]'))
forms_to_add=strrep(forms_to_add,'Light[s]','Light_Compound[s]')
forms_to_add=strrep(forms_to_add,'Light[t]','Light_Compound[t]')


for n=1:length(rxns_to_add)
     model2 = addReaction(model2,rxns_to_add{n},forms_to_add{n},[],0,lb1(n),ub1(n));
end
pos=find(contains(model2.mets, 'NADH-P-OR-NOP[c]'))
[rxnslist,rxnform]=findRxnsFromMets(model2,{'NADH-P-OR-NOP[c]'})
realpos=find(contains(model2.mets,'NADH[c]'))
rxnpos=find(contains(model2.rxns,rxnslist))
model2.S(pos,rxnpos(1))=0
model2.S(pos,rxnpos(2))=0

model2.S(realpos,rxnpos(1))=1
model2.S(realpos,rxnpos(2))=1

pos=find(contains(model2.mets, 'NAD-P-OR-NOP[c]'))
[rxnslist,rxnform]=findRxnsFromMets(model2,{'NAD-P-OR-NOP[c]'})
realpos=find(contains(model2.mets,'NAD[c]'))
rxnpos=find(contains(model2.rxns,rxnslist))
model2.S(pos,rxnpos(1))=0
model2.S(pos,rxnpos(2))=0

model2.S(realpos(1),rxnpos(1))=-1
model2.S(realpos(1),rxnpos(2))=-1

model2=removeMets(model2,{'NAD-P-OR-NOP[c]','NADH-P-OR-NOP[c]'})
%% localizations 
% finding rxns not in rxns_to_add and base1
left_to_check=setdiff(model2.rxns,base_rxns10);
left_to_check=setdiff(left_to_check,blue);
left_to_check=setdiff(left_to_check,rxns_to_add);
% 145 rxns to check location
form_left=printRxnFormula(model2,left_to_check)
jeffloo=cell2table(horzcat(left_to_check,form_left))
% remove rxns
model2=removeRxns(model2,left_to_check)
% add rxns with new localization
newloc=readcell('with_localization.xlsx')
rxns_newloc=newloc(:,1);
forms_newloc=newloc(:,2);

rev1=[];lb1=[];ub1=[];
for n=1:length(forms_newloc)
if contains(forms_newloc(n),'<-->') || contains(forms_newloc(n),'=')
    rev1=[rev1,1];
    lb1=[lb1,-1000];
    ub1=[ub1,1000];
else
    rev1=[rev1,0];
    lb1=[lb1,0];
    ub1=[ub1,1000];
end
end

for n=1:length(rxns_newloc)
     model2 = addReaction(model2,rxns_newloc{n},forms_newloc{n},[],0,lb1(n),ub1(n));
end
model2=changeObjective(model2,'biomass')

pos_bio=(find(contains(model2.mets,'biomass[m]')))
model2.mets{pos_bio}='biomass[c]'
form='biomass[c] ->';
     model2 = addReaction(model2,'EX_Bio',form,[],0,0,1000);

     % bio_mets=findMetsFromRxns(model4,{'biomass'})
     % fix=intersect(bio_mets,newblock)
     % model4=changeRxnBounds(model4,'L-ARABINOKINASE-RXN',-1000,'l')
cel=find(contains(model2.mets,'CELLULOSE'))
     form10='CPD-12575[c]   ->   CELLULOSE_Compound[c] + UDP[c]'
         
     model2 = addReaction(model2,'CELLULOSE-SYNTHASE-UDP-FORMING-RXN',form10,[],0,0,1000);

%% change direction
direct={'2.7.7.14-RXN','L-ARABINOKINASE-RXN','RXN-4308'}
      model2=changeRxnBounds(model2,direct,-1000,'l')
