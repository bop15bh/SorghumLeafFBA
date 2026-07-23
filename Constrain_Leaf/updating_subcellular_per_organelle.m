%% redoing assigning subcellular one organelle at a time
clear 
close ALL

changeCobraSolver('glpk');
%load('Leaf_balanced_FINALSEP25.mat')
 load('Leaf_model_biomass.mat')
 %trans=readcell('old_transport.xlsx');
 %trans=trans(:,1);% nbefore 19 was close 26 closest, 28,29,30,35, 38, 41,42
 % 49 best 
 %  nitrate only: 17 close, 19 bit weird, 22 weirdish, 23 but min mdh, got up to 35
% 3 remo: 19 weirdish, 20 is very good. 
% 3rd try 31 maybe 33 looks reasonable 
%4th try no dups 22 weird, 25 also weird 40 is v close, 42 only1, 49
%mental, 54 close (2 inc), 55 malo 0, 57 wild
 model=removeRxns(model,{'ATR_NITRITE_[cb]_[cm]','RXN-14932_1[M]','RXN-14932_2[M]'})
%model=removeRxns(model,{'ATR_NITRITE_[cb]_[cm]'});
model=removeRxns(model,{'RXN-13697[B]','RXN-13697[M]'});

%  model=removeRxns(model,trans(40));
% deads = model.mets(detectDeadEnds(model));
% [rxns,rxnforms]=findRxnsFromMets(model,deads);
%  model=removeRxns(model,rxns);
%  deads = model.mets(detectDeadEnds(model))
% deads
% 19, 20, 26,30 broken 39 1 rxn, 40 good 1 rxn, 46 sort of broken
% 49 sort of broken, 54 malo broken 1, 55 close 
%model=removeRxns(model,{'TSA-REDUCT-RXN[B]','TSA-REDUCT-RXN[M]','5-METHYLTHIORIBOSE-KINASE-RXN[B]','5-METHYLTHIORIBOSE-KINASE-RXN[M]'})
%40 after dead removal, only 1 inc, 25 close 
dups=readcell('duplicates1.xlsx');
dups2=readcell('dupdup.xlsx');
model=removeRxns(model,dups(:,1))
model=removeRxns(model,dups2(:,1))
emets=model.mets(find(contains(model.mets,'E-[')))
[rxns,rxnforms]=findRxnsFromMets(model,emets)
model=removeRxns(model,rxns)
%%
% setting bounds to all rxns to be 1000 mmol/g/hr
model.ub(model.ub>1000)=1000;
model.lb(model.lb<-1000)=-1000;
%% Setting limit twrio Rubisco
model = changeRxnBounds(model,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]', 6.5, 'u');
%%
%NGAM 
model = changeRxnBounds(model,'ATPASE-RXN[B]', 0.0229, 'l');
model = changeRxnBounds(model,'ATPASE-RXN[M]', 0.0229, 'l');



%% ensuring c4 cycle 
model=removeRxns(model,{'4.1.1.32-RXN[M]','4.1.1.32-RXN[B]'});
%model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.2, 'l');
%model=removeRxns(model,{'ATR_G3P_[cb]_[cm]','ATR_CPD-1777_[cb]_[cm]'})
%model=removeRxns(model,{'ATR_GLYCERATE_[cb]_[cm]'})

%model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.5, 'l');
%model = changeRxnBounds(model,'MALATE-DEH-RXN_2[M]', 0.2, 'u');

%model=removeRxns(model,{'MALATE-DEH-RXN_2[M]','MALATE-DEH-RXN_2[B]'})

model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.3, 'l');
%model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 0.2, 'l');
model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 0.6, 'l');
%model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 10, 'u');
%model = changeRxnBounds(model,'ATR_MAL_[cb]_[cm]', 0.6, 'l');
%model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 0.9, 'l');
%%model = changeRxnBounds(modenl,'D-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN_2[M]', 0, 'l');
%model = changeRxnBounds(model,'ATR_L-ASPARTATE_[cb]_[cm]', 0.2, 'u');
%model = changeRxnBounds(model,'ATR_G3P_[cb]_[cm]', 0.3, 'u');
%model = changeRxnBounds(model,'ATR_CPD-1777_[cb]_[cm]', 0.3, 'u');


%model=removeRxns(model,{'4.1.1.32-RXN[M]','4.1.1.32-RXN[B]'})%,'RXN-13697[M]','RXN-13697[B]'})
%model=removeRxns(model,{'RXN-14932_2[M]','MALATE-DEH-RXN_1[M]','ATR_PLASTOQUINONE-9_[cb]_[cm]','ATR_G3P_[cb]_[cm]'})
%model=removeRxns(model,{'ATR_CPD-12829_[cb]_[cm]','ATR_PLASTOQUINONE-9_[cb]_[cm]'})
%model=removeRxns(model,{'RXN-14932_2[M]','MALATE-DEH-RXN_1[M]','ATR_G3P_[cb]_[cm]'})
%model=removeRxns(model,{'MALATE-DEH-RXN_1[M]','ATR_PLASTOQUINONE-9_[cb]_[cm]','ATR_G3P_[cb]_[cm]'})

 form=' ADP[cm] + PHOSPHO-ENOL-PYRUVATE[cm] + PROTON[cm] -> ATP[cm] + PYRUVATE[cm] ';
 model=addReaction(model,'PEPDEPHOS-RXN_1[M]',form,[],0,0,1000);
 form=' ADP[cb] + PHOSPHO-ENOL-PYRUVATE[cb] + PROTON[cb] -> ATP[cb] + PYRUVATE[cb] ';
 model=addReaction(model,'PEPDEPHOS-RXN_1[B]',form,[],0,0,1000);
form=' ADP[sm] + PHOSPHO-ENOL-PYRUVATE[sm] + PROTON[sm] -> ATP[sm] + PYRUVATE[sm] ';
 model=addReaction(model,'PEPDEPHOS-RXN_2[M]',form,[],0,0,1000);
 form=' ADP[sb] + PHOSPHO-ENOL-PYRUVATE[sb] + PROTON[sb] -> ATP[sb] + PYRUVATE[sb] ';
 model=addReaction(model,'PEPDEPHOS-RXN_2[B]',form,[],0,0,1000);

%% loading rxns to update
rxns_to_change=readcell('Updated_localizedFormulae.xlsx');
%%
deadsnew = model.mets(detectDeadEnds(model));
%% starting with BS chloroplast first  
chloro_rxns1=rxns_to_change(find(contains(rxns_to_change(:,2),'[sb]')),1);
chloro_rxns=rxns_to_change(find(contains(rxns_to_change(:,2),'[sm]')),1);
chloro_listB=chloro_rxns1;

chloro_findB=erase(chloro_listB,'[B]');

chloro_findB=erase(chloro_findB,'_1');
% isolating original model equations that are in cytosol from chloroplast
% list
not_matched=[]; to_rem={};
for n=1:length(chloro_findB)
    pos=find(contains(model.rxns,chloro_findB{n}));
    if length(pos)>1
        options=erase(model.rxns(pos),'[B]');
        options=erase(options,'_1');
        pos1=find(strcmp(options,chloro_findB{n}));
        if length(pos1)==1
            to_rem=[to_rem,model.rxns(pos(pos1))'];
        else
            not_matched=[not_matched,n]
        end
    elseif  length(pos)==1
        to_rem=[to_rem,model.rxns(pos)'];
    else
    end
end


% removing cytosol reactions 
model1 = removeRxns(model, to_rem);
deadsnew1=model1.mets(detectDeadEnds(model1));
og=optimizeCbModel(model1)

% adding chloroplast reactions
chloroB_ids = {}; chloroB_forms = {};
for k = 1:numel(chloro_listB)
    row = find(strcmp(rxns_to_change(:,1), chloro_listB{k}));
    if numel(row)==1
        chloroB_ids{end+1,1}   = chloro_listB{k};
        chloroB_forms{end+1,1} = rxns_to_change{row,2};
    end
end

for n = 1:length(chloroB_ids)
    % unwrap until you get a string
    form = chloroB_forms{n};
    while iscell(form)
        form = form{1};
    end
    
    if contains(form, '<=>')
        model1 = addReaction(model1, chloroB_ids{n}, form, [], 0, -1000, 1000);
    else
        model1 = addReaction(model1, chloroB_ids{n}, form, [], 0, 0, 1000);
    end
end
deadsnew2 = model1.mets(detectDeadEnds(model1));




% 333 deads
%% adding transporters
trans=readcell('mintz_transporters.xlsx')
added={};
for n=1:length(trans)
            model2 = addReaction(model1, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            deads3=model2.mets(detectDeadEnds(model2));
            fixed=setdiff(deadsnew2,deads3);
            new=setdiff(deads3,deadsnew2)
            if ~isempty(fixed) && isempty(new)
                        model1 = addReaction(model1, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            added=[added,trans{n,1}];
            else
            end
          
end
% fixed 20 blocked mets by adding existing transporters
deads4=model1.mets(detectDeadEnds(model1))

%% adding transporters to fix dead ends
c_count=[];deadsnew=deads4;
for n=1:length(deadsnew)
    %  cytosol mets
    if contains(deadsnew{n},'[c')
        
        met=split(deadsnew{n},'[')
       all_mets = model1.mets(startsWith(model1.mets, [met{1} '[']));
       non_c=setdiff(all_mets,all_mets(find(contains(all_mets,'[c'))));
       b_non_c=non_c(find(contains(non_c,'b]')));
       m_non_c=non_c(find(contains(non_c,'m]')));
       if ~isempty(b_non_c)
       for k=1:length(b_non_c)
           c_count=[c_count,n]
           forma=[met{1} '[cb] <=> ' b_non_c{k}];
           model1 = addReaction(model1,[b_non_c{k} '_[c][B]'], forma, [], 0, -1000, 1000);
       end
       else
       end
       if ~isempty(m_non_c)
       for k=1:length(m_non_c)
           formb=[met{1} '[cm] <=> ' m_non_c{k}];
           model1 = addReaction(model1,[m_non_c{k} '_[c][M]'], formb, [], 0, -1000, 1000);
       end
       else
       end
    else
    end
end
deads5=model1.mets(detectDeadEnds(model1))

% 33 left 


%% adding transporters to fix dead ends

deadsnew1_noc=setdiff(deads5,deads5(find(contains(deads5,'[c'))))
B_deads=deadsnew1_noc
for n=1:length(B_deads)  
  c_count=[c_count,n]
        met=split(B_deads{n},'[')
           forma=[met{1} '[cb] <=> ' B_deads{n}];
           model1 = addReaction(model1,[B_deads{n} '_[c][B]'], forma, [], 0, -1000, 1000);
end
      
  deads6 = model1.mets(detectDeadEnds(model1));
%% 
met_list={}
for n=1:length(chloroB_ids)
met_list=[met_list,findMetsFromRxns(model1,chloroB_ids(n))'];
end
met_list=unique(met_list);
model2=model1;sol=[];
for n=[ 132 143 200 228 274 276]  
        met=split(met_list{n},'[');
           forma=[met{1} '[cb] <=> ' met_list{n}];
           model2 = addReaction(model2,[met_list{n} '_[c][B]'], forma, [], 0, -1000, 1000);
           ro=optimizeCbModel(model2);
           sol=[sol,ro.f];

end

%% model grows with updated BS chloroplast 
%% updating mesophyll chloroplast reactions
chloro_listM = chloro_rxns;                 % the [sm] list read earlier
 chloro_findM = erase(chloro_listM,'[M]');
chloro_findM = erase(chloro_findM,{'_1','_2'});   % strip both variants
 
% isolating original model equations that are in cytosol from mesophyll list
not_matchedM = []; to_remM = {};
for n = 1:length(chloro_findM)
    pos = find(strcmp(erase(erase(model2.rxns,'[M]'),{'_1','_2'}), chloro_findM{n}));
    % strcmp on stripped IDs = exact base match, no substring bleed
    if numel(pos) == 1
        to_remM = [to_remM, model2.rxns(pos)'];
    elseif numel(pos) > 1
        not_matchedM = [not_matchedM, n]
    else
    end
end

% remove cytosol mesophyll reactions
model2 = removeRxns(model2, to_remM);
deadsM1 = model2.mets(detectDeadEnds(model2));
og = optimizeCbModel(model2)     % expected 0 mid-surgery
 % add mesophyll chloroplast reactions
chloroM_ids = {}; chloroM_forms = {};
for k = 1:numel(chloro_listM)
    row = find(strcmp(rxns_to_change(:,1), chloro_listM{k}));
    if numel(row) == 1
        chloroM_ids{end+1,1}   = chloro_listM{k};
        chloroM_forms{end+1,1} = rxns_to_change{row,2};
    end
end


for n = 1:length(chloroM_ids)
    form = chloroM_forms{n};
    while iscell(form), form = form{1}; end
    if contains(form,'<=>')
        model2 = addReaction(model2, chloroM_ids{n}, form, [], 0, -1000, 1000);
    else
        model2 = addReaction(model2, chloroM_ids{n}, form, [], 0, 0, 1000);
    end
end
deadsM2 = model2.mets(detectDeadEnds(model2))

%% add existing transporters that fix dead ends without creating new ones
trans = readcell('mintz_transporters.xlsx');
addedM = {};
for n = 1:length(trans)
    model3 = addReaction(model2, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
    deads3 = model3.mets(detectDeadEnds(model3));
    fixed  = setdiff(deadsM2, deads3);
    new    = setdiff(deads3, deadsM2);
    if ~isempty(fixed) && isempty(new)
        model2  = addReaction(model2, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
        addedM  = [addedM, trans{n,1}];
    end
end
deadsM4 = model2.mets(detectDeadEnds(model2))
% fixed 19 rxns


%% dead-end transporters, cytosol mets  (mesophyll side => [cm])
c_count = []; deadsnew = deadsM4;
for n = 1:length(deadsnew)
    if contains(deadsnew{n},'[c')
        met = split(deadsnew{n},'[');
        all_mets = model2.mets(startsWith(model2.mets, [met{1} '[']));
        non_c = setdiff(all_mets, all_mets(find(contains(all_mets,'[c'))));
        m_non_c = non_c(find(contains(non_c,'m]')));
        b_non_c = non_c(find(contains(non_c,'b]')));
        if ~isempty(m_non_c)
            for k = 1:length(m_non_c)
                c_count = [c_count, n];
                formb = [met{1} '[cm] <=> ' m_non_c{k}];
                model2 = addReaction(model2,[m_non_c{k} '_[c][M]'], formb, [], 0, -1000, 1000);
            end
        end
        if ~isempty(b_non_c)
            for k = 1:length(b_non_c)
                forma = [met{1} '[cb] <=> ' b_non_c{k}];
                model2 = addReaction(model2,[b_non_c{k} '_[c][B]'], forma, [], 0, -1000, 1000);
            end
        end
    end
end
deadsM5 = model2.mets(detectDeadEnds(model2))
 % 33 dead ends left

deadsnew1_noc=setdiff(deadsM5,deadsM5(find(contains(deadsM5,'[c'))))
B_deads=deadsnew1_noc
for n=1:length(B_deads)  
  c_count=[c_count,n]
        met=split(B_deads{n},'[')
           forma=[met{1} '[cm] <=> ' B_deads{n}];
           model2 = addReaction(model2,[B_deads{n} '_[c][M]'], forma, [], 0, -1000, 1000);
end
      
  deadsM6 = model2.mets(detectDeadEnds(model2));
  og=optimizeCbModel(model2)
  
  % fixing growth in mesophyll cell
  met_list1={}
for n=1:length(chloroM_ids)
met_list1=[met_list1,findMetsFromRxns(model2,chloroM_ids(n))'];
end
met_list1=unique(met_list1); 
n=[ 132 143 200 228 274 276]
 c_version=strrep(met_list1,'[sm]','[cm]');
 met_list3=intersect(model2.mets,c_version);
 met_list4={};
 for n=1:length(met_list3)
     metty=strrep(met_list3{n},'[cm]','[sm]')
     pos=find(strcmp(met_list1,metty))
 met_list4=[met_list4,met_list1(pos)]
 end

 % first=met_list(n);
 % first=strrep(first,'[sb]','[sm]');
 % shared=intersect(first,met_list1)
model3=model2;sol=[];
 n=[37 48 108 139 151 152]
for n=[37 48 108 139 151 152]
        met=split(met_list4{n},'[');
           forma=[met{1} '[cm] <=> ' met_list4{n}];
           model3 = addReaction(model3,[met_list4{n} '_[c][M]'], forma, [], 0, -1000, 1000);
           ro=optimizeCbModel(model3);
           sol=[sol,ro.f];
end
sol

deads=model3.mets(detectDeadEnds(model3));

save('chloro_updated.mat','model3')
save('model_parameterized.mat','model')

%%

clear
load('chloro_updated.mat');
load('model_parameterized.mat');
%% loading rxns to update
rxns_to_change=readcell('Updated_localizedFormulae.xlsx');
%%
deadsnew = model3.mets(detectDeadEnds(model3));
%% BS mitochondria   

mito_rxns1=rxns_to_change(find(contains(rxns_to_change(:,2),'[mb]')),1);
mito_rxns=rxns_to_change(find(contains(rxns_to_change(:,2),'[mm]')),1);
mito_list=vertcat(mito_rxns1,mito_rxns);
mito_listB=mito_rxns1;

mito_findB=erase(mito_listB,'[B]');
mito_findB=erase(mito_findB,'_1');
mito_findB=erase(mito_findB,'_2');

% isolating original model equations that are in cytosol from chloroplast
% list

not_matched3=[]; to_rem3={};
for n=1:length(mito_findB)
    pos=find(contains(model.rxns,mito_findB{n}));
    if length(pos)>1
        options=erase(model.rxns(pos),'[B]');
        options=erase(options,'_1');
        pos1=find(strcmp(options,mito_findB{n}));
        if length(pos1)==1
            to_rem3=[to_rem3,model.rxns(pos(pos1))'];
        else
            not_matched3=[not_matched3,n]
        end
    elseif  length(pos)==1
        to_rem3=[to_rem3,model.rxns(pos)'];
    else
    end
end


% removing cytosol reactions 
model4 = removeRxns(model3, to_rem3);
deadsnewmb=model4.mets(detectDeadEnds(model4));
og=optimizeCbModel(model4)


% adding mito reactions
mitoB_ids = {}; mitoB_forms = {};
for k = 1:numel(mito_listB)
    row = find(strcmp(rxns_to_change(:,1), mito_listB{k}));
    if numel(row)==1
        mitoB_ids{end+1,1}   = mito_listB{k};
        mitoB_forms{end+1,1} = rxns_to_change{row,2};
    end
end


for n = 1:length(mitoB_ids)
    % unwrap until you get a string
    form = mitoB_forms{n};
    while iscell(form)
        form = form{1};
    end
    
    if contains(form, '<=>')
        model4 = addReaction(model4, mitoB_ids{n}, form, [], 0, -1000, 1000);
    else
        model4 = addReaction(model4, mitoB_ids{n}, form, [], 0, 0, 1000);
    end
end
deadsnew2mb = model4.mets(detectDeadEnds(model4));

% 50 dead ends 

%% adding transporters
trans=readcell('mintz_transporters.xlsx')
addedmb={};
for n=1:length(trans)
            model5 = addReaction(model4, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            deads3mb=model5.mets(detectDeadEnds(model5));
            fixed=setdiff(deadsnew2mb,deads3mb);
            new=setdiff(deads3mb,deadsnew2mb)
            if ~isempty(fixed) && isempty(new)
                        model4 = addReaction(model4, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            addedmb=[addedmb,trans{n,1}];
            else
            end
          
end

% fixed 6 blocked mets by adding existing transporters
deads4mb=model4.mets(detectDeadEnds(model4))

%% adding transporters to fix cytosolic dead ends
c_count3=[];deadsnew=deads4mb;
for n=1:length(deadsnew)
    %  cytosol mets
    if contains(deadsnew{n},'[c')
        
        met=split(deadsnew{n},'[')
       all_mets = model4.mets(startsWith(model4.mets, [met{1} '[']));
       non_c=setdiff(all_mets,all_mets(find(contains(all_mets,'[c'))));
       b_non_c=non_c(find(contains(non_c,'b]')));
       m_non_c=non_c(find(contains(non_c,'m]')));
       if ~isempty(b_non_c)
       for k=1:length(b_non_c)
           c_count3=[c_count3,n]
           forma=[met{1} '[cb] <=> ' b_non_c{k}];
           model4 = addReaction(model4,[b_non_c{k} '_[c][B]'], forma, [], 0, -1000, 1000);
       end
       else
       end
       if ~isempty(m_non_c)
       for k=1:length(m_non_c)
           formb=[met{1} '[cm] <=> ' m_non_c{k}];
           model4 = addReaction(model4,[m_non_c{k} '_[c][M]'], formb, [], 0, -1000, 1000);
       end
       else
       end
    else
    end
end
deads5mb=model4.mets(detectDeadEnds(model4))

% 20 left 
% fixing non-cytosolic dead ends 

deadsnew1_noc=setdiff(deads5mb,deads5mb(find(contains(deads5mb,'[c'))))
mB_deads=deadsnew1_noc; 
for n=1:length(mB_deads)  
  c_count3=[c_count3,n]
        met=split(mB_deads{n},'[')
           forma=[met{1} '[cb] <=> ' mB_deads{n}];
           model4 = addReaction(model4,[mB_deads{n} '_[c][B]'], forma, [], 0, -1000, 1000);
end
      
  deads6mb = model4.mets(detectDeadEnds(model4));

  forma1='3-oxo-octanoyl-ACPs[cb] <=> 3-oxo-octanoyl-ACPs[sb]';
           model4 = addReaction(model4,['3-oxo-octanoyl-ACPs[sb]' '_[c][B]'], forma1, [], 0, -1000, 1000);
    forma2='Acetoacetyl-ACPs[cb] <=> Acetoacetyl-ACPs[sb]';
           model4 = addReaction(model4,['Acetoacetyl-ACPs[sb]' '_[c][B]'], forma2, [], 0, -1000, 1000);
 forma3='MALONYL-ACP[cb] <=> MALONYL-ACP[sb]';
           model4 = addReaction(model4,['MALONYL-ACP[sb]' '_[c][B]'], forma3, [], 0, -1000, 1000);
 
           deads7mb = model4.mets(detectDeadEnds(model4));
           og=optimizeCbModel(model4)

%% % fixing growth in mesophyll cell
  met_list1={}
for n=1:length(mitoB_ids)
met_list1=[met_list1,findMetsFromRxns(model4,mitoB_ids(n))'];
end
met_list1=unique(met_list1); 

 c_version=strrep(met_list1,'[mb]','[cb]');
 met_list3=intersect(model4.mets,c_version);
 met_list4={};
 for n=1:length(met_list3)
     metty=strrep(met_list3{n},'[cb]','[mb]')
     pos=find(strcmp(met_list1,metty))
 met_list4=[met_list4,met_list1(pos)]
 end

 % first=met_list(n);
 % first=strrep(first,'[sb]','[sm]');
 % shared=intersect(first,met_list1)
model5=model4;
 n=52
        met=split(met_list4{n},'[');
           forma=[met{1} '[cb] <=> ' met_list4{n}];
           model5 = addReaction(model5,[met_list4{n} '_[c][B]'], forma, [], 0, -1000, 1000);
           ro=optimizeCbModel(model5);
           
% only needed to add 1 transport rxn to fix growth 

%% mesophyll mitochondria 
mito_listM=mito_rxns;


mito_findM=erase(mito_listM,'[M]');
mito_findM=erase(mito_findM,'_1');
mito_findM=erase(mito_findM,'_2');

% isolating original model equations that are in cytosol from chloroplast
% list

not_matched4=[]; to_rem4={};
for n=1:length(mito_findM)
    pos=find(contains(model.rxns,mito_findM{n}));
    if length(pos)>1
        options=erase(model.rxns(pos),'[M]');
        options=erase(options,'_1');
        pos1=find(strcmp(options,mito_findM{n}));
        if length(pos1)==1
            to_rem4=[to_rem4,model.rxns(pos(pos1))'];
        else
            not_matched4=[not_matched4,n]
        end
    elseif  length(pos)==1
        to_rem4=[to_rem4,model.rxns(pos)'];
    else
    end
end


% removing cytosol reactions 
model6 = removeRxns(model5, to_rem4);
deadsnewmb=model6.mets(detectDeadEnds(model6));
og=optimizeCbModel(model6)


% adding mito reactions
mitoM_ids = {}; mitoM_forms = {};
for k = 1:numel(mito_listM)
    row = find(strcmp(rxns_to_change(:,1), mito_listM{k}));
    if numel(row)==1
        mitoM_ids{end+1,1}   = mito_listM{k};
        mitoM_forms{end+1,1} = rxns_to_change{row,2};
    end
end



for n = 1:length(mitoM_ids)
    % unwrap until you get a string
    form = mitoM_forms{n};
    while iscell(form)
        form = form{1};
    end
    
    if contains(form, '<=>')
        model6 = addReaction(model6, mitoM_ids{n}, form, [], 0, -1000, 1000);
    else
        model6 = addReaction(model6, mitoM_ids{n}, form, [], 0, 0, 1000);
    end
end
deadsnew2mm = model6.mets(detectDeadEnds(model6));
% 47 new dead ends 


%% adding transporters
trans=readcell('mintz_transporters.xlsx')
addedmm={};
for n=1:length(trans)
            model7 = addReaction(model6, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            deads3mm=model7.mets(detectDeadEnds(model7));
            fixed=setdiff(deadsnew2mm,deads3mm);
            new=setdiff(deads3mm,deadsnew2mm)
            if ~isempty(fixed) && isempty(new)
                        model6 = addReaction(model6, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            addedmm=[addedmm,trans{n,1}];
            else
            end
          
end

% fixed 3 blocked mets by adding existing transporters
deads4mm=model6.mets(detectDeadEnds(model6))

%% adding transporters to fix cytosolic dead ends
c_count5=[];deadsnew=deads4mm;
for n=1:length(deadsnew)
    %  cytosol mets
    if contains(deadsnew{n},'[c')
        
        met=split(deadsnew{n},'[')
       all_mets = model6.mets(startsWith(model6.mets, [met{1} '[']));
       non_c=setdiff(all_mets,all_mets(find(contains(all_mets,'[c'))));
       b_non_c=non_c(find(contains(non_c,'b]')));
       m_non_c=non_c(find(contains(non_c,'m]')));
       if ~isempty(b_non_c)
       for k=1:length(b_non_c)
           c_count5=[c_count5,n]
           forma=[met{1} '[cb] <=> ' b_non_c{k}];
           model6 = addReaction(model6,[b_non_c{k} '_[c][B]'], forma, [], 0, -1000, 1000);
       end
       else
       end
       if ~isempty(m_non_c)
       for k=1:length(m_non_c)
            c_count5=[c_count5,n]
           formb=[met{1} '[cm] <=> ' m_non_c{k}];
           model6 = addReaction(model6,[m_non_c{k} '_[c][M]'], formb, [], 0, -1000, 1000);
       end
       else
       end
    else
    end
end
deads5mm=model6.mets(detectDeadEnds(model6))

% 23 left 
% fixing non-cytosolic dead ends 

deadsnew1_noc=setdiff(deads5mm,deads5mm(find(contains(deads5mm,'[c'))))
mm_deads=deadsnew1_noc; 
for n=1:length(mm_deads)  
  c_count5=[c_count5,n]
        met=split(mm_deads{n},'[')
           forma=[met{1} '[cm] <=> ' mm_deads{n}];
           model6 = addReaction(model6,[mm_deads{n} '_[c][M]'], forma, [], 0, -1000, 1000);
end
      
  deads6mm = model6.mets(detectDeadEnds(model6));

   forma4='3-oxo-octanoyl-ACPs[cm] <=> 3-oxo-octanoyl-ACPs[sm]';
           model6 = addReaction(model6,['3-oxo-octanoyl-ACPs[sm]' '_[c][M]'], forma4, [], 0, -1000, 1000);
    forma5='Acetoacetyl-ACPs[cm] <=> Acetoacetyl-ACPs[sm]';
           model6 = addReaction(model6,['Acetoacetyl-ACPs[sm]' '_[c][M]'], forma5, [], 0, -1000, 1000);
 forma6='MALONYL-ACP[cm] <=> MALONYL-ACP[sm]';
           model6 = addReaction(model6,['MALONYL-ACP[sm]' '_[c][M]'], forma6, [], 0, -1000, 1000);
   deads7mm = model6.mets(detectDeadEnds(model6));


%%

%% % fixing growth in mesophyll cell
  met_list1={}
for n=1:length(mitoM_ids)
met_list1=[met_list1,findMetsFromRxns(model6,mitoM_ids(n))'];
end
met_list1=unique(met_list1); 

 c_version=strrep(met_list1,'[mm]','[cm]');
 met_list3=intersect(model6.mets,c_version);
 met_list4={};
 for n=1:length(met_list3)
     metty=strrep(met_list3{n},'[cm]','[mm]')
     pos=find(strcmp(met_list1,metty))
 met_list4=[met_list4,met_list1(pos)]
 end

 % first=met_list(n);
 % first=strrep(first,'[sb]','[sm]');
 % shared=intersect(first,met_list1)
model7=model6;
 n=52
        met=split(met_list4{n},'[');
           forma=[met{1} '[cm] <=> ' met_list4{n}];
           model7 = addReaction(model7,[met_list4{n} '_[c][M]'], forma, [], 0, -1000, 1000);
           ro=optimizeCbModel(model7);
           
% only needed to add 1 transport rxn to fix growth 
%% Peroxisome bundle sheath 

perox_rxns1=rxns_to_change(find(contains(rxns_to_change(:,2),'[xb]')),1);
perox_rxns=rxns_to_change(find(contains(rxns_to_change(:,2),'[xm]')),1);
% 149 peroxisome rxns 

perox_listB=perox_rxns1;

 perox_findB=strrep(perox_listB,'_1[','[');
  perox_findB=strrep(perox_findB,'_2[','[');
perox_findB=erase(perox_findB,'[B]');


perox_findB=unique(perox_findB)
% isolating original model equations that are in cytosol from chloroplast
% list

not_matched5=[]; to_rem5={};
for n=1:length(perox_findB)
    pos=find(contains(model.rxns,perox_findB{n}));
    if length(pos)>1
        options=strrep(model.rxns(pos),'_2[','[');
        options=erase(options,'[B]');
               % options=erase(options,'_1');
                

          pos1=find(strcmp(options,perox_findB{n}));
        if length(pos1)==1
            to_rem5=[to_rem5,model.rxns(pos(pos1))'];
        else
            not_matched5=[not_matched5,n]
        end
    elseif  length(pos)==1
        to_rem5=[to_rem5,model.rxns(pos)'];
    else
    end
end



% removing cytosol reactions 
model8 = removeRxns(model7, to_rem5);
deadsnewxb=model8.mets(detectDeadEnds(model8));
og=optimizeCbModel(model8)


% adding perox reactions
peroxB_ids = {}; peroxB_forms = {};
for k = 1:numel(perox_listB)
    row = find(strcmp(rxns_to_change(:,1), perox_listB{k}));
    if numel(row)==1
        peroxB_ids{end+1,1}   = perox_listB{k};
        peroxB_forms{end+1,1} = rxns_to_change{row,2};
    end
end



for n = 1:length(peroxB_ids)
    % unwrap until you get a string
    form = peroxB_forms{n};
    while iscell(form)
        form = form{1};
    end
    
    if contains(form, '<=>')
        model8 = addReaction(model8, peroxB_ids{n}, form, [], 0, -1000, 1000);
    else
        model8 = addReaction(model8, peroxB_ids{n}, form, [], 0, 0, 1000);
    end
end
deadsnew2xb = model8.mets(detectDeadEnds(model8));
% 47 new dead ends 


%% adding transporters
trans=readcell('mintz_transporters.xlsx')
addedxb={};
for n=1:length(trans)
            model9 = addReaction(model8, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            deads3xb=model9.mets(detectDeadEnds(model9));
            fixed=setdiff(deadsnew2xb,deads3xb);
            new=setdiff(deads3xb,deadsnew2xb)
            if ~isempty(fixed) && isempty(new)
                        model8 = addReaction(model8, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            addedxb=[addedxb,trans{n,1}];
            else
            end
          
end

% fixed 10 blocked mets by adding existing transporters
deads4xb=model8.mets(detectDeadEnds(model8))

%% adding transporters to fix cytosolic dead ends
c_count6=[];deadsnew=deads4xb;
for n=1:length(deadsnew)
    %  cytosol mets
    if contains(deadsnew{n},'[c')
        
        met=split(deadsnew{n},'[')
       all_mets = model8.mets(startsWith(model8.mets, [met{1} '[']));
       non_c=setdiff(all_mets,all_mets(find(contains(all_mets,'[c'))));
       b_non_c=non_c(find(contains(non_c,'b]')));
       m_non_c=non_c(find(contains(non_c,'m]')));
       if ~isempty(b_non_c)
       for k=1:length(b_non_c)
           c_count6=[c_count6,n]
           forma=[met{1} '[cb] <=> ' b_non_c{k}];
           model8 = addReaction(model8,[b_non_c{k} '_[c][B]'], forma, [], 0, -1000, 1000);
       end
       else
       end
       if ~isempty(m_non_c)
       for k=1:length(m_non_c)
            c_count6=[c_count6,n]
           formb=[met{1} '[cm] <=> ' m_non_c{k}];
           model8 = addReaction(model8,[m_non_c{k} '_[c][M]'], formb, [], 0, -1000, 1000);
       end
       else
       end
    else
    end
end
deads5xb=model8.mets(detectDeadEnds(model8))

% 22 left 
% fixing non-cytosolic dead ends 

deadsnew1_noc=setdiff(deads5xb,deads5xb(find(contains(deads5xb,'[c'))))
xb_deads=deadsnew1_noc; 
for n=1:length(xb_deads)  
  c_count6=[c_count6,n]
        met=split(xb_deads{n},'[')
           forma=[met{1} '[cb] <=> ' xb_deads{n}];
           model8 = addReaction(model8,[xb_deads{n} '_[c][B]'], forma, [], 0, -1000, 1000);
end
      
  deads6xb = model8.mets(detectDeadEnds(model8));

  forma7='GLY[cb] <=> GLY[xb]';
           model8 = addReaction(model8,['GLY[xb]' '_[c][B]'], forma7, [], 0, -1000, 1000);
 
  forma8='GLY[cb] <=> GLY[mb]';
           model8 = addReaction(model8,['GLY[mb]' '_[c][B]'], forma8, [], 0, -1000, 1000);
 forma9='GLT[cb] <=> GLT[xb]';
           model8 = addReaction(model8,['GLT[xb]' '_[c][B]'], forma9, [], 0, -1000, 1000);
 forma10='2-KETOGLUTARATE[cb] <=> 2-KETOGLUTARATE[xb]';
           model8 = addReaction(model8,['2-KETOGLUTARATE[xb]' '_[c][B]'], forma10, [], 0, -1000, 1000);
 
           model8=removeRxns(model8, {'CPD-10262[xb]_[c][B]','CPD-14271[xb]_[c][B]','CPD0-2253[xb]_[c][B]','S_CPD-14275[xb]_[c][B]','ATR_S-ALLANTOIN_[cm]_[cb]'})
            
forma11='SUPER-OXIDE[sb] <=> SUPER-OXIDE[cb]';
           model8 = addReaction(model8,['SUPER-OXIDE[sb]' '_[c][B]'], forma11, [], 0, -1000, 1000);
 deads6xb = model8.mets(detectDeadEnds(model8))
           %% % fixing growth in bundle sheath cell
  met_list1={}
for n=1:length(peroxB_ids)
met_list1=[met_list1,findMetsFromRxns(model8,peroxB_ids(n))'];
end
met_list1=unique(met_list1); 

 c_version=strrep(met_list1,'[xb]','[cb]');
 met_list3=intersect(model8.mets,c_version);
 met_list4={};
 for n=1:length(met_list3)
     metty=strrep(met_list3{n},'[cb]','[xb]')
     pos=find(strcmp(met_list1,metty))
 met_list4=[met_list4,met_list1(pos)]
 end

 % first=met_list(n);
 % first=strrep(first,'[sb]','[sm]');
 % shared=intersect(first,met_list1)
model9=model8;sol=[];
 for n=40
        met=split(met_list4{n},'[');
           forma=[met{1} '[cb] <=> ' met_list4{n}];
           model9 = addReaction(model9,[met_list4{n} '_[c][B]'], forma, [], 0, -1000, 1000);
           ro=optimizeCbModel(model9);
           sol=[sol,ro.f];
 end      
 sol
% only needed to add 1 transport rxn to fix growth 

%% peroxisome mesophyll
perox_rxns=rxns_to_change(find(contains(rxns_to_change(:,2),'[xm]')),1);
% 149 peroxisome rxns 

perox_listM=perox_rxns;

 perox_findM=strrep(perox_listM,'_1[','[');
  perox_findM=strrep(perox_findM,'_2[','[');
perox_findM=erase(perox_findM,'[M]');


perox_findM=unique(perox_findM)
% isolating original model equations that are in cytosol from chloroplast
% list

not_matched6=[]; to_rem6={};
for n=1:length(perox_findM)
    pos=find(contains(model.rxns,perox_findM{n}));
    if length(pos)>1
        options=strrep(model.rxns(pos),'_2[','[');
        options=erase(options,'[M]');
               % options=erase(options,'_1');
                

          pos1=find(strcmp(options,perox_findM{n}));
        if length(pos1)==1
            to_rem6=[to_rem6,model.rxns(pos(pos1))'];
        else
            not_matched6=[not_matched6,n]
        end
    elseif  length(pos)==1
        to_rem6=[to_rem6,model.rxns(pos)'];
    else
    end
end



% removing cytosol reactions 
model10 = removeRxns(model9, to_rem6);
deadsnewxm=model10.mets(detectDeadEnds(model10));
og=optimizeCbModel(model10)


% adding perox reactions
peroxM_ids = {}; peroxM_forms = {};
for k = 1:numel(perox_listM)
    row = find(strcmp(rxns_to_change(:,1), perox_listM{k}));
    if numel(row)==1
        peroxM_ids{end+1,1}   = perox_listM{k};
        peroxM_forms{end+1,1} = rxns_to_change{row,2};
    end
end



for n = 1:length(peroxM_ids)
    % unwrap until you get a string
    form = peroxM_forms{n};
    while iscell(form)
        form = form{1};
    end
    
    if contains(form, '<=>')
        model10 = addReaction(model10, peroxM_ids{n}, form, [], 0, -1000, 1000);
    else
        model10 = addReaction(model10, peroxM_ids{n}, form, [], 0, 0, 1000);
    end
end
deadsnew2xm = model10.mets(detectDeadEnds(model10));
% 41 new dead ends 


%% adding transporters
trans=readcell('mintz_transporters.xlsx')
addedxm={};
for n=1:length(trans)
            model11 = addReaction(model10, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            deads3xm=model11.mets(detectDeadEnds(model11));
            fixed=setdiff(deadsnew2xm,deads3xm);
            new=setdiff(deads3xm,deadsnew2xm)
            if ~isempty(fixed) && isempty(new)
                        model10 = addReaction(model10, trans{n,1}, trans{n,2}, [], 0, -1000, 1000);
            addedxm=[addedxm,trans{n,1}];
            else
            end
          
end

model10=removeRxns(model10,{'ATR_S-ALLANTOIN_[cm]_[cb]'});
% fixed 8 blocked mets by adding existing transporters
deads4xm=model10.mets(detectDeadEnds(model10));

%% adding transporters to fix cytosolic dead ends
c_count7=[];deadsnew=deads4xm;
for n=1:length(deadsnew)
    %  cytosol mets
    if contains(deadsnew{n},'[c')
        
        met=split(deadsnew{n},'[')
       all_mets = model10.mets(startsWith(model10.mets, [met{1} '[']));
       non_c=setdiff(all_mets,all_mets(find(contains(all_mets,'[c'))));
       b_non_c=non_c(find(contains(non_c,'b]')));
       m_non_c=non_c(find(contains(non_c,'m]')));
       if ~isempty(b_non_c)
       for k=1:length(b_non_c)
           c_count7=[c_count7,n]
           forma=[met{1} '[cb] <=> ' b_non_c{k}];
           model10 = addReaction(model10,[b_non_c{k} '_[c][B]'], forma, [], 0, -1000, 1000);
       end
       else
       end
       if ~isempty(m_non_c)
       for k=1:length(m_non_c)
            c_count7=[c_count7,n]
           formb=[met{1} '[cm] <=> ' m_non_c{k}];
           model10 = addReaction(model10,[m_non_c{k} '_[c][M]'], formb, [], 0, -1000, 1000);
       end
       else
       end
    else
    end
end
deads5xm=model10.mets(detectDeadEnds(model10))

% 26 left 
% fixing non-cytosolic dead ends 

deadsnew1_noc=setdiff(deads5xm,deads5xm(find(contains(deads5xm,'[c'))))
xm_deads=deadsnew1_noc; 
for n=1:length(xm_deads)  
  c_count7=[c_count7,n]
        met=split(xm_deads{n},'[')
           forma=[met{1} '[cm] <=> ' xm_deads{n}];
           model10 = addReaction(model10,[xm_deads{n} '_[c][M]'], forma, [], 0, -1000, 1000);
end

  deads6xm = model10.mets(detectDeadEnds(model10));
model11=removeRxns(model10, {'CPD-10262[xm]_[c][M]','CPD-10262[xb]_[c][B]','CPD-14271[xm]_[c][M]','CPD-14271[xb]_[c][B]','CPD0-2253[xm]_[c][M]','CPD0-2253[xb]_[c][B]','S_CPD-14275[xm]_[c][M]','S_CPD-14275[xb]_[c][B]', ...
    'TSA-REDUCT-RXN[M]','TSA-REDUCT-RXN[B]'});
            
forma11='SUPER-OXIDE[sm] <=> SUPER-OXIDE[cm]';
           model11 = addReaction(model11,['SUPER-OXIDE[sm]' '_[c][M]'], forma11, [], 0, -1000, 1000);
 deads6xm = model11.mets(detectDeadEnds(model11));


        %% % fixing growth in bundle sheath cell
  met_list1={}
for n=1:length(peroxM_ids)
met_list1=[met_list1,findMetsFromRxns(model11,peroxM_ids(n))'];
end
met_list1=unique(met_list1); 

 c_version=strrep(met_list1,'[xm]','[cm]');
 met_list3=intersect(model11.mets,c_version);
 met_list4={};
 for n=1:length(met_list3)
     metty=strrep(met_list3{n},'[cm]','[xm]')
     pos=find(strcmp(met_list1,metty))
 met_list4=[met_list4,met_list1(pos)]
 end

 % first=met_list(n);
 % first=strrep(first,'[sb]','[sm]');
 % shared=intersect(first,met_list1)
model12=model11;sol=[];

 for n=[3 34]
        met=split(met_list4{n},'[');
           forma=[met{1} '[cm] <=> ' met_list4{n}];
           model12 = addReaction(model12,[met_list4{n} '_[c][M]'], forma, [], 0, -1000, 1000);
           ro=optimizeCbModel(model12);
           sol=[sol,ro.f];
 end      
 sol
  n=[3 34];
  met_list4(n)
% only needed to add 1 transport rxn to fix growth 
model=model12;
%% Adding rxns that   are in multiple organelles
model1=removeRxns(model,{'dihydroorotate dehydrogenase(NAD+)[B]','dihydroorotate dehydrogenase(NAD+)[M]'})
rxns=readcell('orotate.xlsx');
pos1=find(contains(model1.rxns,'HOMOCYSMETB12-RXN'));
model1.rxns(pos1)=strrep(model1.rxns(pos1),'HOMOCYSMETB12-RXN[','HOMOCYSMETB12-RXN_1[');
pos2=find(contains(model1.rxns,'METHYLENETHFDEHYDROG-NADP-RXN['));
model1.rxns(pos2)=strrep(model1.rxns(pos2),'METHYLENETHFDEHYDROG-NADP-RXN[','METHYLENETHFDEHYDROG-NADP-RXN_1[');
pos3=find(contains(model1.rxns,'METHENYLTHFCYCLOHYDRO-RXN['));
model1.rxns(pos3)=strrep(model1.rxns(pos3),'METHENYLTHFCYCLOHYDRO-RXN[','METHENYLTHFCYCLOHYDRO-RXN_1[');
pos4=find(contains(model1.rxns,'METHYLENETHFDEHYDROG-NADP-RXN['));
model1.rxns(pos4)=strrep(model1.rxns(pos4),'METHYLENETHFDEHYDROG-NADP-RXN[','METHYLENETHFDEHYDROG-NADP-RXN_1[');
pos5=find(contains(model1.rxns,'FORMYLTHFDEFORMYL-RXN['));
model1.rxns(pos5)=strrep(model1.rxns(pos5),'FORMYLTHFDEFORMYL-RXN[','FORMYLTHFDEFORMYL-RXN_1[');
pos6=find(contains(model1.rxns,'RXN-2881['));
model1.rxns(pos6)=strrep(model1.rxns(pos6),'RXN-2881[','RXN-2881_1[');

for n = 1:length(rxns(:,1))
    % unwrap until you get a string
    form = rxns(n,2);
    while iscell(form)
        form = form{1};
    end
    
    if contains(form, '<=>')
        model1 = addReaction(model1, rxns{n,1}, form, [], 0, -1000, 1000);
    else
        model1 = addReaction(model1, rxns{n,1}, form, [], 0, 0, 1000);
    end
end
deadsnew2 = model1.mets(detectDeadEnds(model1));

deadsnew1_noc=setdiff(deadsnew2,deadsnew2(find(contains(deadsnew2,'[c'))))
xm_deads=deadsnew1_noc; 
for n=1:length(xm_deads)  
    if contains(xm_deads{n},'m]')
        met=split(xm_deads{n},'[')
           forma=[met{1} '[cm] <=> ' xm_deads{n}];
           model1 = addReaction(model1,[xm_deads{n} '_[c][M]'], forma, [], 0, -1000, 1000);
    else
                met=split(xm_deads{n},'[')
           forma=[met{1} '[cb] <=> ' xm_deads{n}];
           model1 = addReaction(model1,[xm_deads{n} '_[c][B]'], forma, [], 0, -1000, 1000);
    end
    end

  deads6xm = model1.mets(detectDeadEnds(model1));
model=model1;
save('LeafTwoCell_relocalized.mat','model')