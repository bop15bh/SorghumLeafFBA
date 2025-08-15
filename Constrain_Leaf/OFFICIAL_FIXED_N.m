
clear 

%load('constrained_leaf_noredundancies.mat')
load('35deads_withsubsystems.mat')
%load('Constrained_unblocked_MAY25.mat')
deads=model.mets(detectDeadEnds(model))
model1=model;
model1=removeRxns(model1,{'RXN-16398[B]','RXN-16398[M]'});
deads1=model1.mets(detectDeadEnds(model1))

form1='ALPHA-METHYL-5-ALPHA-ERGOSTA[cm] <=> ALPHA-METHYL-5-ALPHA-ERGOSTA[cb]';
model1 = addReaction(model1,'ATR_ALPHA-METHYL-5-ALPHA-ERGOSTA_[cb]_[cm]',form1,[],0,-1000,1000);
form2='CHORISMATE[cm] <=> CHORISMATE[sm]';
model1 = addReaction(model1,'ATR_CHORISMATE_[c]_[s][M]',form2,[],0,-1000,1000);
form3='CPD-2205[cm] <=> CPD-2205[cb]';
model1 = addReaction(model1,'ATR_CPD-2205_[cb]_[cm]',form3,[],0,-1000,1000);
form4='CPD-4081[cm] <=> CPD-4081[cb]';
model1 = addReaction(model1,'ATR_CPD-4081_[cb]_[cm]',form4,[],0,-1000,1000);
form5='MALONATE-S-ALD[cm] <=> MALONATE-S-ALD[cb]';
model1 = addReaction(model1,'ATR_MALONATE-S-ALD_[cb]_[cm]',form5,[],0,-1000,1000);
form6='O-SUCCINYLBENZOATE[cm] <=> O-SUCCINYLBENZOATE[cb]';
model1 = addReaction(model1,'ATR_O-SUCCINYLBENZOATE_[cb]_[cm]',form6,[],0,-1000,1000);
form7='URACIL[cm] <=> URACIL[cb]';
model1 = addReaction(model1,'ATR_URACIL_[cb]_[cm]',form7,[],0,-1000,1000);

model1=removeRxns(model1,{'RXN-10972[B]','2.7.4.24-RXN[B]','RXN-14995[B]','RXN-15386[M]','RXN-15385[M]','RXN-7570[M]','L-GULONOLACTONE-OXIDASE-RXN[M]' ...
    'RXN-7570[B]','RXN-9310[M]','1TRANSKETO-RXN[B]','1TRANSKETO-RXN[M]','SORBITOL-6-PHOSPHATASE-RXN[B]','1.1.99.28-RXN[B]','1.1.99.28-RXN[M]' ...
    'RXN-14515[B]','RXN-14515[M]','GLUCONOLACT-RXN[B]','GLUCONOLACT-RXN[M]','GLUCONOKIN-RXN[B]','GLUCONOKIN-RXN[M]','SERINE-C-PALMITOYLTRANSFERASE-RXN[B]','RXN66-221[B]' ...
    'RXN-8618[M]','4.1.1.80-RXN[B]','RXN-13689[M]','ATR_Light_Compound_PLASTID-STR_THY-LUM[B]','ATR_OXYGEN-MOLECULE_PLASTID-STR_THY-LUM[B]','ATR_WATER_PLASTID-STR_THY-LUM[B]','RXN0-5398_1[M]'})
deads2=model1.mets(detectDeadEnds(model1))
extra=setdiff(deads2,deads1)

%% fixing Red-NADPH-Hemoprotein-Reductases_Compound

[rxns,form]=findRxnsFromMets(model1,{'Red-NADPH-Hemoprotein-Reductases_Compound[cb]','Red-NADPH-Hemoprotein-Reductases_Compound[cm]'});
model2=model1;
for n=1:length(rxns)
    formula1=strrep(form{n},'Red-NADPH-Hemoprotein-Reductases_Compound[cb]','NADPH[cb]')
    formula1=strrep(formula1,'Ox-NADPH-Hemoprotein-Reductases_Compound[cb]','NADP[cb]')
        formula1=strrep(formula1,'Ox-NADPH-Hemoprotein-Reductases_Compound[cm]','NADP[cm]')
    formula1=strrep(formula1,'Red-NADPH-Hemoprotein-Reductases_Compound[cm]','NADPH[cm]')

 if contains(formula1,'<=>')
     model2= addReaction(model2,rxns{n},formula1,[],0,-1000,1000);

 else
    model2= addReaction(model2,rxns{n},formula1,[],0,0,1000);
 end
end 
   
model2=removeMetabolites(model2,{'Red-NADPH-Hemoprotein-Reductases_Compound[cb]','Red-NADPH-Hemoprotein-Reductases_Compound[cm]','Ox-NADPH-Hemoprotein-Reductases_Compound[cb]','Ox-NADPH-Hemoprotein-Reductases_Compound[cm]'});
%%
model=model2;
% results = verifyModel(model,'massBalance', true)
% imbalanced=model.rxns(results.massBalance.imBalancedRxnBool)
ex=model.rxns(find(contains(model.rxns,'EX_')));
% to_fix=setdiff(imbalanced,ex)
% atr=model.rxns(find(contains(model.rxns,'ATR_')));
% missing=[];
% for n=1:length(model.metFormulas)
%     if isnan(model.metCharges{n})
% missing=[missing,n]
%     else
%     end
% end

pos_mal=find(contains(model.mets,'MAL[sm]'))
model.metFormulas(pos_mal)={'C4H4O5'};
model.metCharges(pos_mal)=-2;
pos_oxal=find(contains(model.mets,'OXALACETIC_ACID[sm]'))
model.metFormulas(pos_oxal)={'C4H2O5'};
model.metCharges(pos_oxal)=-2;
pos_112=find(contains(model.mets,'CPD-11232[cb]'))
model.metFormulas(pos_112)={'C18H29NO4'};
model.metCharges(pos_112)=0;
pos_112m=find(contains(model.mets,'CPD-11232[cm]'))
model.metFormulas(pos_112m)={'C18H29NO4'};
model.metCharges(pos_112m)=0;

new=model;

%%
new = changeRxnBounds(new,'EX_Light_Compound_EXTRACELLULAR', 3000, 'u');
new=removeRxns(new,{'ATR_CARBON-DIOXIDE_CYTOSOL_PLASTID-STR[B]'})
new = changeRxnBounds(new,'RXN-15479_1[M]', 2000, 'u');
new = changeRxnBounds(new,'ATPSYN-RXN-Plastid[M]', 2000, 'u');
light=find(contains(new.rxns,'ATR_Light'))
new.lb(light)=-5000;
new.ub(light)=5000;
%new=removeRxns(new,{'EX_CARBON-DIOXIDE_EXTRACELLULAR[M]','EX_CARBON-DIOXIDE_EXTRACELLULAR[B]'})
%new=removeRxns(new,{'RXN-14932_2[M]','MALATE-DEH-RXN_1[M]','ATR_PLASTOQUINONE-9_[cb]_[cm]'})

ex=new.rxns(find(contains(new.rxns,'EX_')));
pos=find(contains(new.rxns,'EX_'))
to_rem=ex(find(contains(ex,'[')));
new=removeRxns(new,to_rem)

new=removeRxns(new,{'TRANS-RXN-395[M]','TRANS-RXN-395[B]','ATR_CARBON-DIOXIDE_CYTOSOL_EXTRACELLULAR[B]'})

 
new=removeRxns(new, {'RXN-11334_1[B]','RXN-11334_1[M]','RXN-11334_3[B]','RXN-11334_3[M]','TRANS-RXN-194[B]','TRANS-RXN-194[M]','TRANS-RXN-206[B]','TRANS-RXN-206[M]','TRANS-RXN-213[B]','TRANS-RXN-213[M]','OXALOACETASE-RXN[B]','OXALOACETASE-RXN[M]'})



form='CPD-4205[cb] + WATER[cb]   ->   CPD-15317[cb] + CPD-4209[cb]';
new=addReaction(new,'RXN-4313_2[B]',form,[],0,0,1000);

form1='CPD-4206[cb] + WATER[cb]   <=>   CPD-15317[cb] + CPD-4210[cb] ';
new=addReaction(new,'RXN-4314_2[B]',form1,[],0,-1000,1000);

form2='CPD-4205[cm] + WATER[cm]   ->   CPD-15317[cm] + CPD-4209[cm]'
new=addReaction(new,'RXN-4313_2[M]',form2,[],0,0,1000);
form3='CPD-4206[cm] + WATER[cm]   <=>   CPD-15317[cm] + CPD-4210[cm]';
new=addReaction(new,'RXN-4314_2[M]',form3,[],0,-1000,1000);
form4='ADENOSINE_DIPHOSPHATE_RIBOSE[cm] + WATER[cm]   <=>   AMP[cm]   + CPD-15317[cm] + 2 PROTON[cm] ';
new=addReaction(new,'RXN0-1441_2[M]',form4,[],0,-1000,1000);
form5='CPD-15317[cm]   <=>   CPD-15895[cm]';
new=addReaction(new,'RXN-15346[M]',form5,[],0,-1000,1000);



form6='MAL[cm]   <=>   MAL[sm]';
new=addReaction(new,'ATR_MAL_[sm]_[cm]',form6,[],0,-1000,1000);
form8=' OXALACETIC_ACID[cm]   <=>    OXALACETIC_ACID[sm]';
new=addReaction(new,'ATR_OXALACETIC_ACID_[sm]_[cm]',form8,[],0,-1000,1000);

form7='MAL[sm] + NADP[sm]   <=>   NADPH[sm] + OXALACETIC_ACID[sm] + PROTON[sm]'
new=addReaction(new,'MALATE-DEHYDROGENASE-NADP+-RXN[M]',form7,[],0,-1000,1000);
new=removeRxns(new, {'MALATE-DEH-RXN_1[M]','MALATE-DEH-RXN_1[B]'})
po=optimizeCbModel(new);
%%
new=removeRxns(new,{'ATR_BUTANAL_[cb]_[cm]','ATR_BUTANOL_[cb]_[cm]','ATR_GAP_[cb]_[cm]','ATR_SUCROSE_[cb]_[cm]','ATR_CAFFEOYLQUINATE_[cb]_[cm]'})

%% adding the extra rxns to old model to find the source of N issue
extra={'NADNUCLEOSID-RXN_1[B]','NADNUCLEOSID-RXN_1[M]','NICONUCADENYLYLTRAN-RXN[M]','NICOTINAMID-RXN[B]','NICOTINAMID-RXN[M]','NITRATE-REDUCTASE-NADH-RXN[B]','RXN0-5289[B]','RXN0-5289[M]','TSA-REDUCT-RXN[M]'};
%extra=setdiff(new.rxns,model.rxns)

og=optimizeCbModel(new);
nit=find(contains(new.rxns,'EX_NITRATE_EXTRACELLULAR'));

%%


%%
new2=new;
sol=[];rgr=[];noo=[];
for n=1:2%length(extra1)
   new2=removeRxns(new2,extra(n))
ro=optimizeCbModel(new2);
    nit=find(contains(new2.rxns,'EX_NITRATE_EXTRACELLULAR'));
   % if ro.v(nit)~=0 
sol=[sol,ro.v(nit)]
rgr=[rgr,ro.f*24/1000]
noo=[noo,n];
end

ro=optimizeCbModel(new2)
co2=find(contains(new2.rxns,'EX_CARBON-DIOXIDE_EXTRACELLULAR'));
CA1=find(contains(new2.rxns,'CARBODEHYDRAT-RXN'))
pepc=find(contains(new2.rxns,'PEPCARBOX'))
mal=find(contains(new2.rxns,'ATR_MAL_[cb]_[cm]'))
rub=find(contains(new2.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
pyr=find(contains(new2.rxns,'ATR_PYRUVATE_[cb]_[cm]'))
    nit=find(contains(new2.rxns,'EX_NITRATE_EXTRACELLULAR'));
model=new2;
save('Constrained_unblocked_Leaf_FINAL0525.mat','model')