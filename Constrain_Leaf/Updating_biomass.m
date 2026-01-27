load('Leaf_balanced_FINALSEP25.mat')
stoich=readcell('biomass_stoich')
%stoich=readcell('test_bio_stoich.xlsx')
vals=cell2mat(stoich(:,1));
mets=stoich(:,2);
metcb=strcat(mets,'[cb]')
metcb{58}='biomass[b]';
metcm=strcat(mets,'[cm]')
metcm{58}='biomass[m]';
pos_m=find(contains(metcm,'Amylose_Compound'))
metcm{pos_m}=strrep(metcm{pos_m},'[cm]','[sm]')
pos_b=find(contains(metcb,'Amylose_Compound'))
metcb{pos_b}=strrep(metcb{pos_b},'[cb]','[sb]')

model = addReaction(model, 'biomass[B]', metcb, vals, 0, 0, 1000);
model = addReaction(model, 'biomass[M]', metcm, vals, 0, 0, 1000);
save('Leaf_model_biomass','model')
