%% unblocking mets 
load("SC_constrainedmodelFINAL.mat");
blocked1=model2.mets(detectDeadEnds(model2))
len=[];pepc=[];
for n=1:length(blocked1)
    pos=find(contains(model2.rxns,'PEPCARBOX-RXN[M]'));
    ro=optimizeCbModel(model2);
    rxns=findRxnsFromMets(model2,blocked1(n));
    pepc=[pepc,ro.v(pos)]
    if contains(rxns,'ATR_')
    len=[len,n]
    else
    end
end
mets=blocked1(len)
rxns=findRxnsFromMets(model2,mets);
model=removeRxns(model2,rxns)
blocked=model.mets(detectDeadEnds(model))

 % removed 6 transport reactions which were obsolete, reducing blocked mets
 % from 81 to 76
 deads=[];    model1=model;
left={};added={}; subcell={};
to_avoid={'D-RIBULOSE-15-P2','UBIQUINONE-9','WATER[tb]'};
match=blocked(find(contains(blocked,to_avoid)));
blocked=setdiff(blocked,match)
for n=1:length(blocked)
     if contains(blocked{n},'[cb]') || contains(blocked{n},'[cm]')
    metminus=erase(blocked{n},'[cb]');
    metminus=erase(metminus,'[cm]');
    form=[metminus '[cb] <=> ' metminus '[cm]'];
    model1=addReaction(model1,['ATR_' metminus '_[cb]_[cm]'],form,[],0,-10000,10000);
        deadnow=detectDeadEnds(model1);
        deads=[deads,length(deadnow)]
        added=[added,['ATR_' metminus '_[cb]_[cm]']];
        subcell=[subcell,blocked{n}]
     elseif contains(blocked{n},'[mm]') || contains(blocked{n},'[sm]') || contains(blocked{n},'[im]')
         metminus=strsplit(blocked{n},'[');    
         form=[blocked{n} ' <=> ' metminus{1} '[cm]'];
    model1=addReaction(model1,['ATR_' blocked{n} '_[cm]'],form,[],0,-10000,10000);
deadnow=detectDeadEnds(model1);
        deads=[deads,length(deadnow)]
                added=[added,['ATR_' blocked{n} '_[cb]_[cm]']];

     elseif contains(blocked{n},'[sb]') || contains(blocked{n},'[mb]') || contains(blocked{n},'[tb]') || contains(blocked{n},'[ib]')
      metminus=strsplit(blocked{n},'[');    
         form=[blocked{n} ' <=> ' metminus{1} '[cb]'];
    model1=addReaction(model1,['ATR_' blocked{n} '_[cb]'],form,[],0,-10000,10000);
            deads=[deads,length(deadnow)]

            added=[added,['ATR_' blocked{n} '_[cb]_[cm]']];
     else
         left=[left,blocked(n)]
     end
end
deadnow=detectDeadEnds(model1);
% 3 blocked mets left
finalmets=model1.mets(deadnow);
intersect(blocked,finalmets)
% 3 new blocked mets so clearly adding transport reaction for these 3 is
% not a fix

% for D-RIBULOSE-15-P2[sb]  i need sb-cb cb-cm cm-sm
% for ubiquinone-9 i need mb-cb cb-ib and same for m
    %         form1='D-RIBULOSE-15-P2[sm] <=> D-RIBULOSE-15-P2[cm]';
    %         form2='D-RIBULOSE-15-P2[cm] <=> D-RIBULOSE-15-P2[cb]';
    %         form3='D-RIBULOSE-15-P2[sb] <=> D-RIBULOSE-15-P2[cb]';
    % 
    % model1=addReaction(model1,'ATR_D-RIBULOSE-15-P2_[cm]_[sm]',form1,[],0,-10000,10000);
    %             model1=addReaction(model1,'ATR_D-RIBULOSE-15-P2_[cm]_[cb]',form2,[],0,-10000,10000);
    % model1=addReaction(model1,'ATR_D-RIBULOSE-15-P2_[cb]_[sb]',form3,[],0,-10000,10000);
% that fixed that met
form4='UBIQUINONE-9[ib] <=> UBIQUINONE-9[cb]';
form5='UBIQUINONE-9[mb] <=> UBIQUINONE-9[cb]';
form6='UBIQUINONE-9[im] <=> UBIQUINONE-9[cm]';
form7='UBIQUINONE-9[mm] <=> UBIQUINONE-9[cm]';

    model1=addReaction(model1,'ATR_D-UBIQUINONE-9_[cb]_[ib]',form4,[],0,-10000,10000);
                model1=addReaction(model1,'ATR_D-UBIQUINONE-9_[cb]_[mb]',form5,[],0,-10000,10000);
   model1=addReaction(model1,'ATR_D-UBIQUINONE-9_[cm]_[im]',form6,[],0,-10000,10000);
                model1=addReaction(model1,'ATR_D-UBIQUINONE-9_[cm]_[mm]',form7,[],0,-10000,10000);
 
                
            %    form8='WATER[tb] <=> WATER[sb]';

  deadnow=detectDeadEnds(model1);

  %% removing two reactions to allow PEPC to function normally 
  % don't need transport of rubp now since it shouldnt have had rubisco in
  % M, now its just   #1579  RXN-18031[B], Bd: -1000 / 1000, , grRules: 
%H2CO3[cb]   <=>   HCO3[cb] + PROTON[cb] 
% RXN0-5224[B]
% its because not_list says bundlesheath shouldnt have CA and it was
% removed from BS


  ro=optimizeCbModel(model1);
      pos=find(contains(model1.rxns,'PEPCARBOX-RXN[M]'));
ro.v(pos)
rxns1={'ATR_H2CO3_[cb]_[cm]' ,'ATR_D-RIBULOSE-15-P2_[cm]_[sm]'}
model1=removeRxns(model1,rxns1)
model1=removeRxns(model1,{'RXN-18031[B]'})
  ro=optimizeCbModel(model1);
      pos=find(contains(model1.rxns,'PEPCARBOX-RXN[M]'));
      ro.f*24/1000
ro.v(pos)
% this works but now 3 mets are blocked
save('SC_constrained_unblocked.mat','model1');
