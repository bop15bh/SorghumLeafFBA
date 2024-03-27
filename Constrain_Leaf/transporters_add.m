%% adding transport rxns for blocked mets
model=final;
form=' LINOLEIC_ACID[cb] <=>  LINOLEIC_ACID[cm] ';
 model = addReaction(model,'ATR_LINOLEIC_ACID[cm]',form,[],0,-1000,1000);
form=' DUMP[cm] <=>  DUMP[cb] ';
 model = addReaction(model,'ATR_DUMP[cm]',form,[],0,-1000,1000);
form=' COUMARALDEHYDE[cm] <=>  COUMARALDEHYDE[cb] ';
 model = addReaction(model,'ATR_COUMARALDEHYDE[cm]',form,[],0,-1000,1000);
form=' TYR[cm] <=>  TYR[cb] ';
 model = addReaction(model,'ATR_TYR[cm]',form,[],0,-1000,1000);
form=' COUMARATE[cb] <=>  COUMARATE[cm] ';
 model = addReaction(model,'ATR_COUMARATE[cb]',form,[],0,-1000,1000);
final=model;
