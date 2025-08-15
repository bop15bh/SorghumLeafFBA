%% pathway comparison 

clear
close all
load('Leaf_balanced_FINAL0625.mat')
model = changeRxnBounds(model,'ATR_PYRUVATE_[cb]_[cm]', 10, 'l');
model = changeRxnBounds(model,'MALATE-DEHYDROGENASE-NADP+-RXN[M]', 40, 'l');
%model = changeRxnBounds(model,'D-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN_2[M]', 0, 'l');
%model = changeRxnBounds(model,'ATR_L-ASPARTATE_[cb]_[cm]', 100, 'u');
model=removeRxns(model,{'4.1.1.32-RXN[M]','4.1.1.32-RXN[B]'})
form=' ADP[cm] + PHOSPHO-ENOL-PYRUVATE[cm] + PROTON[cm] -> ATP[cm] + PYRUVATE[cm] ';
model=addReaction(model,'PEPDEPHOS-RXN_1[M]',form,[],0,0,1000);
form=' ADP[cb] + PHOSPHO-ENOL-PYRUVATE[cb] + PROTON[cb] -> ATP[cb] + PYRUVATE[cb] ';
model=addReaction(model,'PEPDEPHOS-RXN_1[B]',form,[],0,0,1000);
changeCobraSolver('glpk');

%load('sens_controlMay25.mat')
%load('sens_droMay25.mat')
%load('sens_dro2PG.mat')
%load('sens_con2PG.mat')
%load('sens_controlNN.mat')
%load('sens_droNN.mat')
load('sens_controlJune25.mat');
load('sens_droJune25.mat');
%load('sens_control350.mat')
%load('sens_dro350.mat')
% imp_same=readtable('important_same_May.xlsx');
imp_same=readtable('important_same_June.xlsx')
 imp_same=table2cell(imp_same);
%%
model1 = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u'); 
    model1 = changeRxnBounds(model1,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -294.6979, 'l'); 
maxfl=optimizeCbModel(model1);
%pos=find(contains(model1.rxns,lis))
%maxfl=max.v(pos);
 %%
 model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -44, 'l'); 
 model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u'); 
og=optimizeCbModel(model);


%%
for n=1:length(model.subSystems)
   if iscell(model.subSystems{n})
       model.subSystems{n}=' - ';
   else
   end
end

sens_dro_nocell=erase(sens_dro,'[M]')
sens_dro_nocell=erase(sens_dro_nocell,'[B]');
sens_dro_nocell=unique(sens_dro_nocell)
sens_con_nocell=erase(sens_control,'[M]');
sens_con_nocell=erase(sens_con_nocell,'[B]');
sens_con_nocell=unique(sens_con_nocell)
dro_nocell=setdiff(sens_dro_nocell,sens_con_nocell)
con_nocell=setdiff(sens_con_nocell,sens_dro_nocell)
dro_only=setdiff(sens_dro,sens_control);
con_only=setdiff(sens_control,sens_dro);
same=intersect(sens_dro,sens_control)
same_nocell=intersect(sens_con_nocell,sens_dro_nocell);
sames=same(find(contains(same,same_nocell)))
drodro=dro_only(find(contains(dro_only,dro_nocell)))
control=con_only(find(contains(con_only,con_nocell)))
B_con = zeros(length(same_nocell), 1);
M_con = zeros(length(same_nocell), 1);
B_dro = zeros(length(same_nocell), 1);
M_dro = zeros(length(same_nocell), 1);
for n = 1:length(same_nocell)
    pos1 = find(contains(sens_control, same_nocell{n}));
    pos2 = find(contains(sens_dro, same_nocell{n}));
    
    if ~isempty(pos1)
        if sum(contains(sens_control(pos1), '[B]')) > 0
            B_con(n, 1) = 1;
        end
        if sum(contains(sens_control(pos1), '[M]')) > 0
            M_con(n, 1) = 1;
        end
    end
    
    if ~isempty(pos2)
        if sum(contains(sens_dro(pos2), '[B]')) > 0
            B_dro(n, 1) = 1;
        end
        if sum(contains(sens_dro(pos2), '[M]')) > 0
            M_dro(n, 1) = 1;
        end
    end
end
% Initialize arrays to store reactions where cell type swapped
B_swapped = {};  % Initialize as cell array for reactions where B has swapped
M_swapped = {};  % Initialize as cell array for reactions where M has swapped

for n = 1:length(same_nocell)
    % Check for reactions where B has swapped (B in control only, M in drought only)
    if B_con(n) == 1 && M_con(n) == 0 && B_dro(n) == 0 && M_dro(n) == 1
        B_swapped{end+1} = same_nocell{n};  % Append reaction to B_swapped
    end
    
    % Check for reactions where M has swapped (M in control only, B in drought only)
    if M_con(n) == 1 && B_con(n) == 0 && M_dro(n) == 0 && B_dro(n) == 1
        M_swapped{end+1} = same_nocell{n};  % Append reaction to M_swapped
    end
end



B_swapped=B_swapped';
M_swapped=M_swapped';




% data = readtable('77_drought_only_Rxns.xlsx');
% data = data(1:77, :);
% rxns = table2cell(data(:, 1));
dro_only=drodro;
 pos_dro = [];subs_dro={};
for n = 1:length(dro_only)
pos1 = find(strcmp(model.rxns, dro_only{n}));
pos_dro = [pos_dro, pos1];
suby=model.subSystems{pos1};

subs_dro=[subs_dro,model.subSystems{pos1}];

end
con_only=control;
pos_con = [];
subs_con = {};

for n = 1:length(con_only)
    pos1 = find(strcmp(model.rxns, con_only{n}));
    pos_con = [pos_con, pos1];

    suby = model.subSystems{pos1};

    % Skip if empty
    if isempty(suby)
        continue;
    end

    % Add to subs_con depending on type
    if iscell(suby)
        % Filter out empty strings inside the cell
        nonEmptySuby = suby(~cellfun('isempty', suby));
        subs_con = [subs_con, nonEmptySuby];
    elseif ischar(suby)
        % Only add if non-empty string
        if ~isempty(strtrim(suby))
            subs_con = [subs_con, {suby}];
        end
    end
end
same=sames;
pos_same = [];
subs_same = {};

for n = 1:length(same)
    pos1 = find(strcmp(model.rxns, same{n}));
    
    if isempty(pos1)
        disp(['No match at index: ', num2str(n)]);
        continue;
    end

    pos_same = [pos_same, pos1];
    suby = model.subSystems{pos1};

    % Skip if empty
    if isempty(suby)
        continue;
    end

    % Add to subs_same depending on type
    if iscell(suby)
        % Filter out empty strings
        nonEmptySuby = suby(~cellfun('isempty', suby));
        subs_same = [subs_same, nonEmptySuby];
    elseif ischar(suby)
        if ~isempty(strtrim(suby))
            subs_same = [subs_same, {suby}];
        end
    end
end



names_con = model.rxnNames(pos_con);
names_dro=model.rxnNames(pos_dro);
names_same=model.rxnNames(pos_same);
%subs_con = model.subSystems(pos_con);
%subs_dro = model.subSystems(pos_dro);
%subs_same = model.subSystems(pos_same);

dro_subs_ng = setdiff(subs_dro, subs_dro(find(contains(subs_dro, ' - '))));
con_subs_ng = setdiff(subs_con, subs_con(find(contains(subs_con, ' - '))));
same_subs_ng = setdiff(subs_same, subs_same(find(contains(subs_same, ' - '))));

% names = table2cell(data(:, 1));
% atr = setdiff(names, names(find(contains(names, 'ATR_'))));
%% load bulk rna seq
%data=readcell('Corrected.Data.from.v3.1.1.Genome.Counts_BLH.csv')
data2=readcell('DEG_results_Stata.csv')
%titles=data(1,:);
titles=data2(1,:);
%dro=titles(find(contains(titles,'Drought')));
%con=titles(find(contains(titles,'Control')));
load('all_paths.mat')

to_addpaths={'D-serine metabolism','glycerol-3-phosphate shuttle','glycerophosphodiester degradation','serine racemization','D-serine degradation','L-aspartate biosynthesis','L-aspartate degradation I','pyrimidine nucleobases salvage I','L-glutamine biosynthesis I','ammonia assimilation cycle II'};
all_paths=vertcat(all_paths,to_addpaths');

%all_paths = removeDuplicatePathwaysByRxns(model, all_paths,model.rxns);


stata1_day10_FC=(find(contains(titles,'Stata1_Drought_Day10_vs_Control_Day10_log2FC')));
stata1_day10_p=(find(contains(titles,'Stata1_Drought_Day10_vs_Control_Day10_adjP')));
%dro_10=dro(find(contains(dro,'Day10')));
%con_10=con(find(contains(con,'Day10')));
%stat_con=find(contains(con_10,'Stata'));
%stat_con=find(contains(titles,con_10(stat_con)));
% 
% stat_dro=find(contains(dro_10,'Stata'));
% stat_dro=find(contains(titles,dro_10(stat_dro)));
% 
% cheng_con=find(contains(con_10,'Cheng'));
% cheng_con=find(contains(titles,con_10(cheng_con)));
% 
% cheng_dro=find(contains(dro_10,'Cheng'));
% cheng_dro=find(contains(titles,dro_10(cheng_dro)));


%% drought only 
% all_paths=vertcat(paths_dro,paths_con,paths_same);
% all_paths=unique(all_paths);
total_dro_path=[]; pathway={};
for n=1:length(all_paths)
pos=find(contains(model.subSystems,all_paths{n}));
rxns=model.rxns(pos);
epsilon = 1e-6;  % Small value to avoid division by zero
model_Fc = log2((mean(abs(og.v(pos))) + epsilon) / (mean(abs(maxfl.v(pos))) + epsilon));
     [geneList] = findGenesFromRxns(model, rxns);
     genies={};
for j=1:length(geneList)
genies=[genies;geneList{j}]
end
genies=unique(genies)

% control=[];drought=[];
% for k=1:length(genies)
%     pos1=find(contains(data(:,1),genies{k}));
%     dro_pos1=data(pos1,dro)
%     con_pos1=data(pos1,con)
%     mean_dro=mean(cell2mat(dro_pos1))
%     mean_con=mean(cell2mat(con_pos1))
%     control=[control,mean_con];
%     drought=[drought,mean_dro];
% end
% control=control(~isnan(control));
% drought=drought(~isnan(drought));
% total_mcon=mean(control);
% total_mdro=mean(drought);
% FC=log2(total_mdro/total_mcon)
stata_fc = [];
%cheng_fc = [];
stata_pvals = [];
%cheng_pvals = [];

for k = 1:length(genies)
    pos1 = find(contains(data2(:,1), genies{k}));
    if isempty(pos1)
        continue;
    end

    % Extract Day 10 expression values for both genotypes
    % stata_ctrl_vals = cell2mat(data(pos1, stat_con));   % Stata control
    % stata_dro_vals  = cell2mat(data(pos1, stat_dro));   % Stata drought
    % cheng_ctrl_vals = cell2mat(data(pos1, cheng_con));  % Cheng control
    % cheng_dro_vals  = cell2mat(data(pos1, cheng_dro));  % Cheng drought
    % 
    % if isempty(stata_ctrl_vals) || isempty(stata_dro_vals) || ...
    %    isempty(cheng_ctrl_vals) || isempty(cheng_dro_vals)
    %     continue;
    % end
    % 
    % % Log2 Fold Change for each experiment
    % stata_fc(end+1) = log2(mean(stata_dro_vals + 1e-6) / mean(stata_ctrl_vals + 1e-6));
    % cheng_fc(end+1) = log2(mean(cheng_dro_vals + 1e-6) / mean(cheng_ctrl_vals + 1e-6));

    % % T-tests for each genotype
    % [~, p1] = ttest2(stata_ctrl_vals, stata_dro_vals);
    % [~, p2] = ttest2(cheng_ctrl_vals, cheng_dro_vals);
   val1 = data2{pos1, stata1_day10_p};
val2 = data2{pos1, stata1_day10_FC};

if (ischar(val1) && contains(val1, 'NA')) || (ischar(val2) && contains(val2, 'NA'))
    continue;
end

    stata_pvals(end+1) =data2{pos1,stata1_day10_p} ;
     stata_fc(end+1)=data2{pos1, stata1_day10_FC};
   % cheng_pvals(end+1) = p2;
end

% Combine results
%all_fc = [stata_fc, cheng_fc];
%all_pvals = [stata_pvals, cheng_pvals];

% Final pathway-level stats
%mean_log2fc = mean(all_fc, 'omitnan');
%mean_pval = mean(all_pvals, 'omitnan');
mean_log2fc = mean(stata_fc, 'omitnan');
%mean_pval = mean(stata_pvals, 'omitnan');
% Number of p-values
k = length(stata_pvals);

% Fisher's statistic
X2 = -2 * sum(log(stata_pvals));

% Degrees of freedom = 2 * k
df = 2 * k;

% Combined p-value using chi-squared distribution
combined_p = 1 - chi2cdf(X2, df);
mean_pval = combined_p;


% Loop through each ratio and perform a z-test
if mean(abs(og.v(pos)))<1 && mean(abs(maxfl.v(pos)))<1 &&  mean(abs(og.v(pos))) > 0 && mean(abs(maxfl.v(pos))) > 0
count_dro = round(mean(abs(og.v(pos)))*100);
count_og = round(mean(abs(maxfl.v(pos)))*100);
tab = [count_og, 100 - count_og;
count_dro, 100 - count_dro];
[h,p] = fishertest(tab);
p_values_mod = p; % Right-tailed test
elseif xor(mean(abs(og.v(pos))) == 0, mean(abs(maxfl.v(pos))) == 0)
    epsilon = 1;  % fallback count when flux is zero (represents 1%)

flux_dro = mean(abs(og.v(pos)));
flux_og  = mean(abs(maxfl.v(pos)));

% Convert to % scale for contingency table
count_dro = round(flux_dro * 100);
count_og  = round(flux_og * 100);

% Handle zero cases with epsilon
if count_dro == 0
    count_dro = epsilon;
end
if count_og == 0
    count_og = epsilon;
end

% Ensure values do not exceed 100
count_dro = min(count_dro, 100);
count_og  = min(count_og, 100);

% Build contingency table (active vs inactive flux)
tab = [
    count_og,   100 - count_og;
    count_dro,  100 - count_dro
];

% Run Fisher's exact test
[h, p] = fishertest(tab);
                p_values_mod = p; % Right-tailed test

elseif mean(abs(og.v(pos)))==0 && mean(abs(maxfl.v(pos)))==0
            p_values_mod = NaN; % Right-tailed test
else
% Get mean fluxes
flux_control = mean(abs(maxfl.v(pos)));
flux_drought = mean(abs(og.v(pos)));

% Scale to 0–100 based on the maximum flux
max_flux = max([flux_control, flux_drought]);

count_control = round((flux_control / max_flux) * 100);
count_drought = round((flux_drought / max_flux) * 100);

% Build contingency table
tab = [
    count_control, 100 - count_control;
    count_drought, 100 - count_drought
];

% Run Fisher's exact test
[h, p] = fishertest(tab);
p_values_mod = p;

end
% Save pathway-level results
pathway{n,1} = all_paths{n};
pathway{n,2} = mean_log2fc;
pathway{n,3} = mean_pval; 
pathway{n,4} = model_Fc;
pathway{n,5} = p_values_mod;
end

%% FDR Adjustment for Pathway-level Model P-values (ADD HERE)
fprintf('\n=== Applying FDR Correction to Pathway Model P-values ===\n');

% Extract all pathway p-values
pathway_p_values = [];
valid_pathway_indices = [];
for n = 1:size(pathway,1)
    if ~isempty(pathway{n,5}) && ~isnan(pathway{n,5})
        pathway_p_values(end+1) = pathway{n,5};
        valid_pathway_indices(end+1) = n;
    end
end
% Apply FDR correction
if ~isempty(pathway_p_values)
    % Sort p-values and keep track of original positions
    [sorted_p, sort_idx] = sort(pathway_p_values);
    n_tests = length(sorted_p);
    
    % Calculate FDR-adjusted p-values
    adjusted_p = zeros(size(sorted_p));
    for i = 1:n_tests
        adjusted_p(i) = sorted_p(i) * n_tests / i;
    end
    
    % Ensure monotonicity (each p-value >= previous)
    for i = n_tests-1:-1:1
        if adjusted_p(i) > adjusted_p(i+1)
            adjusted_p(i) = adjusted_p(i+1);
        end
    end
    
    % Cap at 1
    adjusted_p(adjusted_p > 1) = 1;
    
    % Map back to original order
    final_adjusted_p = zeros(size(pathway_p_values));
    final_adjusted_p(sort_idx) = adjusted_p;
    
    % Store adjusted p-values in pathway array (new column 6)
    for i = 1:length(valid_pathway_indices)
        pathway{valid_pathway_indices(i),6} = final_adjusted_p(i);
    end
    
    fprintf('Adjusted %d pathway p-values using FDR\n', length(pathway_p_values));
    fprintf('Before FDR: %d significant (p < 0.05)\n', sum(pathway_p_values < 0.05));
    fprintf('After FDR: %d significant (adjusted p < 0.05)\n', sum(final_adjusted_p < 0.05));
end

%% drought only 
paths_dro = {};
for n = 1:length(dro_subs_ng)
pos1 = contains(dro_subs_ng(n), '//');
if ~isempty(pos1)
spit = split(dro_subs_ng{n}, '// ');
for j = 1:length(spit)
paths_dro = [paths_dro, spit{j}];
end
else
paths_dro = [paths_dro, dro_subs_ng{n}];
end
end
paths_dro = unique(paths_dro)';
paths_dro=strtrim(paths_dro);
paths_dro=unique(paths_dro)



count = [];imp_count=[];perc=[];
for n = 1:length(paths_dro)
pos1 = find(contains(subs_dro, paths_dro{n}));
pos2=find(contains(model.subSystems,paths_dro{n}));
total=length(pos2)
count = [count, length(pos1)];
perc=[perc,(length(pos1)/total)*100]
name=pos_dro(pos1);
%match=intersect(imp_dro,model.rxns(name))
%imp_count=[imp_count,(length(match)/total)*100]
end
dro_exp=[];dro_modelfc=[];dro_pexp=[];dro_pmod=[];
for n=1:length(paths_dro)
    pos=find(strcmp(pathway(:,1),paths_dro{n}));
    exp=pathway(pos,2);
    fc=pathway(pos,4);
    pexp=pathway(pos,3);
   % pmod=pathway(pos,5)
   pmod=pathway(pos,6);
    if isempty(exp) 
    dro_exp=[dro_exp,0];
        dro_pexp=[dro_pexp,NaN];
        dro_pmod=[dro_pmod,NaN];
        dro_modelfc=[dro_modelfc,NaN];

    elseif isempty(pathway{pos,6})
                dro_pmod=[dro_pmod,NaN];
                dro_exp=[dro_exp,exp]  ;
                dro_pexp=[dro_pexp,pexp];
                dro_modelfc=[dro_modelfc,fc];

    else
      dro_exp=[dro_exp,exp]  ;
      dro_pexp=[dro_pexp,pexp];
      dro_pmod=[dro_pmod,pmod];
          dro_modelfc=[dro_modelfc,fc];

    end
end
droTAB = cell2table(horzcat(paths_dro, num2cell(perc'),dro_exp',dro_pexp',dro_modelfc',dro_pmod'));
writetable(droTAB, 'Drought_only_pathway_analysis.txt', 'Delimiter', '\t') % tab-delimited
%droTAB = cell2table(horzcat(paths_dro, num2cell(perc'),num2cell(imp_count'),dro_exp'));
%% control only

paths_con = {};
for n = 1:length(con_subs_ng)
pos1 = contains(con_subs_ng(n), '//');
if ~isempty(pos1)
spit = split(con_subs_ng{n}, '// ');
for j = 1:length(spit)
paths_con = [paths_con, spit{j}];
end
else
paths_con = [paths_con, con_subs_ng{n}];
end
end
paths_con = unique(paths_con)';
paths_con=strtrim(paths_con);
paths_con=unique(paths_con)




count = [];imp_count=[];perc=[];
for n = 1:length(paths_con)
pos1 = find(contains(subs_con, paths_con{n}));
pos2=find(contains(model.subSystems,paths_con{n}));

count = [count, length(pos1)];
total=length(pos2)
perc=[perc,(length(pos1)/total)*100]

name=pos_con(pos1);
%match=intersect(imp_con,model.rxns(name))
%imp_count=[imp_count,(length(match)/total)*100]
end
con_exp=[];con_modelfc=[];con_pexp=[];con_pmod=[];
for n=1:length(paths_con)
    pos=find(strcmp(pathway(:,1),paths_con{n}));
    exp=pathway(pos,2);
    fc=pathway(pos,4)
    pexp=pathway(pos,3);
    %pmod=pathway(pos,5);
    pmod=pathway(pos,6);
    if isempty(exp)
    con_exp=[con_exp,0]
    con_pmod=[con_pmod,NaN]
        con_modelfc=[con_modelfc,NaN];
    con_pexp=[con_pexp,NaN]
     elseif isempty(pathway{pos,6})
                dro_pmod=[dro_pmod,NaN];
                dro_exp=[dro_exp,exp]  ;
                dro_pexp=[dro_pexp,pexp];
                dro_modelfc=[dro_modelfc,fc];
    else
      con_exp=[con_exp,exp] 
      con_pmod=[con_pmod,pmod]
      con_pexp=[con_pexp,pexp]
          con_modelfc=[con_modelfc,fc];

    end
end
conTAB = cell2table(horzcat(paths_con, num2cell(perc'),con_exp',con_pexp',con_modelfc',con_pmod'));
writetable(conTAB, 'Control_only_pathway_analysis.txt', 'Delimiter', '\t') % tab-delimited

%conTAB = cell2table(horzcat(paths_con, num2cell(perc'),num2cell(imp_count'),con_exp'));

%% shared 
paths_same = {};
for n = 1:length(same_subs_ng)
pos1 = contains(same_subs_ng(n), '//');
if ~isempty(pos1)
spit = split(same_subs_ng{n}, '// ');
for j = 1:length(spit)
paths_same = [paths_same, spit{j}];
end
else
paths_same = [paths_same, same_subs_ng{n}];
end
end
paths_same = unique(paths_same)';
paths_same=strtrim(paths_same);
paths_same=unique(paths_same)


count = [];imp_count=[];perc=[];imp_freq=[];
for n = 1:length(paths_same)
pos1 = find(contains(subs_same, paths_same{n}));
pos2=find(contains(model.subSystems,paths_same{n}));

total=length(pos2)
perc=[perc,(length(pos1)/total)*100]

count = [count, length(pos1)];
name=pos_same(pos1);
match=intersect(imp_same,model.rxns(name))
imp_count=[imp_count,(length(match)/total)*100]
imp_freq=[imp_freq,length(match)]

end
% adding pathway level expression to tables
same_exp=[];    same_modelfc=[];same_pexp=[];same_pmod=[];
for n=1:length(paths_same)
    pos=find(strcmp(pathway(:,1),paths_same{n}));
    exp=pathway(pos,2);
    fc=pathway(pos,4);
    pexp=pathway(pos,3);
    %pmod=pathway(pos,5);
    pmod=pathway(pos,6);
    if isempty(exp)
    same_exp=[same_exp,0]
    same_pexp=[same_pexp,NaN]
        same_modelfc=[same_modelfc,NaN];
        same_pmod=[same_pmod,NaN]
     elseif isempty(pathway{pos,6})
                dro_pmod=[dro_pmod,NaN];
                dro_exp=[dro_exp,exp]  ;
                dro_pexp=[dro_pexp,pexp];
                dro_modelfc=[dro_modelfc,fc];
    else
      same_exp=[same_exp,exp]  
      same_pmod=[same_pmod,pmod]
      same_pexp=[same_pexp,pexp]
          same_modelfc=[same_modelfc,fc];

    end
end
sameTAB = cell2table(horzcat(paths_same, num2cell(perc'),same_exp',same_pexp',same_modelfc',same_pmod',num2cell(imp_count'),num2cell(imp_freq')));
writetable(sameTAB, 'Shared_pathway_analysis.txt', 'Delimiter', '\t') % tab-delimited

%sameTAB = cell2table(horzcat(paths_same, num2cell(perc'),num2cell(imp_count'),same_exp'));
