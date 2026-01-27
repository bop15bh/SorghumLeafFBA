%% Model validation
% comparing if gene expression per reaction matches model flux change
clear 
close ALL

changeCobraSolver('glpk');
%load('Leaf_balanced_FINALSEP25.mat')
 load('Leaf_model_biomass.mat')
 trans=readcell('old_transport.xlsx');
 trans=trans(:,1);% nbefore 19 was close 26 closest, 28,29,30,35, 38, 41,42
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
%% Setting limit to Rubisco
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

  %form=' 5-HYDROXY-CONIFERALDEHYDE[cb]  <=> 5-HYDROXY-CONIFERALDEHYDE[cm]  ';
 %model=addReaction(model,'ATR_5-HYDROXY-CONIFERALDEHYDE[cb]',form,[],0,-1000,1000);

%load('sens_droSept25.mat')
%load('sens_controlSept25.mat')
load('sens_droNov25.mat')
load('sens_controlNov25.mat')
%load('sens_droJune25.mat')
%load('sens_controlJune25.mat')
%load('sens_droNN.mat')
%load('sens_controlNN.mat')
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

%% dro only 
%model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -44.6512, 'l'); 
model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -0.847, 'l'); 

model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u'); 
og=optimizeCbModel(model);
con_version = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -5.587, 'l'); 
newfl=optimizeCbModel(con_version);
%% 


pos_dro=[];
for n=1:length(drodro)
pos2=find(strcmp(model.rxns,drodro{n}));
pos_dro=[pos_dro,pos2]
end
ogfl_dro=og.v(pos_dro)
lis_dro=model.rxns(pos_dro);
dro_subs=model.subSystems(pos_dro);
con_flux_dro=newfl.v(pos_dro);
list=unique(vertcat(sens_dro,sens_control))
% Pre-filter to get only reactions with genes
reactions_with_genes = {};  % This will store the actual reaction names
for n=1:length(model.rxns)
 % pos=find(strcmp(model.rxns,list{n}));
  %  rxns=model.rxns(pos)
   rxns=model.rxns(n);
    [geneList] = findGenesFromRxns(model, rxns);
    if ~isempty([geneList{:}])  % Check if any genes exist
       reactions_with_genes{end+1} = model.rxns{n};  % Store the reaction name directly
     %  reactions_with_genes{end+1} = model.rxns{pos}; 
    end

end


%%
%%
%% load bulk rna seq
%data=readcell('Corrected.Data.from.v3.1.1.Genome.Counts_BLH.csv')
data2=readcell('DEG_results_Stata.csv')
titles=data2(1,:);
stata1_day10_FC=(find(contains(titles,'Stata1_Drought_Day10_vs_Control_Day10_log2FC')));
stata1_day10_p=(find(contains(titles,'Stata1_Drought_Day10_vs_Control_Day10_adjP')));
%% load bulk rna seq
data2=readcell('DEG_results_Stata.csv')
titles=data2(1,:);
stata1_day10_FC=(find(contains(titles,'Stata1_Drought_Day10_vs_Control_Day10_log2FC')));
stata1_day10_p=(find(contains(titles,'Stata1_Drought_Day10_vs_Control_Day10_adjP')));

% NEW: Find TPM mean columns
stata1_control_mean = find(contains(titles,'Stata1_Drought_Day10') & contains(titles,'mean2'));
stata1_drought_mean = find(contains(titles,'Stata1_Drought_Day10') & contains(titles,'mean1'));

%%
total_dro=[]; pathway={};model_FC1=[];
reactions_failed_filter = {}; % Reactions where no genes passed
reactions_failed_reasons = {}; % Summary of why they failed
gene_failure_details = {}; % Detailed breakdown
% Track reactions with no usable genes
reactions_no_genes = [];
for n=1:length(reactions_with_genes)
pos=find(strcmp(model.rxns,reactions_with_genes{n}));
rxns=model.rxns(pos);
epsilon = 1e-6;  % Small value to avoid division by zero
model_Fc = log2((abs(og.v(pos)) + epsilon) / (abs(newfl.v(pos)) + epsilon));

[geneList] = findGenesFromRxns(model, rxns);
     genies={};
for j=1:length(geneList)
genies=[genies;geneList{j}];
end
genies=unique(genies);
if isempty(genies)
    continue
end
if og.v(pos)==0 && newfl.v(pos)==0
    pathway{n,1} = 0;
pathway{n,2} = 1;

pathway{n,3}=0; % no change
else

% Three-tier scaling approach based on flux magnitude
flux_control = abs(newfl.v(pos));
flux_drought = abs(og.v(pos));

if flux_control == 0 && flux_drought == 0
    p_values_mod = NaN;  % Can't test if both are zero
else
    max_flux = max(flux_control, flux_drought);
    
    if max_flux < 1
    % Small fluxes: ×1000 (max ~999 counts)
    SCALE = 1000;
elseif max_flux < 10
    % Medium-small: ×100 (max ~999 counts)
    SCALE = 100;
elseif max_flux < 100
    % Medium-large: ×10 (max ~999 counts)
    SCALE = 10;
else
    % Large fluxes: ×1 (max ~330 counts)
    SCALE = 1;
end

count_control = max(1, round(flux_control * SCALE));
count_drought = max(1, round(flux_drought * SCALE));
    
    % Create contingency table
    total = count_control + count_drought;
    tab = [count_control, total - count_control;
           count_drought, total - count_drought];
    
    [h, p] = fishertest(tab);
    p_values_mod = p;
end
pathway{n,1} = model_Fc;
pathway{n,2} = p_values_mod;

if  p_values_mod>0.05 || model_Fc==0
pathway{n,3}=0; % no change
elseif p_values_mod<0.05 && model_Fc>0
    pathway{n,3}=1; % increase

elseif p_values_mod<0.05 && model_Fc<0
pathway{n,3}=-1; % decrease

end
end

%% ========== NEW CODE STARTS HERE ==========
    % NEW: Get reaction categorical value from pathway array
    rxn_category = pathway{n,3};
    
    % NEW: Initialize arrays for gene matching
    gene_matches = [];
    gene_categories = [];
    gene_expression_details = {}; % Store details for each gene
    genes_analyzed = 0; % Counter for genes that pass filters  %

    % NEW: Get gene names from RNA-seq data (first column)
    gene_names = data2(:,1); % Assuming first column contains gene names
    
    % NEW: Process each gene for this reaction
    for g = 1:length(genies)
        gene_id = genies{g};
        
        % NEW: Find gene in RNA-seq data
        gene_row = find(contains(gene_names, gene_id));
        
        if ~isempty(gene_row)
             % NEW: Get TPM values and filter
    control_tpm_raw = data2(gene_row, stata1_control_mean);
    drought_tpm_raw = data2(gene_row, stata1_drought_mean);
  
            % NEW: Get fold change and p-value for this gene
fc_value_raw = data2(gene_row, stata1_day10_FC);
p_value_raw = data2(gene_row, stata1_day10_p);

    % Handle potential NA or string values for TPM
    if ischar(control_tpm_raw{1}) || isstring(control_tpm_raw{1})
        if strcmpi(control_tpm_raw{1}, 'NA')
            continue; % Skip this gene
        else
            control_tpm = str2double(control_tpm_raw{1});
        end
    else
        control_tpm = cell2mat(control_tpm_raw);
    end
    
    if ischar(drought_tpm_raw{1}) || isstring(drought_tpm_raw{1})
        if strcmpi(drought_tpm_raw{1}, 'NA')
            continue; % Skip this gene
        else
            drought_tpm = str2double(drought_tpm_raw{1});
        end
    else
        drought_tpm = cell2mat(drought_tpm_raw);
    end
    
    % Apply TPM threshold
    TPM_threshold = 1; % Standard threshold
    if max(control_tpm, drought_tpm) < TPM_threshold
        continue; % Skip lowly expressed genes
    end
% NEW: Check if values are 'NA' strings or already numeric
if ischar(fc_value_raw{1}) || isstring(fc_value_raw{1})
    if strcmpi(fc_value_raw{1}, 'NA')
        continue; % Skip this gene
    else
        fc_value = str2double(fc_value_raw{1});
    end
else
    fc_value = cell2mat(fc_value_raw);
end

if ischar(p_value_raw{1}) || isstring(p_value_raw{1})
    if strcmpi(p_value_raw{1}, 'NA')
        continue; % Skip this gene
    else
        p_value = str2double(p_value_raw{1});
    end
else
    p_value = cell2mat(p_value_raw);
end

% NEW: Additional check for NaN values after conversion
if isnan(p_value) || isnan(fc_value)
    continue; % Skip this gene
end
% If we get here, the gene passed all filters        
genes_analyzed = genes_analyzed + 1;  

% NEW: Categorize gene expression (-1, 0, 1)   % NEW: Categorize gene expression (-1, 0, 1)
            if  p_value > 0.05
                gene_cat = 0; % No significant change
            elseif p_value <= 0.05 && fc_value > 0
                gene_cat = 1; % Significant increase
            elseif p_value <= 0.05 && fc_value < 0
                gene_cat = -1; % Significant decrease
            else
                gene_cat = 0; % Default to no change
            end
            
            % NEW: Store gene category
            gene_categories(end+1) = gene_cat;
            
           % NEW: Check if gene category matches reaction category
            if gene_cat == rxn_category
                gene_matches(end+1) = 1;
            else
                gene_matches(end+1) = 0;
            end
            
            % NEW: Store detailed gene information
            gene_expression_details{end+1,1} = gene_id;
            gene_expression_details{end,2} = fc_value;
            gene_expression_details{end,3} = p_value;
            gene_expression_details{end,4} = gene_cat;
        end
    end
    
    % NOW ADD THIS:
    if genes_analyzed == 0
        % No genes passed filters - mark as not analyzable
        validation_match = NaN; % Use NaN to indicate "not analyzable"
        reactions_no_genes(end+1) = n; % Track this reaction
        
        %% ADD: NEW TRACKING CODE STARTS HERE
        reactions_failed_filter{end+1} = reactions_with_genes{n};
        
        % Count failure types for this reaction
        na_count = 0;
        low_tpm_count = 0;
        not_found_count = 0;
        
        % Re-check each gene to categorize failures
        for g = 1:length(genies)
            gene_id = genies{g};
            gene_row = find(contains(gene_names, gene_id));
            
            if isempty(gene_row)
                not_found_count = not_found_count + 1;
                gene_failure_details{end+1} = sprintf('%s: %s - Not in RNA-seq', reactions_with_genes{n}, gene_id);
            else
                % Check TPM values
                control_tpm_raw = data2(gene_row, stata1_control_mean);
                drought_tpm_raw = data2(gene_row, stata1_drought_mean);
                
                % Check if NA
                if (ischar(control_tpm_raw{1}) && strcmpi(control_tpm_raw{1}, 'NA')) || ...
                   (ischar(drought_tpm_raw{1}) && strcmpi(drought_tpm_raw{1}, 'NA'))
                    na_count = na_count + 1;
                    gene_failure_details{end+1} = sprintf('%s: %s - TPM is NA', reactions_with_genes{n}, gene_id);
                else
                    % Convert and check threshold
                    if ischar(control_tpm_raw{1})
                        control_tpm = str2double(control_tpm_raw{1});
                    else
                        control_tpm = cell2mat(control_tpm_raw);
                    end
                    
                    if ischar(drought_tpm_raw{1})
                        drought_tpm = str2double(drought_tpm_raw{1});
                    else
                        drought_tpm = cell2mat(drought_tpm_raw);
                    end
                    
                    if max(control_tpm, drought_tpm) < 1
                        low_tpm_count = low_tpm_count + 1;
                        gene_failure_details{end+1} = sprintf('%s: %s - Low TPM (%.2f)', reactions_with_genes{n}, gene_id, max(control_tpm, drought_tpm));
                    else
                        % Must be FC or p-value NA
                        na_count = na_count + 1;
                        gene_failure_details{end+1} = sprintf('%s: %s - FC or p-value is NA', reactions_with_genes{n}, gene_id);
                    end
                end
            end
        end
        
        % Store summary
        reactions_failed_reasons{end+1} = sprintf('%s: %d genes total - %d not found, %d NA values, %d low TPM', ...
            reactions_with_genes{n}, length(genies), not_found_count, na_count, low_tpm_count);
        %% ADD: NEW TRACKING CODE ENDS HERE
        
        % Store info indicating no analyzable genes
        pathway{n,4} = NaN;
        pathway{n,5} = 0; % No genes analyzed
        pathway{n,6} = 0; % No matching genes
        pathway{n,7} = {'No genes passed TPM filter'}; % Explanation
    else
        % NEW: Determine if at least one gene matches the reaction
        if ~isempty(gene_matches) && any(gene_matches)
            validation_match = 1;
        else
            validation_match = 0;
        end
        
        % NEW: Add validation result to pathway array (column 4)
        pathway{n,4} = validation_match;
        
        % NEW: Store additional information in pathway array
        pathway{n,5} = length(gene_categories); % Number of genes analyzed
        pathway{n,6} = sum(gene_matches); % Number of matching genes
        pathway{n,7} = gene_expression_details; % Detailed gene info
    end
    
  
end
%% Output key C4 cycle fluxes
fprintf('\n=== Key C4 Cycle Fluxes ===\n');

% Find reaction indices
co2 = find(contains(model.rxns,'EX_CARBON-DIOXIDE_EXTRACELLULAR'));
CA1 = find(contains(model.rxns,'CARBODEHYDRAT-RXN'));
pepc = find(contains(model.rxns,'PEPCARBOX'));
malo = find(contains(model.rxns,'MALATE-DEHYDROGENASE-NADP+-RXN[M]'));
mal = find(contains(model.rxns,'ATR_MAL_[cb]_[cm]'));
me = find(contains(model.rxns,'MALIC-NADP-RXN[B]'));
rub = find(contains(model.rxns,'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN[B]'));
pyr = find(contains(model.rxns,'ATR_PYRUVATE_[cb]_[cm]'));
nit = find(contains(model.rxns,'EX_NITRATE_EXTRACELLULAR'));
asp=find(contains(model.rxns,'ATR_L-ASPARTATE_[cb]_[cm]'));
ala=find(contains(model.rxns,'ATR_L-ALPHA-ALANINE[cb]_[cm]'))

% Display fluxes for drought condition
fprintf('\nDrought condition fluxes:\n');
fprintf('CO2 exchange: %.4f\n', og.v(co2));
if ~isempty(CA1)
    fprintf('Carbonic anhydrase: %.4f\n', og.v(CA1));
end
if ~isempty(pepc)
    fprintf('PEPC: %.4f\n', og.v(pepc));
end
fprintf('Malate dehydrogenase (NADP+): %.4f\n', og.v(malo));
fprintf('Malate transport: %.4f\n', og.v(mal));
if ~isempty(me)
    fprintf('Malic enzyme: %.4f\n', og.v(me));
end
fprintf('Rubisco: %.4f\n', og.v(rub));
fprintf('Pyruvate transport: %.4f\n', og.v(pyr));
if ~isempty(nit)
    fprintf('Nitrate exchange: %.4f\n', og.v(nit));
end
if ~isempty(asp)
    fprintf('Aspartate transport: %.4f\n', og.v(asp));
end
if ~isempty(ala)
    fprintf('Alanine transport: %.4f\n', og.v(ala));
end
% Display fluxes for control condition
fprintf('\nControl condition fluxes:\n');
fprintf('CO2 exchange: %.4f\n', newfl.v(co2));
if ~isempty(CA1)
    fprintf('Carbonic anhydrase: %.4f\n', newfl.v(CA1));
end
if ~isempty(pepc)
    fprintf('PEPC: %.4f\n', newfl.v(pepc));
end
fprintf('Malate dehydrogenase (NADP+): %.4f\n', newfl.v(malo));
fprintf('Malate transport: %.4f\n', newfl.v(mal));
if ~isempty(me)
    fprintf('Malic enzyme: %.4f\n', newfl.v(me));
end
fprintf('Rubisco: %.4f\n', newfl.v(rub));
fprintf('Pyruvate transport: %.4f\n', newfl.v(pyr));
if ~isempty(nit)
    fprintf('Nitrate exchange: %.4f\n', newfl.v(nit));
end

%% FDR Adjustment for Model P-values
fprintf('\n=== Applying FDR Correction to Model P-values ===\n');

% Extract all p-values from pathway array
all_p_values = [];
valid_indices = [];
for n = 1:size(pathway,1)
    if ~isempty(pathway{n,2}) && ~isnan(pathway{n,2})
        all_p_values(end+1) = pathway{n,2};
        valid_indices(end+1) = n;
    end
end

% Apply FDR correction using Benjamini-Hochberg method
if ~isempty(all_p_values)
    % Sort p-values and keep track of original positions
    [sorted_p, sort_idx] = sort(all_p_values);
    n_tests = length(sorted_p);

    % Calculate FDR-adjusted p-values
    adjusted_p = zeros(size(sorted_p));
    for i = 1:n_tests
        adjusted_p(i) = sorted_p(i) * n_tests / i;
    end

    % Ensure monotonicity
    for i = n_tests-1:-1:1
        if adjusted_p(i) > adjusted_p(i+1)
            adjusted_p(i) = adjusted_p(i+1);
        end
    end

    % Cap at 1
    adjusted_p(adjusted_p > 1) = 1;

    % Map back to original order
    final_adjusted_p = zeros(size(all_p_values));
    final_adjusted_p(sort_idx) = adjusted_p;

    % OVERWRITE p-values in column 2 and update categories in column 3
    for i = 1:length(valid_indices)
        n = valid_indices(i);
        pathway{n,2} = final_adjusted_p(i);  % Overwrite with adjusted p-value

        % Recalculate category based on adjusted p-value
        adj_p = pathway{n,2};
        model_Fc = pathway{n,1};

        if adj_p > 0.05 || model_Fc == 0
            pathway{n,3} = 0; % no change
        elseif adj_p <= 0.05 && model_Fc > 0
            pathway{n,3} = 1; % increase
        elseif adj_p <= 0.05 && model_Fc < 0
            pathway{n,3} = -1; % decrease
        end
    end

    fprintf('Adjusted %d p-values using FDR (Benjamini-Hochberg)\n', length(all_p_values));
    fprintf('Before FDR: %d significant (p < 0.05)\n', sum(all_p_values < 0.05));
    fprintf('After FDR: %d significant (adjusted p < 0.05)\n', sum(final_adjusted_p < 0.05));

    fprintf('\n=== Re-validating with FDR-adjusted categories ===\n');
    
    for i = 1:length(valid_indices)
        n = valid_indices(i);
        
        % Skip if no gene data or if reaction had no analyzable genes
        if ~iscell(pathway{n,7}) || isempty(pathway{n,7}) || isnan(pathway{n,4})
            continue;
        end
        
        % Get the UPDATED reaction category (based on FDR-adjusted p-value)
        updated_rxn_category = pathway{n,3};
        
        % Get gene details
        gene_details = pathway{n,7};
        
        % Skip if gene_details is just a message string
        if ~iscell(gene_details) || (size(gene_details,1) == 1 && ischar(gene_details{1}))
            continue;
        end
        
        % Re-check validation with updated category
        gene_matches = [];
        for g = 1:size(gene_details,1)
            gene_cat = gene_details{g,4}; % Gene category (-1, 0, 1)
            
            % Check if gene matches UPDATED reaction category
            if gene_cat == updated_rxn_category
                gene_matches(end+1) = 1;
            else
                gene_matches(end+1) = 0;
            end
        end
        
        % Update validation status
        if any(gene_matches)
            pathway{n,4} = 1; % Validated
        else
            pathway{n,4} = 0; % Not validated
        end
        
        % Update matching gene count
        pathway{n,6} = sum(gene_matches);
    end
    
    fprintf('Re-validation complete with FDR-adjusted categories\n');
end  % This closes the if ~isempty(all_p_values)


% NEW: Calculate summary statistics (excluding non-analyzable reactions)
total_reactions = length(reactions_with_genes);
analyzable_reactions = sum(~isnan(cell2mat(pathway(:,4))));
non_analyzable_reactions = sum(isnan(cell2mat(pathway(:,4))));
validated_reactions = sum(cell2mat(pathway(:,4)) == 1, 'omitnan');
validation_percentage = (validated_reactions / analyzable_reactions) * 100;

fprintf('\n=== Model Validation Summary ===\n');
fprintf('Total reactions with genes: %d\n', total_reactions);
fprintf('Reactions with analyzable genes (TPM >= 1): %d\n', analyzable_reactions);
fprintf('Reactions with no genes passing filters: %d\n', non_analyzable_reactions);
fprintf('Validated reactions (of analyzable): %d\n', validated_reactions);
fprintf('Validation percentage: %.2f%% (of analyzable reactions)\n', validation_percentage);
fprintf('Overall validation: %.2f%% (of all reactions)\n', (validated_reactions/total_reactions)*100);
% NEW: Create summary table for easier viewing
validation_summary = {};
for n = 1:size(pathway,1)
    validation_summary{n,1} = reactions_with_genes{n};
    validation_summary{n,2} = pathway{n,1}; % Model FC
    validation_summary{n,3} = pathway{n,2}; % p-value
    validation_summary{n,4} = pathway{n,3}; % Reaction category
    validation_summary{n,5} = pathway{n,4}; % Validation (0 or 1)
    validation_summary{n,6} = pathway{n,5}; % Number of genes
    validation_summary{n,7} = pathway{n,6}; % Number of matching genes
end

% NEW: Convert to table for display
summary_table = cell2table(validation_summary, ...
    'VariableNames', {'Reaction', 'Model_FC', 'P_value', 'Rxn_Category', ...
    'Validated', 'Total_Genes', 'Matching_Genes'});

% NEW: Display first few rows
disp('Sample validation results:');
disp(summary_table(1:min(10, height(summary_table)), :));

% NEW: Save results
save('model_validation_resultsSEPT.mat', 'pathway', 'validation_summary', 'summary_table');

% NEW: Export to CSV
writetable(summary_table, 'reaction_validation_summarySEPT.csv');
% NEW: Category-wise validation analysis (only for analyzable reactions)
rxn_categories = cell2mat(pathway(:,3));
validation_status = cell2mat(pathway(:,4));
analyzable_idx = ~isnan(validation_status);

categories = unique(rxn_categories);
fprintf('\n=== Validation by Category (Analyzable Reactions Only) ===\n');
for c = 1:length(categories)
    cat_idx = rxn_categories == categories(c) & analyzable_idx;
    cat_validated = sum(validation_status(cat_idx) == 1);
    cat_total = sum(cat_idx);
    if cat_total > 0
        fprintf('Category %d: %d/%d reactions validated (%.1f%%)\n', ...
            categories(c), cat_validated, cat_total, ...
            (cat_validated/cat_total)*100);
    else
        fprintf('Category %d: No analyzable reactions\n', categories(c));
    end
end
%% creating generic gene object per reaction and considering -1+1
%% CREATE GENERIC GENE OBJECT - One gene per reaction
fprintf('\n=== Creating Generic Gene Objects ===\n');

% Initialize storage for generic gene data
generic_gene_data = {};
filtered_reactions = {};
rxn_count = 0;

for n = 1:size(pathway,1)
    % Skip if no pathway data
    if isempty(pathway{n,1})
        continue;
    end
    
    % Get reaction info
    rxn_name = reactions_with_genes{n};
    rxn_category = pathway{n,3};
    
    % Skip if reaction category is 0
    if rxn_category == 0
        continue;
    end
    
    % Get gene details
    if ~iscell(pathway{n,7}) || isempty(pathway{n,7})
        continue;
    end
    
    gene_details = pathway{n,7};
    
    % Check if it's a cell array with actual gene data
    if ~iscell(gene_details) || size(gene_details,2) < 4
        continue;
    end
    
    % Find gene with maximum |log2FC| (or just take the first/only gene)
max_abs_fc = -1;  % Start with -1 to ensure at least one gene is selected
selected_gene_idx = 0;

for g = 1:size(gene_details,1)
    % Extract data from cell array
    gene_fc = gene_details{g,2};
    gene_p = gene_details{g,3};
    gene_cat = gene_details{g,4};
    
    % Consider all genes - use >= instead of > to handle FC=0 cases
    if abs(gene_fc) >= max_abs_fc  % CHANGE: >= instead of >
        max_abs_fc = abs(gene_fc);
        selected_gene_idx = g;
    end
end
    % If no significant gene found, skip this reaction
   % if selected_gene_idx == 0
    %    continue;
    %end
    
    % Extract the selected gene's data directly from gene_details
    selected_gene_id = gene_details{selected_gene_idx, 1};
    selected_gene_fc = gene_details{selected_gene_idx, 2};
    selected_gene_p = gene_details{selected_gene_idx, 3};
    selected_gene_cat = gene_details{selected_gene_idx, 4};
    
    % Double-check gene category is not 0
   % if selected_gene_cat == 0
    %    continue;
    %end
    
    % Now we have a valid reaction-gene pair (both -1 or +1)
    rxn_count = rxn_count + 1;
    
    % Store for generic gene structure
    generic_gene_data{rxn_count} = struct(...
        'reaction_name', rxn_name, ...
        'reaction_category', rxn_category, ...
        'gene_name', selected_gene_id, ...
        'gene_FC', selected_gene_fc, ...
        'gene_pvalue', selected_gene_p, ...
        'gene_category', selected_gene_cat, ...
        'match', rxn_category == selected_gene_cat);
end

% Check if we have any valid pairs
if rxn_count == 0
    fprintf('WARNING: No valid reaction-gene pairs found!\n');
    fprintf('Check that there are reactions and genes with categories -1 or +1\n');
    return;
end
%% Update pathway array with generic gene object validation results
fprintf('\n=== Updating pathway validation with generic gene results ===\n');

% Create lookup for faster searching
for i = 1:length(generic_gene_data)
    rxn_name = generic_gene_data{i}.reaction_name;
    match_status = generic_gene_data{i}.match;
    
    % Find corresponding index in pathway array
    pathway_idx = find(strcmp(reactions_with_genes, rxn_name));
    
    if ~isempty(pathway_idx)
        % Update validation status (column 4) with generic gene match result
        pathway{pathway_idx, 4} = double(match_status);
        
            end
end

fprintf('Updated %d pathway entries with generic gene validation results\n', length(generic_gene_data));
%% CALCULATE COUNTS FOR HYPERGEOMETRIC TEST
fprintf('\n=== Counts for Hypergeometric Test ===\n');

% Extract categories into arrays for easier counting
rxn_cats = cellfun(@(x) x.reaction_category, generic_gene_data);
gene_cats = cellfun(@(x) x.gene_category, generic_gene_data);
% Count categories (now including 0)
rxn_neg = sum(rxn_cats == -1);
rxn_zero = sum(rxn_cats == 0);  % ADD THIS
rxn_pos = sum(rxn_cats == 1);
gene_neg = sum(gene_cats == -1);
gene_zero = sum(gene_cats == 0);  % ADD THIS
gene_pos = sum(gene_cats == 1);
N_total = length(generic_gene_data);

% Count matches (exact matches across all three categories)
matches_neg = sum(rxn_cats == -1 & gene_cats == -1);
matches_zero = sum(rxn_cats == 0 & gene_cats == 0);  % ADD THIS
matches_pos = sum(rxn_cats == 1 & gene_cats == 1);
total_matches = matches_neg + matches_zero + matches_pos;  % UPDATE THIS

% Display counts
fprintf('\nTotal reaction-gene pairs: %d\n', N_total);
fprintf('Reactions: %d negative, %d zero, %d positive\n', rxn_neg, rxn_zero, rxn_pos);
fprintf('Genes: %d negative, %d zero, %d positive\n', gene_neg, gene_zero, gene_pos);
fprintf('\nObserved matches:\n');
fprintf('  Negative matches (-1/-1): %d\n', matches_neg);
fprintf('  Zero matches (0/0): %d\n', matches_zero);  % ADD THIS
fprintf('  Positive matches (+1/+1): %d\n', matches_pos);
fprintf('  Total matches: %d (%.1f%%)\n', total_matches, 100*total_matches/N_total);
%% CALCULATE EXPECTED MATCHES BY CHANCE
fprintf('\n=== Expected Matches by Chance ===\n');

% Calculate expected matches for each category
expected_neg_matches = (rxn_neg * gene_neg) / N_total;
expected_pos_matches = (rxn_pos * gene_pos) / N_total;

% Total expected matches
expected_total = expected_neg_matches + expected_pos_matches;

% Display results
fprintf('Expected -1/-1 matches: %.2f\n', expected_neg_matches);
fprintf('Expected +1/+1 matches: %.2f\n', expected_pos_matches);
fprintf('Total expected matches: %.2f\n', expected_total);

% Also calculate as percentages for reference
expected_percent = (expected_total / N_total) * 100;
observed_percent = (total_matches / N_total) * 100;

fprintf('\nAs percentages:\n');
fprintf('Expected: %.1f%% of pairs\n', expected_percent);
fprintf('Observed: %.1f%% of pairs\n', observed_percent);

%% Test statistical significance
%% Test statistical significance using Hypergeometric Distribution
fprintf('\n=== Hypergeometric Distribution Tests ===\n');

% Test for negative matches
M = N_total;  % Total reaction-gene pairs
K = gene_neg;  % Total genes with -1
N = rxn_neg;   % Total reactions with -1
k = matches_neg;  % Observed matches for -1

% Calculate expected and p-value for negative
expected_neg = (N * K) / M;
p_value_neg = 1 - hygecdf(k-1, M, K, N);
fprintf('\nNegative category:\n');
fprintf('  Expected: %.2f, Observed: %d\n', expected_neg, k);
fprintf('  P-value: %.4f\n', p_value_neg);

% Test for positive matches
K = gene_pos;  % Total genes with +1
N = rxn_pos;   % Total reactions with +1
k = matches_pos;  % Observed matches for +1

% Calculate expected and p-value for positive
expected_pos = (N * K) / M;
p_value_pos = 1 - hygecdf(k-1, M, K, N);
fprintf('\nPositive category:\n');
fprintf('  Expected: %.2f, Observed: %d\n', expected_pos, k);
fprintf('  P-value: %.4f\n', p_value_pos);



% Also report individual significances
if p_value_neg < 0.05
    fprintf('  Negative matches are significantly enriched (p = %.4f)\n', p_value_neg);
else
    fprintf('  Negative matches are not significantly enriched (p = %.4f)\n', p_value_neg);
end

if p_value_pos < 0.05
    fprintf('  Positive matches are significantly enriched (p = %.4f)\n', p_value_pos);
else
    fprintf('  Positive matches are not significantly enriched (p = %.4f)\n', p_value_pos);
end

%% Isolating bottleneck reactions

% Create binary marker for each reaction
in_list_marker = zeros(length(reactions_with_genes), 1);
for n = 1:length(reactions_with_genes)
    if find(strcmp(list,reactions_with_genes{n}))
        in_list_marker(n) = 1;
    else
        in_list_marker(n) = 0;
    end
    % Add to pathway array as column 8
    pathway{n,8} = in_list_marker(n);
end

%% Output all Category 1 reactions (significant increases)
fprintf('\n=== Category 1 Reactions (Significant Increases) ===\n');

cat1_reactions = {};
cat1_count = 0;

for n = 1:length(reactions_with_genes)
    if ~isempty(pathway{n,3}) && pathway{n,3} == 1
        cat1_count = cat1_count + 1;
        cat1_reactions{cat1_count,1} = reactions_with_genes{n};
        cat1_reactions{cat1_count,2} = pathway{n,1}; % Model FC
        cat1_reactions{cat1_count,3} = pathway{n,2}; % Adjusted p-value
        cat1_reactions{cat1_count,4} = pathway{n,4}; % Validation status
        
        % Get subsystem if you want to include it
        pos = find(strcmp(model.rxns, reactions_with_genes{n}));
        if ~isempty(pos)
            cat1_reactions{cat1_count,5} = model.subSystems{pos};
        end
    end
end

fprintf('Total Category 1 reactions: %d\n\n', cat1_count);

% Display the reactions
for i = 1:cat1_count
    fprintf('%d. %s\n', i, cat1_reactions{i,1});
    fprintf('   Model FC: %.4f, Adj p-value: %.4e\n', cat1_reactions{i,2}, cat1_reactions{i,3});
    if ~isnan(cat1_reactions{i,4})
        fprintf('   Validated: %s\n', mat2str(cat1_reactions{i,4}));
    end
    
    fprintf('\n');
end

%% Verify the count with more detail
epsilon_flux = 1e-5;
cat0_indices = find(rxn_categories == 0);

nonzero_both = 0;
zero_drought = 0;
zero_control = 0;
zero_either = 0;

for i = 1:length(cat0_indices)
    rxn_name = reactions_with_genes{cat0_indices(i)};
    pos = find(strcmp(model.rxns, rxn_name));
    
    if ~isempty(pos)
        drought_active = abs(og.v(pos)) > epsilon_flux;
        control_active = abs(newfl.v(pos)) > epsilon_flux;
        
        if drought_active && control_active
            nonzero_both = nonzero_both + 1;
        elseif ~drought_active
            zero_drought = zero_drought + 1;
        elseif ~control_active
            zero_control = zero_control + 1;
        end
        
        if ~drought_active || ~control_active
            zero_either = zero_either + 1;
        end
    end
end

fprintf('\nCategory 0 breakdown (epsilon = %.6f):\n', epsilon_flux);
fprintf('  Active in both: %d\n', nonzero_both);
fprintf('  Inactive in drought: %d\n', zero_drought);
fprintf('  Inactive in control: %d\n', zero_control);
fprintf('  Inactive in at least one: %d\n', zero_either);
fprintf('  Total: %d\n', length(cat0_indices));
