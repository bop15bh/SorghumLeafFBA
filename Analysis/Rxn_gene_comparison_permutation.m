%% Model validation
% comparing if gene expression per reaction matches model flux change
clear 
close ALL

changeCobraSolver('glpk');
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

load('sens_droJune25.mat')
load('sens_controlJune25.mat')
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
model = changeRxnBounds(amodel,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -44, 'l'); 

model = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', 0, 'u'); 
og=optimizeCbModel(model);
con_version = changeRxnBounds(model,'EX_CARBON-DIOXIDE_EXTRACELLULAR', -294.6979, 'l'); 
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
for n=1:length(sens_dro)
  pos=find(strcmp(model.rxns,sens_control{n}));
    rxns=model.rxns(pos)
  % rxns=model.rxns(n);
    [geneList] = findGenesFromRxns(model, rxns);
    if ~isempty([geneList{:}])  % Check if any genes exist
      % reactions_with_genes{end+1} = model.rxns{n};  % Store the reaction name directly
        reactions_with_genes{end+1} = model.rxns{pos}; 
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
genies=[genies;geneList{j}]
end
genies=unique(genies)
if isempty(genies)
    continue
end
if og.v(pos)==0 && newfl.v(pos)==0
    pathway{n,1} = 0;
pathway{n,2} = 1;

pathway{n,3}=0; % no change
else



% Loop through each ratio and perform a z-test
if mean(abs(og.v(pos)))<1 && mean(abs(newfl.v(pos)))<1 &&  mean(abs(og.v(pos))) > 0 && mean(abs(newfl.v(pos))) > 0
count_dro = round(mean(abs(og.v(pos)))*100);
count_og = round(mean(abs(newfl.v(pos)))*100);
tab = [count_og, 100 - count_og;
count_dro, 100 - count_dro];
[h,p] = fishertest(tab);
p_values_mod = p; % Right-tailed test
elseif xor(mean(abs(og.v(pos))) == 0, mean(abs(newfl.v(pos))) == 0)
    epsilon = 1;  % fallback count when flux is zero (represents 1%)

flux_dro = mean(abs(og.v(pos)));
flux_og  = mean(abs(newfl.v(pos)));

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

elseif mean(abs(og.v(pos)))==0 && mean(abs(newfl.v(pos)))==0
            p_values_mod = NaN; % Right-tailed test
else
% Get mean fluxes
flux_control = mean(abs(newfl.v(pos)));
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

    % ADD THIS RE-VALIDATION PART HERE (still inside the if statement)
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

%%
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
save('model_validation_results.mat', 'pathway', 'validation_summary', 'summary_table');

% NEW: Export to CSV
writetable(summary_table, 'reaction_validation_summary.csv');
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
%% Permutation Test for Validation Significance
fprintf('\n=== Permutation Test for Statistical Significance ===\n');

% First, identify which pathway rows actually have data
valid_pathway_rows = [];
for n = 1:size(pathway,1)
    if ~isempty(pathway{n,1}) % Check if row has data
        valid_pathway_rows(end+1) = n;
    end
end

% Get the observed validation percentage (only for analyzable reactions)
observed_validated = 0;
observed_analyzable = 0;
for n = valid_pathway_rows
    if ~isnan(pathway{n,4})
        observed_analyzable = observed_analyzable + 1;
        if pathway{n,4} == 1
            observed_validated = observed_validated + 1;
        end
    end
end
observed_percentage = (observed_validated / observed_analyzable) * 100;

% Set up permutation test
n_permutations = 1000;
permuted_percentages = zeros(n_permutations, 1);

% Create a pool of all gene expression profiles from analyzable reactions
all_gene_data = {};
reaction_gene_counts = [];
reaction_categories = [];  % Store the reaction categories too

for n = valid_pathway_rows
    if ~isnan(pathway{n,4}) && iscell(pathway{n,7})
        gene_details = pathway{n,7};
        if ~isempty(gene_details) && size(gene_details,1) > 0
            % Store gene expression categories
            for g = 1:size(gene_details,1)
                all_gene_data{end+1} = gene_details{g,4}; % Gene category (-1,0,1)
            end
            reaction_gene_counts(end+1) = size(gene_details,1);
            reaction_categories(end+1) = pathway{n,3};  % Store reaction category
        end
    end
end

% Convert to array
all_gene_categories = cell2mat(all_gene_data);

% Perform permutations
fprintf('Running %d permutations...\n', n_permutations);
for perm = 1:n_permutations
    % Shuffle all gene categories
    shuffled_genes = all_gene_categories(randperm(length(all_gene_categories)));
    
    % Reassign shuffled genes to reactions
    gene_idx = 1;
    perm_validated = 0;
    perm_analyzable = 0;
    
    for rxn_idx = 1:length(reaction_gene_counts)
        % Get reaction category
        rxn_cat = reaction_categories(rxn_idx);
        
        % Get number of genes for this reaction
        n_genes = reaction_gene_counts(rxn_idx);
        
        % Extract shuffled genes for this reaction
        reaction_genes = shuffled_genes(gene_idx:gene_idx+n_genes-1);
        
        % Check if any match
        if any(reaction_genes == rxn_cat)
            perm_validated = perm_validated + 1;
        end
        perm_analyzable = perm_analyzable + 1;
        
        gene_idx = gene_idx + n_genes;
    end
    
    % Calculate validation percentage for this permutation
    permuted_percentages(perm) = (perm_validated / perm_analyzable) * 100;
    
    % Progress indicator
    if mod(perm, 100) == 0
        fprintf('  Completed %d permutations\n', perm);
    end
end

% Rest of the code remains the same...

% Calculate p-value
p_value = sum(permuted_percentages >= observed_percentage) / n_permutations;

% Display results
fprintf('\nResults:\n');
fprintf('Observed validation: %.2f%%\n', observed_percentage);
fprintf('Mean random validation: %.2f%% (SD: %.2f%%)\n', ...
    mean(permuted_percentages), std(permuted_percentages));
fprintf('P-value: %.4f\n', p_value);

if p_value < 0.001
    fprintf('Conclusion: Highly significant (p < 0.001)\n');
elseif p_value < 0.01
    fprintf('Conclusion: Very significant (p < 0.01)\n');
elseif p_value < 0.05
    fprintf('Conclusion: Significant (p < 0.05)\n');
else
    fprintf('Conclusion: Not significant (p >= 0.05)\n');
end

% Visualization
figure;
histogram(permuted_percentages, 30, 'FaceColor', [0.7 0.7 0.7]);
hold on;
xline(observed_percentage, 'r', 'LineWidth', 2);
xlabel('Validation Percentage (%)','FontSize', 40);
ylabel('Frequency','FontSize', 40);
title('Permutation Test: Null Distribution vs Observed');
legend('Random permutations', 'Observed validation', 'Location', 'northwest');

% Add percentile info
percentile_rank = (sum(permuted_percentages < observed_percentage) / n_permutations) * 100;
text(observed_percentage + 1, max(ylim)*0.8, ...
    sprintf('%.1f percentile\np = %.4f', percentile_rank, p_value), ...
    'FontSize', 40);
   x_width=30 ;y_width=30;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
              print('permutation','-djpeg','-loose');
              set(gca, 'FontSize', 40)
% Save permutation results
save('permutation_test_results.mat', 'observed_percentage', 'permuted_percentages', ...
    'p_value', 'percentile_rank');

% NEW: Simple visualization
figure;
subplot(1,2,1);
pie([validated_reactions, total_reactions - validated_reactions], ...
    {'Validated', 'Not Validated'});
title('Reaction Validation Overview');

subplot(1,2,2);
bar(categories, [sum(rxn_categories==-1), sum(rxn_categories==0), sum(rxn_categories==1)]);
xlabel('Category','FontSize', 40);
ylabel('Count','FontSize', 40);
title('Distribution of Reaction Categories');
xticklabels({'Decrease (-1)', 'No Change (0)', 'Increase (1)'});
set(gca, 'FontSize', 40)
   x_width=30 ;y_width=30;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
              print('categoricals','-djpeg','-loose');
              
%% ADD: QUICK FILTER FAILURE CHECK
fprintf('\n=== Quick Filter Failure Check ===\n');
fprintf('Reactions that failed filters:\n');
for i = 1:length(reactions_failed_reasons)
    fprintf('%s\n', reactions_failed_reasons{i});
end

% Count failure types across all failed genes
total_na = sum(cellfun(@(x) contains(x, 'NA'), gene_failure_details));
total_low_tpm = sum(cellfun(@(x) contains(x, 'Low TPM'), gene_failure_details));
total_not_found = sum(cellfun(@(x) contains(x, 'Not in RNA-seq'), gene_failure_details));

fprintf('\nTotal gene failures:\n');
fprintf('- Not in RNA-seq: %d\n', total_not_found);
fprintf('- NA values: %d\n', total_na - total_not_found); % Subtract to avoid double counting
fprintf('- Low TPM: %d\n', total_low_tpm);

% Save for quick access
save('quick_failure_check.mat', 'reactions_failed_filter', 'reactions_failed_reasons', 'gene_failure_details');

%% ADD: REACTION-LEVEL FAILURE ANALYSIS
fprintf('\n=== Reaction-Level Failure Analysis ===\n');

% Track which reactions have which types of failures
rxns_with_not_found = 0;
rxns_with_na = 0;
rxns_with_low_tpm = 0;
rxns_with_multiple_issues = 0;

% Detailed tracking for each reaction
reaction_failure_breakdown = {};

% Analyze each failed reaction
for i = 1:length(reactions_failed_filter)
    rxn_name = reactions_failed_filter{i};
    
    % Check what types of failures this reaction has
    rxn_failures = gene_failure_details(contains(gene_failure_details, rxn_name));
    
    has_not_found = any(cellfun(@(x) contains(x, 'Not in RNA-seq'), rxn_failures));
    has_na = any(cellfun(@(x) contains(x, 'NA'), rxn_failures));
    has_low_tpm = any(cellfun(@(x) contains(x, 'Low TPM'), rxn_failures));
    
    % Count reactions by failure type
    if has_not_found
        rxns_with_not_found = rxns_with_not_found + 1;
    end
    if has_na
        rxns_with_na = rxns_with_na + 1;
    end
    if has_low_tpm
        rxns_with_low_tpm = rxns_with_low_tpm + 1;
    end
    
    % Check if reaction has multiple types of failures
    failure_types = has_not_found + has_na + has_low_tpm;
    if failure_types > 1
        rxns_with_multiple_issues = rxns_with_multiple_issues + 1;
    end
    
    % Store detailed breakdown
    reaction_failure_breakdown{i,1} = rxn_name;
    reaction_failure_breakdown{i,2} = has_not_found;
    reaction_failure_breakdown{i,3} = has_na;
    reaction_failure_breakdown{i,4} = has_low_tpm;
    reaction_failure_breakdown{i,5} = failure_types;
end

fprintf('\nReactions by failure type:\n');
fprintf('- Reactions with at least one gene not in RNA-seq: %d\n', rxns_with_not_found);
fprintf('- Reactions with at least one NA value: %d\n', rxns_with_na);
fprintf('- Reactions with at least one low TPM gene: %d\n', rxns_with_low_tpm);
fprintf('- Reactions with multiple types of issues: %d\n', rxns_with_multiple_issues);

%% ADD: EXCLUSIVE CATEGORY ANALYSIS
% Count reactions that ONLY have one type of issue
rxns_only_not_found = 0;
rxns_only_na = 0;
rxns_only_low_tpm = 0;

for i = 1:length(reactions_failed_filter)
    failure_count = reaction_failure_breakdown{i,5};
    if failure_count == 1
        if reaction_failure_breakdown{i,2} % has_not_found
            rxns_only_not_found = rxns_only_not_found + 1;
        elseif reaction_failure_breakdown{i,3} % has_na
            rxns_only_na = rxns_only_na + 1;
        elseif reaction_failure_breakdown{i,4} % has_low_tpm
            rxns_only_low_tpm = rxns_only_low_tpm + 1;
        end
    end
end

fprintf('\nReactions with ONLY one type of issue:\n');
fprintf('- Only "not in RNA-seq" issues: %d\n', rxns_only_not_found);
fprintf('- Only NA value issues: %d\n', rxns_only_na);
fprintf('- Only low TPM issues: %d\n', rxns_only_low_tpm);
fprintf('- Mixed issues: %d\n', rxns_with_multiple_issues);

%% ADD: CHECK FOR POTENTIAL DOUBLE-COUNTING
fprintf('\n=== Checking for overlapping failure types ===\n');

% Re-analyze genes more carefully
genes_low_tpm_only = 0;
genes_na_only = 0;
genes_low_tpm_and_na = 0;

for i = 1:length(reactions_failed_filter)
    rxn_name = reactions_failed_filter{i};
    pos = find(strcmp(model.rxns, rxn_name));
    rxns = model.rxns(pos);
    [geneList] = findGenesFromRxns(model, rxns);
    
    genies = {};
    for j = 1:length(geneList)
        genies = [genies; geneList{j}];
    end
    genies = unique(genies);
    
    for g = 1:length(genies)
        gene_id = genies{g};
        gene_row = find(contains(gene_names, gene_id));
        
        if ~isempty(gene_row)
            % Get all values
            control_tpm_raw = data2(gene_row, stata1_control_mean);
            drought_tpm_raw = data2(gene_row, stata1_drought_mean);
            fc_value_raw = data2(gene_row, stata1_day10_FC);
            p_value_raw = data2(gene_row, stata1_day10_p);
            
            % Check TPM status
            tpm_is_low = false;
            tpm_is_na = false;
            
            % Check control TPM
            if ischar(control_tpm_raw{1}) && strcmpi(control_tpm_raw{1}, 'NA')
                tpm_is_na = true;
            elseif ischar(drought_tpm_raw{1}) && strcmpi(drought_tpm_raw{1}, 'NA')
                tpm_is_na = true;
            else
                % Convert values
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
                    tpm_is_low = true;
                end
            end
            
            % Check FC/p-value status
            fc_p_is_na = false;
            if (ischar(fc_value_raw{1}) && strcmpi(fc_value_raw{1}, 'NA')) || ...
               (ischar(p_value_raw{1}) && strcmpi(p_value_raw{1}, 'NA'))
                fc_p_is_na = true;
            end
            
            % Categorize
            if tpm_is_low && fc_p_is_na
                genes_low_tpm_and_na = genes_low_tpm_and_na + 1;
            elseif tpm_is_low
                genes_low_tpm_only = genes_low_tpm_only + 1;
            elseif tpm_is_na || fc_p_is_na
                genes_na_only = genes_na_only + 1;
            end
        end
    end
end

fprintf('\nGene-level overlap analysis:\n');
fprintf('- Genes with ONLY low TPM: %d\n', genes_low_tpm_only);
fprintf('- Genes with ONLY NA values: %d\n', genes_na_only);
fprintf('- Genes with BOTH low TPM AND NA p-value/FC: %d\n', genes_low_tpm_and_na);

% Update the save to include new analysis
save('quick_failure_check.mat', 'reactions_failed_filter', 'reactions_failed_reasons', ...
    'gene_failure_details', 'reaction_failure_breakdown');

fprintf('\nTotal failed reactions: %d out of %d\n', length(reactions_failed_filter), total_reactions);