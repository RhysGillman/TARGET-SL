function [] = main_SCS( network, cancer_type, run_mode)

[original,~,~]=importdata('SCS/EXPR_internate.txt');
EXPR_data=original.textdata;
% Get the sample names
sample_names = string(EXPR_data(1,2:end));
% Get the sample numbers
sample_number = (1:size(EXPR_data,2)-1);
samples = struct('Sample_Number', [], 'Sample_Name', []);
for i = 1:size(EXPR_data,2)-1
    samples(i).Sample_Number = i;
    samples(i).Sample_Name = EXPR_data(1,i+1);
end
sample_table = struct2table(samples);

% Write the table to a new file
writetable(sample_table, strcat('../results/',run_mode,'/network_',network,'/SCS/',cancer_type,'/sample_names.txt'));

end
