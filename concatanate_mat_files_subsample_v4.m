tic
clear all
close all


files_to_load = {'10X79_1_05-May-2017.mat'
    '10X79_2_05-May-2017.mat'
    '10X80_1_05-May-2017.mat'
    '10X80_2_05-May-2017.mat'
    '10X82_1_05-May-2017.mat'
    '10X83_2_05-May-2017.mat'
    '10X83_3_05-May-2017.mat'
    '10X83_4_05-May-2017.mat'
    '10X84_2_05-May-2017.mat'
    '10X84_3_05-May-2017.mat'
    '10X84_4_18-Jun-2017.mat'};
        
        
    
ncells_files = zeros(length(files_to_load),1);
for i=1:length(files_to_load)
    fprintf(['load ',num2str(i),'\n']);
%     load(['/data2/C1_stuff/DentGyr/loom_files_10X/',files_to_load{i}]);
    load(['/mnt/sanger-data2/C1_stuff/DentGyr/loom_files_10X/',files_to_load{i}]);
    eval(['ncells_files(i) = length(cellid_',files_to_load{i}(1:end-16),');']);
end
% geneid = cellfun(@(x) x(3:end-1), genes_sym_10X79_1,'UniformOutput',false);
geneid = genes_sym_10X79_1;
data = zeros(length(geneid),sum(ncells_files));
cellid = cell(sum(ncells_files),1);
source = cell(sum(ncells_files),1);
sex = cell(sum(ncells_files),1);
age = cell(sum(ncells_files),1);
strain = cell(sum(ncells_files),1);

k = 0;
for i=1:length(files_to_load)
    fprintf(['concatnate ',num2str(i),'\n']);
    eval(['data(:,k+1:k+ncells_files(i)) = moldata_',files_to_load{i}(1:end-16),''';']);
    eval(['cellid(k+1:k+ncells_files(i)) = cellid_',files_to_load{i}(1:end-16),';']);
    
    eval(['source(k+1:k+ncells_files(i)) = cellannot_',files_to_load{i}(1:end-16),'.SampleID;']);
    eval(['sex(k+1:k+ncells_files(i)) = cellannot_',files_to_load{i}(1:end-16),'.Sex;']);
    eval(['age(k+1:k+ncells_files(i)) = cellannot_',files_to_load{i}(1:end-16),'.Age;']);
    eval(['strain(k+1:k+ncells_files(i)) = cellannot_',files_to_load{i}(1:end-16),'.Strain;']);
    
    eval(['k = k+ncells_files(i);']);
end

[geneid,ia] = unique(geneid );
data = data(ia,:);

load /mnt/sanger-data2/C1_stuff/DentGyr/10X43_1_v1.mat
[geneid_43_1,ia] = unique(geneid_43_1);
data_43_1 = data_43_1(ia,:);
data = [data,data_43_1];
source = [source;repmat({'10X43_1'},length(cellid_43_1),1)];
sex = [sex;repmat({'?'},length(cellid_43_1),1)];
% age = [age;repmat({'1000'},length(cellid_43_1),1)];
strain = [strain;repmat({'CD1'},length(cellid_43_1),1)];
cellid = [cellid;cellfun(@(x) ['10X43_1_',x], cellid_43_1,'uniformoutput',0) ];
% % % % % % % % % % % % % % % % % % % % % % % % fix sex by expression
sum_m_g = sum([data_43_1(strcmpi(geneid_43_1,'ddx3y'),:);data_43_1(strcmpi(geneid_43_1,'uty'),:);data_43_1(strcmpi(geneid_43_1,'eif2s3y'),:)]);
sum_f_g = sum([data_43_1(strcmpi(geneid_43_1,'xist'),:);data_43_1(strcmpi(geneid_43_1,'tsix'),:)]);
in_f = sum_f_g > sum_m_g;
in_m = sum_f_g < sum_m_g;
in_uk = sum_f_g==0 & sum_m_g==0;%probably male
sexcalc = zeros(length(cellid_43_1),1);
sexcalc(in_f) = 1;
sexcalc(in_m) = -1;
sexcalc(in_uk) = 0;
agetmp = zeros(size(sexcalc));
agetmp(sexcalc==1) = 12;
agetmp(sexcalc==-1) = 35;
agetmp(sexcalc==0) = 30;%probably, not sure
age = [age;cellfun(@(x) ['p',num2str(x)], m2c(agetmp),'uniformoutput',0)];

% % % % % % % % % % % % % % % % % % % % % % % % % 

load /mnt/sanger-data2/C1_stuff/DentGyr/10X46_1_v1.mat
[geneid_46_1,ia] = unique(geneid_46_1);
data_46_1 = data_46_1(ia,:);
data = [data,data_46_1];
source = [source;repmat({'10X46_1'},length(cellid_46_1),1)];
sex = [sex;repmat({'?'},length(cellid_46_1),1)];
% age = [age;repmat({'1000'},length(cellid_46_1),1)];
strain = [strain;repmat({'CD1'},length(cellid_46_1),1)];
cellid = [cellid;cellfun(@(x) ['10X46_1_',x], cellid_46_1,'uniformoutput',0) ];
% % % % % % % % % % % % % % % % % % % % % % % % fix sex by expression
sum_m_g = sum([data_46_1(strcmpi(geneid_46_1,'ddx3y'),:);data_46_1(strcmpi(geneid_46_1,'uty'),:);data_46_1(strcmpi(geneid_46_1,'eif2s3y'),:)]);
sum_f_g = sum([data_46_1(strcmpi(geneid_46_1,'xist'),:);data_46_1(strcmpi(geneid_46_1,'tsix'),:)]);
in_f = sum_f_g > sum_m_g;
in_m = sum_f_g < sum_m_g;
in_uk = sum_f_g==0 & sum_m_g==0;%probably male
sexcalc = zeros(length(cellid_46_1),1);
sexcalc(in_f) = 1;
sexcalc(in_m) = -1;
sexcalc(in_uk) = 0;
agetmp = zeros(size(sexcalc));
agetmp(sexcalc==1) = 16;
agetmp(sexcalc==-1) = 24;
agetmp(sexcalc==0) = 22;%probably, not sure
age = [age;cellfun(@(x) ['p',num2str(x)], m2c(agetmp),'uniformoutput',0)];
% % % % % % % % % % % % % % % % % % % % % % % % % 

tot_mol = sum(data)';
ind_subsam = [];

in = find(strcmpi(source,'10X79_1') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:end)];
in = find(strcmpi(source,'10X80_1') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:3500)];
in = find(strcmpi(source,'10X80_2') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:3500)];

in = find(strcmpi(source,'10X79_2') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:end)];
in = find(strcmpi(source,'10X83_4') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:3500)];
in = find(strcmpi(source,'10X84_4') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:3000)];

in = find(strcmpi(source,'10X84_2') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:2500)];
in = find(strcmpi(source,'10X84_3') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:2500)];

in = find(strcmpi(source,'10X83_2') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:2500)];
in = find(strcmpi(source,'10X83_3') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:2500)];

in = find(strcmpi(source,'10X82_1') & tot_mol>1000);
in = in(randperm(length(in)));
ind_subsam = [ind_subsam; in(1:2500)];

% in = find(strcmpi(source,'10X43_1') & tot_mol>1000);
% in = in(randperm(length(in)));
% ind_subsam = [ind_subsam; in(1:1000)];
% in = find(strcmpi(source,'10X46_1') & tot_mol>1000);
% in = in(randperm(length(in)));
% ind_subsam = [ind_subsam; in(1:1000)];

data = data(:,ind_subsam);
sex = sex(ind_subsam);
age = age(ind_subsam);
strain = strain(ind_subsam);
source = source(ind_subsam);
cellid = cellid(ind_subsam);



save(['DG_V2kit_subsample_samples_merged_',date],'data','sex','age','strain','source','cellid','geneid','-v7.3')








toc