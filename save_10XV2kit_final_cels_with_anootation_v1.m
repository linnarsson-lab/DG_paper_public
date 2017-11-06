
tic
clear all
close all

load /mnt/sanger-data2/C1_stuff/DentGyr/DG_V2kit_subsample_samples_merged_18-Jun-2017.mat


tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
tot_genes = sum(data>0);

validcells = (tot_mol>800 & (tot_mol./tot_genes)>1.2 & tot_mol<3e4 & tot_genes>800);
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
source = source(validcells);
age = age(validcells);
sex = sex(validcells);
strain = strain(validcells);

age = regexprep(age,'p','');
age = regexprep(age,'P','');
age = regexprep(age,'E','');
age = cellfun(@str2double, age);
age(age==16.5) = -16.5;

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);


marker_genes = {'Stmn2','Mog','Aldoc','C1qc','Cldn5'};
% figure('color','w','position',[100,100,900,900]); 
inrmv = false(length(cellid),1);
for i=1:length(marker_genes)
    for j=i+1:length(marker_genes)
%         subplot(4,4,(i-1)*4+j-1);
        tmp1 = data(strcmpi(geneid,marker_genes{i}),:);
        tmp2 = data(strcmpi(geneid,marker_genes{j}),:);
        thtmp1 = 1;%prctile(tmp1(tmp1>0),20);
        thtmp2 = 1;%prctile(tmp2(tmp2>0),20);
        inrmv(tmp2>thtmp2 & tmp1>thtmp1) = true;
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % 
validcells = ~inrmv;
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
source = source(validcells);
age = age(validcells);

% data = round(data./repmat(sum(data),length(data(:,1)),1)*5e3);

% % % % % % % % % % % % % % % % % % % % % % % % % fix sex by expression
% sum_m_g = sum([data(strcmpi(geneid,'ddx3y'),:);data(strcmpi(geneid,'uty'),:);data(strcmpi(geneid,'eif2s3y'),:)]);
% sum_f_g = sum([data(strcmpi(geneid,'xist'),:);data(strcmpi(geneid,'tsix'),:)]);
% in_f = sum_f_g > sum_m_g;
% in_m = sum_f_g < sum_m_g;
% in_uk = sum_f_g==0 & sum_m_g==0;%probably male
% sexcalc = zeros(length(cellid),1);
% sexcalc(in_f) = 1;
% sexcalc(in_m) = -1;
% sexcalc(in_uk) = 0;
% age = zeros(size(sexcalc));
% age(source==1 & sexcalc==1) = 12;
% age(source==1 & sexcalc==-1) = 35;
% age(source==1 & sexcalc==0) = 30;%probably, not sure
% age(source==2 & sexcalc==1) = 16;
% age(source==2 & sexcalc==-1) = 24;
% age(source==2 & sexcalc==0) = 22;%probably, not sure
% % % % % % % % % % % % % % % % % % % % % % % % % % 

% in = find(sum(data>0,2)>20 & sum(data>0,2)<length(data(1,:))*0.6);
% m_v = mean(data(in,:),2);
% cv_v = std(data(in,:),[],2)./m_v;
% log2_m = log2(m_v);
% log2_cv = log2(cv_v);
% x0 = [-0.5,1];
% % [param_fit,fval] =  run_min_cvfit(log2_m,log2_cv,x0);
% param_fit = robustfit(log2_m,log2_cv);
% param_fit = [param_fit(2), param_fit(1)];
% param_fit = round(param_fit*100)/100;
% log2_cv_fit = log2( (2.^log2_m).^param_fit(1)) + param_fit(2);
% tmp = log2_cv - log2_cv_fit;
% [~,xi1] = sort(tmp,'descend');
% corr_filt = xi1(1:12000);
% figure; plot(log2_m,log2_cv,'.','markersize',3); hold on;
% plot(log2_m(corr_filt),log2_cv(corr_filt),'.r','markersize',3); hold on;
% 
% 
% sex_genes = upper({'EHD2', 'ESPL1', 'JARID1D', 'PNPLA4', 'RPS4Y1', 'XIST','tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'});
% [~,loc] = ismember(sex_genes,upper(geneid(in)));
% corr_filt = setdiff(corr_filt,loc);
% stress_genes = upper({'Rpl26','Gstp1','Rpl35a','Erh','Slc25a5','Pgk1','Eno1','Tubb2a','Emc4','Scg5','Btg2','Fos','Klf6'....
%     'Nr4a1','Dusp1','Egr1','Jun','Fosb','Dnajb1','Ier2','Junb','Csrnp1',});
% [~,loc] = ismember(stress_genes,upper(geneid(in)));
% corr_filt = setdiff(corr_filt,loc);
% 
% data = data(in(corr_filt),:);
% geneid = geneid(in(corr_filt));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % permute the cells order
% rng(1);
% cell_perm = randperm(length(cellid));
% gene_perm = randperm(length(geneid));
% geneid = geneid(gene_perm);
% cellid = cellid(cell_perm);
% source = source(cell_perm);
% age = age(cell_perm);
% moldata = data(gene_perm,cell_perm);
moldata = data;


cellid_cluster = loadCellFile('mnt/sanger-data2/C1_stuff/DentGyr/DG_10X_V2_markertable_knnMutal_MCL_1p4_18-Jun-2017_manualinspection.txt');
cellid_cluster = cellid_cluster(2:end,:);
[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid = cellid(loc);
source = source(loc);
age = age(loc);
moldata = moldata(:,loc);


grlev2 = cellid_cluster(:,3);
grlev2 = regexprep(grlev2,'nb2','immature_granul');
% lev2uni = unique(lev2gr);
grlev2_uni = {'rgl','astro_old','astro_adol','immature_astro','rgl_young','ependymal'....
    ,'pvm','mgl','endo','vlmc','mol','nfol','opc','nipc_young','nipc','nb1'.....
    ,'immature_granul','immature_pyr','cr','gaba_young','gaba_old','pyr_ca3','granul_adol','granul_old'};%,'pyr_ca1'
% grlev2_uni = unique(grlev2);
% grlev2_uni = setdiff(grlev2_uni,'ex');

dend_order_names = regexprep(grlev2_uni,'_','-');

T_cells_tmp = zeros(length(grlev2),1);
for i=1:length(grlev2_uni)
    T_cells_tmp( strcmpi(grlev2, grlev2_uni{i}) ) = i;
end
[T_cells_tmp,xi] = sort(T_cells_tmp);
xi(T_cells_tmp==0) = [];
T_cells_tmp(T_cells_tmp==0) = [];
T_cells_tmp_uni = unique(T_cells_tmp);
cellid = cellid(xi);
source = source(xi);
age = age(xi);
moldata = moldata(:,xi);
grlev2 = grlev2(xi);





table1 = [cellid';m2c(age)';grlev2'];
table1 = [{'cellid';'age(days)';'cluster name'},table1];

saveCellFile(table1,['cell_annotation_10XV2kit_data_DG_',date,'.txt']);

clear table
T = table(moldata,'RowNames',geneid);
writetable(T,['moldata_10XV2kit_data_DG_no_annot_',date,'.txt'],'WriteRowNames',true);


toc





