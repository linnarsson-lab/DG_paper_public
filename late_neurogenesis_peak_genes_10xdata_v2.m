
tic
clear all
close all


files_to_load = {'10X43_1_DG_01-Nov-2016.mat'
    '10X46_1_DG_01-Nov-2016.mat'};
    
ncells_files = zeros(length(files_to_load),1);
for i=1:length(files_to_load)
    fprintf(['load ',num2str(i),'\n']);
    load(['/data2/10X_adol_mice/loom_data_files/mat_data_files/',files_to_load{i}]);
    eval(['ncells_files(i) = length(cellid_',files_to_load{i}(4:end-16),');']);
end
geneid = cellfun(@(x) x(3:end-1), genes_sym_43_1_DG,'UniformOutput',false);
data = zeros(length(geneid),sum(ncells_files));
cellid = cell(sum(ncells_files),1);
source = zeros(sum(ncells_files),1);
k = 0;
for i=1:length(files_to_load)
    fprintf(['concatnate ',num2str(i),'\n']);
    eval(['data(:,k+1:k+ncells_files(i)) = moldata_',files_to_load{i}(4:end-16),''';']);
    eval(['cellid(k+1:k+ncells_files(i)) = cellid_',files_to_load{i}(4:end-16),';']);
    eval(['source(k+1:k+ncells_files(i)) = i;']);
    eval(['k = k+ncells_files(i);']);
end


tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
tot_genes = sum(data>0);

validcells = (tot_mol>800 & (tot_mol./tot_genes)>1.2 & tot_mol<2e4 & tot_genes>600);
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
source = source(validcells);

tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
tot_genes = sum(data>0);
% figure('color','w','position',[100,100,325,700]); 
% subplot(3,1,1)
% hist(tot_mol,50);
% ylabel('#detected molecules');
% axis tight
% subplot(3,1,2)
% hist(tot_genes,50);
% ylabel('#detected genes');
% axis tight
% subplot(3,1,3)
% plot(tot_mol,tot_genes,'or'); hold on;
% axis tight
% xlabel('#detected molecules');
% ylabel('#detected genes');

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
%         plot(tmp1,tmp2,'or','markerfacecolor','r');hold on;
%         plot(tmp1(tmp2>thtmp2 & tmp1>thtmp1),tmp2(tmp2>thtmp2 & tmp1>thtmp1),'sg','markerfacecolor','g');hold on;
        inrmv(tmp2>thtmp2 & tmp1>thtmp1) = true;
%         plot(thtmp1*[1,1],[0,max(tmp2)],'k')
%         plot([0,max(tmp1)],thtmp2*[1,1],'k')
%         axis tight
%         xlabel(marker_genes{i})
%         ylabel(marker_genes{j})
%         title([marker_genes{i},'(pos)',num2str(sum(tmp1>thtmp1)),', ',marker_genes{j},'(pos)',num2str(sum(tmp2>thtmp2)),', both=',num2str(sum(tmp2>thtmp2 & tmp1>thtmp1))]);
%         set(gca,'fontsize',6)
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % 
validcells = ~inrmv;
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
source = source(validcells);

% % % % % % % % % % % % % % % % % % % % % % % % % fix sex by expression
sum_m_g = sum([data(strcmpi(geneid,'ddx3y'),:);data(strcmpi(geneid,'uty'),:);data(strcmpi(geneid,'eif2s3y'),:)]);
sum_f_g = sum([data(strcmpi(geneid,'xist'),:);data(strcmpi(geneid,'tsix'),:)]);
in_f = sum_f_g > sum_m_g;
in_m = sum_f_g < sum_m_g;
in_uk = sum_f_g==0 & sum_m_g==0;%probably male
sexcalc = zeros(length(cellid),1);
sexcalc(in_f) = 1;
sexcalc(in_m) = -1;
sexcalc(in_uk) = 0;
age = zeros(size(sexcalc));
age(source==1 & sexcalc==1) = 12;
age(source==1 & sexcalc==-1) = 35;
age(source==1 & sexcalc==0) = 30;%probably, not sure
age(source==2 & sexcalc==1) = 16;
age(source==2 & sexcalc==-1) = 24;
age(source==2 & sexcalc==0) = 22;%probably, not sure
% % % % % % % % % % % % % % % % % % % % % % % % % % 

in = find(sum(data>0,2)>20 & sum(data>0,2)<length(data(1,:))*0.6);
m_v = mean(data(in,:),2);
cv_v = std(data(in,:),[],2)./m_v;
log2_m = log2(m_v);
log2_cv = log2(cv_v);
x0 = [-0.5,1];
% [param_fit,fval] =  run_min_cvfit(log2_m,log2_cv,x0);
param_fit = robustfit(log2_m,log2_cv);
param_fit = [param_fit(2), param_fit(1)];
param_fit = round(param_fit*100)/100;
log2_cv_fit = log2( (2.^log2_m).^param_fit(1)) + param_fit(2);
tmp = log2_cv - log2_cv_fit;
[~,xi1] = sort(tmp,'descend');
corr_filt = xi1(1:12000);
figure; plot(log2_m,log2_cv,'.','markersize',3); hold on;
plot(log2_m(corr_filt),log2_cv(corr_filt),'.r','markersize',3); hold on;


sex_genes = upper({'EHD2', 'ESPL1', 'JARID1D', 'PNPLA4', 'RPS4Y1', 'XIST','tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'});
[~,loc] = ismember(sex_genes,upper(geneid(in)));
corr_filt = setdiff(corr_filt,loc);
stress_genes = upper({'Rpl26','Gstp1','Rpl35a','Erh','Slc25a5','Pgk1','Eno1','Tubb2a','Emc4','Scg5'});
[~,loc] = ismember(stress_genes,upper(geneid(in)));
corr_filt = setdiff(corr_filt,loc);

data = data(in(corr_filt),:);
geneid = geneid(in(corr_filt));
geneid_orig = geneid;
cellid_orig = cellid;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% permute the cells order
rng(1);
cell_perm = randperm(length(cellid));
gene_perm = randperm(length(geneid));
geneid = geneid(gene_perm);
cellid = cellid(cell_perm);
source = source(cell_perm);
age = age(cell_perm);
moldata = data(gene_perm,cell_perm);
moldata = moldata./repmat(sum(moldata),length(data(:,1)),1)*1e4;


cellid_cluster = loadCellFile('/data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
cellid_cluster = cellid_cluster(2:end,:);
[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid = cellid(loc);
source = source(loc);
age = age(loc);
moldata = moldata(:,loc);


grlev2 = cellid_cluster(:,3);
% lev2uni = unique(lev2gr);
grlev2_uni = {'gran_sox11','gran_fxyd7','gran_mature'};
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

age_uni = unique(age);
age_bar = zeros(length(age_uni),length(age));
for i=1:length(age_uni)
    age_bar(i,age==age_uni(i)) = 1;
end


in = find(sum(moldata>0,2)>20 & sum(moldata>0,2)<length(moldata(1,:))*0.6);
m_v = mean(moldata(in,:),2);
cv_v = std(moldata(in,:),[],2)./m_v;
log2_m = log2(m_v);
log2_cv = log2(cv_v);
x0 = [-0.5,1];
% [param_fit,fval] =  run_min_cvfit(log2_m,log2_cv,x0);
param_fit = robustfit(log2_m,log2_cv);
param_fit = [param_fit(2), param_fit(1)];
param_fit = round(param_fit*100)/100;
log2_cv_fit = log2( (2.^log2_m).^param_fit(1)) + param_fit(2);
tmp = log2_cv - log2_cv_fit;
[~,xi1] = sort(tmp,'descend');
corr_filt = xi1(1:5000);
figure; plot(log2_m,log2_cv,'.','markersize',3); hold on;
plot(log2_m(corr_filt),log2_cv(corr_filt),'.r','markersize',3); hold on;
moldata = moldata(in(corr_filt),:);
geneid = geneid(in(corr_filt));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
TF_list = loadCellFile('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/list_of_TF_mm');
tf = ismember(geneid,TF_list);
load /data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/mm_unigo_genes_cell.mat

% % % cell cycle
% k = find(~cellfun(@isempty ,strfind(lower(unigo_genes_cell(:,2)),'cell cycle')) | ...
%     ~cellfun(@isempty ,strfind(lower(unigo_genes_cell(:,2)),'mitosis')) | ...
%     ~cellfun(@isempty ,strfind(lower(unigo_genes_cell(:,2)),'cell division')));
proc_names = {'G1/S transition of mitotic cell cycle'
    'cell cycle arrest'
    'regulation of cell cycle'};
tmp = [];
for i=1:length(proc_names)
    tmp = [tmp; unigo_genes_cell{ strcmpi(unigo_genes_cell(:,2),proc_names{i}), 3}];
end
cc_genes = tmp;

cc_genes = loadCellFile('/data2/C1_stuff/DentGyr/goatools_guide_cellcycle_genes_mouse.txt');
cc_genes = unique([cc_genes(:,1);tmp]);
tf_cc = ismember(geneid,cc_genes);

log_m_gr_lev2 = zeros(length(geneid),length(T_cells_tmp_uni));
se_log_m_gr_lev2 = zeros(length(geneid),length(T_cells_tmp_uni));
log_m_gr_lev2_norm = zeros(length(geneid),length(T_cells_tmp_uni));
se_log_m_gr_lev2_norm = zeros(length(geneid),length(T_cells_tmp_uni));
present_mean = zeros(length(geneid),length(T_cells_tmp_uni));
moldata_1_norm = moldata./repmat(sum(moldata),length(geneid),1)*1e4;
for i=1:length(T_cells_tmp_uni)
    i
    in = T_cells_tmp==T_cells_tmp_uni(i);
    log_m_gr_lev2(:,i) = mean(log2(moldata(:,in)+1),2);
    se_log_m_gr_lev2(:,i) = std(log2(moldata(:,in)+1),[],2)/sqrt(sum(in));
    log_m_gr_lev2_norm(:,i) = mean(log2(moldata_1_norm(:,in)+1),2);
    se_log_m_gr_lev2_norm(:,i) = std(log2(moldata_1_norm(:,in)+1),[],2)/sqrt(sum(in));
    present_mean(:,i) = mean(moldata(:,in)>0,2);
end

clust_to_comp = {'gran_sox11','gran_fxyd7','gran_mature'};
pval_mat = zeros(length(geneid),length(clust_to_comp)*(length(clust_to_comp)-1));
qval_mat = zeros(length(geneid),length(clust_to_comp)*(length(clust_to_comp)-1));
comp_title = cell(1,length(clust_to_comp)*(length(clust_to_comp)-1));
p_rs1 = ones(length(geneid),length(clust_to_comp)-1);
for i=2:length(clust_to_comp)
    i
    gr1 = strcmpi(grlev2,clust_to_comp{i-1});
    gr2 = strcmpi(grlev2,clust_to_comp{i});
    for k=1:length(geneid);
        p_rs1(k,i-1) = qval_from_pval(ranksum(moldata(k,gr1),moldata(k,gr2)));
    end
end

fold_enrich = log_m_gr_lev2./repmat(mean(log_m_gr_lev2,2),1,length(log_m_gr_lev2(1,:)));


insig = find( min(p_rs1,[],2)<0.01 & max(fold_enrich,[],2)>1.5 & max(present_mean,[],2)>0.05 & (fold_enrich(:,1)<0.4 | fold_enrich(:,3)<0.4) & tf );

[~,itmax] = max(fold_enrich(insig,:),[],2);
fcatmax = zeros(length(insig),1);
for i=1:length(insig)
    fcatmax(i) = fold_enrich(insig(i),itmax(i));
end

[~,xi] = sortrows([itmax,-fcatmax]);
figure;
set(gcf,'position',[100,100,800,1000],'color','w')
axes('position',[0.1,0.05,0.88,0.95]);
imagesc(fold_enrich(insig(xi),:));
set(gca,'ytick',[1:length(insig)], 'yticklabel',geneid(insig(xi)),'xtick',[1:length(grlev2_uni)],'xticklabel',regexprep(grlev2_uni,'_','-'),'fontsize',8)
colormap('summer');
freezeColors(gca);
eval(['export_fig Late_Ngenesis_peakTime_TF_',date,'.pdf']);
table_TFs = [geneid(insig(xi)),grlev2_uni(itmax(xi))'];
saveCellFile(table_TFs,['Late_Ngenesis_peakTime_TF_',date,'.txt'])


insig = find( min(p_rs1,[],2)<0.01 & max(fold_enrich,[],2)>1.5 & max(present_mean,[],2)>0.05 & (fold_enrich(:,1)<0.4 | fold_enrich(:,3)<0.4) & ~tf & ~tf_cc);

[~,itmax] = max(fold_enrich(insig,:),[],2);
fcatmax = zeros(length(insig),1);
for i=1:length(insig)
    fcatmax(i) = fold_enrich(insig(i),itmax(i));
end

[~,xi] = sortrows([itmax,-fcatmax]);
figure;
set(gcf,'position',[100,100,800,1000],'color','w')
axes('position',[0.1,0.05,0.88,0.95]);
imagesc(fold_enrich(insig(xi),:));
set(gca,'ytick',[1:length(insig)], 'yticklabel',geneid(insig(xi)),'xtick',[1:length(grlev2_uni)],'xticklabel',regexprep(grlev2_uni,'_','-'),'fontsize',8)
colormap('summer');
freezeColors(gca);
eval(['export_fig Late_Ngenesis_peakTime_nonTFnonCC_',date,'.pdf']);
table_nonTFnonCC = [geneid(insig(xi)),grlev2_uni(itmax(xi))'];
saveCellFile(table_nonTFnonCC,['Late_Ngenesis_peakTime_nonTFnonCC_',date,'.txt'])


insig = find( min(p_rs1,[],2)<0.01 & max(fold_enrich,[],2)>1.5 & max(present_mean,[],2)>0.05 & (fold_enrich(:,1)<0.4 | fold_enrich(:,3)<0.4) & tf_cc & ~tf);

[~,itmax] = max(fold_enrich(insig,:),[],2);
fcatmax = zeros(length(insig),1);
for i=1:length(insig)
    fcatmax(i) = fold_enrich(insig(i),itmax(i));
end

[~,xi] = sortrows([itmax,-fcatmax]);
figure;
set(gcf,'position',[100,100,800,1000],'color','w')
axes('position',[0.1,0.05,0.88,0.95]);
imagesc(fold_enrich(insig(xi),:));
set(gca,'ytick',[1:length(insig)], 'yticklabel',geneid(insig(xi)),'xtick',[1:length(grlev2_uni)],'xticklabel',regexprep(grlev2_uni,'_','-'),'fontsize',8)
colormap('summer');
freezeColors(gca);
eval(['export_fig Late_Ngenesis_peakTime_CCgenes_',date,'.pdf']);
table_CCgenes = [geneid(insig(xi)),grlev2_uni(itmax(xi))'];
saveCellFile(table_CCgenes,['Late_Ngenesis_peakTime_CCgenes_',date,'.txt'])


insig = find( min(p_rs1,[],2)<0.01 & max(fold_enrich,[],2)>1.5 & max(present_mean,[],2)>0.05 & (fold_enrich(:,1)<0.4 | fold_enrich(:,3)<0.4));

[~,itmax] = max(fold_enrich(insig,:),[],2);
fcatmax = zeros(length(insig),1);
for i=1:length(insig)
    fcatmax(i) = fold_enrich(insig(i),itmax(i));
end

[~,xi] = sortrows([itmax,-fcatmax]);
xi = xi([find(itmax(xi)==1,20, 'first');find(itmax(xi)==2,20, 'first');find(itmax(xi)==3,20, 'first')]);
figure;
set(gcf,'position',[100,100,800,1000],'color','w')
axes('position',[0.1,0.05,0.88,0.95]);
imagesc(fold_enrich(insig(xi),:));
set(gca,'ytick',[1:length(insig)], 'yticklabel',geneid(insig(xi)),'xtick',[1:length(grlev2_uni)],'xticklabel',regexprep(grlev2_uni,'_','-'),'fontsize',8)
colormap('summer');
freezeColors(gca);
eval(['export_fig Late_Ngenesis_peakTime_Top20_genes_',date,'.pdf']);
table_ALL= [geneid(insig(xi)),grlev2_uni(itmax(xi))'];
saveCellFile(table_ALL,['Late_Ngenesis_peakTime_Top20_genes_',date,'.txt'])







% 
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% moldata = data;
% geneid = geneid_orig;
% cellid = cellid_orig;
% % % % % % % % % % % % % % 
% 
% cellid_cluster = loadCellFile('/data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
% cellid_cluster = cellid_cluster(2:end,:);
% [~,loc] = ismember(cellid_cluster(:,1), cellid);
% cellid = cellid(loc);
% moldata = moldata(:,loc);
% 
% 
% grlev2 = cellid_cluster(:,3);
% % lev2uni = unique(lev2gr);
% grlev2_uni = {'endo','peric','vlmc','mgl','pvm','OL','nfol','opc','astro','nsc1','gran_cycling','gran_eomes',...
%     'gran_sox11','gran_fxyd7','gran_mature','cck_tox','mossy_adcyap1','mossy_cy26b1','mossy_klk8','gaba_cnr1','gaba_lhx6','reln' };
% endo_gr = {'endo','peric','vlmc'};
% micro_gr = {'mgl','pvm'};
% oligo_gr = {'OL','nfol','opc'};
% granul_gr = {'gran_fxyd7','gran_mature'};
% mossy_gr = {'mossy_adcyap1','mossy_cy26b1','mossy_klk8','cck_tox'};
% gaba_gr = {'gaba_cnr1','gaba_lhx6'};
% grlev2(ismember(grlev2,endo_gr)) = repmat({'endo'},sum(ismember(grlev2,endo_gr)),1);
% grlev2(ismember(grlev2,micro_gr)) = repmat({'micro'},sum(ismember(grlev2,micro_gr)),1);
% grlev2(ismember(grlev2,oligo_gr)) = repmat({'oligo'},sum(ismember(grlev2,oligo_gr)),1);
% grlev2(ismember(grlev2,granul_gr)) = repmat({'granul'},sum(ismember(grlev2,granul_gr)),1);
% grlev2(ismember(grlev2,mossy_gr)) = repmat({'mossy'},sum(ismember(grlev2,mossy_gr)),1);
% grlev2(ismember(grlev2,gaba_gr)) = repmat({'gaba'},sum(ismember(grlev2,gaba_gr)),1);
% 
% macro_gr = {'endo','micro','oligo','granul','mossy','gaba'};
% 
% grlev2_uni = {'endo','micro','oligo','granul','mossy','gaba','reln','astro','nsc1','gran_cycling','gran_eomes','gran_sox11' };
% 
% dend_order_names = regexprep(grlev2_uni,'_','-');
% T_cells_tmp = zeros(length(grlev2),1);
% for i=1:length(grlev2_uni)
%     T_cells_tmp( strcmpi(grlev2, grlev2_uni{i}) ) = i;
% end
% [T_cells_tmp,xi] = sort(T_cells_tmp);
% xi(T_cells_tmp==0) = [];
% T_cells_tmp(T_cells_tmp==0) = [];
% T_cells_tmp_uni = unique(T_cells_tmp);
% cellid = cellid(xi);
% moldata = moldata(:,xi);
% grlev2 = grlev2(xi);
% 
% 
% 
% log_m_gr_lev2 = zeros(length(geneid),length(T_cells_tmp_uni));
% se_log_m_gr_lev2 = zeros(length(geneid),length(T_cells_tmp_uni));
% log_m_gr_lev2_norm = zeros(length(geneid),length(T_cells_tmp_uni));
% se_log_m_gr_lev2_norm = zeros(length(geneid),length(T_cells_tmp_uni));
% present_mean = zeros(length(geneid),length(T_cells_tmp_uni));
% moldata_1_norm = moldata./repmat(sum(moldata),length(geneid),1)*1e4;
% for i=1:length(T_cells_tmp_uni)
%     i
%     in = T_cells_tmp==T_cells_tmp_uni(i);
%     log_m_gr_lev2(:,i) = mean(log2(moldata(:,in)+1),2);
%     se_log_m_gr_lev2(:,i) = std(log2(moldata(:,in)+1),[],2)/sqrt(sum(in));
%     log_m_gr_lev2_norm(:,i) = mean(log2(moldata_1_norm(:,in)+1),2);
%     se_log_m_gr_lev2_norm(:,i) = std(log2(moldata_1_norm(:,in)+1),[],2)/sqrt(sum(in));
%     present_mean(:,i) = mean(moldata(:,in)>0,2);
% end
% 
% [~,inds] = ismember(table_TFs(:,1),geneid);
% figure;
% set(gcf,'position',[100,100,400,1000],'color','w')
% axes('position',[0.2,0.05,0.75,0.95]);
% imagesc(present_mean(inds,1:7));
% set(gca,'ytick',[1:length(inds)], 'yticklabel',table_TFs(:,1),'xtick',[1:length(grlev2_uni(1:7))],'xticklabel',regexprep(grlev2_uni(1:7),'_','-'),'fontsize',8,'XTickLabelRotation',45)
% set(gcf,'colormap',flipud(colormap('gray')))
% freezeColors(gca);
% colorbar
% eval(['export_fig Late_Ngenesis_fractioMacroclusters_TFs_',date,'.pdf']);
% 
% [~,inds] = ismember(table_nonTFnonCC(:,1),geneid);
% figure;
% set(gcf,'position',[100,100,400,1000],'color','w')
% axes('position',[0.2,0.05,0.75,0.95]);
% imagesc(present_mean(inds,1:7));
% set(gca,'ytick',[1:length(inds)], 'yticklabel',table_nonTFnonCC(:,1),'xtick',[1:length(grlev2_uni(1:7))],'xticklabel',regexprep(grlev2_uni(1:7),'_','-'),'fontsize',8,'XTickLabelRotation',45)
% set(gcf,'colormap',flipud(colormap('gray')))
% freezeColors(gca);
% colorbar
% eval(['export_fig Late_Ngenesis_fractioMacroclusters_nonTFnonCC_',date,'.pdf']);
% 
% 
% 
% [~,inds] = ismember(table_CCgenes(:,1),geneid);
% figure;
% set(gcf,'position',[100,100,400,1000],'color','w')
% axes('position',[0.2,0.05,0.75,0.95]);
% imagesc(present_mean(inds,1:7));
% set(gca,'ytick',[1:length(inds)], 'yticklabel',table_CCgenes(:,1),'xtick',[1:length(grlev2_uni(1:7))],'xticklabel',regexprep(grlev2_uni(1:7),'_','-'),'fontsize',8,'XTickLabelRotation',45)
% set(gcf,'colormap',flipud(colormap('gray')))
% freezeColors(gca);
% colorbar
% eval(['export_fig Late_Ngenesis_fractioMacroclusters_CC_',date,'.pdf']);
% 
% 












toc