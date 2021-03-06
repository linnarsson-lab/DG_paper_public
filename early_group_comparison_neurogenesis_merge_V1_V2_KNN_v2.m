tic
clear all
close all


load 10X43_1_v1.mat
load 10X46_1_v1.mat

cellid_43_1 = cellfun(@(x) ['10X43_1_',x,'-1'], cellid_43_1,'uniformoutput',0);
cellid_46_1 = cellfun(@(x) ['10X46_1_',x,'-1'], cellid_46_1,'uniformoutput',0);




load /mnt/sanger-data2/C1_stuff/DentGyr/DG_V2kit_subsample_samples_merged_18-Jun-2017.mat
% load /data2/C1_stuff/DentGyr/DG_V2kit_subsample_samples_merged_18-Jun-2017.mat

[geneid,ia,ib] = intersect(geneid,geneid_46_1);
data = [data(ia,:), data_43_1(ib,:), data_46_1(ib,:)];
data1 = data(ia,:); data2= [ data_43_1(ib,:), data_46_1(ib,:)]; 
% dm = mean(log2(data1+1),2) - mean(log2(data2+1),2);
% data(abs(dm)>0.5,:) = [];
% geneid(abs(dm)>0.5) = [];

source = [source; repmat({'10X43_1'},length(cellid_43_1),1); repmat({'10X46_1'},length(cellid_46_1),1)];
cellid = [cellid;cellid_43_1;cellid_46_1];
age = [age;cell(length([cellid_43_1;cellid_46_1]),1)];
sex = [sex;cell(length([cellid_43_1;cellid_46_1]),1)];
strain = [strain;repmat({'CD-1'},length([cellid_43_1;cellid_46_1]),1)];
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

age(strcmpi(source,'10X43_1') & sexcalc==1) = repmat({'p12'}, sum(strcmpi(source,'10X43_1') & sexcalc==1),1);
age(strcmpi(source,'10X43_1') & sexcalc==-1) = repmat({'p35'}, sum(strcmpi(source,'10X43_1') & sexcalc==-1),1);
age(strcmpi(source,'10X43_1') & sexcalc==0) = repmat({'p30'}, sum(strcmpi(source,'10X43_1') & sexcalc==0),1);%probably, not sure
age(strcmpi(source,'10X46_1') & sexcalc==1) = repmat({'p16'}, sum(strcmpi(source,'10X46_1') & sexcalc==1),1);
age(strcmpi(source,'10X46_1') & sexcalc==-1) = repmat({'p24'}, sum(strcmpi(source,'10X46_1') & sexcalc==-1),1);;
age(strcmpi(source,'10X46_1') & sexcalc==0) = repmat({'p22'}, sum(strcmpi(source,'10X46_1') & sexcalc==0),1);%probably, not sure

sex(strcmpi(source,'10X43_1') & sexcalc==1) = repmat({'F'}, sum(strcmpi(source,'10X43_1') & sexcalc==1),1);
sex(strcmpi(source,'10X43_1') & sexcalc==-1) = repmat({'M'}, sum(strcmpi(source,'10X43_1') & sexcalc==-1),1);
sex(strcmpi(source,'10X43_1') & sexcalc==0) = repmat({'M'}, sum(strcmpi(source,'10X43_1') & sexcalc==0),1);%probably, not sure
sex(strcmpi(source,'10X46_1') & sexcalc==1) = repmat({'F'}, sum(strcmpi(source,'10X46_1') & sexcalc==1),1);
sex(strcmpi(source,'10X46_1') & sexcalc==-1) = repmat({'M'}, sum(strcmpi(source,'10X46_1') & sexcalc==-1),1);;
sex(strcmpi(source,'10X46_1') & sexcalc==0) = repmat({'M'}, sum(strcmpi(source,'10X46_1') & sexcalc==0),1);%probably, not sure

% % % % % % % % % % % % % % % % % % % % % % % % % % 
tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
tot_genes = sum(data>0);

% validcells = (tot_mol>800 & (tot_mol./tot_genes)>1.2 & tot_mol<3e4 & tot_genes>800);
% sum(validcells)
% data = data(:,validcells);
% cellid = cellid(validcells);
% source = source(validcells);
% age = age(validcells);
% sex = sex(validcells);
% strain = strain(validcells);

age = regexprep(age,'p','');
age = regexprep(age,'P','');
age = regexprep(age,'E','');
age = cellfun(@str2double, age);
age(age==16.5) = -16.5;

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);

marker_genes = {'Stmn2','Mog','Aldoc','C1qc','Cldn5'};

% inrmv = false(length(cellid),1);
% for i=1:length(marker_genes)
%     for j=i+1:length(marker_genes)
% %         subplot(4,4,(i-1)*4+j-1);
%         tmp1 = data(strcmpi(geneid,marker_genes{i}),:);
%         tmp2 = data(strcmpi(geneid,marker_genes{j}),:);
%         thtmp1 = 1;%prctile(tmp1(tmp1>0),20);
%         thtmp2 = 1;%prctile(tmp2(tmp2>0),20);
%         inrmv(tmp2>thtmp2 & tmp1>thtmp1) = true;
%     end
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% validcells = ~inrmv;
% sum(validcells)
% data = data(:,validcells);
% cellid = cellid(validcells);
% source = source(validcells);
% age = age(validcells);

data = round(data./repmat(sum(data),length(data(:,1)),1)*5e3);
moldata = data;


cellid_cluster = loadCellFile('mnt/sanger-data2/C1_stuff/DentGyr/DG_10X_V2_markertable_knnMutal_MCL_1p4_18-Jun-2017_manualinspection.txt');
% cellid_cluster = loadCellFile('/data2/C1_stuff/DentGyr/DG_10X_V2_markertable_knnMutal_MCL_1p4_18-Jun-2017_manualinspection.txt');
cellid_cluster = cellid_cluster(2:end,:);
cellid_cluster2 = loadCellFile('mnt/sanger-data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
% cellid_cluster2 = loadCellFile('/data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
cellid_cluster2 = cellid_cluster2(2:end,:);
cellid_cluster = [cellid_cluster;cellid_cluster2];

[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid_cluster(loc==0,:) = [];
loc(loc==0) = [];
cellid = cellid(loc);
source = source(loc);
age = age(loc);
moldata = moldata(:,loc);
data_full = data(:,loc);

grlev2 = cellid_cluster(:,3);
% lev2uni = unique(lev2gr);
% grlev2_uni = {'rgl','astro_old','astro_adol','immature_astro','rgl_young','nipc_young','nipc','nb1','nb2'.....
%     'astro','nsc1','gran_cycling','gran_eomes','gran_sox11'};
% grlev2_uni = {'rgl','astro_old','astro_adol','immature_astro','rgl_young','nipc_young','nipc'.....
%     'astro','nsc1','gran_cycling'};
grlev2_uni = {'astro','astro_adol','astro_old','rgl','nsc1','immature_astro',....
    'rgl_young','nipc_young','nipc','gran_cycling','gran_eomes','nb1','gran_sox11','nb2'};

tf = ismember(grlev2,grlev2_uni);
[~,p] = ttest2(log2(moldata(:,tf)+1)', log2(moldata(:,~tf)+1)','tail','right');
in_leave = fdr_proc(p,0.2);
% % % % % % % % % % % % % % % % % % % % % % % % 



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
data_full = data_full(:,xi);
geneid_full = geneid;

age_uni = unique(age);
age_bar = zeros(length(age_uni),length(age));
for i=1:length(age_uni)
    age_bar(i,age==age_uni(i)) = 1;
end

moldata = moldata(in_leave,:);
geneid = geneid(in_leave);

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
corr_filt = xi1(1:3000);
figure; plot(log2_m,log2_cv,'.','markersize',3); hold on;
plot(log2_m(corr_filt),log2_cv(corr_filt),'.r','markersize',3); hold on;
moldata = moldata(in(corr_filt),:);
geneid = geneid(in(corr_filt));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% sumpergene = sum(dataout_sorted_all_nonn,2);
meanpergene = mean(moldata,2);%
molenrich_mat = zeros(length(moldata(:,1)),length(T_cells_tmp_uni));
% meangr_mat = zeros(length(moldata(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(moldata(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(moldata(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(moldata(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
[~,grord] = sort(gr_center);
meangrpos_mat = meangrpos_mat(:,grord);
molenrich_mat = molenrich_mat(:,grord);
[~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0.00001,'descend');
[~,xi0p5] = sort(molenrich_mat.*meangrpos_mat.^0.5,'descend');
[~,xi1] = sort(molenrich_mat.*meangrpos_mat.^1,'descend');


ind_gr_tmp_mark = [xi0(1:20,:);xi0p5(1:20,:);xi1(1:20,:)];
ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
rmv = false(size(ind_gr_tmp_mark));
for i=2:length(ind_gr_tmp_mark)
    if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
        rmv(i) = true;
    end
end
ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = geneid(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));

molenrich_mat_mark = zeros(length(ind_gr_tmp_mark),length(T_cells_tmp_uni));
meangrpos_mat_mark = zeros(length(ind_gr_tmp_mark),length(T_cells_tmp_uni));
order_gr = [T_cells_tmp(diff(T_cells_tmp)~=0);T_cells_tmp(end)];
for jjj=1:length(T_cells_tmp_uni)
    jjj
    molenrich_mat_mark(:,jjj) = mean(moldata(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj)),2)./meanpergene(ind_gr_tmp_mark);
    meangrpos_mat_mark(:,jjj) = mean(moldata(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj))>0,2);
end
[~,imax] = max(molenrich_mat_mark.*meangrpos_mat_mark.^0.5,[],2);
[~,xi] = sort(imax);
ind_gr_tmp_mark = ind_gr_tmp_mark(xi);

datamarkers = moldata(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
gr_tmp_mark = gr_tmp_mark(xi);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

gr_to_comp = [{'astro_adol','astro_old'}; {'astro','astro_adol'}; {'rgl','nsc1'}; ...
    {'nipc_young','nipc'}; {'rgl_young','rgl'}; {'immature_astro','astro_adol'};....
    {'rgl','astro_adol'}; {'nsc1','astro'}];
qval_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
fc_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(length(gr_to_comp),1);
for i=1:length(gr_to_comp)
    i
    in1 = find(strcmpi(grlev2,gr_to_comp{i,1}));
    in2 = find(strcmpi(grlev2,gr_to_comp{i,2}));
    compname{i} = [gr_to_comp{i,1},' vs ',gr_to_comp{i,2}];
    for j=1:length(ind_gr_tmp_mark)
        qval_table(j,i) = ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
        fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
    end
    qval_table(:,i) = qval_from_pval( qval_table(:,i));
    axes(ha(i))
    plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
    [~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
    top_g = 40;
    tmp = xi(1:top_g/2);
    [~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
    xi = [tmp;xi(1:top_g/2)];
    
    plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
    text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
    title(regexprep(compname{i},'_','-'));
    axis tight
end
eval(['export_fig group_comp1_volcano_ranksum_v1_v2_early_Ngenesis_',date,'.pdf']);
table1 = [qval_table, fc_table];
tmp = [1:2*length(gr_to_comp)];
tmp = reshape(tmp,length(gr_to_comp),2)';
tmp = tmp(:);
table1 = table1(:,tmp);
table1 = [cell(1,1+2*length(gr_to_comp)) ; [gr_tmp_mark, m2c(table1)] ];
tmp = [compname';compname'];
tmp = tmp(:);
table1(1,2:end) = tmp';

saveCellFile(table1,['group_comp1_ranksum_v1_v2_early_Ngenesis_',date,'.txt']);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(8,1);

qval_table = zeros(length(ind_gr_tmp_mark),8);
fc_table = zeros(length(ind_gr_tmp_mark),8);
i=1;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl_young') & age==5);
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl-young E16 vs p5')
axis tight
compname{i} = 'rgl-young E16 vs p5';
i=2;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl-young E16 vs rgl p18-23')
axis tight
compname{i} = 'rgl-young E16 vs rgl p18-23';
i=3;
in1 = find(strcmpi(grlev2,'rgl_young') & (age==-16.5 | age==0 | age==5));
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl-young E16-p5 vs rgl p18-23')
axis tight
compname{i} = 'rgl-young E16-p5 vs rgl p18-23';
i=4;
in1 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'rgl') & (age==120 | age==132));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl p18-23 vs rgl p120-132')
axis tight
compname{i} = 'rgl rgl p18-23 vs p120-132';

i=5;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc_young') & (age==5));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) = ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc-young E16.5 vs p5')
axis tight
compname{i} = 'nipc-young E16.5 vs p5';
i=6;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc-young E16.5 vs nipc p18-23')
axis tight
compname{i} = 'nipc-young E16.5 vs nipc p18-23';
i=7;
in1 = find(strcmpi(grlev2,'nipc_young') & (age<=5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) = ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc-young <p5 vs nipc p18-23')
axis tight
compname{i} = 'nipc-young <p5 vs nipc p18-23';
i=8;
in1 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'nipc') & (age>=120));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc p18-23 vs nipc p120-132')
axis tight
compname{i} = 'nipc p18-23 vs nipc p120-132';
eval(['export_fig group_comp2_volcano_ranksum_v1_v2_early_Ngenesis_',date,'.pdf']);
table1 = [qval_table, fc_table];
tmp = [1:2*length(gr_to_comp)];
tmp = reshape(tmp,length(gr_to_comp),2)';
tmp = tmp(:);
table1 = table1(:,tmp);
table1 = [cell(1,1+2*length(gr_to_comp)) ; [gr_tmp_mark, m2c(table1)] ];
tmp = [compname';compname'];
tmp = tmp(:);
table1(1,2:end) = tmp';

saveCellFile(table1,['group_comp2_ranksum_v1_v2_early_Ngenesis_',date,'.txt']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % scatter plots
gr_to_comp = [{'astro_adol','astro_old'}; {'astro','astro_adol'}; {'rgl','nsc1'}; ...
    {'nipc_young','nipc'}; {'rgl_young','rgl'}; {'immature_astro','astro_adol'};....
    {'rgl','astro_adol'}; {'nsc1','astro'}];
qval_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
fc_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(length(gr_to_comp),1);
for i=1:length(gr_to_comp)
    i
    in1 = find(strcmpi(grlev2,gr_to_comp{i,1}));
    in2 = find(strcmpi(grlev2,gr_to_comp{i,2}));
    compname{i} = [gr_to_comp{i,1},' vs ',gr_to_comp{i,2}];
    for j=1:length(ind_gr_tmp_mark)
        qval_table(j,i) = ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
        fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
    end
    qval_table(:,i) = qval_from_pval( qval_table(:,i));
    axes(ha(i))
    plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
%     plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
    [~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
    top_g = 40;
    tmp = xi(1:top_g/2);
    [~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
    xi = [tmp;xi(1:top_g/2)];
    
    plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
    text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
    sig = qval_table(:,i)<0.001;
    plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
    axis tight
    xl = get(gca,'xlim'); yl = get(gca,'ylim');
    plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
    set(gca,'xlim',xl,'ylim',yl)    
    title(regexprep(compname{i},'_','-'));
end
eval(['export_fig group_comp1_scatter_ranksum_v1_v2_early_Ngenesis_',date,'.pdf']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(8,1);

qval_table = zeros(length(ind_gr_tmp_mark),8);
fc_table = zeros(length(ind_gr_tmp_mark),8);
i=1;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl_young') & age==5);
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('rgl-young E16 vs p5')
compname{i} = 'rgl-young E16 vs p5';
i=2;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('rgl-young E16 vs rgl p18-23')
compname{i} = 'rgl-young E16 vs rgl p18-23';
i=3;
in1 = find(strcmpi(grlev2,'rgl_young') & (age==-16.5 | age==0 | age==5));
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('rgl-young E16-p5 vs rgl p18-23')
compname{i} = 'rgl-young E16-p5 vs rgl p18-23';
i=4;
in1 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'rgl') & (age==120 | age==132));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('rgl p18-23 vs rgl p120-132')
compname{i} = 'rgl p18-23 vs rgl p120-132';

i=5;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc_young') & (age==5));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) = ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('nipc-young E16.5 vs p5')
compname{i} = 'nipc-young E16.5 vs p5';
i=6;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('nipc-young E16.5 vs nipc p18-23')
compname{i} = 'nipc-young E16.5 vs nipc p18-23';
i=7;
in1 = find(strcmpi(grlev2,'nipc_young') & (age<=5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) = ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('nipc-young <p5 vs nipc p18-23')
compname{i} = 'nipc-young <p5 vs nipc p18-23';
i=8;
in1 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'nipc') & (age>=120));
for j=1:length(ind_gr_tmp_mark)
    qval_table(j,i) =  ranksum(moldata(ind_gr_tmp_mark(j),in1), moldata(ind_gr_tmp_mark(j),in2)) ;
    fc_table(j,i) = mean(log2(moldata(ind_gr_tmp_mark(j),in1)+1)) - mean(log2(moldata(ind_gr_tmp_mark(j),in2)+1));
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(log2(moldata(ind_gr_tmp_mark,in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark,in2)+1),2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),'.r'); hold on;
text(mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)+1),2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(log2(moldata(ind_gr_tmp_mark(sig),in1)+1),2), mean(log2(moldata(ind_gr_tmp_mark(sig),in2)+1),2),'.g'); hold on;
axis tight
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([0,1]*xl(2),[0,1]*xl(2),'k-');plot([0,1]*xl(2),1+[0,1]*xl(2),'k--');plot([0,1]*xl(2),-1+[0,1]*xl(2),'k--');
set(gca,'xlim',xl,'ylim',yl)
title('nipc p18-23 vs nipc p120-132')
compname{i} = 'nipc p18-23 vs nipc p120-132';
eval(['export_fig group_comp2_scatter_ranksum_v1_v2_early_Ngenesis_',date,'.pdf']);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % fraction of cells
% test with binomial difference normal approximation
gr_to_comp = [{'astro_adol','astro_old'}; {'astro','astro_adol'}; {'rgl','nsc1'}; ...
    {'nipc_young','nipc'}; {'rgl_young','rgl'}; {'immature_astro','astro_adol'};....
    {'rgl','astro_adol'}; {'nsc1','astro'}];
qval_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
fc_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(length(gr_to_comp),1);
for i=1:length(gr_to_comp)
    i
    in1 = find(strcmpi(grlev2,gr_to_comp{i,1}));
    in2 = find(strcmpi(grlev2,gr_to_comp{i,2}));
    compname{i} = [gr_to_comp{i,1},' vs ',gr_to_comp{i,2}];
    for j=1:length(ind_gr_tmp_mark)
        x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
        x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
        n1 = length(in1);
        n2 = length(in2);
        if (x1+x2)<0.01*min([length(in1), length(in2)])
            qval_table(j,i) = 1;
        else
            p = (x1+x2) / (n1+n2);
            SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
            z = (x1/n1-x2/n2)/SE;
            qval_table(j,i) = 2*normcdf(-abs(z));
        end
        fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
    end
    qval_table(:,i) = qval_from_pval(qval_table(:,i));
    axes(ha(i))
    plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
    [~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
    top_g = 40;
    tmp = xi(1:top_g/2);
    [~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
    xi = [tmp;xi(1:top_g/2)];
    plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
    text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
    title(regexprep(compname{i},'_','-'))
    axis tight
end
eval(['export_fig group_comp1_volcano_fracbino_v1_v2_early_Ngenesis_',date,'.pdf']);
table1 = [qval_table, fc_table];
tmp = [1:2*length(gr_to_comp)];
tmp = reshape(tmp,length(gr_to_comp),2)';
tmp = tmp(:);
table1 = table1(:,tmp);
table1 = [cell(1,1+2*length(gr_to_comp)) ; [gr_tmp_mark, m2c(table1)] ];
tmp = [compname';compname'];
tmp = tmp(:);
table1(1,2:end) = tmp';

saveCellFile(table1,['group_comp1_fracbino_v1_v2_early_Ngenesis_',date,'.txt']);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
gr_to_comp = [{'astro_adol','astro_old'}; {'astro','astro_adol'}; {'rgl','nsc1'}; ...
    {'nipc_young','nipc'}; {'rgl_young','rgl'}; {'immature_astro','astro_adol'};....
    {'rgl','astro_adol'}; {'nsc1','astro'}];
qval_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
fc_table = zeros(length(ind_gr_tmp_mark),length(gr_to_comp));
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(length(gr_to_comp),1);
for i=1:length(gr_to_comp)
    i
    in1 = find(strcmpi(grlev2,gr_to_comp{i,1}));
    in2 = find(strcmpi(grlev2,gr_to_comp{i,2}));
    compname{i} = [gr_to_comp{i,1},' vs ',gr_to_comp{i,2}];
    for j=1:length(ind_gr_tmp_mark)
        x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
        x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
        n1 = length(in1);
        n2 = length(in2);
        if (x1+x2)<0.01*min([length(in1), length(in2)])
            qval_table(j,i) = 1;
        else
            p = (x1+x2) / (n1+n2);
            SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
            z = (x1/n1-x2/n2)/SE;
            qval_table(j,i) = 2*normcdf(-abs(z));
        end
        fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
    end
    qval_table(:,i) = qval_from_pval(qval_table(:,i));
    axes(ha(i))
    plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
    [~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
    top_g = 40;
    tmp = xi(1:top_g/2);
    [~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
    xi = [tmp;xi(1:top_g/2)];
    
    plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
    text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
    sig = qval_table(:,i)<0.001;
    plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
    title(regexprep(compname{i},'_','-'))
    axis tight
end
eval(['export_fig group_comp1_scatter_fracbino_v1_v2_early_Ngenesis_',date,'.pdf']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(8,1);

qval_table = zeros(length(ind_gr_tmp_mark),8);
fc_table = zeros(length(ind_gr_tmp_mark),8);
i=1;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl_young') & age==5);
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl-young E16 vs p5')
axis tight
compname{i} = 'rgl-young E16 vs p5';
i=2;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl-young E16 vs rgl p18-23')
axis tight
compname{i} = 'rgl-young E16 vs rgl p18-23';
i=3;
in1 = find(strcmpi(grlev2,'rgl_young') & (age==-16.5 | age==0 | age==5));
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl-young E16-p5 vs rgl p18-23')
axis tight
compname{i} = 'rgl-young E16-p5 vs rgl p18-23';
i=4;
in1 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'rgl') & (age==120 | age==132));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('rgl p18-23 vs rgl p120-132')
axis tight
compname{i} = 'rgl p18-23 vs rgl p120-132';

i=5;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc_young') & (age==5));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc-young E16.5 vs p5')
axis tight
compname{i} = 'nipc-young E16.5 vs p5';
i=6;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc-young E16.5 vs nipc p18-23')
axis tight
compname{i} = 'nipc-young E16.5 vs nipc p18-23';
i=7;
in1 = find(strcmpi(grlev2,'nipc_young') & (age<=5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc-young <p5 vs nipc p18-23')
axis tight
compname{i} = 'nipc-young <p5 vs nipc p18-23';
i=8;
in1 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'nipc') & (age>=120));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(fc_table(:,i),-log10(qval_table(:,i)),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),'.r'); hold on;
text(fc_table(xi(1:top_g),i),-log10(qval_table(xi(1:top_g),i)),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
title('nipc p18-23vs nipc p120-132 ')
axis tight
compname{i} = 'nipc p18-23 vs nipc p120-132';
eval(['export_fig group_comp2_volcano_fracbino_v1_v2_early_Ngenesis_',date,'.pdf']);
table1 = [qval_table, fc_table];
tmp = [1:2*length(gr_to_comp)];
tmp = reshape(tmp,length(gr_to_comp),2)';
tmp = tmp(:);
table1 = table1(:,tmp);
table1 = [cell(1,1+2*length(gr_to_comp)) ; [gr_tmp_mark, m2c(table1)] ];
tmp = [compname';compname'];
tmp = tmp(:);
table1(1,2:end) = tmp';

saveCellFile(table1,['group_comp2_fracbino_v1_v2_early_Ngenesis_',date,'.txt']);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure('position',[100,100,1600,800],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.05], [0.07,0.07], [0.03,0.03]);
compname = cell(8,1);

qval_table = zeros(length(ind_gr_tmp_mark),8);
fc_table = zeros(length(ind_gr_tmp_mark),8);
i=1;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl_young') & age==5);
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('rgl-young E16 vs p5')
axis tight
compname{i} = 'rgl-young E16 vs p5';
i=2;
in1 = find(strcmpi(grlev2,'rgl_young') & age==-16.5);
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('rgl-young E16 vs rgl p18-23')
axis tight
compname{i} = 'rgl-young E16 vs rgl p18-23';
i=3;
in1 = find(strcmpi(grlev2,'rgl_young') & (age==-16.5 | age==0 | age==5));
in2 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('rgl-young E16-p5 vs rgl p18-23')
axis tight
compname{i} = 'rgl-young E16-p5 vs rgl p18-23';
i=4;
in1 = find(strcmpi(grlev2,'rgl') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'rgl') & (age==120 | age==132));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('rgl p18-23 vs rgl p120-132')
axis tight
compname{i} = 'rgl p18-23 vs rgl p120-132';

i=5;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc_young') & (age==5));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('nipc-young E16.5 vs p5')
axis tight
compname{i} = 'nipc-young E16.5 vs p5';
i=6;
in1 = find(strcmpi(grlev2,'nipc_young') & (age==-16.5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('nipc-young E16.5 vs nipc p18-23')
axis tight
compname{i} = 'nipc-young E16.5 vs nipc p18-23';
i=7;
in1 = find(strcmpi(grlev2,'nipc_young') & (age<=5));
in2 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('nipc-young <p5 vs nipc p18-23')
axis tight
compname{i} = 'nipc-young <p5 vs nipc p18-23';
i=8;
in1 = find(strcmpi(grlev2,'nipc') & (age==18 | age==19 | age==23));
in2 = find(strcmpi(grlev2,'nipc') & (age>=120));
for j=1:length(ind_gr_tmp_mark)
    x1 = sum(moldata(ind_gr_tmp_mark(j),in1)>0);
    x2 = sum(moldata(ind_gr_tmp_mark(j),in2)>0);
    n1 = length(in1);
    n2 = length(in2);
    if (x1+x2)<0.01*min([length(in1), length(in2)])
        qval_table(j,i) = 1;
    else
        p = (x1+x2) / (n1+n2);
        SE = sqrt( p * ( 1 - p ) * ((1/n1) + (1/n2) ) );
        z = (x1/n1-x2/n2)/SE;
        qval_table(j,i) = 2*normcdf(-abs(z));
    end
    fc_table(j,i) = mean(moldata(ind_gr_tmp_mark(j),in1)>0) - mean(moldata(ind_gr_tmp_mark(j),in2)>0);
end
qval_table(:,i) = qval_from_pval( qval_table(:,i));
axes(ha(i))
plot(mean(moldata(ind_gr_tmp_mark,in1)>0,2), mean(moldata(ind_gr_tmp_mark,in2)>0,2),'.'); hold on;
[~,xi] = sort(fc_table(:,i).*log10(qval_table(:,i)));
top_g = 40;
tmp = xi(1:top_g/2);
[~,xi] = sort(-fc_table(:,i).*log10(qval_table(:,i)));
xi = [tmp;xi(1:top_g/2)];
plot(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),'.r'); hold on;
text(mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in1)>0,2), mean(moldata(ind_gr_tmp_mark(xi(1:top_g)),in2)>0,2),gr_tmp_mark(xi(1:top_g)),'fontsize',6);
sig = qval_table(:,i)<0.001;
plot(mean(moldata(ind_gr_tmp_mark(sig),in1)>0,2), mean(moldata(ind_gr_tmp_mark(sig),in2)>0,2),'.g'); hold on;
title('nipc p18-23 vs nipc p120-132')
axis tight
compname{i} = 'nipc p18-23 vs nipc p120-132';
eval(['export_fig group_comp2_scatter_fracbino_v1_v2_early_Ngenesis_',date,'.pdf']);








