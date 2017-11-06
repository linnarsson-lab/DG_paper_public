
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
age = age(validcells);

data = round(data./repmat(sum(data),length(data(:,1)),1)*5e3);

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
stress_genes = upper({'Rpl26','Gstp1','Rpl35a','Erh','Slc25a5','Pgk1','Eno1','Tubb2a','Emc4','Scg5','Btg2','Fos','Klf6'....
    'Nr4a1','Dusp1','Egr1','Jun','Fosb','Dnajb1','Ier2','Junb','Csrnp1',});
[~,loc] = ismember(stress_genes,upper(geneid(in)));
corr_filt = setdiff(corr_filt,loc);

data = data(in(corr_filt),:);
geneid = geneid(in(corr_filt));

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

age_uni = unique(age);
age_bar = zeros(length(age_uni),length(age));
for i=1:length(age_uni)
    age_bar(i,age==age_uni(i)) = 1;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % %
% T_cells_tmp_uni = unique(T_cells_tmp);
% T_cells_tmp_uni(T_cells_tmp_uni==0) = [];


age_cluster_mat = zeros(length(age_uni),length(grlev2_uni));

for i=1:length(T_cells_tmp_uni)
    i
    in = T_cells_tmp==T_cells_tmp_uni(i);
    age_cluster_mat(:,i) = sum(age_bar(:,in),2);
end

age_cluster_mat_table = [ [{'age/cluster'};m2c(age_uni)],[dend_order_names;m2c(age_cluster_mat)] ];
tmp = age_cluster_mat./repmat(sum(age_cluster_mat),length(age_uni),1)*100;
tmp(tmp==0) = 0.001;
for i=1:length(T_cells_tmp_uni)
    dend_order_names{i} = [dend_order_names{i},'(',num2str(sum(age_cluster_mat(:,i))),')'];
end

figure;
set(gcf,'position',[100,100,900,600],'color','w')
axes('position',[0.2,0.13,0.75,0.85]);
x = repmat([1:length(grlev2_uni)],length(age_uni),1);
y = repmat([1:length(age_uni)]',1,length(grlev2_uni));
scatter(x(:), y(:), 8*tmp(:),'facecolor','b');
axis tight
set(gca,'xtick',[1:length(grlev2_uni)],'XTickLabel',dend_order_names,'ytick',[1:length(age_uni)],'YTickLabel',age_uni,'XTickLabelRotation',45,'ydir','reverse')
eval(['export_fig blobs_10X_V2kit_age_vs_clusters_',date,'.pdf']);

saveCellFile(age_cluster_mat_table,['age_cluster_mat_table_10X_V2kit_age_vs_clusters_',date,'txt'])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

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


mean_exp_gr_lev2 = zeros(length(geneid),length(T_cells_tmp_uni));
mean_exp_gr_lev2_norm = zeros(length(geneid),length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
    mean_exp_gr_lev2(:,i) = mean(moldata(:,in),2);
    mean_exp_gr_lev2_norm(:,i) = mean(moldata_1_norm(:,in),2);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % %

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

figure;
set(gcf,'position',[100,100,1000,900],'color','w')
ax1 = axes('position',[0.1,0.08,0.88,0.80]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
cells_bor = find(diff(T_cells_tmp)>0);
for jj=1:length(cells_bor)
    plot(cells_bor(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',gr_center,'xticklabel',dend_order_names,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 6,'XTickLabelRotation',45)
colormap('summer');
freezeColors(gca);

ax2 = axes('position',[0.1,0.88,0.88,0.12]);
imagesc(~age_bar); hold on;
for i=1:length(age_uni)-1
    plot([0,length(age_bar(1,:))],(i+0.5)*[1,1],'-k');
end
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1:length(age_uni)],'yticklabel',age_uni);
linkaxes([ax1,ax2],'x');

eval(['export_fig DG_V210X_markertable_finalClsuter_',date,'.pdf'])

table1 = [[{'cellid','cluster'},gr_tmp_mark'] ;[cellid,m2c([T_cells_tmp,datamarkers'])] ];
saveCellFile(table1,['DG_V210X_markertable_finalClsuter_',date,'.txt']);




% % % % % % % % % % % % % % % % % % % % % % % % % 

log_m_mat_link = log_m_gr_lev2_norm(ind_gr_tmp_mark,:);%log2(mean_exp_gr_lev2_norm(ind_gr_tmp_mark,:)+1);%
gene_link = gr_tmp_mark;
% Z = linkage(log_m_mat_link','ward','correlation');
% Rcelltypes = corr_mat(log_m_mat_link);
% D = pdist(log_m_mat_link','correlation');
% leaforder = optimalleaforder(Z,D);[1:length(grlev2_uni)]; %
% 
% Z = linkage(log_m_mat_link(:,leaforder)','ward','correlation');

Z = linkage(log_m_mat_link','ward','Euclidean');
Rcelltypes = corr_mat(log_m_mat_link);
D = pdist(log_m_mat_link','Euclidean');
leaforder = [1:length(grlev2_uni)]; %optimalleaforder(Z,D);

Z = linkage(log_m_mat_link(:,leaforder)','ward','Euclidean');

figure('color','w','position',[200,200,1450,800],'visible','on');
axes('position',[0.03,0.05,0.4,0.9])
[H,T,outperm] = dendrogram(Z,length(T_cells_tmp_uni),'Orientation','left','reorder',...
    [1:length(T_cells_tmp_uni)],'labels',regexprep(grlev2_uni(leaforder),'_','-'),'colorthreshold',20);
set(H,'linewidth',2)
xl = get(gca,'xlim');
set(gca,'ylim',[0.5,length(T_cells_tmp_uni)+0.5],'xtick',[],'ydir','normal','fontsize',8)%,'xlim',[0,1.1]
% set(gca,'xlim',[0.5,length(T_cells_tmp_uni)+1],'ytick',[],'fontsize',8)
% rotateticklabel(gca,45);

axes('position',[0.58,0.05,0.35,0.9])
imagesc(Rcelltypes(leaforder,leaforder),[0.1,1]) 
set(gca,'xtick',[],'xticklabel',leaforder,'xdir','normal','ydir','normal','ytick',[])
colormap('summer')
hcm = colorbar;
set(hcm,'position',[0.96,0.7,0.02,0.25],'ytick',[-0.4:0.1:1],'ydir','normal')
% eval(['export_fig oligos_dendogram_by_lev2_',date,'.pdf'])
save2pdf(['DG_V210X_dendogram_by_finalClsuter_',date,'.pdf'],gcf,300)

figure('color','w','position',[200,200,800,1000],'visible','on');
axes('position',[0.03,0.03,0.85,0.96])
[H,T,outperm] = dendrogram(Z,length(T_cells_tmp_uni),'Orientation','left','reorder',...
    [1:length(T_cells_tmp_uni)],'labels',regexprep(grlev2_uni(leaforder),'_','-'),'colorthreshold',0.2);
set(H,'linewidth',1)
for i=1:length(H)
    xd = get(H(i),'XData');
    yd = get(H(i),'YData');
    text(xd(2),mean(yd(2:3)),num2str(i),'color','r');
end
eval(['export_fig -r 2000 DG_V210X_dendogram_finalClsuter_withnumbers_',date,'.pdf'])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     

data_tsne = log2_cent(datamarkers);
% data_tsne = cent_norm(log2(datamarkers+1));

no_dims = 2;
initial_dims = 50;%length(data(:,1));%20;
perplexity = 100;
epsilon = 100;
dist_flag = 2;
max_iter = 1000;
theta = 0.2;

% tic
% mappedX_cell = fast_tsne((data_tsne)', no_dims, initial_dims, perplexity,theta);
% toc
% 
% mappedX_cell(:,2) = -mappedX_cell(:,2);
% saveCellFile([cellid,m2c(mappedX_cell)]...
%     ,['tsne_coordinates_V2kit_finalClsuter_perplexity_',num2str(perplexity),'_Ngenes=',...
%     num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.txt']);
% a = loadCellFile('mnt/sanger-data2/C1_stuff/DentGyr/tsne_coordinates_finalClsuter_perplexity_100_Ngenes=802_NPCA=50_24-Jun-2017.txt');
a = loadCellFile('mnt/sanger-data2/C1_stuff/DentGyr/tsne_coordinates_V2kit_finalClsuter_perplexity_100_Ngenes=756_NPCA=50_04-Aug-2017.txt');
[~,loc] = ismember(cellid,a(:,1));
a = a(loc,:);
mappedX_cell = cell2mat(a(:,2:3));                           



figure;
set(gcf,'color','w','position',[20,20,900,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(T_cells_tmp_uni));
cluster_bound = cell(length(T_cells_tmp_uni),1);
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
    if mod(i,8)+1==1
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    else
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    end
    x = mappedX_cell(in,1); y=mappedX_cell(in,2);
    %     d = sqrt( (x-mean(x)).^2 + (y-mean(y)).^2 );
    %     rmv = d>sqrt(std(x).^2 + std(y)^2);
    rmv = abs(x-median(x))>2*std(x) | abs(y-median(y))>2*std(y);
    x(rmv) = [];
    y(rmv) = [];
%     k = boundary(x,y,0.8);
    
    for h=1:3;
        k = boundary(x,y,0.8);
        rmv = k;
        x(rmv) = [];
        y(rmv) = [];
    end
    k = boundary(x,y,0.8);
    cluster_bound{i} = [smoothn(x(k),5),smoothn(y(k),5)];
    cluster_bound{i} = [cluster_bound{i};cluster_bound{i}(1,:)];
end
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
%     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
    ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(i),'_','-'));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
    plot(cluster_bound{i}(:,1),cluster_bound{i}(:,2));
end

axis tight
axis off
title(['tsne perplexity = ',num2str(perplexity),';Ngenes=',...
    num2str(length(data_tsne(:,1))),';NPCA=',num2str(initial_dims),';eps=',num2str(epsilon),';Niter=',num2str(max_iter)]);
eval(['export_fig tsne_DG_V210X_by_finalClsuter_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


figure;
set(gcf,'color','w','position',[20,20,900,900])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(age_uni));
for i=1:length(age_uni)
    in = age==age_uni(i);
    subplot(3,3,i);
    plot(mappedX_cell(:,1),mappedX_cell(:,2),'.','color',0.6*[1,1,1],'markersize',3);hold on;
    plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',5);hold on;
    axis tight
    title(['age=p',num2str(age_uni(i))]);
    axis off
end

eval(['export_fig tsne_DG_V210X_finalClsuter_byAge_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

figure;
set(gcf,'color','w','position',[20,20,900,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(age_uni));
for i=1:length(age_uni)
    in = age==age_uni(i);
%     if i==3 | i==5
%         plot(mappedX_cell(in,1),mappedX_cell(in,2),'.','color',0.6*[1,1,1],'markersize',3);hold on;
%     else
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',8);hold on;
%     end
    axis tight
%     title(['age=p',num2str(age_uni(i))]);
    axis off
end
legend(cellfun(@num2str,m2c(age_uni),'uniformoutput',false))

for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
%     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
    ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(i),'_','-'));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
    plot(cluster_bound{i}(:,1),cluster_bound{i}(:,2));
end
eval(['export_fig tsne_DG_V210X_finalClsuter_byAge4Colors_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


figure;
set(gcf,'color','w','position',[20,20,900,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
%     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
    ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(i),'_','-')); hold on;
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
    plot(cluster_bound{i}(:,1),cluster_bound{i}(:,2),'color',colors(i,:));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% reelin cells
list = {'Lhx1','Lhx5','Reln','Ndnf','Nhlh2','Lhx1os','Trp73','Diablo','Cacna2d2'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Reln_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% neurogenesis cells
list = {'Ascl1','Cdk1','Neurog2','Eomes','Calb2','Igfbpl1','Sox11','Fxyd7','Dcx'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Neurogenesis_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% astrocytes cells
list = {'Aldoc','Gfap','Rhcg','Fabp7','Ednrb','Itih3','Htra1','Frzb','Tfap2c'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Astrocytes_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% granule cells
list = {'Snap25','Prox1','C1ql2','Ntng1','Neurod6','Stmn2','Emx1','Fxyd7','Icam5'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Granule_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% mossy cells
list = {'Snap25','Slc17a7','Slc32a1','Gad1','Neurod6','Stmn2','Emx1','Sst','Sv2b'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Mossy_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% interneurons cells
list = {'Gad1','Slc32a1','Cnr1','Lhx6','Sst','Arx','Igf1','Tox2','Nos1'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Interneurons_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% oligo cells
list = {'Pdgfra','Cspg4','Tmem100','Neu4','Bmp4','Gpr17','Mog','Plp1','Sox10'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Oligo_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% endothelial microglia cells
list = {'Cldn5','Vtn','Pdgfrb','Igf2','Dcn','Aif1','C1qc','Pf4','Mrc1'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_Endo_micro_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% endothelial microglia cells
list = {'Cldn5','Snap25','Aldoc','Sox10','Aif1','Slc32a1','Slc17a7','Cdk1','Sox11'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),50,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),50,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_MainClass_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);



% additional neuronal
list = {'Tfap2c','Rhcg','Wnt8b','Jag1','Dner','Notch1','Notch2','Lpar1','Dlk2'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_additionGenes1_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% % % % % additional genes
list = {'Prom1', 'Mki67','Top2a','Aurkb','Neurod4','Tac2','Calb2','Gal'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_additionGenes2_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% % % % % additional genes
list = {'Aqp4','Gfap','Rhcg','Vim','Hes1','Hes5','Wnt8b','Tfap2c','Sox2'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
%     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
       interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,3,i);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    for ii=1:length(T_cells_tmp_uni)
        in = T_cells_tmp==T_cells_tmp_uni(ii);
        %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
        plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
    end
    axis tight
    axis off
end
eval(['export_fig tsne_DG_V210X_additionGenes3_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% calc reln cluster stats
ind_reln = strcmpi(grlev2,'cr');
lhx1 = moldata(strcmpi(geneid,'Lhx1'),ind_reln);
lhx1os = moldata(strcmpi(geneid,'Lhx1os'),ind_reln);
lhx5 = moldata(strcmpi(geneid,'Lhx5'),ind_reln);
ndnf = moldata(strcmpi(geneid,'Ndnf'),ind_reln);
reln = moldata(strcmpi(geneid,'Reln'),ind_reln);
trp73 = moldata(strcmpi(geneid,'Trp73'),ind_reln);
for i=1:length(age_uni)
    age_reln(i) = sum(age(ind_reln)==age_uni(i));
    lhx1_pos(i) = sum(lhx1(age(ind_reln)==age_uni(i))>0);
    lhx1os_pos(i) = sum(lhx1os(age(ind_reln)==age_uni(i))>0);
    lhx5_pos(i) = sum(lhx5(age(ind_reln)==age_uni(i))>0);
    ndnf_pos(i) = sum(ndnf(age(ind_reln)==age_uni(i))>0);
    reln_pos(i) = sum(reln(age(ind_reln)==age_uni(i))>0);
    trp73_pos(i) = sum(trp73(age(ind_reln)==age_uni(i))>0);
end
table_reln_stats = [age_uni';age_reln;lhx1_pos;lhx1os_pos;lhx5_pos;ndnf_pos;reln_pos;trp73_pos];
table_reln_stats = [{'Age';'CR_cells';'Lhx1_pos';'Lhx1os_pos';'Lhx5_pos';'Ndnf_pos';'Reln_pos';'Trp73_pos'},m2c(table_reln_stats)];
saveCellFile(table_reln_stats,['Reln_cluster_stats_age_',date,'.txt'])





toc














