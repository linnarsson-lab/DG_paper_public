
tic
clear all
close all

load /data2/C1_stuff/DentGyr/afterLoading_analysis_DG_fromdatabase_23-Feb-2017.mat

% % % % % % % % % % % % % 
geneuni = genes_DG;

total_mol = sum(moldata_DG);

uni_platevec = unique(chip_DG);
samples_all = cell(length(chip_DG),1);
for i=1:length(chip_DG)
    samples_all{i} = [chip_DG{i},'_',well_DG{i}];
end

% % % % % exclude low wells
goodwells = total_mol>1.5e3 & total_mol<70e3;
moldata_1 = moldata_DG(:,goodwells);
chip_1 = chip_DG(:,goodwells);
well_1 = well_DG(:,goodwells);
age_1 = age_DG(:,goodwells);
sex_1 = sex_DG(:,goodwells);%% female=1, male =-1, mix=0;
diam_1 = diam_DG(:,goodwells);
tissue_all_1 = tissue_DG(goodwells);
total_mol_1 = total_mol(:,goodwells);
image_all_1 = image_all_DG(goodwells);
moldata_spikes = moldata_DG_ercc(:,goodwells);
green_1 = green_DG(goodwells);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
moldata = moldata_1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
tissue_uni = unique(tissue_all_1);

k = strfind(geneuni(:,1),'mt-');
ind_mito = ~cellfun(@isempty,k);
indgood = sum(moldata_1,2)>10 & ~ind_mito;
moldata_1 = moldata_1(indgood,:);
geneuni_1 = geneuni(indgood,1);
cellid_1 = cell(length(chip_1),1);
for i=1:length(chip_1)
    cellid_1{i} = [chip_1{i},'_',well_1{i}];
end

data = moldata_1;
cellid = cellid_1;
geneid = geneuni_1;

tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
tot_genes = sum(data>0);

validcells = (tot_mol>600 & (tot_mol./tot_genes)>1.2 & tot_mol<2e4 & tot_genes>500);
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
age = age_1(validcells);
% source = source(validcells);

tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
tot_genes = sum(data>0);

marker_genes = {'C1ql3','Mog','Aldoc','C1qc','Cldn5'};
% figure('color','w','position',[100,100,900,900]); 
inrmv = false(length(cellid),1);
for i=1:length(marker_genes)
    for j=i+1:length(marker_genes)
%         subplot(4,4,(i-1)*4+j-1);
        tmp1 = data(strcmpi(geneid,marker_genes{i}),:);
        tmp2 = data(strcmpi(geneid,marker_genes{j}),:);
        thtmp1 = 0;%prctile(tmp1(tmp1>0),20);
        thtmp2 = 0;%prctile(tmp2(tmp2>0),20);
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
marker_genes = {'Stmn2','Mog','Aldoc','C1qc','Cldn5'};
for i=1:length(marker_genes)
    for j=i+1:length(marker_genes)
        tmp1 = data(strcmpi(geneid,marker_genes{i}),:);
        tmp2 = data(strcmpi(geneid,marker_genes{j}),:);
        thtmp1 = 0;%prctile(tmp1(tmp1>0),20);
        thtmp2 = 0;%prctile(tmp2(tmp2>0),20);
        inrmv(tmp2>thtmp2 & tmp1>thtmp1) = true;
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % 
validcells = ~inrmv;
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
age = age(validcells);
% source = source(validcells);

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
corr_filt = xi1(1:5000);
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% permute the cells order
rng(1);
cell_perm = randperm(length(cellid));
gene_perm = randperm(length(geneid));
geneid = geneid(gene_perm);
cellid = cellid(cell_perm);
% source = source(cell_perm);
age = age(cell_perm);
moldata = data(gene_perm,cell_perm);
moldata = moldata./repmat(sum(moldata),length(data(:,1)),1)*1e4;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

n_pca = 100;
datalog_tmp = cent_norm([(log2(moldata+1))]);
% datalog_tmp = log2(moldata+1); datalog_tmp = datalog_tmp - repmat(mean(datalog_tmp,2),1,length(datalog_tmp(1,:)));
[U]=pcasecon(datalog_tmp,min([length(moldata(1,:)),n_pca]));
prj = [U'*datalog_tmp]';
% prj = pca_wis(datalog_tmp',n_pca);
% pks = zeros(length(prj(1,:)),1);
% for jj=2:length(prj(1,:))
%     [~,pks(jj)] = kstest2(prj(:,jj-1),prj(:,jj));
% end
% prj = prj(:,pks<0.1);

Rcorr = corr_mat(prj');
% D = squareform(pdist(prj));
k = 10;
[~,XI] = sort(Rcorr,'descend');
% [~,XI] = sort(D);
knnmat = zeros(size(Rcorr));
for i=1:length(knnmat)
    knnmat(i,XI(1:k,i)) = 1;
%     knnmat(XI(1:k,i),i) = 1;
end
knnmat = knnmat.*knnmat';
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
jacknn = squareform(pdist(knnmat,'jaccard'));
jacknn = jacknn.*knnmat;
% jacknn(jacknn<0.8) = 0;

% in = find(sum(knnmat)>1);
[S, C] = graphconncomp(sparse(knnmat),'Directed',false);
f = hist(C,max(C));
in = true(length(knnmat),1);
for i=1:length(f)
    if f(i)<10
        in(C==i) = false;
    end
end
in = find(in);
included_idx = in;


% in = find(sum(jacknn>0)>1);
G = graph(jacknn(in,in),'OmitSelfLoops');
figure('position',[100,100,1000,1000],'color','w');
P = plot(G,'layout','force');


tic
[Z] = cluster_mat_with_mcl(jacknn(in,in),0,[5,1000],1.25);
toc

clust_edge_con = eye(max(Z));
for i=1:max(Z)
    intmp1 = find(Z==i);
    for j=i+1:max(Z)
        intmp2 = find(Z==j);
        clust_edge_con(i,j) = mean(mean(jacknn(in(intmp1),in(intmp2))));
        clust_edge_con(j,i) = clust_edge_con(i,j);
    end
end
% ind_sort = sort_mat_by_neighborhood(1-clust_edge_con,1,5,2);

z = linkage(1-clust_edge_con,'single');
lefaorder = optimalleaforder(z, 1-clust_edge_con);
% lefaorder = ind_sort;
newZ = zeros(size(Z));
for i=1:length(lefaorder);
    newZ(Z==lefaorder(i)) = i;
end
   

xdata = get(P,'xdata');
ydata = get(P,'ydata');
figure('position',[100,100,1000,1000],'color','w');
Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
color_vec = distinguishable_colors(min([max(Z),100]));
for i=1:min([max(newZ),100])
    highlight(Ptmp,find(newZ==i),'NodeColor',color_vec(i,:),'MarkerSize',3) ;
    ht = text(mean(xdata(newZ==i)), mean(ydata(newZ==i)), num2str(i));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
end
axis tight
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',4)
eval(['export_fig DG_C1_knnMutal_MCL_1p28_',date,'.pdf'])
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
[T_cells_tmp,xi] = sort(newZ);
xi(T_cells_tmp==0) = [];
T_cells_tmp(T_cells_tmp==0) = [];
T_cells_tmp = T_cells_tmp;
T_cells_tmp_uni = unique(T_cells_tmp);
% T_cells_tmp_uni = [37,38,8:31,3:4,39:40,36,32,33,34,1,2,5:7,35]
gr_center = zeros(length(T_cells_tmp_uni),1);
data_sorted_all = moldata(:,included_idx(xi));
cells_bor = find(diff(T_cells_tmp)>0);
cellid_sorted = cellid(included_idx(xi));
age_sorted = age(included_idx(xi));


age_uni = unique(age_sorted);
age_bar = zeros(length(age_uni),length(age));
for i=1:length(age_uni)
    age_bar(i,age_sorted==age_uni(i)) = 1;
end



sumpergene = sum(data_sorted_all,2);
meanpergene = mean(data_sorted_all,2);%
molenrich_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
meangr_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
[~,grord] = sort(gr_center);
meangrpos_mat(meangrpos_mat<0.2) = 0;
molenrich_mat(meanpergene==0 | isnan(meanpergene),:) = 0;
% meangrpos_mat = meangrpos_mat(:,grord);
% molenrich_mat = molenrich_mat(:,grord);
[~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0.001,'descend');
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
    molenrich_mat_mark(:,jjj) = mean(data_sorted_all(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj)),2)./meanpergene(ind_gr_tmp_mark);
    meangrpos_mat_mark(:,jjj) = mean(data_sorted_all(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj))>0,2);
end
[~,imax] = max(molenrich_mat_mark.*meangrpos_mat_mark.^0.5,[],2);
[~,xi] = sort(imax);
ind_gr_tmp_mark = ind_gr_tmp_mark(xi);

datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
gr_tmp_mark = gr_tmp_mark(xi);

figure;
set(gcf,'position',[100,100,450,750],'color','w')
ax1 = axes('position',[0.1,0.02,0.88,0.86]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor)
    plot(cells_bor(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',gr_center,'xticklabel',[1:length(gr_center)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 8)
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

eval(['export_fig DG_C1_markertable_knnMutal_MCL_1p28_',date,'.pdf'])

table1 = [[{'cellid','cluster'},gr_tmp_mark'] ;[cellid_sorted,m2c([T_cells_tmp,datamarkers'])] ];
saveCellFile(table1,['DG_C1_markertable_knnMutal_MCL_1p28_',date,'.txt']);


% % % % % % % % % % % % % % % % % % % % % % % % % 
log_m_gr_lev2 = zeros(length(geneid),length(T_cells_tmp_uni));
se_log_m_gr_lev2 = zeros(length(geneid),length(T_cells_tmp_uni));
log_m_gr_lev2_norm = zeros(length(geneid),length(T_cells_tmp_uni));
se_log_m_gr_lev2_norm = zeros(length(geneid),length(T_cells_tmp_uni));
present_mean = zeros(length(geneid),length(T_cells_tmp_uni));
moldata_1_norm = data_sorted_all./repmat(sum(data_sorted_all),length(geneid),1)*1e4;
for i=1:length(T_cells_tmp_uni)
    i
    in = T_cells_tmp==T_cells_tmp_uni(i);
    log_m_gr_lev2(:,i) = mean(log2(data_sorted_all(:,in)+1),2);
    se_log_m_gr_lev2(:,i) = std(log2(data_sorted_all(:,in)+1),[],2)/sqrt(sum(in));
    log_m_gr_lev2_norm(:,i) = mean(log2(moldata_1_norm(:,in)+1),2);
    se_log_m_gr_lev2_norm(:,i) = std(log2(moldata_1_norm(:,in)+1),[],2)/sqrt(sum(in));
    present_mean(:,i) = mean(data_sorted_all(:,in)>0,2);
end


log_m_mat_link = log_m_gr_lev2_norm(ind_gr_tmp_mark,:);%log2(mean_exp_gr_lev2_norm(ind_gr_tmp_mark,:)+1);%
gene_link = gr_tmp_mark;
Z = linkage(log_m_mat_link','ward','correlation');
Rcelltypes = corr_mat(log_m_mat_link);
D = pdist(log_m_mat_link','correlation');
leaforder = optimalleaforder(Z,D);%[1:length(grlev2_uni)]; %

Z = linkage(log_m_mat_link(:,leaforder)','ward','correlation');


figure('color','w','position',[200,200,1450,800],'visible','on');
axes('position',[0.03,0.05,0.4,0.9])
[H,T,outperm] = dendrogram(Z,length(T_cells_tmp_uni),'Orientation','left','reorder',...
    [1:length(T_cells_tmp_uni)],'labels',cellfun(@num2str,m2c(leaforder),'UniformOutput', false),'colorthreshold',0.8);%
set(H,'linewidth',2)
xl = get(gca,'xlim');
set(gca,'ylim',[0.5,length(T_cells_tmp_uni)+0.5],'xtick',[],'ydir','normal','fontsize',8)%,'xlim',[0,1.1]
% set(gca,'xlim',[0.5,length(T_cells_tmp_uni)+1],'ytick',[],'fontsize',8)
% rotateticklabel(gca,45);

axes('position',[0.58,0.05,0.35,0.9])
imagesc(Rcelltypes(leaforder,leaforder),[0.3,1]) 
set(gca,'xtick',[],'xticklabel',leaforder,'xdir','normal','ydir','normal','ytick',[])
colormap('summer')
hcm = colorbar;
set(hcm,'position',[0.96,0.7,0.02,0.25],'ytick',[-0.4:0.1:1],'ydir','normal')
% eval(['export_fig oligos_dendogram_by_lev2_',date,'.pdf'])
save2pdf(['DG_C1_dendogram_by_NotFinalClsuter_',date,'.pdf'],gcf,300)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
ind_gr_tmp_mark = [xi0(1:2,1:end)];
% ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
% rmv = false(size(ind_gr_tmp_mark));
% for i=2:length(ind_gr_tmp_mark)
%     if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
%         rmv(i) = true;
%     end
% end
% ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = geneid(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));
figure('position',[100,100,1000,1000],'color','w');
Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
color_vec = distinguishable_colors(min([max(newZ),100]));
for i=1:min([max(newZ),100])
    highlight(Ptmp,find(newZ==i),'NodeColor',color_vec(i,:),'MarkerSize',3) ;
    ht = text(mean(xdata(newZ==i)), mean(ydata(newZ==i)), gr_tmp_mark((i-1)*2+1:i*2));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
end
axis tight
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',4)
eval(['export_fig DG_C1_bestMarker_knnMutal_MCL_1p28_',date,'.pdf'])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

figure('position',[100,100,1000,1000],'color','w');
age_uni = unique(age);
Ptmp = plot(G,'xdata',xdata,'ydata',ydata); hold on;
color_vec = distinguishable_colors(length(age_uni));
for i=1:length(age_uni)
    plot(xdata(find(age(included_idx)==age_uni(i))),ydata(find(age(included_idx)==age_uni(i))),'.','color',color_vec(i,:),'MarkerSize',6); hold on
%     highlight(Ptmp,find(age(in)==age_uni(i)),'NodeColor',color_vec(i,:),'MarkerSize',2) 
end
axis tight

set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','none','markersize',0.1)
legend( [{'graph'},cellfun(@num2str,m2c(age_uni),'UniformOutput', false)] )
eval(['export_fig DG_C1_byAge_knnMutal_MCL_1p28_',date,'.pdf'])









