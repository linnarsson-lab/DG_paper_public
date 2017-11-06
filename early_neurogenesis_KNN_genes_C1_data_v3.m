
tic
clear all
close all

load /data2/C1_stuff/DentGyr/afterLoading_analysis_DG_fromdatabase_23-Dec-2015.mat

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
green_facs = green_1(validcells);
% source = source(validcells);

tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
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
age = age(validcells);
green_facs = green_facs(validcells);
% source = source(validcells);

% % % % % % % % % % % % % % % % % % % % % % % % % % 

% in = find(sum(data>0,2)>20 & sum(data>0,2)<length(data(1,:))*0.6);

data_full = data;
geneid_full = geneid;
% data = data(in,:);
% geneid = geneid(in);
moldata = data;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

cellid_cluster = loadCellFile_turbo('/data2/C1_stuff/DentGyr/DG_C1_revised_markertable_knnMutal_MCL_1p28_04-Feb-2017_manule_inspection.txt',1);
cellid_cluster = cellid_cluster(2:end,:);
% cellid_cluster(strcmpi(cellid_cluster(:,2),'28'),3) = repmat({'rmv'},sum(strcmpi(cellid_cluster(:,2),'28')),1);
[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid = cellid(loc);
green_facs = green_facs(loc);
% source = source(loc);
age = age(loc);
moldata = moldata(:,loc);
data_full = data_full(:,loc);


grlev2 = cellid_cluster(:,3);
% lev2uni = unique(grlev2);
grlev2_uni = {'astro','rgl','cycling','granul_sox11'};


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
cellid_cluster_1 = cellid_cluster(xi,:);
% source = source(xi);
age = age(xi);
green_facs = green_facs(xi);
moldata = moldata(:,xi);
grlev2 = grlev2(xi);
data_full = data_full(:,xi);

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
corr_filt = xi1(1:1000);
figure; plot(log2_m,log2_cv,'.','markersize',3); hold on;
plot(log2_m(corr_filt),log2_cv(corr_filt),'.r','markersize',3); hold on;
moldata = moldata(in(corr_filt),:);
geneid = geneid(in(corr_filt));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

moldata = moldata./repmat(sum(moldata),length(moldata(:,1)),1)*1e4;
n_pca = 100;
datalog_tmp = cent_norm([(log2(moldata+1))]);
[U]=pcasecon(datalog_tmp,min([length(moldata(1,:)),n_pca]));
prj = [U'*datalog_tmp]';


Rcorr = corr_mat(prj');
% D = squareform(pdist(prj));
k_knn = 40;
[~,XI] = sort(Rcorr,'descend');
% [~,XI] = sort(D);
knnmat = zeros(size(Rcorr));
for i=1:length(knnmat)
    knnmat(i,XI(1:k_knn,i)) = 1;
%     knnmat(XI(1:k,i),i) = 1;
end
knnmat = knnmat.*knnmat';

jacknn = squareform(pdist(knnmat,'jaccard'));
jacknn = jacknn.*knnmat;

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

G = graph(jacknn(in,in),'OmitSelfLoops');
figure('position',[100,100,1000,1000],'color','w');
P = plot(G,'layout','force');

xdata = get(P,'xdata');
ydata = get(P,'ydata');
figure('position',[100,100,1000,1000],'color','w');
Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
color_vec = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    highlight(Ptmp,find(T_cells_tmp(in)==T_cells_tmp_uni(i)),'NodeColor',color_vec(i,:),'MarkerSize',3) ;
    ht = text(mean(xdata(T_cells_tmp(in)==T_cells_tmp_uni(i))), ...
        mean(ydata(T_cells_tmp(in)==T_cells_tmp_uni(i))), dend_order_names{i});
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
end
axis tight
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',4)
title(['early neurogenesis mutal KNN, k=',num2str(k_knn)])
eval(['export_fig Early_C1_Ngenesis_MKNN40_',date,'.pdf']);

% % % % % % % % % % % % % % % % % % % % % % % 
xdata = get(P,'xdata');
ydata = get(P,'ydata');
figure('position',[100,100,1000,1000],'color','w');
Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
color_vec = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    highlight(Ptmp,find(T_cells_tmp(in)==T_cells_tmp_uni(i)),'NodeColor',color_vec(i,:),'MarkerSize',3) ;
    ht = text(mean(xdata(T_cells_tmp(in)==T_cells_tmp_uni(i))), ...
        mean(ydata(T_cells_tmp(in)==T_cells_tmp_uni(i))), dend_order_names{i});
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
end
hold on;
plot(xdata(age(in)>=50), ydata(age(in)>=50), 'sb');
plot(xdata(age(in)<=20), ydata(age(in)<=20), 'dy');
plot(xdata(green_facs(in)>0), ydata(green_facs(in)>0), 'og');


axis tight
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',4)
title(['early neurogenesis mutal KNN, k=',num2str(k_knn)])
eval(['export_fig Early_C1_Age_Ngenesis_MKNN40_',date,'.pdf']);
% % % % % % % % % % % % % % % % % % % % 
clust_edge_con = zeros(length(T_cells_tmp_uni));
clust_edge_con_sum = zeros(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    intmp1 = find(T_cells_tmp(in)==T_cells_tmp_uni(i));
    for j=i:length(T_cells_tmp_uni)
        intmp2 = find(T_cells_tmp(in)==T_cells_tmp_uni(j));
        
        
        clust_edge_con_sum(i,j) = sum(sum(knnmat(in(intmp1),in(intmp2))));
        clust_edge_con_sum(j,i) = clust_edge_con_sum(i,j);   
        
        clust_edge_con(i,j) = clust_edge_con_sum(i,j) / (k_knn*max([length(intmp2),length(intmp1)]));
        clust_edge_con(j,i) = clust_edge_con(i,j);
    end
end

figure('position',[100,100,500,500],'color','w');
imagesc(clust_edge_con*100,[0,10]);
set(gca,'xtick',[1:length(grlev2_uni)],'xticklabel',dend_order_names,'ytick'...
    ,[1:length(grlev2_uni)],'yticklabel',dend_order_names,'XTickLabelRotation',45)
colormap('summer');
freezeColors(gca);
title(['transition frequency mutal KNN, k=',num2str(k_knn)]);
% eval(['export_fig Early_Ngenesis_transition_frequency_MKNN40_',date,'.pdf']);
save2pdf(['Early_C1_Ngenesis_transition_frequency_MKNN40_',date,'.pdf'],gcf,300);

saveCellFile(m2c(clust_edge_con_sum),['Early_C1_Ngenesis_transition_counts_MKNN40_',date,'.txt'])


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% xdata = get(Ptmp,'xdata');
% ydata = get(Ptmp,'ydata');
% figure('position',[100,100,1000,1000],'color','w');
% genePlot = 'cdkn1a';
% markergene = (moldata(strcmpi(geneid,genePlot),in));
% tmpthlow = prctile(markergene(markergene>0),1);
% tmpthhigh = prctile(markergene(markergene>0),95);
% markergene(markergene>tmpthhigh) = tmpthhigh;
% markergene(markergene<tmpthlow) = tmpthlow;
% c_rgb = [1,0,0];rand([1,3]);
% 
% markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
%     interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
%     ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
% h_graph = plot(G,'xdata',xdata,'ydata',ydata);
% set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',8)
% set(h_graph,'NodeColor',markergene_color);
% title(genePlot);
% axis tight
figure('position',[100,100,1000,1000],'color','w');
genePlot = 'c1ql3';
markergene = (data_full(strcmpi(geneid_full,genePlot),in));
tmpthlow = 0;%prctile(markergene(markergene>0),1);
tmpthhigh = prctile(markergene(markergene>0),90);
markergene(markergene>tmpthhigh) = tmpthhigh;
markergene(markergene<tmpthlow) = tmpthlow;
c_rgb = [1,0,0];rand([1,3]);

markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
h_graph = plot(G,'xdata',xdata,'ydata',ydata);
set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',8)
set(h_graph,'NodeColor',markergene_color);
title(genePlot);
axis tight


list_to_plot = {'Cdk1','Cdk2','Cdk4','Cdk6','Mcm2','Mcm3','Mcm4','Mcm5','Mcm6','Mcm7','Mcm8','Mki67','Aurkb','cdkn2c','Cenpa','Cenpe','Cenpf','Cenph','Cenpm',...
    'Gmnn','Top2a','Chek1','Chek2','Ccnb1','Ccna2','Ccnd1','Ccne1',...
    'Neurog2','Neurod2','Neurod1','Neurod4','Eomes','Tbr1','Igfbpl1','Ascl1','Gfap','Sox2','Sox9','Sox11','Aqp4','Aldoc',...
    'Hes1','Hes5','Nes','Lpar1','Dcx','Tfap2c'};
    

xdata = get(Ptmp,'xdata');
ydata = get(Ptmp,'ydata');

for kk=1:length(list_to_plot)
    kk
    hf = figure('position',[100,100,1000,1000],'color','w','visible','off');
    % genePlot = 'dcx';
    genePlot = list_to_plot{kk};
    markergene = (data_full(strcmpi(geneid_full,genePlot),in));
    tmpthlow = 0;%prctile(markergene(markergene>0),1);
    tmpthhigh = 1;%prctile(markergene(markergene>0),90);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    h_graph = plot(G,'xdata',xdata,'ydata',ydata);
    set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',8)
    set(h_graph,'NodeColor',markergene_color);
    title(genePlot);
    axis tight
    
    eval(['export_fig ',genePlot,'_C1_Early_Ngenesis_MKNN40__',date,'.pdf']);
    close(hf)
end


toc