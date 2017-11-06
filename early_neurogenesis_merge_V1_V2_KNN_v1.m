tic
clear all
close all


load 10X43_1_v1.mat
load 10X46_1_v1.mat

cellid_43_1 = cellfun(@(x) ['10X43_1_',x,'-1'], cellid_43_1,'uniformoutput',0);
cellid_46_1 = cellfun(@(x) ['10X46_1_',x,'-1'], cellid_46_1,'uniformoutput',0);




load /mnt/sanger-data2/C1_stuff/DentGyr/DG_V2kit_subsample_samples_merged_18-Jun-2017.mat

[geneid,ia,ib] = intersect(geneid,geneid_46_1);
data = [data(ia,:), data_43_1(ib,:), data_46_1(ib,:)];
data1 = data(ia,:); data2= [ data_43_1(ib,:), data_46_1(ib,:)]; 
dm = mean(log2(data1+1),2) - mean(log2(data2+1),2);
data(abs(dm)>0.8,:) = [];
geneid(abs(dm)>0.8) = [];


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
cellid_cluster = cellid_cluster(2:end,:);
cellid_cluster2 = loadCellFile('mnt/sanger-data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
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
grlev2_uni = {'rgl','astro_old','astro_adol','immature_astro','rgl_young','nipc_young','nipc','nb1',.....
    'astro','nsc1','gran_cycling','gran_eomes','gran_sox11'};

tf = ismember(grlev2,grlev2_uni);
[~,p] = ttest2(log2(moldata(:,tf)+1)', log2(moldata(:,~tf)+1)','tail','right');
in_leave = fdr_proc(p,0.1);


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

moldata = moldata./repmat(sum(moldata),length(moldata(:,1)),1)*5e3;

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
corr_filt = xi1(1:2000);
figure; plot(log2_m,log2_cv,'.','markersize',3); hold on;
plot(log2_m(corr_filt),log2_cv(corr_filt),'.r','markersize',3); hold on;
moldata = moldata(in(corr_filt),:);
geneid = geneid(in(corr_filt));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


n_pca = 30;
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
Ptmp = plot(G,'layout','force');

xdata = get(Ptmp,'xdata');
ydata = get(Ptmp,'ydata');
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
axis off
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',6)
title(['early neurogenesis mutal KNN, k=',num2str(k_knn)])
eval(['export_fig Early_Ngenesis_V1V2_MKNN40_',date,'.pdf']);


% clust_edge_con = zeros(length(T_cells_tmp_uni));
% clust_edge_con_sum = zeros(length(T_cells_tmp_uni));
% for i=1:length(T_cells_tmp_uni)
%     intmp1 = find(T_cells_tmp(in)==T_cells_tmp_uni(i));
%     for j=i+1:length(T_cells_tmp_uni)
%         intmp2 = find(T_cells_tmp(in)==T_cells_tmp_uni(j));
%         clust_edge_con(i,j) = mean(mean(knnmat(in(intmp1),in(intmp2))));
%         clust_edge_con(j,i) = clust_edge_con(i,j);
%         
%         clust_edge_con_sum(i,j) = sum(sum(knnmat(in(intmp1),in(intmp2))));
%         clust_edge_con_sum(j,i) = clust_edge_con_sum(i,j);        
%     end
% end
% 
% figure('position',[100,100,500,500],'color','w');
% imagesc(clust_edge_con*100,[0,1]);
% set(gca,'xtick',[1:length(grlev2_uni)],'xticklabel',dend_order_names,'ytick'...
%     ,[1:length(grlev2_uni)],'yticklabel',dend_order_names,'XTickLabelRotation',45)
% colormap('summer');
% freezeColors(gca);
% title(['transition frequency mutal KNN, k=',num2str(k_knn)]);
% eval(['export_fig Early_Ngenesis_transition_frequency_MKNN40_',date,'.pdf']);
% save2pdf(['Early_Ngenesis_transition_frequency_MKNN40_',date,'.pdf'],gcf,300);



% % % % % % % % % % % % % % % % % % % 
% xdata = get(Ptmp,'xdata');
% ydata = get(Ptmp,'ydata');
color_vec = distinguishable_colors(length(age_uni));
figure('position',[100,100,1000,1000],'color','w');
Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','none'); hold on;

for i=1:length(age_uni)
    plot(xdata(age(in)==age_uni(i)),ydata(age(in)==age_uni(i)),'.','Color',color_vec(i,:),'MarkerSize',6) ; hold on;
end

axis tight
axis off
title(['early neurogenesis mutal KNN, k=',num2str(k_knn)])
legend([{'graph'};cellfun(@num2str,m2c(age_uni),'uniformoutput',false)])
eval(['export_fig Early_Ngenesis_age_V1V2_MKNN40_',date,'.pdf']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
color_vec = distinguishable_colors(length(age_uni));
figure('position',[100,100,1600,600],'color','w');
[ha, pos] = tight_subplot(2, 7, [0.01,0.01], [0.01,0.05], [0.01,0.01]);
for i=1:length(age_uni)
%     subplot(2,4,i);
    axes(ha(i));
    Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
    set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','none'); hold on;
    plot(xdata(age(in)==age_uni(i)),ydata(age(in)==age_uni(i)),'.','Color',color_vec(i,:),'MarkerSize',6) ; hold on;
    axis tight
    axis off
    title(['age = ',num2str(age_uni(i))])
end
eval(['export_fig Early_Ngenesis_age_subplots_V1V2_MKNN40_',date,'.pdf']);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

list_to_plot = {'Cdk1','Top2a','Neurog2','Neurod6','Eomes','Tbr1','Emx1','Igfbpl1','Ascl1'};


xdata = get(Ptmp,'xdata');
ydata = get(Ptmp,'ydata');
figure('position',[100,100,1000,1000],'color','w');
[ha, pos] = tight_subplot(3, 3, [0.01,0.01], [0.01,0.05], [0.01,0.01]);
for kk=1:length(list_to_plot)
%     hf = figure('position',[100,100,1000,1000],'color','w');
    axes(ha(kk));
    % genePlot = 'dcx';
    genePlot = list_to_plot{kk};
    markergene = (data_full(strcmpi(geneid_full,genePlot),in));
    tmpthlow = 0;%prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),90);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    inpos = markergene>0;
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    h_graph = plot(G,'xdata',xdata,'ydata',ydata);
%     set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',8)
    set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','none'); hold on;
    scatter(xdata(~inpos),ydata(~inpos),30,markergene_color(~inpos,:),'.'); hold on;
    scatter(xdata(inpos),ydata(inpos),30,markergene_color(inpos,:),'.'); hold on;
%     set(h_graph,'NodeColor',markergene_color);
    title(genePlot);
    axis tight
    axis off    
%     eval(['export_fig ',genePlot,'_Early_Ngenesis_MKNN40__',date,'.pdf']);
%     close(hf)
end
eval(['export_fig Early_Ngenesis_markers1_V1V2_MKNN40_',date,'.pdf']);

list_to_plot = {'Gfap','Aqp4','Aldoc','Stmn2','Tfap2c','Rhcg','Prox1','Neurod2','Neurod1'};

xdata = get(Ptmp,'xdata');
ydata = get(Ptmp,'ydata');
figure('position',[100,100,1000,1000],'color','w');
[ha, pos] = tight_subplot(3, 3, [0.01,0.01], [0.01,0.05], [0.01,0.01]);
for kk=1:length(list_to_plot)
%     hf = figure('position',[100,100,1000,1000],'color','w');
    axes(ha(kk));
    % genePlot = 'dcx';
    genePlot = list_to_plot{kk};
    markergene = (data_full(strcmpi(geneid_full,genePlot),in));
    tmpthlow = 0;%prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),90);
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    inpos = markergene>0;
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    h_graph = plot(G,'xdata',xdata,'ydata',ydata);
%     set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',8)
    set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','none'); hold on;
    scatter(xdata(~inpos),ydata(~inpos),30,markergene_color(~inpos,:),'.'); hold on;
    scatter(xdata(inpos),ydata(inpos),30,markergene_color(inpos,:),'.'); hold on;
%     set(h_graph,'NodeColor',markergene_color);
    title(genePlot);
    axis tight
    axis off    
%     eval(['export_fig ',genePlot,'_Early_Ngenesis_MKNN40__',date,'.pdf']);
%     close(hf)
end

eval(['export_fig Early_Ngenesis_markers2_V1V2_MKNN40_',date,'.pdf']);




% list_to_plot = {'Cdk1','Cdk2','Cdk4','Cdk6','Mcm2','Mcm3','Mcm4','Mcm5','Mcm6','Mcm7','Mcm8','Mki67','Aurkb','cdkn2c','Cenpa','Cenpe','Cenpf','Cenph','Cenpm',...
%     'Gmnn','Top2a','Chek1','Chek2','Ccnb1','Ccna2','Ccnd1','Ccne1',...
%     'Neurog2','Neurod2','Neurod1','Neurod4','Eomes','Tbr1','Igfbpl1','Ascl1','Gfap','Sox2','Sox9','Sox11','Aqp4','Aldoc',...
%     'Hes1','Hes5','Nes','Lpar1','Dcx'};
%     
% 
% list_to_plot = {'Cdk1','Top2a',...
%     'Neurog2','Neurod6','Eomes','Tbr1','Emx1','Igfbpl1','Ascl1','Gfap','Aqp4','Aldoc',...
%     'Stmn2','Tfap2c','Rhcg'};
% 
% xdata = get(Ptmp,'xdata');
% ydata = get(Ptmp,'ydata');
% 
% for kk=1:length(list_to_plot)
%     hf = figure('position',[100,100,1000,1000],'color','w');
%     % genePlot = 'dcx';
%     genePlot = list_to_plot{kk};
%     markergene = (data_full(strcmpi(geneid_full,genePlot),in));
%     tmpthlow = 0;%prctile(markergene(markergene>0),1);
%     tmpthhigh = 1;%prctile(markergene(markergene>0),90);
%     markergene(markergene>tmpthhigh) = tmpthhigh;
%     markergene(markergene<tmpthlow) = tmpthlow;
%     c_rgb = [1,0,0];rand([1,3]);
%     
%     markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
%         interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
%         ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
%     h_graph = plot(G,'xdata',xdata,'ydata',ydata);
%     set(h_graph,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',8)
%     set(h_graph,'NodeColor',markergene_color);
%     title(genePlot);
%     axis tight
%     
% %     eval(['export_fig ',genePlot,'_Early_Ngenesis_MKNN40__',date,'.pdf']);
% %     close(hf)
% end
% 

