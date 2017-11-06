
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

data_full = data;
geneid_full = geneid;
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
moldata = moldata./repmat(sum(moldata),length(data(:,1)),1)*1e4;
data_full = data_full(:,cell_perm);


cellid_cluster = loadCellFile('/data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
cellid_cluster = cellid_cluster(2:end,:);
[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid = cellid(loc);
source = source(loc);
age = age(loc);
moldata = moldata(:,loc);
data_full = data_full(:,loc);

grlev2 = cellid_cluster(:,3);
% lev2uni = unique(lev2gr);
grlev2_uni = {'astro','nsc1','gran_cycling','gran_eomes','gran_sox11'};
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

% to remove gene specific to cycling nIPC cluster, test with ranksum and
% remove all genes with q<0.05
gr1 = T_cells_tmp==3;
gr2 = T_cells_tmp~=3;
for k=1:length(geneid);
    p_rs1(k) = ranksum(moldata(k,gr1),moldata(k,gr2),'tail','both');
end
q1 = qval_from_pval(p_rs1)';
moldata = moldata(q1>0.05,:);
geneid = geneid(q1>0.05);
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
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',4)
title(['early neurogenesis mutal KNN, k=',num2str(k_knn)])
eval(['export_fig Early_Ngenesis_NO_CCgenes_MKNN40__',date,'.pdf']);


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
% % eval(['export_fig Early_Ngenesis_transition_frequency_MKNN40_',date,'.pdf']);
% save2pdf(['Early_Ngenesis_transition_frequency_MKNN40_',date,'.pdf'],gcf,300);
% 
% saveCellFile(m2c(clust_edge_con_sum),['Early_Ngenesis_transition_counts_MKNN40_',date,'.txt'])


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

xdata = get(Ptmp,'xdata');
ydata = get(Ptmp,'ydata');
figure('position',[100,100,1000,1000],'color','w');
genePlot = 'cdk1';
markergene = (data_full(strcmpi(geneid_full,genePlot),in));
tmpthlow = 0;prctile(markergene(markergene>0),1);
tmpthhigh = prctile(markergene(markergene>0),95);
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


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

list_to_plot = {'Cdk1','Cdk2','Cdk4','Cdk6','Mcm2','Mcm3','Mcm4','Mcm5','Mcm6','Mcm7','Mcm8','Mki67','Aurkb','cdkn2c','Cenpa',...
    'Cenpe','Cenpf','Cenph','Cenpm',...
    'Gmnn','Top2a','Chek1','Chek2','Ccnb1','Ccna2','Ccnd1','Ccne1',...
    'Neurog2','Neurod2','Neurod1','Neurod4','Eomes','Tbr1','Igfbpl1','Ascl1','Gfap','Sox2','Sox9','Sox11','Aqp4','Aldoc',...
    'Hes1','Hes5','Nes','Lpar1','Dcx','Rgs5', 'Rhcg', 'Vnn1','Tfap2c','Etv4','Id3','Olig1','Olig2','Tshz2','Hopx',....
    'Ckap2l','Nhlh2','Tac2','Gal','Bhlhe22','Hs3st1','Egfr','Emx1','Emx2','Calb2','Calb1','Foxg1'};
    

xdata = get(Ptmp,'xdata');
ydata = get(Ptmp,'ydata');
marker_posfrac_nipc1 = zeros(length(list_to_plot),1);
marker_posfrac_nipc2 = zeros(length(list_to_plot),1);
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
    hold on;
    nipc_x = xdata(T_cells_tmp(in)==3);
    nipc_y = ydata(T_cells_tmp(in)==3);    
%     highlight(h_graph,find(T_cells_tmp(in)==3),'NodeColor',color_vec(i,:),'MarkerSize',3,'marker','s') ;
    title(genePlot);
    axis tight
    
    markergene_nipc = markergene(T_cells_tmp==3);
    line_p_x = [-1.5,2.5];%[0.402 0.63];
    line_p_y = [2.3,-0.8];%[0.76 0.53];
    m = diff(line_p_y)/diff(line_p_x);
    y_proj_nipc = m*(nipc_x - line_p_x(1)) + line_p_y(1);
    nipc_g1 = find( (nipc_y-y_proj_nipc)>0 & (nipc_y-y_proj_nipc)<3 );
    nipc_g2 = find( (nipc_y-y_proj_nipc)<0 );
    plot(nipc_x(nipc_g1),nipc_y(nipc_g1),'sb');
    plot(nipc_x(nipc_g2),nipc_y(nipc_g2),'oc');
    plot(line_p_x,line_p_y,'k','linewidth',2)
%     annotation(hf,'line',line_p_x,line_p_y);%annotation(hf,'line',[0.402 0.63],[0.76 0.53]);
    marker_posfrac_nipc1(kk) = sum(markergene_nipc(nipc_g1))/length(nipc_g1);
    marker_posfrac_nipc2(kk) = sum(markergene_nipc(nipc_g2))/length(nipc_g2);


    
    eval(['export_fig ',genePlot,'_nIPC_dissect_noCC_Early_Ngenesis_MKNN40__',date,'.pdf']);
    close(hf)
end

table1 = [{'genesym','fraction_upper_nIPC','fraction_lower_nIPC'};[list_to_plot',m2c([marker_posfrac_nipc1,marker_posfrac_nipc2])]];
saveCellFile(table1,['nIPC_dissect_noCC_Early_Ngenesis_MKNN40_fracpos_',date,'.txt']);




toc