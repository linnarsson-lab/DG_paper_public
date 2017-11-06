
tic
clear all
close all

load /mnt/sanger-data2/C1_stuff/DentGyr/DG_V2kit_subsample_samples_merged_18-Jun-2017.mat


tot_mol = sum(data);
tot_mol(tot_mol>2e4) = 3e4;
tot_genes = sum(data>0);

validcells = (tot_mol>1000 & (tot_mol./tot_genes)>1.2 & tot_mol<3e4 & tot_genes>800);
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
corr_filt = xi1(1:5000);
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
% moldata = round(moldata./repmat(sum(moldata),length(data(:,1)),1)*1e4);

n_pca = 20;
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
k = 30;
[~,XI] = sort(Rcorr,'descend');
% [~,XI] = sort(D);
knnmat = zeros(size(Rcorr));
for i=1:length(knnmat)
    knnmat(i,XI(1:k,i)) = 1;
    knnmat(XI(1:k,i),i) = 1;
end
knnmat = knnmat.*knnmat';
% % % % % % % % % % % % % % % % % % % 
% knnmat2 = knnmat;
% for i=1:length(knnmat)
%     if sum(knnmat2(i,:))>10
%         tmpind = find(knnmat2(i,:)>0);
%         [~,cord] = sort(Rcorr(i,tmpind),'descend');
%         knnmat2(i,tmpind(cord(11:end))) = 0;
%         knnmat2(tmpind(cord(11:end)),i) = 0;
%     end
% end
% knnmat = knnmat2;
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
% ind1 = find(C(in)==1,1,'first');
% Cuni = unique(C(in));
% for i=2:length(Cuni)
%     indtmp(i) = find(C(in)==Cuni(i),1,'first');
%     jacknn(in(ind1),in(indtmp(i))) = 1;
%     jacknn(in(indtmp(i)),in(ind1)) = 1;
% end
    

% in = find(sum(jacknn>0)>1);
G = graph(jacknn(in,in),'OmitSelfLoops');
figure('position',[100,100,1000,1000],'color','w');
P = plot(G,'layout','force');
% layout(P,'force','Iterations',20);
% color_vec = distinguishable_colors(100);
% [COMTY ending] = cluster_jl(jacknn(in,in),1,0);
% for i=1:100
%     highlight(P,find(COMTY.COM{1}==i),'NodeColor',color_vec(i,:),'MarkerSize',6) 
% end
% axis tight


% G = graph(knnmat(in,in),'OmitSelfLoops');
% figure('position',[100,100,1000,1000],'color','w');
% P = plot(G,'layout','force');
% color_vec = distinguishable_colors(100);
% [COMTY ending] = cluster_jl(knnmat(in,in),0,0);
% for i=1:100
%     highlight(P,find(COMTY.COM{1}==i),'NodeColor',color_vec(i,:),'MarkerSize',6) 
% end

tic
[Z] = cluster_mat_with_mcl(jacknn(included_idx,included_idx),0,[5,3000],1.4);
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
% eval(['export_fig DG_10X_V1V2_knnMutal_MCL_1p25_',date,'.pdf'])
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


ind_gr_tmp_mark = [xi0(1:5,:);xi0p5(1:5,:);xi1(1:5,:)];
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
set(gca,'xtick',gr_center,'xticklabel',[1:length(gr_center)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 6)
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

eval(['export_fig DG_10X_V2_markertable_knnMutal_MCL_1p4_',date,'.pdf'])

% table1 = [[{'cellid','cluster'},gr_tmp_mark'] ;[cellid_sorted,m2c([T_cells_tmp,datamarkers'])] ];
% saveCellFile(table1,['DG_10X_V1V2_markertable_knnMutal_MCL_1p4_',date,'.txt']);

table1 = [[{'cellid','cluster'}] ;[cellid_sorted,m2c([T_cells_tmp])] ];
saveCellFile(table1,['DG_10X_V1V2_markertable_knnMutal_MCL_1p4_',date,'.txt']);
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
color_vec = distinguishable_colors(min([max(Z),100]));
for i=1:min([max(newZ),100])
    highlight(Ptmp,find(newZ==i),'NodeColor',color_vec(i,:),'MarkerSize',3) ;
    ht = text(mean(xdata(newZ==i)), mean(ydata(newZ==i)), gr_tmp_mark((i-1)*2+1:i*2));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
end
axis tight
set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',4)
eval(['export_fig DG_10X_V2_bestMarker_knnMutal_MCL_1p25_',date,'.pdf'])
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
legend( [{'graph'};cellfun(@num2str,m2c(age_uni),'UniformOutput', false)] )
% eval(['export_fig DG_10X_V1V2_byAge_knnMutal_MCL_1p25_',date,'.pdf'])


% 
% xdata = get(P,'xdata');
% ydata = get(P,'ydata');
% figure('position',[100,100,1000,1000],'color','w');
% Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
% color_vec = distinguishable_colors(100);
% for i=1:100
%     highlight(Ptmp,find(COMTY.COM{1}==i),'NodeColor',color_vec(i,:),'MarkerSize',6) 
% end
% axis tight
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% xdata = get(P,'xdata');
% ydata = get(P,'ydata');
% figure('position',[100,100,1000,1000],'color','w');
% genePlot = 'fxyd7';
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
% Ptmp = plot(G,'xdata',xdata,'ydata',ydata);
% set(Ptmp,'linewidth',0.1,'edgecolor',0.8*[1,1,1],'edgealpha',0.8,'marker','.','markersize',3)
% set(Ptmp,'NodeColor',markergene_color);
% title(genePlot);
% axis tight









% % % % % % % % % % % % % % % % % % % % % % % % % % % % 



data_tsne = log2_cent(datamarkers);
% data_tsne = cent_norm(log2(datamarkers+1));

no_dims = 2;
initial_dims = 20;%length(data(:,1));%20;
perplexity = 80;
epsilon = 100;
dist_flag = 2;
max_iter = 1000;
theta = 0.5;
% 
% tic
% mappedX_cell = tsne((data_tsne)', [], no_dims, initial_dims, perplexity,epsilon,dist_flag,max_iter);
% toc
% tic
tic
mappedX_cell = fast_tsne((data_tsne)', no_dims, initial_dims, perplexity,theta);
toc

% saveCellFile([cellid(included_idx),m2c(mappedX_cell)]...
%     ,['tsne_coordinates_MCLcluster_perplexity_',num2str(perplexity),'_Ngenes=',...
%     num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.txt']);


figure;
set(gcf,'color','w','position',[20,20,900,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
    if mod(i,8)+1==1
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    else
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    end
end
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
    ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
end
axis tight
title(['tsne perplexity = ',num2str(perplexity),';Ngenes=',...
    num2str(length(data_tsne(:,1))),';NPCA=',num2str(initial_dims),';eps=',num2str(epsilon),';Niter=',num2str(max_iter)]);
eval(['export_fig tsne_DG_10X_V2_byMCLcluster_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


figure;
set(gcf,'color','w','position',[20,20,900,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(age_uni));
for i=1:length(age_uni)
    in = age_sorted==age_uni(i);
    if mod(i,8)+1==1
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    else
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    end
end
axis tight
legend(cellfun(@num2str,m2c(age_uni),'UniformOutput',false))
eval(['export_fig tsne_DG_10X_V2_byAge_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% % % % % % % % % % % % % % % % % % % % % % /

% reelin cells
list = {'Lhx1','Lhx5','Reln','Ndnf','Nhlh2','Lhx1os','Trp73','Diablo','Cacna2d2'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_CR_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% neurogenesis cells
list = {'Ascl1','Cdk1','Neurog2','Eomes','Calb2','Igfbpl1','Sox11','Fxyd7','Dcx'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_Neurogenesis1_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);
% astrocytes cells
list = {'Aldoc','Gfap','Rhcg','Fabp7','Ednrb','Itih3','Htra1','Frzb','Tfap2c'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end

eval(['export_fig tsne_DG_10X_V2_astro1_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% oligos cells
list = {'Pdgfra','Cspg4','Tmem100','Neu4','Bmp4','Gpr17','Mog','Plp1','Sox10'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_oligos_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);



% granule cells
list = {'Snap25','Prox1','C1ql2','Ntng1','Neurod6','Stmn2','Emx1','Fxyd7','Icam5'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_Neuros1_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);
% astrocytes cells
list = {'Aqp4','Gfap','Rhcg','Vim','Hes1','Hes5','Wnt8b','Tfap2c','Sox2'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_astro2_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% astrocytes cells
list = {'Tnc','Slc1a3','Gm10561','Vim','Ptprz1','Egfr','Padi2','Notch2','Sox1'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_astro3_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);



% granule cells
list = {'Snap25','Slc17a7','Slc32a1','Gad1','Neurod6','Stmn2','Emx1','Sst','Sv2b'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_neurons2_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);



% astrocytes cells
list = {'Tfap2c','Prom1','Lpar1','Top2a','Pbk','Neurod4','Neurog2','Neurod1','Neurod2'  };
figure;
set(gcf,'color','w','position',[20,20,1100,960])
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_sorted_all(strcmpi(geneid,genePlot),:));
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
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10X_V2_astro4_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);



