
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
% moldata = moldata./repmat(sum(moldata),length(data(:,1)),1)*1e4;


cellid_cluster = loadCellFile('/data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
cellid_cluster = cellid_cluster(2:end,:);
[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid = cellid(loc);
source = source(loc);
age = age(loc);
moldata = moldata(:,loc);


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

age_uni = unique(age);
age_bar = zeros(length(age_uni),length(age));
for i=1:length(age_uni)
    age_bar(i,age==age_uni(i)) = 1;
end


% in = find(sum(moldata>0,2)>20 & sum(moldata>0,2)<length(moldata(1,:))*0.6);
% m_v = mean(moldata(in,:),2);
% cv_v = std(moldata(in,:),[],2)./m_v;
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
% corr_filt = xi1(1:10000);
% figure; plot(log2_m,log2_cv,'.','markersize',3); hold on;
% plot(log2_m(corr_filt),log2_cv(corr_filt),'.r','markersize',3); hold on;
% moldata = moldata(in(corr_filt),:);
% geneid = geneid(in(corr_filt));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
TF_list = loadCellFile('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/list_of_TF_mm');
tf = ismember(geneid,TF_list);



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
moldata_validcells = moldata;
dend_order_names = grlev2_uni;%regexprep(grlev2_uni,'_','-');
% Tcelluni = unique(lev2_gr_num_treeorder);
space_gr = 20;
colorvec = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    gr_cent(i) = mean(ind + (i-1)*space_gr);
end


listtoplot = loadCellFile('/data2/C1_stuff/DentGyr/Early_Ngenesis_peakTime_TF_20-Jan-2017_shortlist.txt');
list2 = loadCellFile('/data2/C1_stuff/DentGyr/Early_Ngenesis_peakTime_nonTFnonCC_19-Jan-2017.txt');
intotake = [];
for i=1:length(grlev2_uni)
    tmp1 = find(strcmpi(list2(:,2),grlev2_uni{i}));
    tmp2 = find(strcmpi(listtoplot(:,2),grlev2_uni{i}));
    listtoplot = [listtoplot(1:tmp2(1)-1,:);listtoplot(tmp2,:);list2(tmp1(1:3),:);listtoplot(tmp2(end)+1:end,:)];
end


listtoplot = flipud(listtoplot(:,1));

yax_lim = [0,space_gr*(length(T_cells_tmp_uni)-1) + length(moldata_validcells(1,:))];
figure;
set(gcf,'position',[100,100,450,750],'color','w')
barhight = 0.90/length(listtoplot);
for k=1:length(listtoplot)
    k
    geneind = find(strcmpi(geneid,listtoplot{k}));
    %     hf = figure('color','w','position',[200,200,1000,800],'visible','on');
    %     axes('position',[0.395,0.05,0.58,0.92])
    ax_vec(k) = axes('position',[0.1,0.02+barhight*(k-1),0.84,barhight]);
    hold on;
    
    maxval = sort(moldata_validcells(geneind,:),'descend');
    maxval = maxval(5);
    x_vec = [];
    if maxval>0
        for i=1:length(T_cells_tmp_uni)
            ind = find(T_cells_tmp==i);
            tmp = moldata_validcells(geneind,ind);
            inpos = find(tmp>0);
%             inneg = find(tmp==0);
%             if length(inpos)<10 & length(inneg)>9
%                 inneg = inneg(randperm(length(inneg)));
%                 inpos = [inpos,inneg(1:10)];
%             end
            plot((ind(1) + (i-1)*space_gr-1)*[1,1],[0,maxval], '-','color',[0,0,0]);%220/255*[1,1,1]
            plot((ind(end) + (i-1)*space_gr+0.5)*[1,1],[0,maxval], '-','color',[0,0,0]);
            x_vec = [x_vec;ind(inpos) + (i-1)*space_gr];
%             if ~isempty(inpos)
%                 h = bar(ind(inpos) + (i-1)*space_gr,tmp(inpos));
%                 set(h,'facecolor',colorvec(i,:),'edgecolor',colorvec(i,:),'barwidth',1);%colorvec(i,:)     get_RGB('NavyBlue')
%             end
        end
        tmp = moldata_validcells(geneind,:);
        h = bar(x_vec,tmp(tmp>0));
        set(h,'facecolor',get_RGB('NavyBlue'),'edgecolor',get_RGB('NavyBlue'),'barwidth',1);%colorvec(i,:)     get_RGB('NavyBlue')
    else
        maxval = 1;
        for i=1:length(T_cells_tmp_uni)
            ind = find(T_cells_tmp==i);
            plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
            plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
        end
    end
    if mod(k,2)==0
        set(gca,'color',0.8*[1,1,1])
    end
    if k==1
        set(gca,'xtick',gr_cent,'xticklabel',dend_order_names,'fontsize',5,'ylim',[0,maxval],'xlim',yax_lim,'ydir','normal','ytick',[])
        %         ylabel(geneid{geneind},'fontsize',6);
    else
        set(gca,'xtick',[],'fontsize',7,'ylim',[0,maxval],'xlim',yax_lim,'ydir','normal','ytick',[])
        %         ylabel(geneid{geneind},'fontsize',6);
    end
    text(yax_lim(2),maxval/2,[num2str(round(maxval))],'fontsize',5,'horizontalalignment','left')
    text(-10,maxval/2,[geneid{geneind}],'fontsize',6,'horizontalalignment','right')
    
end

eval(['export_fig Early_Ngenesis_SCbarplots_TFs_',date,'.pdf']);











toc

