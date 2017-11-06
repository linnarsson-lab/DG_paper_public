
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
green_facs = green_facs(cell_perm);
% source = source(cell_perm);
age = age(cell_perm);
moldata = data(gene_perm,cell_perm);
moldata = moldata./repmat(sum(moldata),length(data(:,1)),1)*1e4;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


cellid_cluster = loadCellFile_turbo('/data2/C1_stuff/DentGyr/DG_C1_markertable_knnMutal_MCL_1p28_23-Feb-2017_manual_inspection.txt',1);
cellid_cluster = cellid_cluster(2:end,:);
cellid_cluster(strcmpi(cellid_cluster(:,3),'rmv'),:) = [];
cellid_cluster(:,3) = regexprep(cellid_cluster(:,3),'astro_young','astro');
cellid_cluster(:,3) = regexprep(cellid_cluster(:,3),'astro_old','astro');
[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid = cellid(loc);
green_facs = green_facs(loc);
% source = source(loc);
age = age(loc);
moldata = moldata(:,loc);


grlev2 = cellid_cluster(:,3);
grlev2_uni = unique(grlev2);
grlev2_uni = {'mgl','pvm','endo','peric','ol','astro','rgl','opc','nipc','nb','granul_fxyd7','mossy_calb2',...
    'reln','gaba_cnr1','gaba_lhx6','granul_mature','mossy_sv2b'};


% grlev2_uni = {'endo','peric','vlmc','mgl','pvm','OL','nfol','opc','nsc2','astro','nsc1','gran_cycling','gran_eomes',...
%     'gran_sox11','gran_fxyd7','gran_mature','cck_tox','mossy_adcyap1','mossy_cy26b1','mossy_klk8','gaba_cnr1','gaba_lhx6','reln' };
% grlev2_uni = {'mgl','pvm','endo','vlmc','ol','opc','astro','nsc','cycling','gran_early_cd24a','gran_early_sepw1'...
%     ,'gran_intermi','gran_mature','mossy_calb2','mosssy_sv2b','gaba','reln'};


    

    
    

T_cells_tmp = zeros(length(grlev2),1);
for i=1:length(grlev2_uni)
    T_cells_tmp( strcmpi(grlev2, grlev2_uni{i}) ) = i;
end
[T_cells_tmp,xi] = sort(T_cells_tmp);
xi(T_cells_tmp==0) = [];
T_cells_tmp(T_cells_tmp==0) = [];
T_cells_tmp_uni = unique(T_cells_tmp);
cellid = cellid(xi);
green_facs = green_facs(xi);
% source = source(xi);
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
ax1 = axes('position',[0.1,0.02,0.88,0.86]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
cells_bor = find(diff(T_cells_tmp)>0);
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

eval(['export_fig DG_C1_markertable_finalClsuter_',date,'.pdf'])

table1 = [[{'cellid','cluster'},gr_tmp_mark'] ;[cellid,m2c([T_cells_tmp,datamarkers'])] ];
saveCellFile(table1,['DG_C1_markertable_finalClsuter_',date,'.txt']);




% % % % % % % % % % % % % % % % % % % % % % % % % 

log_m_mat_link = log_m_gr_lev2_norm(ind_gr_tmp_mark,:);%log2(mean_exp_gr_lev2_norm(ind_gr_tmp_mark,:)+1);%
gene_link = gr_tmp_mark;
Z = linkage(log_m_mat_link','ward','correlation');
Rcelltypes = corr_mat(log_m_mat_link);
D = pdist(log_m_mat_link','correlation');
leaforder = [1:length(grlev2_uni)]; %optimalleaforder(Z,D);%

Z = linkage(log_m_mat_link(:,leaforder)','ward','correlation');


figure('color','w','position',[200,200,1450,800],'visible','on');
axes('position',[0.03,0.05,0.4,0.9])
[H,T,outperm] = dendrogram(Z,length(T_cells_tmp_uni),'Orientation','left','reorder',...
    [1:length(T_cells_tmp_uni)],'labels',regexprep(grlev2_uni(leaforder),'_','-'),'colorthreshold',0.8);
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
save2pdf(['DG_C1_dendogram_by_finalClsuter_',date,'.pdf'],gcf,300)

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
eval(['export_fig -r 2000 DG_C1_dendogram_finalClsuter_withnumbers_',date,'.pdf'])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     

data_tsne = log2_cent(datamarkers);
% data_tsne = cent_norm(log2(datamarkers+1));

no_dims = 2;
initial_dims = 40;%length(data(:,1));%20;
perplexity = 60;
epsilon = 100;
dist_flag = 2;
max_iter = 1000;
theta = 0.1;
% 
tic
mappedX_cell = tsne((data_tsne)', [], no_dims, initial_dims, perplexity,epsilon,dist_flag,max_iter);
toc
% mappedX_cell = fliplr(mappedX_cell);
% 
saveCellFile([cellid,m2c(mappedX_cell)]...
    ,['tsne_coordinates_C1_finalClsuter_perplexity_',num2str(perplexity),'_Ngenes=',...
    num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.txt']);
% a = loadCellFile('/data2/C1_stuff/DentGyr/tsne_coordinates_C1_finalClsuter_perplexity_80_Ngenes=535_NPCA=80_01-Dec-2016.txt');
% mappedX_cell = cell2mat(a(:,2:3));                           
gr_label_xy = [-26 -27
    -34 -30
    -34 -10
    -26 -5
    -32 10
    -16 47
    -28 36
    -23 15
    -15 29
    14 30
    8 5
    20 11
    27 15
    32 0
    32 10
    10 -20
    23 -1];   
% names_as_in_paper = regexprep(grlev2_uni,'_','-');
names_as_in_paper = flipud({'Mossy-Sv2b'
    'Granule-mature'
    'GABA-Lhx6'
    'GABA-Cnr1'
    'Cajal-Retzius'
    'Mossy-Calb2'
    'Granule-immature'
    'Neuroblast'
    'nIPC'        
    'OPC'
    'Radial Glia-like'
    'Astrocytes'
    'OL'    
    'VLMC'
    'Endothelial'
    'PVM'
    'Microglia'});

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
%     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
%     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),names_as_in_paper{i});
    ht = text(gr_label_xy(i,1),gr_label_xy(i,2),names_as_in_paper{i},'fontsize',8)%,'backgroundcolor',0.8*[1,1,1]
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6);
end
set(gca,'ydir','reverse');
axis tight
title(['tsne perplexity = ',num2str(perplexity),';Ngenes=',...
    num2str(length(data_tsne(:,1))),';NPCA=',num2str(initial_dims),';eps=',num2str(epsilon),';Niter=',num2str(max_iter)]);
eval(['export_fig tsne_DG_C1_by_finalClsuter_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


figure;
set(gcf,'color','w','position',[20,20,1300,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(age_uni));
for i=1:length(age_uni)
    in = age==age_uni(i);
    subplot(3,6,i);
    plot(mappedX_cell(:,1),mappedX_cell(:,2),'.','color',0.6*[1,1,1],'markersize',3);hold on;
    plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',5);hold on;
    axis tight
    set(gca,'ydir','reverse');
    title(['age=p',num2str(age_uni(i))]);
    axis off
end
eval(['export_fig tsne_DG_C1_finalClsuter_byAge_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);




figure;
set(gcf,'color','w','position',[20,20,900,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(T_cells_tmp_uni));
plot(mappedX_cell(:,1),mappedX_cell(:,2),'.','color',0.8*[1,1,1],'markersize',6);hold on;
markergene = (moldata(strcmpi(geneid,'Gfap'),:));
inpos = markergene>0;
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'p','color',0.8*[1,0,0],'markersize',10);hold on;
markergene = (moldata(strcmpi(geneid,'Cdk1'),:));
inpos = markergene>0;
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'s','color',0.8*[0,0,0],'markersize',10);hold on;
inpos = green_facs>0;
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'.','color',0.8*[0,1,0],'markersize',10);hold on;
inpos = age<=18;
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'o','color',0.8*[0,0,1],'markersize',6);hold on;
inpos = age>=50;
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'o','color',0.8*[1,0,1],'markersize',6);hold on;
set(gca,'ydir','reverse');
legend({'all','Gfap>0','Cdk1>0','FACS-GFP','age<=p18','age>=p50'})
title(['tsne perplexity = ',num2str(perplexity),';Ngenes=',...
    num2str(length(data_tsne(:,1))),';NPCA=',num2str(initial_dims),';eps=',num2str(epsilon),';Niter=',num2str(max_iter)]);
eval(['export_fig tsne_DG_C1_finalClsuter_FACS_GFAP_CDK1_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


figure;
set(gcf,'color','w','position',[20,20,900,800])
markervec = '...............so*+><p';
colors = distinguishable_colors(length(T_cells_tmp_uni));
plot(mappedX_cell(:,1),mappedX_cell(:,2),'.','color',0.8*[1,1,1],'markersize',6);hold on;

inpos = green_facs>0;
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'.','color',0.8*[0,1,0],'markersize',14);hold on;
inpos = age<=20;
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'.','color',0.8*[0,0,1],'markersize',14);hold on;
inpos = age>=50 & ~(green_facs>0);
plot(mappedX_cell(inpos,1),mappedX_cell(inpos,2),'.','color',0.8*[1,0,1],'markersize',14);hold on;
set(gca,'ydir','reverse');
legend({'all','FACS-GFP','age<=p20','age>=p50'})
title(['tsne perplexity = ',num2str(perplexity),';Ngenes=',...
    num2str(length(data_tsne(:,1))),';NPCA=',num2str(initial_dims),';eps=',num2str(epsilon),';Niter=',num2str(max_iter)]);
for i=1:length(T_cells_tmp_uni)
    text(gr_label_xy(i,1),gr_label_xy(i,2),names_as_in_paper{i},'fontsize',8)%,'backgroundcolor',0.8*[1,1,1]
end
eval(['export_fig tsne_DG_C1_finalClsuter_FACS_age_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

set(gca,'xlim',[-35 23],'ylim',[20,65])
eval(['export_fig tsne_DG_C1_zoominAstro_finalClsuter_FACS_age_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% reelin cells
list = {'Lhx1','Lhx5','Reln','Ndnf','Nhlh2','Trp73','Diablo','Cacna2d2'  };
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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Reln_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Neurogenesis_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% astrocytes cells
list = {'Aldoc','Gfap','Veph1','Fabp7','Ednrb','Itih3','Htra1','Frzb','Tfap2c'  };
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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Astrocytes_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% astrocytes cells 2
list = {'Aqp4','Nes','Hes5','Lpar1','Sox2','Hes1','Top2a','Bub1','Cks2' };
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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Astrocytes2_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% granule cells
list = {'Stmn2','Prox1','C1ql3','Ntng1','Fst','C1ql2','Plk5','Slc17a6','Icam5'  };
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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Granule_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);

% mossy cells
list = {'Neurod6','Sv2b','Klk8','Kcnc2','Slc17a7','Dkk3','Bok','Hs3st4','Cyp26b1'  };
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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Mossy_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% interneurons cells
list = {'Gad1','Slc32a1','Cnr1','Lhx6','Sst','Arx','Igf1','Dlx1','Npy'  };
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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Interneurons_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


% oligo cells
list = {'Pdgfra','Cspg4','Tmem100','Neu4','Bmp4','Gpr17','Mog','Plp1','Ermn'  };
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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Oligo_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);


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
    set(gca,'ydir','reverse');
    title(genePlot);
    axis tight
    axis off
end
eval(['export_fig tsne_DG_C1_Endo_micro_cells_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% extract stats of GFP Gfap Cdk1
markergene = (moldata(strcmpi(geneid,'Gfap'),:));
inpos_gfap = markergene>0;
markergene = (moldata(strcmpi(geneid,'Cdk1'),:));
inpos_cdk1 = markergene>0;
inpos_gfp = green_facs>0;
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
    grsum(i) = sum(in);
    sum_gfp(i) = sum(inpos_gfp(in));
    sum_gfap(i) = sum(inpos_gfap(in));
    sum_cdk1(i) = sum(inpos_cdk1(in));
    sum_gfp_gfap(i) = sum(inpos_gfp(in) & inpos_gfap(in));
end
table1 = [grlev2_uni;m2c([grsum;sum_gfp;sum_gfap;sum_gfp_gfap;sum_cdk1])];
table1 = [ {'cluster_name';'sum_cells';'sum_GFP';'sum_Gfap';'sum_GFP&Gfap';'sum_Cdk1'},table1];
saveCellFile(table1,['table_GFP_Gfap_Cdk1_per_cluster_',date,'.txt']);


toc










