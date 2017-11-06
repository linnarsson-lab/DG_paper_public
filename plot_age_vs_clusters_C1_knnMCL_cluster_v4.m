
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


data = data(in,:);
geneid = geneid(in);
moldata = data;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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
% lev2uni = unique(grlev2);
% grlev2_uni = {'endo','peric','vlmc','mgl','pvm','OL','nfol','opc','nsc2','astro','nsc1','gran_cycling','gran_eomes',...
%     'gran_sox11','gran_fxyd7','gran_mature','cck_tox','mossy_adcyap1','mossy_cy26b1','mossy_klk8','gaba_cnr1','gaba_lhx6','reln' };
grlev2_uni = {'mgl','pvm','endo','peric','ol','astro','rgl','opc','nipc','nb','granul_fxyd7','mossy_calb2',...
    'reln','gaba_cnr1','gaba_lhx6','granul_mature','mossy_sv2b'};
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


age_cluster_mat = zeros(length(age_uni),length(grlev2_uni));

for i=1:length(T_cells_tmp_uni)
    i
    in = T_cells_tmp==T_cells_tmp_uni(i);
    age_cluster_mat(:,i) = sum(age_bar(:,in),2);
end

tmp = age_cluster_mat./repmat(sum(age_cluster_mat),length(age_uni),1)*100;
tmp(tmp==0) = 0.001;

figure;
set(gcf,'position',[100,100,900,600],'color','w')
axes('position',[0.2,0.13,0.75,0.85]);
x = repmat([1:length(grlev2_uni)],length(age_uni),1);
y = repmat([1:length(age_uni)]',1,length(grlev2_uni));
scatter(x(:), y(:), 3*tmp(:),'facecolor','b');
axis tight
set(gca,'xtick',[1:length(grlev2_uni)],'XTickLabel',names_as_in_paper,'ytick',[1:length(age_uni)],'YTickLabel',age_uni,'XTickLabelRotation',45,'ydir','reverse')
eval(['export_fig blobs_C1_age_vs_clusters_',date,'.pdf']);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% colorvec = distinguishable_colors(length(T_cells_tmp_uni));
% figure;
% set(gcf,'position',[100,100,450,750],'color','w')
% barhight = 0.90/length(T_cells_tmp_uni);
% for k=1:length(T_cells_tmp_uni)
%     k
%     %     hf = figure('color','w','position',[200,200,1000,800],'visible','on');
%     %     axes('position',[0.395,0.05,0.58,0.92])
%     ax_vec(k) = axes('position',[0.4,0.07+barhight*(k-1),0.5,barhight]);
%     hold on;
%     plot(tmp(:,k),'--o','color',colorvec(k,:),'markerfacecolor',colorvec(k,:));
%     maxval = max(tmp(:,k));
%     text(-1,maxval/2,names_as_in_paper{k},'horizontalalignment','right','fontsize',8);
%     if k>1
%         set(gca,'xtick',[],'ylim',[0,maxval*1.2],'fontsize',8)
%     else
%         set(gca,'xtick',[1:length(age_uni)],'xTickLabel',age_uni,'ylim',[0,maxval*1.2],'fontsize',8)
%     end
%     axis tight
% end
% 
% eval(['export_fig lines_C1_age_vs_clusters_',date,'.pdf']);

toc