
tic
clear all
close all


files_to_load = {'10X43_1_DG_01-Nov-2016.mat'
    '10X46_1_DG_01-Nov-2016.mat'};
    
ncells_files = zeros(length(files_to_load),1);
for i=1:length(files_to_load)
    fprintf(['load ',num2str(i),'\n']);
    load(['/mnt/sanger-data2/10X_adol_mice/loom_data_files/mat_data_files/',files_to_load{i}]);
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

in = find(sum(data>0,2)>10 & sum(data>0,2)<length(data(1,:))*1);


data = data(in,:);
geneid = geneid(in);
moldata = data;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 
% % permute the cells order
% rng(1);
% cell_perm = randperm(length(cellid));
% gene_perm = randperm(length(geneid));
% geneid = geneid(gene_perm);
% cellid = cellid(cell_perm);
% source = source(cell_perm);
% age = age(cell_perm);
% moldata = data(gene_perm,cell_perm);
% moldata = moldata./repmat(sum(moldata),length(data(:,1)),1)*1e4;


cellid_cluster = loadCellFile('/mnt/sanger-data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
cellid_cluster = cellid_cluster(2:end,:);
[~,loc] = ismember(cellid_cluster(:,1), cellid);
cellid = cellid(loc);
source = source(loc);
age = age(loc);
moldata = moldata(:,loc);


grlev2 = cellid_cluster(:,3);
% lev2uni = unique(lev2gr);
% grlev2_uni = {'endo','peric','vlmc','mgl','pvm','OL','nfol','opc','nsc2','astro','nsc1','gran_cycling','gran_eomes',...
%     'gran_sox11','gran_fxyd7','gran_mature','cck_tox','mossy_adcyap1','mossy_cy26b1','mossy_klk8','gaba_cnr1','gaba_lhx6','reln' };
grlev2_uni = flipud(regexprep(loadCellFile('/mnt/sanger-data2/C1_stuff/DentGyr/DG_10X_dendrogram_order.txt'),'-','_'));
names_as_in_paper = flipud({'Cajal-Retzius'
    'GABA-Lhx6'
    'GABA-Cnr1'
    'Mossy-Klk8'
    'Mossy-Cyp26b1'
    'Mossy-Adcyap1'
    'Cck-Tox'
    'Granule-mature'
    'Granule-immature'
    'Neuroblast 2'
    'Neuroblast 1'
    'nIPC'
    'Radial Glia-like'
    'Astrocytes'
    'OPC'
    'NFOL'
    'OL'
    'PVM'
    'Microglia'
    'VLMC'
    'Pericytes'
    'Endothelial'});


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


a = loadCellFile('/mnt/sanger-data2/C1_stuff/DentGyr/tsne_coordinates_finalClsuter_perplexity_50_Ngenes=648_NPCA=80_19-Jan-2017.txt');
[~,loc] = ismember(cellid,a(:,1));
a = a(loc,:);
mappedX_cell = cell2mat(a(:,2:3));
gr_label_xy = [-52.38,50.4
    -43,56
    -37,59
    43,49
    33,57
    -70,6.5
    -65,-6
    -58,-13
    -19,-45
    -1,-60
    30,-50
    33,-46
    46,-21
    20,8
    5,30
    5,48
    10,58
    7,61
    3,56
    44,20
    44,30
    -52,-37];                           

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
fname = dir('10X_level2_singlecell_barplots_all_genes_jan20_2017/*png');
fname = {fname.name}';
fname = regexprep(fname,'_level2_group.png','');
genesleft = setdiff(geneid,fname);
[~,loc] = ismember(genesleft,geneid);
moldata = moldata(loc,:);
geneid = geneid(loc);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
moldata_validcells = moldata;
dend_order_names = names_as_in_paper;%regexprep(grlev2_uni,'_','-');
% Tcelluni = unique(lev2_gr_num_treeorder);
space_gr = 20;
colorvec = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    gr_cent(i) = mean(ind + (i-1)*space_gr);
end

yax_lim = [0,space_gr*length(T_cells_tmp_uni) + length(moldata_validcells(1,:))];
poolobj = parpool(8);
tic
parfor k=1:length(moldata_validcells(:,1))
    k
    hf = figure('color','w','position',[200,200,1000,800],'visible','off');
    %     axes('position',[0.395,0.05,0.58,0.92])
    subplot(1,2,1);
    hold on;
    
    maxval = 0.8*max(moldata_validcells(k,:));
    if maxval>0
        for i=1:length(T_cells_tmp_uni)
            ind = find(T_cells_tmp==i);
            plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
            plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
            h = barh(ind + (i-1)*space_gr,moldata_validcells(k,ind));
            set(h,'facecolor',colorvec(i,:),'edgecolor',colorvec(i,:),'barwidth',1);
        end
    else
        maxval = 1;
        for i=1:length(T_cells_tmp_uni)
            ind = find(T_cells_tmp==i);
            plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
            plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
        end
    end
    set(gca,'position',[0.1,0.05,0.2,0.92])
    set(gca,'ytick',gr_cent,'yticklabel',dend_order_names,'fontsize',7,'xlim',[0,maxval],'ylim',yax_lim,'ydir','reverse')
    xlabel(geneid{k},'fontsize',10);
    %     yl = get(gca,'ylim');
    %     disp(num2str(yl))
    set(gca,'ylim',yax_lim);
    
    
    % % % % % % % % % % % % % % % % % % % % % % % %     tsne part
    genePlot = geneid{k};
    markergene = moldata_validcells(k,:);%(moldata(strcmpi(geneid,genePlot),:));
    inpos = markergene>0;
    tmpthlow = 0;%prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),95);
%     if tmpthlow==tmpthhigh
%         tmpthlow==0;
%     end
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    if length(unique(markergene))>1
        %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
        %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
        markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
            interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
            ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
        subplot(1,2,2);
        scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),70,markergene_color(~inpos,:),'.'); hold on;
        scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),70,markergene_color(inpos,:),'.'); hold on;
    end
    for i=1:length(T_cells_tmp_uni)
        text(gr_label_xy(i,1),gr_label_xy(i,2),dend_order_names{i},'fontsize',8)%,'backgroundcolor',0.8*[1,1,1]
    end
    set(gca,'position',[0.33,0.05,0.65,0.92])
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    axis tight
    axis off
    % % % % % % % % % % % % % % % % % % % % % % % % % %
    set(gcf,'Renderer','OpenGL')
    %     eval(['export_fig -r200 level2_barplots_all_genes_dec25_2015/',geneuni_1{k},'_level2_group.png']);
    %     export_fig(gcf,'-r200',['level2_barplots_all_genes_dec25_2015_v2/',geneuni_1{k},'_level2_group.png']);
    set(gcf,'PaperPositionMode','auto')
    %     print(['level2_barplots_all_genes_dec25_2015_v2/',geneuni_1{k},'_level2_group'],'-dpng','-r0')
    %     save2png(['level2_singlecell_barplots_all_genes_jan14_2017_v1/',geneuni_1{k},'_level2_group'],hf,200)
    save2png(['10X_level2_singlecell_barplots_all_genes_jan20_2017/',geneid{k},'_level2_group'],hf,200)
    close(hf)
end
toc
delete(poolobj)




% % % % % % % % % % % % % % % % % % % 
hf = figure('color','w','position',[200,200,1000,800],'visible','on');
subplot(1,2,1);
hold on;
maxval = 1;
for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
    plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
end


set(gca,'position',[0.1,0.05,0.2,0.92])
set(gca,'ytick',gr_cent,'yticklabel',dend_order_names,'fontsize',7,'xlim',[0,maxval],'ylim',yax_lim,'ydir','reverse')
set(gca,'ylim',yax_lim);
text(1,yax_lim(2)/4,'Gene not found','fontsize',40,'color','r')

% % % % % % % % % % % % % % % % % % % % % % % %     tsne part
subplot(1,2,2);
scatter(mappedX_cell(:,1),mappedX_cell(:,2),70,[0.7,0.7,0.7],'.'); hold on;

for i=1:length(T_cells_tmp_uni)
    text(gr_label_xy(i,1),gr_label_xy(i,2),dend_order_names{i},'fontsize',8)%,'backgroundcolor',0.8*[1,1,1]
end
set(gca,'position',[0.33,0.05,0.65,0.92])
set(gca,'xlim',[-150,150],'ylim',[-150,150])
axis tight
axis off
set(gcf,'Renderer','OpenGL')
set(gcf,'PaperPositionMode','auto')
save2png(['10X_level2_singlecell_barplots_all_genes_jan20_2017/','Empty','_level2_group'],hf,200)




% % % % % % % % % % % % % % % % % % % 
hf = figure('color','w','position',[198,101,1198,901],'visible','on');
subplot(1,2,1);
hold on;
maxval = 1;
for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
    plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
end


set(gca,'position',[0.1,0.05,0.2,0.92])
set(gca,'ytick',gr_cent,'yticklabel',dend_order_names,'fontsize',7,'xlim',[0,maxval],'ylim',yax_lim,'ydir','reverse')
set(gca,'ylim',yax_lim);
text(1,yax_lim(2)/4,'Gene not found','fontsize',40,'color','r')

% % % % % % % % % % % % % % % % % % % % % % % %     tsne part
subplot(1,2,2);
scatter(mappedX_cell(:,1),mappedX_cell(:,2),70,[0.7,0.7,0.7],'.'); hold on;

for i=1:length(T_cells_tmp_uni)
    text(gr_label_xy(i,1),gr_label_xy(i,2),dend_order_names{i},'fontsize',8)%,'backgroundcolor',0.8*[1,1,1]
end
set(gca,'position',[0.33,0.05,0.65,0.92])
set(gca,'xlim',[-150,150],'ylim',[-150,150])
axis tight
axis off
set(gcf,'Renderer','OpenGL')
set(gcf,'PaperPositionMode','auto')
save2png(['10XV2kit_level2_singlecell_barplots_all_genes_sept5_2017/','Empty','_level2_group'],hf,200)





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % make barplots for molecules and genes
total_mol = sum(moldata_validcells);
total_genes = sum(moldata_validcells>0);
% % % % % % % % % % % % % % % % 

hf = figure('color','w','position',[198,101,1198,901],'visible','on');
%     axes('position',[0.395,0.05,0.58,0.92])
subplot(1,2,1);
hold on;
gr_ave = zeros(length(T_cells_tmp_uni),1);
gr_std = zeros(length(T_cells_tmp_uni),1);
maxval = 0.9*max(total_mol);
if maxval>0
    for i=1:length(T_cells_tmp_uni)
        ind = find(T_cells_tmp==i);
        plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
        plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
        h = barh(ind + (i-1)*space_gr,total_mol(ind));
        set(h,'facecolor',colorvec(i,:),'edgecolor',colorvec(i,:),'barwidth',1);
        gr_ave(i) = mean(total_mol(ind));
        gr_std(i) = std(total_mol(ind));
    end
else
    maxval = 1;
    for i=1:length(T_cells_tmp_uni)
        ind = find(T_cells_tmp==i);
        plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
        plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
    end
end
h = herrorbar(gr_ave, gr_cent,gr_std,'.');
set(h,'color',[0.1,0.9,0.05],'linewidth',3);
set(gca,'position',[0.1,0.05,0.2,0.92])
set(gca,'ytick',gr_cent,'yticklabel',dend_order_names,'fontsize',7,'xlim',[0,maxval],'ylim',yax_lim,'ydir','reverse')
xlabel('total-mol','fontsize',10);
%     yl = get(gca,'ylim');
%     disp(num2str(yl))
set(gca,'ylim',yax_lim);


% % % % % % % % % % % % % % % % % % % % % % % %     tsne part

markergene = total_mol;%(moldata(strcmpi(geneid,genePlot),:));
inpos = markergene>0;
tmpthlow = 0;%prctile(markergene(markergene>0),1);
tmpthhigh = prctile(markergene(markergene>0),98);
%     if tmpthlow==tmpthhigh
%         tmpthlow==0;
%     end
markergene(markergene>tmpthhigh) = tmpthhigh;
markergene(markergene<tmpthlow) = tmpthlow;
c_rgb = [1,0,0];rand([1,3]);
if length(unique(markergene))>1
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(1,2,2);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),70,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),70,markergene_color(inpos,:),'.'); hold on;
end
for i=1:length(T_cells_tmp_uni)
    text(gr_label_xy(i,1),gr_label_xy(i,2),dend_order_names{i},'fontsize',8);%,'backgroundcolor',0.8*[1,1,1]htext(i) =
    %         in = T_cells_tmp == T_cells_tmp_uni(i);
    %         htext(i) = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)),'fontsize',8)%,'backgroundcolor',0.8*[1,1,1]
end
set(gca,'position',[0.33,0.05,0.65,0.92])
set(gca,'xlim',[-150,150],'ylim',[-150,150])
title('total-mol');
axis tight
axis off
% % % % % % % % % % % % % % % % % % % % % % % % % %
set(gcf,'Renderer','OpenGL')
set(gcf,'PaperPositionMode','auto')

eval(['export_fig 10X_level2_singlecell_barplots_all_genes_jan20_2017/barplot_DG_10XV1_totmol_percell_',date,'.pdf'])



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
hf = figure('color','w','position',[198,101,1198,901],'visible','on');
%     axes('position',[0.395,0.05,0.58,0.92])
subplot(1,2,1);
hold on;
gr_ave = zeros(length(T_cells_tmp_uni),1);
gr_std = zeros(length(T_cells_tmp_uni),1);
maxval = 0.9*max(total_genes);
if maxval>0
    for i=1:length(T_cells_tmp_uni)
        ind = find(T_cells_tmp==i);
        plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
        plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
        h = barh(ind + (i-1)*space_gr,total_genes(ind));
        set(h,'facecolor',colorvec(i,:),'edgecolor',colorvec(i,:),'barwidth',1);
        gr_ave(i) = mean(total_genes(ind));
        gr_std(i) = std(total_genes(ind));
    end
else
    maxval = 1;
    for i=1:length(T_cells_tmp_uni)
        ind = find(T_cells_tmp==i);
        plot([0,maxval],(ind(1) + (i-1)*space_gr-1)*[1,1], '--','color',220/255*[1,1,1]);
        plot([0,maxval],(ind(end) + (i-1)*space_gr+0.5)*[1,1], '--','color',220/255*[1,1,1]);
    end
end
h = herrorbar(gr_ave, gr_cent,gr_std,'.');
set(h,'color',[0.1,0.9,0.05],'linewidth',3);
set(gca,'position',[0.1,0.05,0.2,0.92])
set(gca,'ytick',gr_cent,'yticklabel',dend_order_names,'fontsize',7,'xlim',[0,maxval],'ylim',yax_lim,'ydir','reverse')
xlabel('total-genes','fontsize',10);
%     yl = get(gca,'ylim');
%     disp(num2str(yl))
set(gca,'ylim',yax_lim);


% % % % % % % % % % % % % % % % % % % % % % % %     tsne part

markergene = total_genes;%(moldata(strcmpi(geneid,genePlot),:));
inpos = markergene>0;
tmpthlow = 0;%prctile(markergene(markergene>0),1);
tmpthhigh = prctile(markergene(markergene>0),98);
%     if tmpthlow==tmpthhigh
%         tmpthlow==0;
%     end
markergene(markergene>tmpthhigh) = tmpthhigh;
markergene(markergene<tmpthlow) = tmpthlow;
c_rgb = [1,0,0];rand([1,3]);
if length(unique(markergene))>1
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(1,2,2);
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),70,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),70,markergene_color(inpos,:),'.'); hold on;
end
for i=1:length(T_cells_tmp_uni)
    text(gr_label_xy(i,1),gr_label_xy(i,2),dend_order_names{i},'fontsize',8);%,'backgroundcolor',0.8*[1,1,1]htext(i) =
    %         in = T_cells_tmp == T_cells_tmp_uni(i);
    %         htext(i) = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)),'fontsize',8)%,'backgroundcolor',0.8*[1,1,1]
end
set(gca,'position',[0.33,0.05,0.65,0.92])
set(gca,'xlim',[-150,150],'ylim',[-150,150])
title('total-genes');
axis tight
axis off
% % % % % % % % % % % % % % % % % % % % % % % % % %
set(gcf,'Renderer','OpenGL')
set(gcf,'PaperPositionMode','auto')

eval(['export_fig 10X_level2_singlecell_barplots_all_genes_jan20_2017/barplot_DG_10XV1_totgenes_percell_',date,'.pdf'])





% % % % % % % % % % % % % % % % % % plot excluded genes regarding revision
% minor comment
list = {'Rpl26', 'Rpl35a', 'Erh','Gstp1','Slc25a5','Pgk1','Eno1' ,'Tubb2a','Emc4','Scg5',};
figure('visible','on','position',[100,73,950,520],'color','w');
    [ha, pos] = tight_subplot(2, 5, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
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
%     subplot(3,3,i);
    set(gcf,'CurrentAxes',ha(i))
    scatter(mappedX_cell(~inpos,1),mappedX_cell(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mappedX_cell(inpos,1),mappedX_cell(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
%     for ii=1:length(T_cells_tmp_uni)
%         in = T_cells_tmp==T_cells_tmp_uni(ii);
%         %     ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),num2str(T_cells_tmp_uni(i)));
% %         ht = text(median(mappedX_cell(in,1)),median(mappedX_cell(in,2)),regexprep(grlev2_uni(ii),'_','-'));
% %         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
%         plot(cluster_bound{ii}(:,1),cluster_bound{ii}(:,2));
%     end
    axis equal
    axis tight
    axis off
end
eval(['export_fig tsne_DG_10XV1_RemovedGenes1_cells','_',date,'.pdf']);

eval(['export_fig tsne_DG_10XV1_RemovedGenes1_cells','_',date,'.png']);

toc










