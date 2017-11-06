tic
clear all
close all


load 10X43_1_v1.mat
load 10X46_1_v1.mat

cellid_43_1 = cellfun(@(x) ['10X43_1_',x,'-1'], cellid_43_1,'uniformoutput',0);
cellid_46_1 = cellfun(@(x) ['10X46_1_',x,'-1'], cellid_46_1,'uniformoutput',0);




load /mnt/sanger-data2/C1_stuff/DentGyr/DG_V2kit_subsample_samples_merged_18-Jun-2017.mat
% load /data2/C1_stuff/DentGyr/DG_V2kit_subsample_samples_merged_18-Jun-2017.mat

[geneid,ia,ib] = intersect(geneid,geneid_46_1);
data = [data(ia,:), data_43_1(ib,:), data_46_1(ib,:)];
data1 = data(ia,:); data2= [ data_43_1(ib,:), data_46_1(ib,:)]; 
% dm = mean(log2(data1+1),2) - mean(log2(data2+1),2);
% data(abs(dm)>0.5,:) = [];
% geneid(abs(dm)>0.5) = [];

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
% cellid_cluster = loadCellFile('/data2/C1_stuff/DentGyr/DG_10X_V2_markertable_knnMutal_MCL_1p4_18-Jun-2017_manualinspection.txt');
cellid_cluster = cellid_cluster(2:end,:);
cellid_cluster2 = loadCellFile('mnt/sanger-data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
% cellid_cluster2 = loadCellFile('/data2/C1_stuff/DentGyr/cellid_cluster_knnMutal_k20_MCL_1p25_manual_inspection_nov28_2016.txt');
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
% grlev2_uni = {'rgl','astro_old','astro_adol','immature_astro','rgl_young','nipc_young','nipc'.....
%     'astro','nsc1','gran_cycling'};
% grlev2_uni = {'astro','astro_adol','astro_old','rgl','nsc1','immature_astro',....
%     'rgl_young','nipc_young','nipc','gran_cycling','gran_eomes','nb1','gran_sox11','nb2'};

grlev2_uni = {'immature_astro','astro_adol','astro_old','rgl_young','rgl',....
    'nipc_young','nipc','nb1'};

% tf = ismember(grlev2,grlev2_uni);
% [~,p] = ttest2(log2(moldata(:,tf)+1)', log2(moldata(:,~tf)+1)','tail','right');
% in_leave = fdr_proc(p,0.2);
% % % % % % % % % % % % % % % % % % % % % % % % 

n_gr_subsam = 200;
rng(100)
dend_order_names = regexprep(grlev2_uni,'_','-');
T_cells_tmp = zeros(length(grlev2),1);
for i=1:length(grlev2_uni)
    in = find(strcmpi(grlev2, grlev2_uni{i}));
    in = in(randperm(length(in)));
    in = in(1:min([n_gr_subsam,length(in)]));
%     T_cells_tmp( strcmpi(grlev2, grlev2_uni{i}) ) = i;
    T_cells_tmp( in ) = i;
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

% moldata = moldata(in_leave,:);
% geneid = geneid(in_leave);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
moldata_validcells = moldata;
dend_order_names = regexprep(grlev2_uni,'_','-');%grlev2_uni;%
% % Tcelluni = unique(lev2_gr_num_treeorder);
space_gr = 20;
colorvec = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    gr_cent(i) = mean(ind + (i-1)*space_gr);
end


listtoplot = loadCellFile('/mnt/sanger-data2/C1_stuff/DentGyr/list_to_plot_Early_Ngenesis_V2kit_july7_2017.txt');
listtoplot = flipud(listtoplot(:,1));
listtoplot(strcmpi(listtoplot,'-')) = [];

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
    tmp_ave = [];
    if maxval>0
        for i=1:length(T_cells_tmp_uni)
            ind = find(T_cells_tmp==i);
            tmp = moldata_validcells(geneind,ind);
            z = [ [(ind(1) + (i-1)*space_gr-1),(ind(end) + (i-1)*space_gr+0.5)]',repmat(mean(tmp),2,1)];
            z = [z(1,1),0;z;z(2,1),0];
            tmp_ave = [tmp_ave;  z ];
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
        plot(tmp_ave(:,1),tmp_ave(:,2),'m');
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

ax2 = axes('position',[0.1,0.92,0.84,0.08]);
imagesc(~age_bar); hold on;
for i=1:length(age_uni)-1
    plot([0,length(age_bar(1,:))],(i+0.5)*[1,1],'-k');
end
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1:length(age_uni)],'yticklabel',age_uni,'fontsize',6);


eval(['export_fig Early_Ngenesis_SCbarplots_V2kit_',date,'.pdf']);



toc






