tic
clear all
close all

a = loadCellFile('/data2/C1_stuff/DentGyr/example_nIPCfocus_genes_NOCC_KNN40_early_Ngenesis_feb5_2017/nIPC_dissect_noCC_Early_Ngenesis_MKNN40_fracpos_06-Feb-2017.txt');


geneorder = {'Cdk1','Cdk2','Cdk4','Cdk6','Mcm2','Mcm3','Mcm4','Mcm5','Mcm6','Mcm7','Mki67','Aurkb','Cdkn2c','Cenpa','Cenpe','Cenpf',....
    'Cenph','Cenpm','Gmnn','Top2a','Chek1','Chek2','Ccnb1','Ccna2','Ccnd1','Ccne1',....
    'Neurod1','Neurod2','Neurod4','Neurog2','Eomes','Tbr1','Igfbpl1','Sox11','Dcx','Emx1','Emx2',....
    'Calb1','Calb2','Tac2','Gal','Bhlhe22','Hs3st1','Nhlh2',.....
    'Hes5','Ascl1','Tfap2c','Ckap2l','Foxg1','Egfr'...
    'Aldoc','Gfap','Sox2','Sox9','Hes1','Nes','Lpar1','Vnn1',...
    'Etv4','Id3','Olig1','Olig2','Tshz2','Hopx'};
[~,loc] = ismember(upper(geneorder),upper(a(2:end,1)));
a = a(loc+1,:);

genesym = a(2:end,1);
frac_nb_like = cell2mat(a(2:end,2));
frac_rgl_like = cell2mat(a(2:end,3));
nb_like_enrich = frac_nb_like./(frac_nb_like+frac_rgl_like);
rgl_like_enrich = frac_rgl_like./(frac_nb_like+frac_rgl_like);

genesym(isnan(nb_like_enrich)) = [];
nb_like_enrich(isnan(nb_like_enrich)) = [];
rgl_like_enrich(isnan(rgl_like_enrich)) = [];
genesym(isnan(nb_like_enrich)) = [];

figure('position',[100,100,1200,300],'color','w');
axes('position',[0.05,0.2,0.93,0.75])
hb = bar(100*[nb_like_enrich,rgl_like_enrich],'stacked');
set(hb,'barwidth',0.5)
axis tight
set(gca,'xtick',[1:length(genesym)],'xticklabel',genesym,'ytick',[0,25,50,75,100],'XTickLabelRotation',90)
ylabel('% positive cell distribution')

eval(['export_fig nIPC_pos_cells_dist_Early_Ngenesis_NO_CCgenes_MKNN40_',date,'.pdf']);



toc