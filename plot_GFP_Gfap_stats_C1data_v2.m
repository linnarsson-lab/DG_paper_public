tic
clear all
close all


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

a = loadCellFile('/data2/C1_stuff/DentGyr/table_GFP_Gfap_Cdk1_per_cluster_23-Feb-2017.txt');


grname = a(1,2:end);
sumCells = cell2mat(a(2,2:end));
sumGFP = cell2mat(a(3,2:end));
sumGfap = cell2mat(a(4,2:end));
sumGFP_Gfap = cell2mat(a(5,2:end));
sumCdk1 = cell2mat(a(6,2:end));

colorvec = distinguishable_colors(length(grname));
hf = figure('color','w','position',[200,200,1000,400],'visible','on');
subplot(1,4,1);
hold on;
tmp = sumGFP./sumCells*100;
for i=1:length(grname)
    plot(rand(1), tmp(i), '.','color',colorvec(i,:),'markersize',20);
end
set(gca,'xlim',[-1,2],'xtick',[]);
ylabel('%hGfap-GFP','fontsize',12);

subplot(1,4,2);
hold on;
tmp = sumGfap./sumCells*100;
for i=1:length(grname)
    plot(rand(1), tmp(i), '.','color',colorvec(i,:),'markersize',20);
end
set(gca,'xlim',[-1,2],'xtick',[]);
ylabel('%Gfap','fontsize',12);

subplot(1,4,3);
hold on;
tmp = sumGFP_Gfap./sumGFP*100;
for i=1:length(grname)
    plot(rand(1), tmp(i), '.','color',colorvec(i,:),'markersize',20);
end
set(gca,'xlim',[-1,2],'xtick',[]);
ylabel('%hGfap-GFP & Gfap (of GFP)','fontsize',12);

subplot(1,4,4);
hold on;
tmp = sumCdk1./sumCells*100;
for i=1:length(grname)
    plot(rand(1), tmp(i), '.','color',colorvec(i,:),'markersize',20);
end
set(gca,'xlim',[-1,2],'xtick',[]);
ylabel('%Cdk1','fontsize',12);
legend(names_as_in_paper)

eval(['export_fig C1data_Gfap_GFP_stats_',date,'.pdf']);


toc