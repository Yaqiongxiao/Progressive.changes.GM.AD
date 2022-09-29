%%%calculate coefficience GCA
clear
clc

%%%%计算结构因果协变
roits = load('first_C1_C2_C3_continue_overlap.txt');
[Result_X2Y,Result_Y2X,ROI_sequence] = restgca_CROI_Bivariate(roits,1,[]);
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\';
outname1 = [outpath,'First_X2Y'];
outname2 = [outpath,'First_Y2X'];
save(outname1,'Result_X2Y');
save(outname2,'Result_Y2X');

clear
clc

roits = load('second_C1_C2_C3_continue_overlap.txt');
[Result_X2Y,Result_Y2X,ROI_sequence] = restgca_CROI_Bivariate(roits,1,[]);
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\';
outname1 = [outpath,'Second_X2Y'];
outname2 = [outpath,'Second_Y2X'];
save(outname1,'Result_X2Y');
save(outname2,'Result_Y2X');


clear
clc

roits = load('hc_C1_C2_C3_continue_overlap.txt');
[Result_X2Y,Result_Y2X,ROI_sequence] = restgca_CROI_Bivariate(roits,1,[]);
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\';
outname1 = [outpath,'HC_X2Y'];
outname2 = [outpath,'HC_Y2X'];
save(outname1,'Result_X2Y');
save(outname2,'Result_Y2X');

clear;clc;
hc=load('hc_C1_C2_C3_continue_overlap.txt');
ad1=load('first_C1_C2_C3_continue_overlap.txt');
all_data = hc;all_data(46:85,:)=ad1;
rdiff_X2Y = zeros(10,10000);
rdiff_Y2X = zeros(10,10000);
for i = 1:10000
    idx = randperm(85);
    idhc = idx(1,1:45)';
    idad = idx(1,46:85)';
    hcdata = all_data(idhc,:);addata = all_data(idad,:);
    [hc_rx2y,hc_ry2x] = restgca_CROI_Bivariate(hcdata,1,[]); 
    [ad_rx2y,ad_ry2x] = restgca_CROI_Bivariate(addata,1,[]); 
    rdiff_X2Y(:,i)=(ad_rx2y(:,1)')-(hc_rx2y(:,1)');
    rdiff_Y2X(:,i)=(ad_ry2x(:,1)')-(hc_ry2x(:,1)');
end
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\';
outname1 = [outpath,'AD1-HC_X2Y'];
outname2 = [outpath,'AD1-HC_Y2X'];
save(outname1,'rdiff_X2Y');
save(outname2,'rdiff_Y2X');


clear;clc;
hc=load('hc_C1_C2_C3_continue_overlap.txt');
ad1=load('second_C1_C2_C3_continue_overlap.txt');
all_data = hc;all_data(46:85,:)=ad1;
rdiff_X2Y = zeros(10,10000);
rdiff_Y2X = zeros(10,10000);
for i = 1:10000
    idx = randperm(85);
    idhc = idx(1,1:45)';
    idad = idx(1,46:85)';
    hcdata = all_data(idhc,:);addata = all_data(idad,:);
    [hc_rx2y,hc_ry2x] = restgca_CROI_Bivariate(hcdata,1,[]); 
    [ad_rx2y,ad_ry2x] = restgca_CROI_Bivariate(addata,1,[]); 
    rdiff_X2Y(:,i)=(ad_rx2y(:,1)')-(hc_rx2y(:,1)');
    rdiff_Y2X(:,i)=(ad_ry2x(:,1)')-(hc_ry2x(:,1)');
end
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\';
outname1 = [outpath,'AD2-HC_X2Y'];
outname2 = [outpath,'AD2-HC_Y2X'];
save(outname1,'rdiff_X2Y');
save(outname2,'rdiff_Y2X');



clear;clc;
hc=load('first_C1_C2_C3_continue_overlap.txt');
ad1=load('second_C1_C2_C3_continue_overlap.txt');
all_data = hc;all_data(41:80,:)=ad1;
rdiff_X2Y = zeros(10,10000);
rdiff_Y2X = zeros(10,10000);
for i = 1:10000
    idx = randperm(80);
    idhc = idx(1,1:40)';
    idad = idx(1,41:80)';
    hcdata = all_data(idhc,:);addata = all_data(idad,:);
    [hc_rx2y,hc_ry2x] = restgca_CROI_Bivariate(hcdata,1,[]); 
    [ad_rx2y,ad_ry2x] = restgca_CROI_Bivariate(addata,1,[]); 
    rdiff_X2Y(:,i)=(ad_rx2y(:,1)')-(hc_rx2y(:,1)');
    rdiff_Y2X(:,i)=(ad_ry2x(:,1)')-(hc_ry2x(:,1)');
end
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\';
outname1 = [outpath,'AD2-AD1_X2Y'];
outname2 = [outpath,'AD2-AD1_Y2X'];
save(outname1,'rdiff_X2Y');
save(outname2,'rdiff_Y2X');




%%%%计算结构因果协变以及输出T值
clear
clc

roits = load('first_second_overlapin1st_C123_continue2nd.txt');
[Result_X2Y,Result_Y2X,ROI_sequence,F_X2Y,F_Y2X] = restgca_CROI_Bivariate(roits22,1);
PX2Y = 1-fcdf(F_X2Y.*F_X2Y,1,37);
PY2X = 1-fcdf(F_Y2X.*F_Y2X,1,37);
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\First_2nd_overlap_to_C123continue\';
outname1 = [outpath,'First2ndoverlap_2_seconddecrease'];
save(outname1,'Result_X2Y');
outname2 = [outpath,'Pvalues_First2ndoverlap_2_seconddecrease'];
save(outname2,'PX2Y');

%%%%只计算结构因果协变通过置换确定平均作用
clear;clc;
roits = load('first_second_overlapin1st_C123_continue2nd.txt');
rdiff_X2Y = zeros(10,10000);
rdiff_Y2X = zeros(10,10000);
for i = 1:10000
    idx = randperm(40);
    newdata = roits(idx',:);
    [hc_rx2y,hc_ry2x] = restgca_CROI_Bivariate(newdata,1); 
    rdiff_X2Y(:,i)=hc_rx2y(:,1);
    rdiff_Y2X(:,i)=hc_ry2x(:,1);
end
outpath='E:\ADNI_longitudinal_data\GMV_skeleton\CaSCN\';
outname1 = [outpath,'First_2nd_overlap_to_C123continue_AD2-AD1_X2Y'];
outname2 = [outpath,'First_2nd_overlap_to_C123continue_AD2-AD1_Y2X'];
save(outname1,'rdiff_X2Y');
save(outname2,'rdiff_Y2X');




% hist(a(:,4),100)
% 
% 
% x=[b(4),b(4)]
% y=[0,330]
% 
% hold on
% plot(x,y,'r')
% ylim([0 350])
% 
