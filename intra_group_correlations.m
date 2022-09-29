clear;clc;
hc=load('hc_C1_C2_C3_continue_overlap.txt');
ad1=load('first_C1_C2_C3_continue_overlap.txt');
ad2=load('second_C1_C2_C3_continue_overlap.txt');

rhcad1ad2 = zeros(3,10);
zid = 1;
for i = 1:4
    for j = (i+1):5
        rhcad1ad2(1,zid) = corr(hc(:,i),hc(:,j));
        rhcad1ad2(2,zid) = corr(ad1(:,i),ad1(:,j));
        rhcad1ad2(3,zid) = corr(ad2(:,i),ad2(:,j));
        zid = zid + 1;
    end
end

save intra_group_correlations_hc_ad1_ad2 rhcad1ad2