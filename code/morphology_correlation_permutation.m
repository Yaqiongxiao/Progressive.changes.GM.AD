clear;clc;
hc=load('hc_C1_C2_C3_continue_overlap.txt');
ad1=load('first_C1_C2_C3_continue_overlap.txt');
all_data = hc;
all_data(46:85,:)=ad1;

rdiff = zeros(10,10000);
for i = 1:10000
    idx = randperm(85);
    idhc = idx(1,1:45)';
    idmdd = idx(1,46:85)';
    hcdata = all_data(idhc,:);mdddata = all_data(idmdd,:);
    zid = 1;
    for j = 1:5
        for k = (j+1):5
            if j <=4
                rhc = corr(hcdata(:,j),hcdata(:,k));
                rmdd = corr(mdddata(:,j),mdddata(:,k));
                rdiff(zid,i)=rmdd-rhc;
                zid = zid + 1;
            end
        end
    end
end
save hc_first_differences rdiff


clear;clc;
hc=load('hc_C1_C2_C3_continue_overlap.txt');
ad1=load('second_C1_C2_C3_continue_overlap.txt');
all_data = hc;all_data(46:85,:)=ad1;
rdiff = zeros(10,10000);
for i = 1:10000
    idx = randperm(85);
    idhc = idx(1,1:45)';
    idmdd = idx(1,46:85)';
    hcdata = all_data(idhc,:);mdddata = all_data(idmdd,:);
    zid = 1;
    for j = 1:5
        for k = (j+1):5
            if j <=4
                rhc = corr(hcdata(:,j),hcdata(:,k));
                rmdd = corr(mdddata(:,j),mdddata(:,k));
                rdiff(zid,i)=rmdd-rhc;
                zid = zid + 1;
            end
        end
    end
end
save hc_second_differences rdiff

clear;clc;
hc=load('first_C1_C2_C3_continue_overlap.txt');
ad1=load('second_C1_C2_C3_continue_overlap.txt');
all_data = hc;all_data(41:80,:)=ad1;
rdiff = zeros(10,10000);
for i = 1:10000
    idx = randperm(80);
    idhc = idx(1,1:40)';
    idmdd = idx(1,41:80)';
    hcdata = all_data(idhc,:);mdddata = all_data(idmdd,:);
    zid = 1;
    for j = 1:5
        for k = (j+1):5
            if j <=4
                rhc = corr(hcdata(:,j),hcdata(:,k));
                rmdd = corr(mdddata(:,j),mdddata(:,k));
                rdiff(zid,i)=rmdd-rhc;
                zid = zid + 1;
            end
        end
    end
end
save first_second_differences rdiff
    
