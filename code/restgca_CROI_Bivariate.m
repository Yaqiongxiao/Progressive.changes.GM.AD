function [Result_X2Y,Result_Y2X,ROI_sequence,F_X2Y,F_Y2X] = restgca_CROI_Bivariate(theROITimeCourses,Order)
nDim4=length(theROITimeCourses);
numROIs=size(theROITimeCourses,2);
ROI_sequence=combntns(1:numROIs,2);
Past_1=zeros(nDim4-Order,Order);
Past_2=zeros(nDim4-Order,Order);
Result_X2Y=zeros(Order*2,size(ROI_sequence,1))';
Result_Y2X=zeros(Order*2,size(ROI_sequence,1))';
F_X2Y=zeros(1,size(ROI_sequence,1))';
F_Y2X=zeros(1,size(ROI_sequence,1))';
for i=1:size(ROI_sequence,1),
    ROI_used=theROITimeCourses(:,ROI_sequence(i,:));
    Now=ROI_used(Order+1:end,:);
    for j=1:Order,
        Past_1(:,j)=ROI_used(j:nDim4-Order+j-1,1);
        Past_2(:,j)=ROI_used(j:nDim4-Order+j-1,2);
        Regressors1=[Past_1,Past_2];
        Regressors2=[Past_2,Past_1];
    end
    [b_1,~,~,~, T1]=rest_regress_ss(Now(:,2),Regressors1);
    [b_2,~,~,~, T2]=rest_regress_ss(Now(:,1),Regressors2);
    Result_X2Y(i,:)=b_1(1:Order*2);
    Result_Y2X(i,:)=b_2(1:Order*2);
    F_X2Y(i,:)=T1(1);
    F_Y2X(i,:)=T2(1);
end
