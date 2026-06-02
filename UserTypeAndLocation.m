%%
clc; clear; close all;
load('PGEData.mat')
%%
PeakEnerCon = zeros(NumUsers,NumMonths);
for month = 1:NumMonths
    PeakEnerCon(:,month) = max(EnerCon(:,MonthStartMet(month):MonthEndMet(month)),[],2);
end
%%
PeakPowCon = 4*PeakEnerCon;
MeanPeakPowCon = mean(PeakPowCon,2);
UserType = zeros(NumUsers,1);
for user = 1:NumUsers
    if MeanPeakPowCon(user) <= 75
        UserType(user) = 1;
    elseif MeanPeakPowCon(user) > 499
        UserType(user) = 3;
    else
        UserType(user) = 2;
    end
end

%%

%%
sum(double(UserType==1))
sum(double(UserType==2))
sum(double(UserType==3))
%%
load('SFMap.mat')
NumUsersPerZip = [25,49,12,35,63,32,19,21,53,7,5,8,9,4,13,5,16,22,36,1,2,2,26,5,0,0,0];
NumUsers = sum(NumUsersPerZip);
rng('default')
MeshLen = 100;
NumUsersPerZipCumSum = [0,cumsum(NumUsersPerZip)];
BuildingLoc = zeros(NumUsers,2);
for i = 1:24
    lat = SFGeometry{i,1}(:,1);
    lon = SFGeometry{i,1}(:,2);
    polyin = polyshape({lat},{lon});

    latMin = min(lat); latMax = max(lat);
    lonMin = min(lon); lonMax = max(lon);
    xp = linspace(latMin,latMax,MeshLen);
    yp = linspace(lonMin,lonMax,MeshLen);

    in = inpolygon(repmat(xp',1,MeshLen),repmat(yp,MeshLen,1),lat,lon);
    x = find(double(in) == 1);
    y = randperm(length(x),NumUsersPerZip(i));
    yMat = ceil(x(y)/MeshLen);
    xMat = x(y) - (yMat-1)*MeshLen;

    BuildingLoc(NumUsersPerZipCumSum(i)+1:NumUsersPerZipCumSum(i+1),1) = xp(xMat)';
    BuildingLoc(NumUsersPerZipCumSum(i)+1:NumUsersPerZipCumSum(i+1),2) = yp(yMat)';
end
%%
z = find(UserType==1);
EnerConType2_3 = EnerCon;
EnerConType2_3(z,:) = [];
%%
close all;
set(0,'defaultpatchlinewidth',2);
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',14);
fig = figure(1);
fig.Position = [0 0 400 300];
for i = 25:27
    lat = SFGeometry{i,1}(:,1);
    lon = SFGeometry{i,1}(:,2);
    fig = figure(1);
    fig.Position = [0 0 400 300];
    fill(lat,lon,[0.5 0.5 0.5],'FaceAlpha','0.3','LineWidth',1)
    hold on;
end
for i = 1:24
    lat = SFGeometry{i,1}(:,1);
    lon = SFGeometry{i,1}(:,2);
    fig = figure(1);
    fig.Position = [0 0 400 300];
    plot(lat,lon,'-k','LineWidth',1)
    hold on;
    xlim([-122.53 -122.35])
    ylim([37.7 37.85])
    xlabel('Latitude')
    ylabel('Longitude')
end
for usertype = 1:3
    z = find(UserType == usertype);
    if usertype == 1
        h1 = plot(BuildingLoc(z,1),BuildingLoc(z,2),'x','MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'c');
    elseif usertype == 2
        h2 = plot(BuildingLoc(z,1),BuildingLoc(z,2),'x','MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'b');
    else
        h3 = plot(BuildingLoc(z,1),BuildingLoc(z,2),'x','MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    end
    hold on
end
% legend([h1,h2,h3],{'Small (B-1)','Medium (B-10)','Large (B-19)'},'Orientation','horizontal','NumColumns',2,'Location','northwest')
%%
UserZip = repelem(AllZIP,NumUsersPerZip);
UserInfo(:,1) = UserZip';
UserInfo(:,2:3) = BuildingLoc;
UserInfo(:,4) = UserType;
UserInfo = array2table(UserInfo,'VariableNames',{'ZIP','x','y','Type'});
%%
% save('BuildingDis','xRelCordinates','yRelCordinates');

%% Information about medium and large users only
UserInfoType2_3 = UserInfo(UserInfo.Type ==2 | UserInfo.Type == 3,:);
CSx = mean(UserInfoType2_3.x);
CSy = mean(UserInfoType2_3.y);
UserInfoType2 = UserInfo(UserInfo.Type ==2,:);
UserInfoType3 = UserInfo(UserInfo.Type ==3,:);
%%
close all;
set(0,'defaultpatchlinewidth',2);
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',14);
fig = figure(1);
fig.Position = [0 0 400 300];
for i = 25:27
    lat = SFGeometry{i,1}(:,1);
    lon = SFGeometry{i,1}(:,2);
    fill(lat,lon,[0.5 0.5 0.5],'FaceAlpha','0.3','LineWidth',1)
    hold on;
end
for i = 1:24
    lat = SFGeometry{i,1}(:,1);
    lon = SFGeometry{i,1}(:,2);
    plot(lat,lon,'-k','LineWidth',1)
    hold on;
    xlim([-122.53 -122.35])
    ylim([37.7 37.85])
    xlabel('Latitude')
    ylabel('Longitude')
end
h2 = plot(UserInfoType2.x,UserInfoType2.y,'x','MarkerSize', 4, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
hold on
h3 = plot(UserInfoType3.x,UserInfoType3.y,'x','MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
hold on
char = plot(CSx,CSy,'+','MarkerSize', 4, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
hold off;
legend([h2,h3,char],{'B-10 users','B-19 users','Charging station'},'Orientation','horizontal','NumColumns',1,'Location','northwest')
%% Haversine Distance
NumUsersType2_3 = size(UserInfoType2_3,1);
UserInfoType2_3.xRel = haversine(CSx*ones(NumUsersType2_3,1),CSy*ones(NumUsersType2_3,1),...
                        UserInfoType2_3.x,CSy*ones(NumUsersType2_3,1));
UserInfoType2_3.yRel = haversine(CSx*ones(NumUsersType2_3,1),CSy*ones(NumUsersType2_3,1),...
                        CSx*ones(NumUsersType2_3,1),UserInfoType2_3.y);
%%
% save('PGEDataMod','EnerConType2_3','UserInfoType2_3', 'DaysInMonth','MonthStartMet','MonthEndMet','NumMonths')