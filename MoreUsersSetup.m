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
EnerCon = [EnerConType2_3;EnerCon(z,:)];

%%
UserZip = repelem(AllZIP,NumUsersPerZip);
UserInfo(:,1) = UserZip';
UserInfo(:,2:3) = BuildingLoc;
UserInfo(:,4) = UserType;
UserInfo = array2table(UserInfo,'VariableNames',{'ZIP','x','y','Type'});

%% Information about medium and large users only
UserInfoType2_3 = UserInfo(UserInfo.Type ==2 | UserInfo.Type == 3,:);
CSx = mean(UserInfoType2_3.x);
CSy = mean(UserInfoType2_3.y);
UserInfo = [UserInfoType2_3;UserInfo(UserInfo.Type == 1,:)];

%% Haversine Distance
UserInfo.xRel = haversine(CSx*ones(NumUsers,1),CSy*ones(NumUsers,1),...
                        UserInfo.x,CSy*ones(NumUsers,1));
UserInfo.yRel = haversine(CSx*ones(NumUsers,1),CSy*ones(NumUsers,1),...
                        CSx*ones(NumUsers,1),UserInfo.y);
%% Rescale EnerCon of users who are of type 1
ModUserTypeUser1 = 2+ binornd(1,sum(double(UserType==3))/(sum(double(UserType==3))+sum(double(UserType==2))),sum(double(UserType==1)),1);
UserInfo.Type(UserInfo.Type == 1,:) = ModUserTypeUser1;
%%
MeanPeakPowCon = [MeanPeakPowCon(find(UserType==2 | UserType==3));MeanPeakPowCon(z,:)];
PeakPowCon = [PeakPowCon(find(UserType==2 | UserType==3),:);PeakPowCon(z,:)];
RescaleFactor = zeros(length(z),1);
for n = 1:length(z)
    if ModUserTypeUser1(n) == 2
        RescaleFactor(n) = 75 + (MeanPeakPowCon(136+n)/75)*(500-75);
    else
        RescaleFactor(n) = 500 + (MeanPeakPowCon(136+n)/75)*(1000-500);
    end
end
EnerCon(137:end,:) = EnerCon(137:end,:).*(RescaleFactor./MeanPeakPowCon(137:end,:));

%%
save('PGEModMoreUsers','EnerCon','UserInfo','DaysInMonth','MonthStartMet','MonthEndMet','NumMonths')