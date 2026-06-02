%% load data
clc; clear all; close all;
% StartDate: May 01st 2014 (Thursday)
% EndDate: July 31st 2016 (Sunday)
ZIP = [94102,94103,94104,94105,94107,94108,94109,94110,94111,94112,94114,94115,...,
       94116,94117,94118,94121,94122,94123,94124,94127,94131,94132,94133,...,
       94134];
PowerData = [];
NumUsersPerZip = zeros(1,length(ZIP));
for i = 1:length(ZIP)
    load([num2str(ZIP(i)),'.mat']);
    NumUsersPerZip(1,i) = size(File1,1);
    PowerData = [PowerData;File1];
end
%%
NumUsers = size(PowerData,1);
NumMonths = 27;
DaysInMonth = [31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,29,31,30,31,30,31];
CumDays = [0,cumsum(DaysInMonth)];
PowCon = PowerData(:,1:sum(DaysInMonth)*96);
PowCon = [NaN*ones(NumUsers,96*3),PowCon];
yWeekAvg = zeros(NumUsers,96*7);
for user = 1:NumUsers
    yReshape = reshape(PowCon(user,:),96*7,118)';
    yWeekAvg(user,:) = mean(yReshape,'omitnan');
    for w = 1:118
        idx = find(isnan(yReshape(w,:)));
        yReshape(w,idx) = yWeekAvg(user,idx);
    end
    yReshape = yReshape';
    PowCon(user,:) = yReshape(:);
end
PowCon(:,1:96*3) = [];
MonthPeak = zeros(NumMonths,NumUsers);
MonthSecondPeak = zeros(NumMonths,NumUsers);
MonthPeakMet = zeros(NumMonths,NumUsers);
for month = 1:NumMonths
    MonthStartMet(month) = (CumDays(month)*96) + 1;
    MonthEndMet(month) = CumDays(month+1)*96;
    [MonthPeak(month,:),MonthPeakMet(month,:)] = max(PowCon(:,MonthStartMet(month):MonthEndMet(month))');
    x = sort(PowCon(:,MonthStartMet(month):MonthEndMet(month))','descend');
    MonthSecondPeak(month,:) = x(2,:);
end
MonthPeak = MonthPeak';
MonthSecondPeak = MonthSecondPeak';
MonthPeakMet = MonthPeakMet';

%% Identifying users with multiple peaks in any month or all NaN values in any month
mult_user = [];
a = cell(NumUsers,NumMonths);
for month = 1:NumMonths
    for user = 1:NumUsers
        a{user,month}= find(PowCon(user,MonthStartMet(month):MonthEndMet(month)) == MonthPeak(user,month));
        if size(a{user,month},2) > 10
            mult_user = [mult_user,user];
        end
    end
end
mult_user = unique(mult_user);
PowCon(mult_user,:) = [];
MonthPeak(mult_user,:) = [];
MonthPeakMet(mult_user,:) = [];
yWeekAvg(mult_user,:) = [];
NumUsers = size(PowCon,1);


