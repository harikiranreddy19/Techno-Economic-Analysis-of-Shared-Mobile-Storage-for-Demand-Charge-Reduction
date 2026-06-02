function [e_dis,e_tr,x,y] = DischargeTransitEnergy(ucParam,ciParam,evParam,M,u)

    T = ucParam.T;
    NumUsers = ciParam.NumUsers;
    NumDrivers = evParam.NumDrivers;

    e_dis = zeros(NumDrivers,T);
    e_tr = zeros(NumDrivers,T);
    for j = 1:NumDrivers
        temp = reshape(M(:,j,:),NumUsers,T);
        e_dis(j,:) = sum(ciParam.UserCharLim.*u.* temp.* [zeros(NumUsers,1), temp(:,1:T-1)],1);
        x(j,:) = ciParam.XDes'*temp;
        y(j,:) = ciParam.YDes'*temp;
    end
    e_tr(:,2:T) = evParam.EnerPerMile* (abs(x(:,2:T)-x(:,1:T-1)) + abs(y(:,2:T)-y(:,1:T-1)));


end