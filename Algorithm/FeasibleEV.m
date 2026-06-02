function [Jf] = FeasibleEV(ucParam,ciParam,evParam,M,i_star,T_ser,T_tr)

    NumDrivers = evParam.NumDrivers;
    NumUsers = ciParam.NumUsers;
    T = ucParam.T;
    Jf = [];
    ServicePlusTransitInt = unique([T_ser,T_tr]);
    x = zeros(NumDrivers,T);
    y = zeros(NumDrivers,T);
    for j = 1:NumDrivers
        x(j,:) = ciParam.XDes'*reshape(M(:,j,:),NumUsers,T);
        y(j,:) = ciParam.YDes'*reshape(M(:,j,:),NumUsers,T);

        temp1 = sum(double(ismember(x(j,ServicePlusTransitInt), [0,ciParam.XDes(i_star)])));
        temp2 = sum(double(ismember(y(j,ServicePlusTransitInt), [0,ciParam.YDes(i_star)])));

        if temp1 == length(ServicePlusTransitInt) && temp2 == length(ServicePlusTransitInt)
            Jf = [Jf,j];
        end
    end

end