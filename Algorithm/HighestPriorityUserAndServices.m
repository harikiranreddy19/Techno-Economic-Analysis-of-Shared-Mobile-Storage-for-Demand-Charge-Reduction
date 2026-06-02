function [i_star,s_star,T_ser,T_tr] = HighestPriorityUserAndServices(ciParam,u,M,P,evParam,ucParam)

    P_max = max(max(P));
    [i_star_pot,s_star_pot] = find(P == P_max);
    s_star = max(s_star_pot);
    s_star_ind = find(s_star_pot == s_star);
    
    if length(s_star_ind) == 1
        i_star = i_star_pot(s_star_ind);
    else
        i_star_pot_can = i_star_pot(s_star_ind); 
        [~,i_star_ind_min_dis] = min(abs(ciParam.XDes(i_star_pot_can))+abs(ciParam.YDes(i_star_pot_can)));
        i_star = i_star_pot_can(i_star_ind_min_dis);
    end

    [~,T_ser] = NewDC(i_star,s_star,ciParam,ucParam,evParam,u,M);
    T_tr = T_ser-1;

end