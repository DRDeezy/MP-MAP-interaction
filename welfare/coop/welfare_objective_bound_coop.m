function [fval]=welfare_objective_bound_coop(xopt)

global M_ options_ oo_


switch options_.policy_regime


    case 0

        if xopt(1)>=1 || xopt(1)<=0 || xopt(2)>=5 || xopt(2)<=1
            fval=150+sum(xopt.^2);
            return
        end

        set_param_value('rho_ib',xopt(1));
        set_param_value('phi_pie',xopt(2));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=150+sum([xopt(1),xopt(2)].^2); %return with penalty
        else
            fval=-oo_.mean;
        end

    case 1

        if xopt(1)>=1 || xopt(1)<=0 || xopt(2)>=5 || xopt(2)<=0 || xopt(3)>=1 || xopt(3)<=0 || xopt(4)>=5 || xopt(4)<=1
            fval=150+sum(xopt.^2);
            return
        end

        set_param_value('rho_vi',xopt(1));
        set_param_value('phi_vi',xopt(2));
        set_param_value('rho_ib',xopt(3));
        set_param_value('phi_pie',xopt(4));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=150+sum([xopt(1),xopt(2),xopt(3),xopt(4)].^2); %return with penalty
        else
            fval=-oo_.mean;
        end

    case 2

        if xopt(1)>=0 || xopt(1)<=-5 || xopt(2)>=1 || xopt(2)<=0 || xopt(3)>=5 || xopt(3)<=1
            fval=150+sum(xopt.^2);
            return
        end

        set_param_value('rho_m_e',xopt(1));
        set_param_value('rho_ib',xopt(2));
        set_param_value('phi_pie',xopt(3));


        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=150+sum([xopt(1),xopt(2),xopt(3)].^2); %return with penalty
        else
            fval=-oo_.mean;
        end

    case 3

        if xopt(1)>=1 || xopt(1)<=0 || xopt(2)>=0 || xopt(2)<=-20 || xopt(3)>=1 || xopt(3)<=0 || xopt(4)>=0 || xopt(4)<=-20 || xopt(5)>=1 || xopt(5)<=0 || xopt(6)>=5 || xopt(6)<=1
            fval=150+sum(xopt.^2);
            return
        end

        set_param_value('rho_w_bu',xopt(1));
        set_param_value('chi_w_bu',xopt(2));
        set_param_value('rho_w_b',xopt(3));
        set_param_value('chi_w_b',xopt(4));
        set_param_value('rho_ib',xopt(5));
        set_param_value('phi_pie',xopt(6));

        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=150+sum([xopt(1),xopt(2),xopt(3),xopt(4),xopt(4),xopt(5)].^2); %return with penalty
        else
            fval=-oo_.mean;
        end

    case 4

        if xopt(1)>=1 || xopt(1)<=0 || xopt(2)>=5 || xopt(2)<=1 || xopt(3)>=0 
            fval=150+sum(xopt.^2);
            return
        end

        set_param_value('rho_ib',xopt(1));
        set_param_value('phi_pie',xopt(2));
        set_param_value('phi_B',xopt(3));


        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=150+sum([xopt(1),xopt(2),xopt(3)].^2); %return with penalty
        else
            fval=-oo_.mean;
        end

    case 5

        if xopt(1)>=1 || xopt(1)<=0 || xopt(2)>=5 || xopt(2)<=1 || xopt(3)>=10 || xopt(3)<=0 
            fval=150+sum(xopt.^2);
            return
        end

        set_param_value('rho_ib',xopt(1));
        set_param_value('phi_pie',xopt(2));
        set_param_value('phi_AP',xopt(3));


        var_list_={'W_EMU'};
        %info=stoch_simul(var_list_);
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=150+sum([xopt(1),xopt(2),xopt(3)].^2); %return with penalty
        else
            fval=-oo_.mean;
        end


end