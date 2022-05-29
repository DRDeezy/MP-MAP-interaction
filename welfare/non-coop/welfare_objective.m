function [fval]=welfare_objective(xopt)

global M_ options_ oo_


switch options_.policy_regime


    case 0 % Monetary policy only
        set_param_value('rho_ib',xopt(1));
        set_param_value('phi_pie',xopt(2));
        %set_param_value('phi_y',xopt(3));
        %for ii=1:size(x_opt_name,1)
        %   set_param_value(x_opt_name{ii,1},x_opt(ii));
        %end

        var_list_={'W_EMU'};
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=10e6+sum([xopt(1),xopt(2)].^2); %return with penalty
        else
            fval=-oo_.mean(strmatch('W_EMU',var_list_,'exact')); %extract Welfare gap measure;
        end


    case 1 % Active MAP
        set_param_value('rho_vi',xopt(1));
        set_param_value('phi_vi',xopt(2));

        var_list_={'W_EMU'};

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=10e6+sum([xopt(1),xopt(2)].^2); %return with penalty
        else
            fval=-oo_.mean(strmatch('W_EMU',var_list_,'exact')); %extract Welfare gap measure;
        end

    case 2 % LTV cap
        set_param_value('rho_m_e',xopt(1));

        var_list_={'W_EMU'};

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=10e6+sum([xopt(1)].^2); %return with penalty
        else
            fval=-oo_.mean(strmatch('W_EMU',var_list_,'exact')); %extract Welfare gap measure;
        end

    case 3 % Weights

        set_param_value('rho_w_bu',xopt(1));
        set_param_value('chi_w_bu',xopt(2));
        set_param_value('rho_w_b',xopt(3));
        set_param_value('chi_w_b',xopt(4));

        var_list_={'W_EMU'};

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=10e6+sum([xopt(1),xopt(2),xopt(3),xopt(4)].^2); %return with penalty
        else
            fval=-oo_.mean(strmatch('W_EMU',var_list_,'exact')); %extract Welfare gap measure;
        end

    case 4 %LAW: augmented with borrowing

        set_param_value('rho_ib',xopt(1));
        set_param_value('phi_pie',xopt(2));
        set_param_value('phi_B',xopt(3));
        %set_param_value('phi_y',xopt(4));

        var_list_={'W_EMU'};

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=10e6+sum([xopt(1),xopt(2),xopt(3)].^2); %return with penalty
        else
            fval=-oo_.mean(strmatch('W_EMU',var_list_,'exact')); %extract Welfare gap measure;
        end


    case 5 %LAW: augmented with borrowing

        set_param_value('rho_ib',xopt(1));
        set_param_value('phi_pie',xopt(2));
        set_param_value('phi_AP',xopt(3));


        var_list_={'W_EMU'};

        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
        if info(1) %solution was not successful
            fval=10e6+sum([xopt(1),xopt(2),xopt(3)].^2); %return with penalty
        else
            fval=-oo_.mean(strmatch('W_EMU',var_list_,'exact')); %extract Welfare gap measure;
        end


end