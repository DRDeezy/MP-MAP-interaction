% rho = 0.55:0.01:0.99;
% phi_p = 4.5:0.1:10;
format longG
% 
% %Define loss-function
% 
% colNames = {'Welfare','rho_ib','phi_pie'};
% 
% opt_mat = zeros(length(rho)*length(phi_p),3);
% numit = 0;
% first_time = 1;
% 
% for j=1:length(phi_p)
%     for i=1:length(rho)
%         if first_time == 1 %Launches for the first time to store all the information
%             dynare final noclearall;
%             numit = numit + 1;
%             first_time = 0;
%             mean_s = mean(oo_.occbin.simul.piecewise(1:end,2));
%             opt_mat(numit,1:3) = [mean_s,rho(1),phi_p(1)];
%         else
%             numit = numit + 1;
%             set_param_value('rho_ib',rho(i));
%             set_param_value('phi_pie',phi_p(j));
%             [oo_, out]= occbin.solver(M_, oo_, options_);
%             if ~out.error_flag
%                 mean_s = mean(oo_.occbin.simul.piecewise(1:end,2));
%                 opt_mat(numit,1:3) = [mean_s,rho(i),phi_p(j)];
%             else
%                 mean_s = NaN;
%             end
%         end
%     end
% end
% 
% max_table = array2table(opt_mat,'VariableNames',colNames);
% [~,maxidx] = max(max_table.Welfare);
% optimal_rho = max_table(maxidx,:);
% 


dynare final;
loss = ((std(oo_.occbin.simul.piecewise(1:end,3)))^2 + (std(oo_.occbin.simul.piecewise(1:end,2)))^2 + (std(oo_.occbin.simul.piecewise(1:end,4)))^2 + ...
   (std(oo_.occbin.simul.piecewise(1:end,5)))^2)/5 ;
loss


