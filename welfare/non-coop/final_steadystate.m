function [ys,params,check] = final_steadystate(ys,exo,M_,options_)
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] 0 if steady state computation worked and to
%                        1 of not (allows to impose restrictions on parameters)


%18.05 - added banks profits to households to get rid of -b_eeu*deltae term

%% Step 0: initialize indicator and set options for numerical solver
check = 0;
%Shows output if optimset('Display','Final','TolX',1e-10,'TolFun',1e-10);
options = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
flagScreenOutput = 0;
%% Step 1: read out parameters to access them with their name
for ii = 1:M_.param_nbr
  eval([ M_.param_names{ii} ' = M_.params(' int2str(ii) ');']);
end

%% Step 2: Enter model equations here
%New block with monopolistic banks

%By now this block stays the same.
eps_K_b = 1;
pie = 1;
A_e  = 1;
%r_ib = (1/beta_p) - 1;
q_k  = 1;
mc_E = 1 + (1 - mk_y_ss)/mk_y_ss;
x    = 1/mc_E;
%R_b  = r_ib;
%r_b  = r_ib + mcspread*1.5;
%r_b  = r_ib + mcspread;
vi = vi_ss;
deltae = 0.5;
%r_bu = (r_ib + mcspread)*1.25;
%spread = 1.5*mcspread;
%spread = mcspread;
%r_bu = (1/beta_e) - 1;
m_e = m_e_ss;
%----------------------------------------------------------------------------------------------------%
eps_d        = -1.46025;                          
%eps_b_ee       = 2.932806;
eps_b_ee       = 3.154542134;
eps_b_eeu = 2.117626240126008;

r_d = (1/beta_p - 1);
r_ib = (1/beta_p - 1) * (1-eps_d)/(-eps_d);
r_bu = (1/beta_e) - 1;
r_b = r_ib * eps_b_ee/(eps_b_ee-1); 
R_b  = r_ib;
R_bu  = r_ib;
R_d = r_ib;
weight_bu = 1;
weight_b  = 1;

%Equation 1: stays the same
%Solves for rental rate of capital
r_k_0 = 0.5;
[r_k,~,exitflag] = fsolve(@(r_k) ...
    (1 - beta_e*(1+r_b)) * (m_e_ss*q_k*(1-deltak))/(1+r_b) +  beta_e * (q_k*(1-deltak) + r_k) - q_k...
    , r_k_0,options);
if exitflag <= 0
    check = 1; % set failure indicator
    return     % return without updating steady states
end


%b_multiplier = (vi*delta_b - r_b + r_ib - r_ib*vi)/(r_bu*(1-deltae) - vi*delta_b - r_ib + r_ib * vi);
b_multiplier = (vi*delta_b - r_b + r_ib - r_ib*vi)/(r_bu*(1-deltae) - vi*delta_b*(1-deltae) - r_ib*(1-deltae) + r_ib * vi*(1-deltae));

%Equation 2
%l stands for labor demand, i.e. l_pd; k stands for k_e
k_l = ( (r_k * x) / (ksi * A_e) )^(1/(ksi-1));

%Equation 3
w_p = (k_l)^ksi * ((1-ksi)*A_e)/(x);

%Equation 4
%y stands for y_e; l stands for l_pd
y_l = A_e * (k_l)^ksi;

%Equation 5
%b_ee_ stands for b_ee; l stands for l_pd
b_ee_l = ( (m_e_ss*q_k*(1-deltak)) / (1+r_b) ) * k_l;

%Equation 11
%b_eeu_l stands for b_eeu; l stands for l_pd; solves for b_ee
b_eeu_l_0 = 0.5; %0.5
[b_eeu_l,~,exitflag] = fsolve(@(b_eeu_l) ...
    vi * delta_b * (b_ee_l + b_eeu_l ) ...
    - (r_b*b_ee_l + r_bu*b_eeu_l*(1-deltae) - r_d * (b_ee_l + b_eeu_l - vi*(b_ee_l + b_eeu_l) ) ) ...
    , b_eeu_l_0,options);
if exitflag <= 0
    check = 1; % set failure indicator
    return     % return without updating steady states
end

if b_eeu_l <= 0
    warning('Issue with b_eeu_l')
    check = 1; % set failure indicator
    return     % return without updating steady states
end


%%Equation 12
%b_eeu_l stands for b_eeu; l stands for l_pd; solves for d_p_l
d_p_l_0 = 0.5;
[d_p_l,~,exitflag] = fsolve(@(d_p_l) ...
    vi * delta_b * (b_ee_l + b_eeu_l) - r_b * b_ee_l - r_bu * b_eeu_l*(1-deltae) + r_d * d_p_l...
    , d_p_l_0,options);
if exitflag <= 0
    warning('Issue with d_p_l')
    check = 1; % set failure indicator
    return     % return without updating steady states
end

%%Equation 8
%%d_p stands for d_p; l stands for l_pd
%d_p_l = b_ee_l * (delta_b - r_b)/(delta_b - r_ib);  

%Equation 9.1
%J_R stands for J_R; l stands for l_pd
J_R_l = y_l * (1 - (1/x));

%Equation (taken form bank's) section
J_B_l  = (r_b*(b_ee_l) + r_bu*(b_eeu_l)*(1-deltae) - r_d*d_p_l);

%Equation 9
%c_p stands for c_p; l stands for l_pd
c_p_l = w_p + r_d * d_p_l + J_R_l + J_B_l;


%Equation 10
%Solves for optimal labor supply
l_pd_0 = 1/3;
[l_pd,~,exitflag] = fsolve(@(l_pd) ...
    c_p_l - (w_p)/((l_pd)^(1 + phi))...
    , l_pd_0,options);
if exitflag <= 0
    check = 1; % set failure indicator
    return     % return without updating steady states
end

%Equation 2
k_e = k_l * l_pd;
%Equation 4
y_e = y_l * l_pd;
%Equation 5
b_ee = b_ee_l * l_pd;
%Equation 8
d_p = d_p_l * l_pd;
%Equation 9.1
J_R = J_R_l * l_pd;
%Equation 9
c_p = c_p_l * l_pd;
%Equation 11
b_eeu = b_eeu_l * l_pd; 
%Consumption of entrepreneurs
c_e = - (1-deltae) * (1 + r_bu) * (b_eeu) - (1+r_b) * (b_ee)  -  (w_p)*(l_pd)  - (q_k) * (k_e) + (y_e) / (x) + (b_ee) + ...
(b_eeu) + (q_k) * (1-deltak) * (k_e);

if c_e<=0
    warning('Issue with c_e ALARM')
    check =1; 
    return
end

%Lagrange multipliers
lam_p = (c_p)^(-1);
lam_e = (c_e)^(-1); 
s_e = (lam_e) - beta_e*(lam_e)*(1+r_b);


%Market clearing
B              = gamma_e*(b_ee) + gamma_e*(b_eeu); 
D              = gamma_p*(d_p) ; 
K_b = B - D;


%Aggregation
C              = gamma_p*(c_p) + gamma_e*(c_e); 
l_p = gamma_e * l_pd/gamma_p; 
K              = gamma_e*(k_e); 
%B              = (D) + (K_b);  
Y              = gamma_e*(y_e); 
%Y              = (C) + (q_k)*((K) - (1-deltak)*(K)) + delta_b*(K_b) - deltae*b_eeu; 
%Y              = (C) + (q_k)*((K) - (1-deltak)*(K)); %for case with banks
%profits to patient household

%Investment
I = deltak * k_e;
mk_y = mk_y_ss;

%Banks
Y1 = C+I;
lev = (B)/(K_b);
rr = r_b - (pie);
%J_B  = (r_b*(b_ee) + r_bu*(b_eeu)*(1-deltae) - r_ib*(D));
J_B = J_B_l * l_pd;

%Omega
omega = lam_e/((b_eeu) * (deltae) * (1 + r_bu));

%Welfare
V_S = (log((c_p)) - ((l_p)^(1+phi)/(1+phi)))/(1-beta_p);
V_B = (log(c_e) - 0.5*omega*((1+r_bu)*b_eeu*deltae)^2)/(1-beta_e);
W_S = 0;
W_B = 0;
W_EMU = 0;

%Risk-weighted assets
B_ee = gamma_e * b_ee;
B_eeu = gamma_e * b_eeu;
B_rw = weight_bu * B_eeu + weight_b * B_ee;

log_b_ee  = log(b_ee);
log_q_k = 0;
log_k_e = log(k_e);
one_plus_r_b = 1 + r_b;
log_one_plus_r_b = log(one_plus_r_b);
q_k_exp = 1;
log_q_k_exp = 0;
log_c_e = 0;

%--Log-deviations--%
log_r_ib = 0;

%% Step 4: Update parameters and variables
params=NaN(M_.param_nbr,1);
for iter = 1:M_.param_nbr %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

for ii = 1:M_.orig_endo_nbr %auxiliary variables are set automatically
  eval(['ys(' int2str(ii) ') = ' M_.endo_names{ii} ';']);
end