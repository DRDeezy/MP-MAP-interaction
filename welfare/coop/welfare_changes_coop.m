%Loads alternative welfare improving rules

load 'MAP_Chat_1_coop.mat'
load 'LTV_Chat_2_coop.mat'
load 'weights_Chat_3_coop.mat'
load 'law_borr_Chat_4_coop.mat'
load 'law_asset_prices_Chat_5_coop.mat'

%Loads initial values for the baseline model (basic rule)

load 'MP_only_Chat_0_coop.mat'


%Calculates changes in consumption

change_1_coop = (MAP_Chat_1_coop-MP_only_Chat_0_coop)./MP_only_Chat_0_coop;
change_2_coop = (LTV_Chat_2_coop-MP_only_Chat_0_coop)./MP_only_Chat_0_coop;
change_3_coop = (weights_Chat_3_coop-MP_only_Chat_0_coop)./MP_only_Chat_0_coop;
change_4_coop = (law_borr_Chat_4_coop-MP_only_Chat_0_coop)./MP_only_Chat_0_coop;
change_5_coop = (law_asset_prices_Chat_5_coop-MP_only_Chat_0_coop)./MP_only_Chat_0_coop;

n = 1;
lambd = 0.5;

%Calculates changes in welfare

change_EMU_1_coop = n*(lambd*change_1_coop(1) + (1-lambd)*change_1_coop(2));
change_EMU_2_coop = n*(lambd*change_2_coop(1) + (1-lambd)*change_2_coop(2));
change_EMU_3_coop = n*(lambd*change_3_coop(1) + (1-lambd)*change_3_coop(2));
change_EMU_4_coop = n*(lambd*change_4_coop(1) + (1-lambd)*change_4_coop(2));
change_EMU_5_coop = n*(lambd*change_5_coop(1) + (1-lambd)*change_5_coop(2));

%Stores pairs of consumption change / changes in welfare 
%From zero to 1
welfare_result_1_coop = [change_1_coop, change_EMU_1_coop];
welfare_result_2_coop = [change_2_coop, change_EMU_2_coop];
welfare_result_3_coop = [change_3_coop, change_EMU_3_coop];
welfare_result_4_coop = [change_4_coop, change_EMU_4_coop];
welfare_result_5_coop = [change_5_coop, change_EMU_5_coop];


%Creates and saves table
%1000 is a scaling factor
welfare_results_coop = 1000 * [welfare_result_1_coop;welfare_result_2_coop;welfare_result_3_coop;welfare_result_4_coop;welfare_result_5_coop];
rounded_results_coop = round(welfare_results_coop,5);
welfare_results_table_coop = array2table(rounded_results_coop, 'VariableNames',{'Savers','Borrowers','Weighted'}, ...
    'RowNames',{'MAP+','LTV+', 'Weights+', 'LAW_B+', 'LAW_P+'} );
%save('welfare_result','welfare_result');
save('welfare_results_table_coop', 'welfare_results_table_coop')
welfare_results_table_coop