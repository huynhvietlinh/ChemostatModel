function Chemostat_Model
    plot_diff_ratio;
    plot_best_Hill_model    
    simulate_and_export;
    %search_parameter_Hill_function_model
end

function simulate_and_export()
    % Load experimental data
	[D1 D2 D3] = get_exp_data();
    % Herber model
	global mu_max C_max K_C V_max v_i v_o;
	assign_Herbert_model_parameter;
	period = 300;
	[T1 Y1] = ode45(@Herbert_chemostat_ode,[0, period], [0.1, V_max, 0.1, 0.15]);
	[T2 Y2] = ode45(@Herbert_chemostat_ode,[0, period], [0.2, V_max, 0.2, 0.15]);
	[T3 Y3] = ode45(@Herbert_chemostat_ode,[0, period], [0.4, V_max, 0.4, 0.35]);
	[T4 Y4] = ode45(@Herbert_chemostat_ode,[0, period], [0.8, V_max, 0.8, 0.15]);
	extract_and_print(T1, Y1, 4, 'Data/Herbert_glycerol_1.dat');
	extract_and_print(T2, Y2, 4, 'Data/Herbert_glycerol_2.dat');
	extract_and_print(T3, Y3, 4, 'Data/Herbert_glycerol_4.dat');
	extract_and_print(T4, Y4, 4, 'Data/Herbert_glycerol_8.dat');
	figure;
	p = plot(T1, Y1(:,4), '-g', T2, Y2(:,4), '-b', T3, Y3(:,4), '-r', T4, Y4(:,4), '-m');
	hold('on');
	scatter(D1(:,1), D1(:,2), 30, 'g', 'filled');
	hold('on');
	scatter(D1(:,1), D1(:,6), 30, 'b', 'filled');
	hold('on');
	scatter(D1(:,1), D1(:,10), 30, 'm', 'filled');
	hold('on');
	scatter(D2(:,1), D2(:,2), 30, 'r', 'filled');
	% Extended model
	[T1 Y1] = simulate_extended_model(0.1, 0.15, 2, 150);
	[T2 Y2] = simulate_extended_model(0.2, 0.15, 3, 155);
	[T3 Y3] = simulate_extended_model(0.4, 0.35, 4, 160);
	[T4 Y4] = simulate_extended_model(0.8, 0.15, 5, 155);
	extract_and_print(T1,Y1, 5, 'Data/Extended_glycerol_1.dat');
	extract_and_print(T2,Y2, 5, 'Data/Extended_glycerol_2.dat');
	extract_and_print(T3,Y3, 5, 'Data/Extended_glycerol_4.dat');
	extract_and_print(T4,Y4, 5, 'Data/Extended_glycerol_8.dat');
	figure;
	p = plot(T1, Y1(:,5), '-g', T2, Y2(:,5), '-b', T3, Y3(:,5), '-r', T4, Y4(:,5), '-m');
    hold('on');
	scatter(D1(:,1), D1(:,2), 30, 'g', 'filled');
	hold('on');
	scatter(D1(:,1), D1(:,6), 30, 'b', 'filled');
	hold('on');
	scatter(D1(:,1), D1(:,10), 30, 'm', 'filled');
	hold('on');
	scatter(D2(:,1), D2(:,2), 30, 'r', 'filled');
end

function plot_best_Hill_model
    global alpha_H K_H n_H;
    alpha_H = 20;
    K_H = 20;
    n_H = 1.5;
    plot_Hill_function_trial(true);
end

function search_parameter_Hill_function_model
    global alpha_H K_H n_H;
    alpha_H_list = [15,16,17,18,19,20,21,22,23,24,25];
    K_H_list = [15,16,17,18,19,20,21,22,23,24,25];
    n_H_list = [1.5,2];
    disp([2.1, 3.15, 4.4, 5.5]);
    alpha_H_best = -1;
    K_H_best = -1;
    n_H_best = -1;
    min_diff = 1e6;
    for i = 1:length(alpha_H_list)
        disp(i);
        for j = 1:length(K_H_list)
            for k = 1:length(n_H_list)
                 alpha_H = alpha_H_list(i);
                 K_H = K_H_list(j);
                 n_H = n_H_list(k);                
                 OD = plot_Hill_function_trial(false);
                 diff = abs(OD(1) - 2.1) + abs(OD(2) - 3.15) + abs(OD(3) - 4.4) + abs(OD(4) - 5.5);
                 if (diff < min_diff)
                     alpha_H_best = alpha_H;
                     K_H_best = K_H;
                     n_H_best = n_H;
                    min_diff = diff;
                 end
            end
        end        
    end
    disp([alpha_H_best, K_H_best, n_H_best]);
    alpha_H = alpha_H_best;
    K_H = K_H_best;
    n_H = n_H_best;                
    plot_Hill_function_trial(true);
end

function OD = plot_Hill_function_trial(is_plotted)
    % Extended model but with Hill function for the glycerol consumption
    global mu_max C_max K_C V_max v_i v_o;
	assign_extended_model_parameter;
    wt_phase_time = 80;
	[T1 Y1] = ode45(@extended_model_ode_with_Hill_function,[0, wt_phase_time], [0.1, V_max, 0.1, 0.15]);
    [T2 Y2] = ode45(@extended_model_ode_with_Hill_function,[0, wt_phase_time], [0.2, V_max, 0.2, 0.15]);
    [T3 Y3] = ode45(@extended_model_ode_with_Hill_function,[0, wt_phase_time], [0.4, V_max, 0.4, 0.35]);
    [T4 Y4] = ode45(@extended_model_ode_with_Hill_function,[0, wt_phase_time], [0.8, V_max, 0.8, 0.15]);
    OD = [Y1(length(Y1),4), Y2(length(Y2),4), Y3(length(Y3),4), Y4(length(Y4),4)];    
    %disp(OD);
    if (is_plotted)
        figure;
        p = plot(T1, Y1(:,4), '-g', T2, Y2(:,4), '-b', T3, Y3(:,4), '-r', T4, Y4(:,4), '-m');
        [D1 D2 D3] = get_exp_data();
        hold('on');
        scatter(D1(:,1), D1(:,2), 30, 'g', 'filled');
        hold('on');
        scatter(D1(:,1), D1(:,6), 30, 'b', 'filled');
        hold('on');
        scatter(D1(:,1), D1(:,10), 30, 'm', 'filled');
        hold('on');
        scatter(D2(:,1), D2(:,2), 30, 'r', 'filled');
        %
        extract_and_print(T1, Y1, 4, 'Data/Hill_model_glycerol_1.dat');
        extract_and_print(T2, Y2, 4, 'Data/Hill_model_glycerol_2.dat');
        extract_and_print(T3, Y3, 4, 'Data/Hill_model_glycerol_4.dat');
        extract_and_print(T4, Y4, 4, 'Data/Hill_model_glycerol_8.dat');
    end
end

function plot_diff_ratio
	% Herbert model with different ratio
	global mu_max C_max K_C V_max v_i v_o K_ratio;
	assign_Herbert_model_parameter;
	period = 300;
	K_ratio = 1;
	[T1 Y1] = ode45(@Herbert_chemostat_ode,[0, period], [0.1, V_max, 0.1, 0.15]);
	K_ratio = 2;
	[T2 Y2] = ode45(@Herbert_chemostat_ode,[0, period], [0.1, V_max, 0.1, 0.15]);
	K_ratio = 0.5;
	[T3 Y3] = ode45(@Herbert_chemostat_ode,[0, period], [0.1, V_max, 0.1, 0.15]);
	extract_and_print(T1,Y1, 4, 'Data/Herbert_K_C_1.dat');
	extract_and_print(T2,Y2, 4, 'Data/Herbert_K_C_2.dat');
	extract_and_print(T3,Y3, 4, 'Data/Herbert_K_C_0_5.dat');
	figure;
	%p = plot(T1, log(Y1(:,4)), '-g', T2, log(Y2(:,4)), '-b', T3, log(Y3(:,4)), '-r');
	%hold on;
	p = plot(T1, Y1(:,4), '-g', T2, Y2(:,4), '-b', T3, Y3(:,4), '-r');
	legend('K_C = K_M', 'K_C = 2K_M', 'K_C = 0.5K_M');
	xlabel('Time [h]');
	ylabel('OD');
	%plot2svg('Diff_K_ratio.svg');
end

function [T, Y] = simulate_extended_model(init_conc, init_OD, mutated_strain, wt_phase_time)
	global mu_max C_max K_C V_max v_i v_o;
	assign_extended_model_parameter;
	% WT: strain = 1
	%wt_phase_time = 100;
	[T1 Y1] = ode45(@extended_model_ode,[0, wt_phase_time], [1, init_conc, V_max, init_conc, init_OD]);
	% Mutant: strain = 2
	[T2 Y2] = ode45(@extended_model_ode,[0, 300 - wt_phase_time], [mutated_strain, Y1(size(Y1,1),2), Y1(size(Y1,1),3), Y1(size(Y1,1),4), Y1(size(Y1,1),5)]);
	T = cat(1,T1,T2 + wt_phase_time);
	Y = cat(1,Y1,Y2);	
end

function dydt = extended_model_ode(t,y)
	global V_max v_i v_o mu_max C_max K_C;
	% Parameters
	% V_max	
	% v_o
	% v_i
	% mu_max	maximum growth rate
	% v_C_max
	% K_C
	
	% -------------------------------------------
	% Variables
	Strain = y(1);	%
	C_0 = y(2);	% In-flux glycerol concentration	
	V = y(3);	% Culture volume
	C = y(4);	% Extracellular carbon source
	B = y(5);	% Biomass
	
	dydt = zeros(5,1);
	% V - Volume
	if (V > V_max)
		rho_V = v_o;
	else
		rho_V = 0;
	end
	dydt(3) = v_i - rho_V;
	mu = mu_max;
	switch Strain
		case 1
			Stress = 1/(1 + exp(0.09*(115 - t)));
			mu = 0.9*mu_max*(C/(C + 2*K_C))*(1 - Stress); % Fit
			dydt(4) = (v_i/V)*(C_0 - C) - (C/(C + K_C))*(B/V)*0.8*C_max*(1 + 0.055*B^2);
		case 2
			mu = 0.32*mu_max*(C/(C + K_C));
			dydt(4) = (v_i/V)*(C_0 - C) - (C/(C + K_C))*(B/V)*0.45*C_max*(1 + 0.055*B^2);
		case 3
			mu = 0.35*mu_max*(C/(C + K_C));
			dydt(4) = (v_i/V)*(C_0 - C) - (C/(C + K_C))*(B/V)*0.45*C_max*(1 + 0.055*B^2);
		case 4
			mu = 0.315*mu_max*(C/(C + K_C));
			dydt(4) = (v_i/V)*(C_0 - C) - (C/(C + K_C))*(B/V)*0.55*C_max*(1 + 0.055*B^2);
		case 5
			mu = 0.31*mu_max*(C/(C + K_C));
			dydt(4) = (v_i/V)*(C_0 - C) - (C/(C + K_C))*(B/V)*0.5*C_max*(1 + 0.055*B^2);
	end
	% Biomass
	dydt(5)  = mu*B - rho_V*(B/V);
end
function dydt = extended_model_ode_with_Hill_function(t,y)
	global V_max v_i v_o mu_max C_max K_C;
    global alpha_H K_H n_H;
	% Parameters
	% V_max	
	% v_o
	% v_i
	% mu_max	maximum growth rate
	% v_C_max
	% K_C
	
	% -------------------------------------------
	% Variables
	C_0 = y(1);	% In-flux glycerol concentration	
	V = y(2);	% Culture volume	
	C = y(3);	% Effective carbon source
	B = y(4);	% Biomass
	
	dydt = zeros(4,1);
	% V - Volume
	if (V > V_max)
		rho_V = v_o;
	else
		rho_V = 0;
	end
	dydt(2) = v_i - rho_V;	
    Stress = 1/(1 + exp(0.09*(115 - t)));
	mu = 0.9*mu_max*(C/(C + 2*K_C))*(1 - Stress); % Fit		
	%dydt(3) = C_0*v_i/V - C_eff*v_i/V - 0.8*(1 + 0.055*B^2)*C_max*(C_eff/(C_eff + K_C))*(B/V);
    %dydt(3) = (v_i/V)*(C_0 - C) - C_max*(C/(C + K_C))*(B/V)*0.8*(1 + 0.055*B^2);
    dydt(3) = (v_i/V)*(C_0 - C) - C_max*(C/(C + K_C))*(B/V)*alpha_H*(1/(1 + (K_H/B)^n_H));
	% Biomass
	dydt(4)  = mu*B - rho_V*(B/V);
end


function dydt = Herbert_chemostat_ode(t,y)
	global mu_max C_max K_C V_max v_i v_o K_ratio;
	% -------------------------------------------
	% Variables
	C_0 = y(1);	% In-flux glycerol concentration	
	V = y(2);	% Culture volume
	C = y(3);	% Extracellular glycerol
	B = y(4);	% Biomass
	
	dydt = zeros(4,1);
	% V - Volume
	if (V > V_max)
		rho_V = v_o;
	else
		rho_V = 0;
	end
	dydt(2) = v_i - rho_V;	
	mu = mu_max*C/(C + K_C);
	% Glycerol
	dydt(3) = C_0*v_i/V - C*v_i/V - C_max*(C/(C + K_ratio*K_C))*(B/V);
	% Biomass
	dydt(4)  = mu*B - rho_V*(B/V);
end

function assign_extended_model_parameter
	global mu_max C_max K_C V_max v_i v_o;
	mu_max = 0.41;
	%C_max = 0.0096; original
	C_max = 1.3*0.0096;
	K_C = 0.000975;
	V_max = 1;
	v_i = 0.00167*60;
	v_o = 0.0025*60;
end

function assign_Herbert_model_parameter
	global mu_max C_max K_C V_max v_i v_o K_ratio;
	% Palsson Nat Comm 0.23 - 0.45 h-1
	% Palsson Nature 2002 0.2 - 0.5 h-1
	% Herbert 1956: mu_max = 0.85 h-1, K_C = 12.3 ug/ml = 0.0123 g/l
	mu_max = 0.2;
	% Palsson Nat Comm 0.075 g l-1 hr-1 OD -1 = 0.006%, 0.12 = 0.0096
	C_max = 0.018;
	% K_C = 12.3 ug/ml = 0.0123 g/l = 0.000975%
	K_C = 0.001;
	% Chemostat parameters
	V_max = 1;
	v_i = 0.00167*60;
	v_o = 0.0025*60;
	K_ratio = 1;
end

function extract_and_print(T, Y, index, filename)
	fileID = fopen(filename, 'w');
	current_time = 0;
	resolution = 0.5;
	for t = 1:length(T)
		t_opt = -1;
		if (t > 1)
			t1 = T(t - 1);
			t2 = T(t);
			if (t1 <= current_time && t2 >= current_time)
				if (abs(t1 - current_time) < abs(t2 - current_time))
					t_opt = t - 1;
				else 
					t_opt = t;
				end
			end
		else
			t_opt = t;
		end
		if (t_opt > 0)
			fprintf(fileID, '%f\t%f\n', T(t_opt), Y(t_opt,index));
			current_time = current_time + resolution;
		end
	end
	fclose(fileID);
end

function y = piecewise_estimation(x, input_list, output_list)
	if (x <= input_list(1))
		y = output_list(1);
	end
	if (x > input_list(length(input_list)))
		y = output_list(length(output_list));
	end
	for i = 2:length(input_list)
		if (x > input_list(i - 1) && x <= input_list(i))
			y = output_list(i);
		end
	end
end

function [Data_1 Data_2 Data_3] = get_exp_data
	Data_1 = [0	0.154266667	0.001451819	0.1422	0.003187998	0.155566667	0.006013965	0.146466667	0.001121507	0.151233333	0.000996	0.152233333	0.005480977;
	4	0.451966667	0.005375975	0.4349	0.006301058	0.545433333	0.015984611	0.476966667	0.001560271	0.528466667	0.004664166	0.4714	0.00412351;
	6	1.0341	0.019679177	0.962733333	0.015002926	1.119666667	0.017618299	1.095133333	0.004861184	1.110333333	0.017105685	1.074466667	0.002182761;
	8	1.657833333	0.007055337	1.757	0.017820868	1.739666667	0.016340985	1.630333333	0.006495725	1.679	0.026581635	1.555166667	0.020415136;
	12	2.584166667	0.0120289	2.936666667	0.020083852	2.927333333	0.005494947	2.790666667	0.009426617	2.621	0.009013878	2.491166667	0.022520978;
	14	3.1135	0.01682508	3.278666667	0.019929738	3.457166667	0.01713022	3.298833333	0.002185813	3.182166667	0.020413095	3.114666667	0.011417579;
	22.5	3.283333333	0.011215069	3.426666667	0.019641226	4.818333333	0.024360715	5.061	0.030022214	7.006666667	0.036130935	6.664333333	0.013531855;
	24	3.556333333	0.034372146	3.325333333	0.027063711	4.7	0.009073772	4.904666667	0.04062156	7.826666667	0.04428067	7.387666667	0.074310908;
	26	3.217333333	0.017246578	3.338666667	0.039725447	4.539333333	0.031178696	5.152666667	0.052663502	8.053	0.044230457	8.277	0.053116225;
	28	2.646333333	0.044468466	2.953333333	0.025575596	4.127333333	0.007512952	4.588	0.004932883	8.501	0.03508561	8.695333333	0.013220355;
	30	2.15	0.027465129	2.809666667	0.005206833	4.080333333	0.023411773	4.539333333	0.018853234	7.985333333	0.083275313	7.908	0.013650397;
	32	2.161666667	0.028880982	2.576333333	0.034700304	3.890333333	0.012454361	4.136	0.02251666	6.787666667	0.037176755	7.169	0.044377171;
	36	2.015333333	0.036933875	2.086666667	0.02136456	3.493666667	0.015677301	3.631666667	0.051656988	6.370333333	0.064317787	6.140666667	0.067720832;
	47	2.195333333	0.060938038	2.013666667	0.019029217	3.294333333	0.009614803	3.567	0.01040833	5.82	0.091689694	5.909666667	0.020512869;
	48	2.182	0.040153871	2.103	0.023437861	3.277	0.031262331	3.439	0.019857828	5.818	0.061220911	5.812666667	0.035516819;
	52	1.975	0.040079088	2.133333333	0.023168465	3.226	0.029484459	3.536	0.011789826	5.892666667	0.027290617	5.749333333	0.029790006;
	61	2.073333333	0.045498474	2.099	0.026274195	3.214333333	0.071405727	3.287	0.029143324	5.886666667	0.058913307	5.834333333	0.033378303;
	72	2.087333333	0.016271994	2.168	0.020305993	3.147333333	0.043731504	3.197666667	0.037315472	5.501333333	0.112512715	5.514333333	0.057027284;
	77	2.09	0.060434538	2.173333333	0.042981908	3.221333333	0.098014171	3.348666667	0.046333333	5.776666667	0.027551971	5.614333333	0.13269807;
	96	1.970666667	0.074117774	2.163666667	0.035309741	3.119	0.058386642	3.263333333	0.084422219	5.742	0.032516662	5.669333333	0.070314373;
	118	1.883	0.113799531	2.114333333	0.013860415	3.085	0.058386642	2.993666667	0.082335358	4.880666667	0.078859228	5.364666667	0.035068188;
	129	1.703	0.089968513	1.772666667	0.120666667	2.451666667	0.034891897	2.763	0.133166562	3.899666667	0.125153683	4.079	0.048013887;
	144	1.085	0.074848736	1.076666667	0.076116869	1.447	0.058197938	1.616133333	0.011830093	1.949333333	0.038254992	1.843666667	0.139833631;
	149	0.662833333	0.059941648	0.409666667	0.041074458	0.805333333	0.040834354	0.891	0.025238859	1.259333333	0.06385748	1.118	0.087789521;
	166.5	0.466966667	0.034553356	0.453666667	0.012706079	0.544366667	0.036528087	0.5638	0.064212953	1.290233333	0.030788328	1.342333333	0.039252742;
	180	0.654666667	0.049626617	0.805333333	0.040834354	0.927666667	0.047837689	0.952066667	0.043537966	1.591666667	0.067907617	1.718666667	0.035713365;
	192	1.019333333	0.038181947	1.424166667	0.028852113	1.639	0.035157977	1.776666667	0.015624056	2.0585	0.035792224	2.1265	0.033752778;
	215	1.794	0.095379942	2.26	0.018583146	2.965666667	0.115404121	2.178333333	0.109748703	3.189333333	0.038110949	3.755666667	0.014836142;
	223	1.82	0.031895663	2.211666667	0.027702788	3.119666667	0.112595638	2.651	0.052048055	4.248333333	0.010974718	4.269333333	0.040883303;
	238	1.847666667	0.072778965	2.281333333	0.022666667	3.329666667	0.044945646	2.747333333	0.050406128	5.821	0.02953529	4.795	0.02116601;
	245	1.694	0.035019042	2.261333333	0.03006844	3.405	0.024131584	2.831666667	0.094014774	5.748666667	0.00788106	4.788	0.127735404;
	271	1.863333333	0.065077731	2.160333333	0.048015044	3.378333333	0.030333333	2.837	0.035557465	5.865333333	0.063572356	4.980333333	0.073984983];
	%Data_2 	= [0	0.35657	0.00622;
	%	24	6.19733	0.08145;
	%	48	5.29933	0.02157;
	%	72	5.27767	0.03182;
	%	81	5.24667	0.06;
	%	103	3.30467	0.08113;
	%	125	0.52557	0.00614;
	%	144	0.4204	0.00157;
	%	152	0.4619	0.00663;
	%	167	1.799	0.0448;
	%	175.5	2.79033	0.06034;
	%	192	3.258	0.03568;
	%	215	4.19033	0.05662;
	%	236	5.02833	0.03912;
	%	248	4.95	0.10318;
	%	266	4.973	0.05478;
	%	288	5.026	0.05511;
	%	296	4.91833	0.03001;
	%	312	5.03133	0.02871;
	%	336	4.99367	0.0189;
	%	344	5.01733	0.06885
	%];
	%Data_2 = [0	0.3529	0.00315
	%	24	5.28167	0.05183;
	%	42	4.52533	0.02173;
	%	48	4.38933	0.03053;
	%	66	4.35033	0.01206;
	%	72	4.35567	0.00306;
	%	96	3.83433	0.06058;
	%	114	2.99467	0.0125;
	%	122	2.465	0.02805;
	%	138	1.08283	0.01109;
	%	144	0.7708	0.00541;
	%	170	0.60177	0.00332;
	%	186	0.79293	0.00838;
	%	192	0.89027	0.0031;
	%	210	1.73933	0.01366;
	%	216.25	1.93767	0.01332;
	%	234.5	3.11367	0.04692;
	%	242	3.29133	0.02854;
	%	257.5	3.288	0.03251;
	%	264	3.23267	0.05237;
	%	286	3.963	0.026;
	%	306	3.95167	0.00219;
	%	330	3.93367	0.02133;
	%	354	3.95267	0.01519
	%];
	Data_2 = [0	0.2975	0.00377;
		24	5.54633	0.02747;
		48	4.452	0.01552;
		55	4.42367	0.02237;
		78	4.41467	0.034;
		102	3.46167	0.04119;
		119.5	2.59067	0.01955;
		128	1.87	0.0361;
		144	0.62293	0.00451;
		156	0.41247	0.00125;
		168	0.43707	0.00393;
		176	0.46507	0.00885;
		191	0.94067	0.00643;
		200	1.12883	0.0058;
		215	1.66567	0.0165;
		221.75	1.861	0.01997;
		226	2.013	0.02563;
		241	3.39467	0.04631;
		262.5	3.905	0.05624;
		272	4.07133	0.01387;
		287	4.10333	0.06757;
	%	311	4.06167	0.03301;
	%	335	4.087	0.04015;
	%	360	4.00867	0.09292;
	%	383	4.08567	0.06929
	];
	Data_3 = [0	0.3529	0.00182;
		24	5.28167	0.02992;
		42	4.52533	0.01255;
		48	4.38933	0.01763;
		66	4.35033	0.00696;
		72	4.35567	0.00176;
		96	3.83433	0.03498;
		114	2.99467	0.00722;
		122	2.465	0.0162;
		138	1.08283	0.00641;
		144	0.7708	0.00312;
		170	0.60177	0.00192;
		186	0.79293	0.00484;
		192	0.89027	0.00179;
		210	1.73933	0.00789;
		216.25	1.93767	0.00769;
		234.5	3.247	0.11086;
		242	4.02133	0.01648;
		257.5	4.018	0.01877;
		264	3.96267	0.03023;
		286	3.963	0.026;
		306	3.95167	0.00219;
		330	3.93367	0.02133;
		354	3.95267	0.01519
	];
end

