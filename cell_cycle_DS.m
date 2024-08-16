clc; clear; close all; format compact;
scale = 5; % to increase the size of the system

%% Random Number Generators
fix_RNG = 0;
if fix_RNG
    rng_init_state = rng(595); % fix random seed
    % Create 5 Random Number Streams, 1 for each slow reaction.
    GEN = RandStream.create('mrg32k3a','NumStreams',5,'Seed',79,'CellOutput',1)';
else % randomize seed
    rng_init_state = rng('shuffle');
    GEN = RandStream.create('mrg32k3a','NumStreams',5,'Seed',randi(2^32),'CellOutput',1)';
end
% to save the RNG states before entering a downscaling window
GEN_states = cell(size(GEN,1),1);

%% control variables
on_PC = 1; % run on PC (use Figures) or cluster (no Figures, write files)? 
DS = 1; % perform Downscale, 1=true, 0=false

plot_f_UN = true; % plot all future UN  -scaled trajectories
plot_f_DS = true; % plot all future Down-Scaled trajectories
plot_all_species = false; % plot all species = T, or the independent ones = F
plot_QE = false;  % plot quasi-equilibrated quantities

% Positions for figures (laptop screen)
% from left, from bottom, width, height
p2 = [575  50   560 420]; % bottom left
p4 = [10   50   560 420]; % bottom middle
p5 = [800  250  560 420]; % top    right - 2
p6 = [800  50   560 420]; % bottom right - 1
p7 = [810  55   560 420]; % bottom right - 2
p8 = [250  200  560 420];

%% user-defined parameters
t_max = 500;       % mins
sampling_dt= 0.01; % mins

% N: the number of KMC steps after which the DS algorithm is invoked.
steps_DS = 100000;

% w = tf – ts: the width, in KMC time units, of the window over which
% the algorithm performs the downscaling evaluations.
min_t_ahead = 10; % MIN time width of a downscaling window
max_t_ahead = 15; % MAX time width of a downscaling window

% dfb: the base downscale factor used in the generation of multiple downscaled trajectories.
down_factor = 5;

% the total number of downscaled trajectories generated inside the downscaling window 
multi_DS = 9; 

% the minimum timescale separation, in orders of magnitude 
TSSmin = 2.0;

% Parameters for objective function. Only some value pairs are used:
% Pair 1: alpha=1, beta =  0.5
% Pair 2: alpha=1, beta =  1.0
% Pair 3: alpha=1, beta =  2.0
% Pair 4: alpha=1, beta = 10.0
alpha_beta_pair = 3;

% emax: maximum absolute error threshold for any chosen downscale factor
error_threshold = 0.05;

% einc: maximum allowed absolute increase in orders of magnitude of the overall error norm
e_inc = 2;

% how many DS attempts to do
N_sc_tot=5;

%% reaction constants and reaction matrix
NA  = 6.0221367e14; % 1/nano-mol
Vb  = 17e-16;       % L. Base volume that corresponds to O=2, Yi=207
Vol = Vb*scale;     % L. Volume of the system
NV  = NA * Vol;     % L/nano-mol

zeta     = 100;     % controls the Time-Scale separation
k0       = 0.15;
k1       = 50;
k2       = 9.635e-3 / NV;
chi_init = 50 * zeta/ NV; chi = chi_init;  % 50
phi_init = 20 * zeta/ NV; phi = phi_init;  % 20
alfa     = 31.135   * NV;
beta     = 1.038    * NV;
lamda1   = 0.048    / NV;
lamda2   = 0.01;
% rates vector instead of individual variables
rates_init = [k0 k1 k2 chi phi alfa beta lamda1 lamda2]';

%         1  2  3  4  5  6  7  8  9 reactions, 6 species
r_mat = [ 0  0  0  0 -1  1  0  0  0; ...  % O
          1  1 -2  2  0  0  0 -1  0; ...  % X
          0  0  1 -1 -1  1  0  0  0; ...  % X2
          0  0  0  0  1 -1  0  0  0; ...  % OX2
          0  0  0  0  0  0 -1  0  1; ...  % Yi
          0  0  0  0  0  0  1  0 -1];     % Y

%% initial population
P_O   = 0;
P_X   = 0;
P_X2  = 0;
P_OX2 = 2   * scale; % 2
P_Yi  = 207 * scale; % 207
P_Y   = 0;
% population vector instead of individual variables
popul_init = [P_O P_X P_X2 P_OX2 P_Yi P_Y]';

%% initialize arrays and other necessary variables
N_sc = 0; N_sc_idx=1; 
Accepted=0; Rejected=0; Upscaled=0;
Decision = NaN(1,N_sc_tot);
t_DS_pos = NaN(1,N_sc_tot);
t_ahead_width = NaN(1,N_sc_tot); % width of each DS time window
accur_idx = 100;

steps = ceil(t_max/sampling_dt);
data = NaN(steps+1, 10); % t, dt, 6 species, chi, phi. 1+1+6+1+1
steps_vs_t = NaN(steps+1,3);
steps_vs_t(1,:) = 0; % 1st sample at t=0, steps=0
map = [1 2 0 0 0 0 3 4 5]; % used in saving the occurrence times of slow reactions
% What is the index of the i-th slow reaction ?
map_inv = [1  2  7  8  9];
cost_DS = NaN(multi_DS+1,N_sc_tot); % actual cost during DS window

d_factors = NaN(multi_DS,1);
df = ones(N_sc_tot,1);
IAT_dif = cell(5,multi_DS); % 1 row per slow reaction. 1 col per DS traj
CE_data = cell(N_sc_tot,3);% to save the cost and error "lines"

a  = zeros(9, 1);
P  = zeros(9, 1);
Dt = zeros(9, 1);

%% initialize/reset variables
t = 0;
T = zeros(9, 1);
sampling_t = sampling_dt;
samples = 1;

popul = popul_init;
rates = rates_init;
counters = zeros(9, 2);
j_total = 0; j_chunk = 0;
f_traj = 1;
% counters_f_?? keep track of extra work done by the algoritmh while running future trajectories.
% Apart from a SINGLE trajectory (among the UNscaled-DownScaled) that will be used to
% fill the "Down-scaling window" time interval, the other f_traj are "wasted" and are not actively
% used anywhere (apart from the t-tests that decide which SET={UN, DS} is going to be accepted)
% Double-counting KMC steps might occur, and is taken care of when a future trajectory is accepted
counters_f_UN = zeros(11, N_sc_tot); % 9 species, 1 blank, 1 sums
counters_f_DS = zeros(11, f_traj+multi_DS, N_sc_tot);

% 4 lines for population error,  1 blank, 4 for error norm
popul_error_Un = zeros(9,N_sc_tot);   % for pair  Un-DS

% Occurrence times of slow reactions 1, 2, 7, 8, 9
ocT    = NaN(25100100,5); % for overall error calculations
ocT_Un = NaN(15000,5);
ocT_DS = NaN(15000,5);

IAT_err_Un = NaN(5,multi_DS,N_sc_tot);
N_r = NaN(5, multi_DS);

data(1,1  ) = t;
data(1,2  ) = NaN; %dt
data(1,3:8) = popul;
data(1,9  ) = rates(4); % chi
data(1,10 ) = rates(5); % phi

% propensities
a(1) = popul(1)                   * rates(1);          % O      * k0
a(2) = popul(4)                   * rates(2);          % OX2    * k1
a(3) = popul(2) * abs(popul(2)-1) * rates(4);          % X(X-1) * chi
a(4) = popul(3)                   * rates(7)*rates(4); % X2     * beta*chi
a(5) = popul(1) * popul(3)        * rates(5);          % O * X2 * phi
a(6) = popul(4)                   * rates(6)*rates(5); % OX2    * alfa*phi
a(7) = popul(5) * popul(3)        * rates(3);          % Yi*X2  * k2
a(8) = popul(2) * popul(6)        * rates(8);          % X * Y  * lamda1
a(9) = popul(6)                   * rates(9);          % Y      * lamda2

P(1) = -log(rand(GEN{1}));
P(2) = -log(rand(GEN{2}));
P(3) = -log(rand());
P(4) = -log(rand());
P(5) = -log(rand());
P(6) = -log(rand());
P(7) = -log(rand(GEN{3}));
P(8) = -log(rand(GEN{4}));
P(9) = -log(rand(GEN{5}));

Dt = ( P - T ) ./ a; % P, T, a, Dt are all vectors
[D, m] = min(Dt);

%% Main KMC Loop
criterion_cur_val = j_chunk;
crit_checkpnt_val = steps_DS;

disp(datetime('now'));
tic; jj=0; % j for reporting, jj for Global step counter
N_UN = 0; N_DS = 0; % count attempts (for verification)

while t <= t_max
    %% Modified Next Reaction Method
    t = t + D;

    %% sample (before updating population)
    while (t > sampling_t && sampling_t <= t_max)
        samples = samples + 1;
        data(samples,1   ) = sampling_t;
        data(samples,2   ) = D;
        data(samples,3:8 ) = popul;
        data(samples,9:10)= rates(4:5);
        steps_vs_t(samples,1) = sampling_t;
        steps_vs_t(samples,2) = jj; % golbal counter, includes DS window
        steps_vs_t(samples,3) = j_total; % KMC steps outside DS window
        sampling_t = sampling_t + sampling_dt;
    end

    %% update population and counters
%   popul = popul + r_mat(:,m); % slower than explicit loop
    for i=1:6 % same performance as the switch-case construct
        popul(i) = popul(i) + r_mat(i,m);
    end
    counters(m) = counters(m) + 1;

    % update Modified_Next_Reaction_Method quantities
    [D,m,T,P,a,Dt] = Mod_NRM(D,Dt,T,P,a,m,popul,rates,map,GEN);
    j_total = j_total + 1; jj=jj+1;
    j_chunk = j_chunk + 1;

    if DS
        criterion_cur_val = j_chunk;
    end

    %% check for downscaling
    if ( DS && N_sc==(N_sc_idx-1) && criterion_cur_val >= crit_checkpnt_val && t>1)
        %% variables, counter and Figure set-up
        N_sc = N_sc + 1;
        N_sc_idx = min(N_sc_idx+1, N_sc_tot);
        % erase occurrence times from previous downscalings
        ocT_Un = NaN(size(ocT_Un));
        % Update DS criterion_current_value and checkpoint accordingly
        t_DS_pos(N_sc) = t;
        j_chunk = 0; % reset KMC chunk counter

        % decide on width of time window (t_ahead)
        if (N_sc==1)
            t_ahead_i = max(min_t_ahead, 1e6 * t / jj);
            t_ahead_i = min(max_t_ahead, t_ahead_i);
        else
            % Time Window in Previous Attempt
            TW_PA = t_ahead_width(N_sc-1);
            % what was the outcome of previous attempt?
            OC_PA = Decision(N_sc-1); % 1-2 DS, 3 Rej
            % desired steps on each trajectory
            deS = 1e6;
            % correction factor
            cf = 1;
            if (OC_PA==1 || OC_PA==2)
                % average steps in previous attempt
                AS_PA = counters_f_UN(end, N_UN)/f_traj;
                % correct estimation, rates have been downscaled
                cf = down_factor;
            elseif (OC_PA==3 )
                % average steps in previous attempt
                AS_PA = counters_f_UN(end, N_UN)/f_traj;
                % no correction needed.
            end
            sps = AS_PA / TW_PA; % steps per KMC second
            ETW = deS/sps * cf;  % Estimated Time Window

            % MAX of  {max_t_ahead, ETW (can be tiny)}
            t_ahead_i = max(min_t_ahead, ETW);
            % MIN of  {min_t_ahead, ETW (can be huge)}
            t_ahead_i = min(max_t_ahead, t_ahead_i);
        end
%       t_ahead_i = 25.1; % for testing
        t_ahead_width(N_sc) = t_ahead_i;
        % Number of points in f_traj, depends on time window width
        N_points = floor( t_ahead_width(N_sc) / sampling_dt);
        f_unscaled   = zeros(N_points,10); 
        f_DownScaled = zeros(N_points,10,multi_DS+1);

        if on_PC
            fig1 = figure('Position', p4);
            title(['Downscaling attempt: ', num2str(N_sc)]);
            xlabel('KMC Time (min)'); ylabel('Population');
            hold on; grid on; box on; pause(1e-9);
        end

        %% calculate UNscaled trajectories
        % save the state of the "slow" random number generators
        for i=1:size(GEN,1)
            GEN_states{i} = GEN{i}.State;
        end
        N_UN = N_UN + 1;

        sampling_t_f = sampling_t; samples_f = 0;
        popul_fun=popul; t_f=t; % future unscaled population
        a_f=a; P_f=P; T_f=T; Dt_f=Dt; [D_f, m_f]=min(Dt_f);
        while(t_f < t + t_ahead_width(N_sc))
            t_f = t_f + D_f;

            % the 2nd condition prevents sampling 1 more point beyond the intented interval
            while(t_f > sampling_t_f && samples_f < N_points)
                samples_f = samples_f + 1;
                f_unscaled(samples_f,  1) = sampling_t_f;
                f_unscaled(samples_f,  2) = D_f;
                f_unscaled(samples_f,3:8) = popul_fun;
                f_unscaled(samples_f,9:10)= rates(4:5);
                sampling_t_f = sampling_t_f + sampling_dt;
            end

            % update population & counters
            for i=1:6
                popul_fun(i) = popul_fun(i) + r_mat(i,m_f);
            end
            counters_f_UN(m_f, N_UN) = counters_f_UN(m_f, N_UN) + 1;

            % save occurrence time for slow reactions
            if ~(m_f>2 && m_f<7) % m = 1, 2, 7, 8, 9
                ocT_Un(counters_f_UN(m_f,N_UN), map(m_f)) = t_f;
            end

            % update Modified_Next_Reaction_Method quantities
            [D_f, m_f, T_f, P_f, a_f, Dt_f] = ...
                Mod_NRM(D_f, Dt_f, T_f, P_f, a_f, m_f, popul_fun,rates,map,GEN);

            jj=jj+1;
        end
        counters_f_UN(11, N_UN) = sum( counters_f_UN(1:9, N_UN) );
        % plot each trajectory
        if on_PC && plot_f_UN
            % active figure changes randomly from time to time so...
            figure(fig1); % retrieve figure I want to plot on
            tt  = f_unscaled(:,1); % time
            pO  = f_unscaled(:,3);
            pX  = f_unscaled(:,4);
            pX2 = f_unscaled(:,5);
            pOX2= f_unscaled(:,6); % dependent species
            pYi = f_unscaled(:,7); % dependent species
            pY  = f_unscaled(:,8);
            if plot_all_species
                p_un=plot(tt,pO,':k', tt,pX,'--k', tt,pX2,'-.k',...
                    tt,pOX2,'--k',tt,pYi,':k',tt,pY,'-.k');
            else
                p_un=plot(tt,pO,':k', tt,pX,'--k', tt,pX2,'-.k',...
                    tt,pY,'-.k');
            end
            pause(1e-9);
        end
        cost_DS(1,N_sc) = counters_f_UN(11, N_UN)/t_ahead_i;
        pause(1e-9);

        %% calculate the Inter-Arrival Times of the Unscaled trajectory
        R1_IAT_Un = diff(ocT_Un(:,1));
        R2_IAT_Un = diff(ocT_Un(:,2));
        R7_IAT_Un = diff(ocT_Un(:,3));
        R8_IAT_Un = diff(ocT_Un(:,4));
        R9_IAT_Un = diff(ocT_Un(:,5));

        %% calculate reaction frequencies (of the UNscaled)
        freq_UN(1:9,1) = 1:9; % reaction number
        freq_UN(1:9,2) = counters_f_UN(1:9,N_sc) / t_ahead_width(N_sc);
        freq_UN_sorted = sortrows(freq_UN,2,'descend');
        if on_PC
            figure();
            barh(flipud(freq_UN(:,2)));
            title(['Downscaling attempt: ', num2str(N_sc)]);
            xlabel('Frequency (min^{-1})');
            ylabel('Reaction Number');
            set(gca,'XScale','log'); grid on;
            lbls = flipud(yticklabels); % retrieve and flip labels
            yticklabels(lbls);          % set labels to the updated values
        end

        % Over the stiff region, the fast reactions are 3-4, 5-6.
        % The 1st slow reaction must appear 5th in sorted frequencies.
        % Over the non-stiff region, ORDERING MAY CHANGE
        % Check if reactions 3-4, 5-6 are indeed the fastest
        freq_ord = ismember(freq_UN_sorted(1:4,1), [3 4 5 6]);
        if nnz(freq_ord) < 4 % another reaction appears as fast
            % do NOT proceed to downscaling
            early_stop = 1;
        else
            % Compare the slowest fast (4th) against the fastest slow (5th)
            % Take log10(frequencies) to compare orders of magitude
            freq_dif = log10(freq_UN_sorted(4,2)) - ...
                       log10(freq_UN_sorted(5,2));
            if freq_dif < TSSmin
                early_stop = 1; % the fast are not fast enough
            else
                early_stop = 0;
            end
        end
        if on_PC
            % line for fastest slow reaction
            line(freq_UN_sorted(5,2)*[1 1],[0.5 9.5],'Color','red','LineWidth',2,'LineStyle','--');
            % line for the slowest fast
            line(freq_UN_sorted(4,2)*[1 1],[0.5 9.5],'Color','red','LineWidth',2,'LineStyle','--');
            % print a text with the difference of order of magn.
            text(freq_UN(8,2)*1.5, 9.5, ['TSS = ' num2str(freq_dif,4)])
            % print decision in text
            if early_stop
                text(freq_UN(8,2)*1.5, 9,"FAIL",'Color','red');
            else
                text(freq_UN(8,2)*1.5, 9,"PASS",'Color',[0 0.5 0]);
            end
            pause(1e-9);
        end

        if early_stop % accept the UNscaled trajectory
            Reject % calls external file (for brevity & clarity)
        else
            %% calculate Down-Scaled trajectories
            N_DS = N_DS + 1;
            f_traj = multi_DS; % how many DS trajectories to generate
            for f_tr=1:f_traj
                % erase occurrence times from previous attempts
                ocT_DS = NaN(size(ocT_DS));
                % the downscaled rates are dfferent in every 'iteration'
                d_factors(f_tr) = down_factor ^ (f_tr-1);
                chi_ds(f_tr) = rates(4)/down_factor ^ (f_tr-1);
                phi_ds(f_tr) = rates(5)/down_factor ^ (f_tr-1);
                rates_ds = [k0 k1 k2 chi_ds(f_tr) phi_ds(f_tr) alfa beta lamda1 lamda2]';
                % RESTORE the state of the "slow" RNGs
                for i=1:size(GEN,1)
                    GEN{i}.State = GEN_states{i};
                end
                sampling_t_f = sampling_t; samples_f = 0;
                popul_fds=popul; t_f=t; % future downscaled population
                a_f=a; P_f=P; T_f=T; Dt_f=Dt; [D_f, m_f]=min(Dt_f);
                while(t_f < t + t_ahead_width(N_sc))
                    t_f = t_f + D_f;

                    % the 2nd condition prevents sampling 1 more point beyond the intented interval
                    while(t_f > sampling_t_f && samples_f+1 <= N_points)
                        samples_f = samples_f + 1;
                        f_DownScaled(samples_f,  1, f_tr) = sampling_t_f;
                        f_DownScaled(samples_f,  2, f_tr) = D_f;
                        f_DownScaled(samples_f,3:8, f_tr) = popul_fds;
                        f_DownScaled(samples_f,9:10, f_tr)= rates_ds(4:5);
                        sampling_t_f = sampling_t_f + sampling_dt;
                    end

                    % update population & counters
                    for i=1:6
                        popul_fds(i) = popul_fds(i) + r_mat(i,m_f);
                    end
                    counters_f_DS(m_f, f_tr, N_DS) = counters_f_DS(m_f, f_tr, N_DS) + 1;

                    % save occurrence time for slow reactions
                    if ~(m_f>2 && m_f<7) % m = 1, 2, 7, 8, 9
                        ocT_DS(counters_f_DS(m_f,f_tr,N_DS), map(m_f)) = t_f;
                    end

                    % update Modified_Next_Reaction_Method quantities
                    [D_f, m_f, T_f, P_f, a_f, Dt_f] = ...
                        Mod_NRM(D_f, Dt_f, T_f, P_f, a_f, m_f, popul_fds,rates_ds,map,GEN);

                    jj=jj+1;
                end
                counters_f_DS(11, f_tr, N_DS) = sum( counters_f_DS(1:9, f_tr, N_DS) );

                % plot each DownScaled trajectory
                if on_PC && plot_f_DS
                    % active figure changes randomly from time to time so...
                    figure(fig1); % retrieve figure I want to plot on
                    tt = f_DownScaled(:,1,f_tr);
                    pO  = f_DownScaled(:,3,f_tr);
                    pX  = f_DownScaled(:,4,f_tr);
                    pX2 = f_DownScaled(:,5,f_tr);
                    pOX2= f_DownScaled(:,6,f_tr); % dependent species
                    pYi = f_DownScaled(:,7,f_tr); % dependent species
                    pY  = f_DownScaled(:,8,f_tr);
                    if plot_all_species
                        p_un=plot(tt,pO,'--r', tt,pX,':r', tt,pX2,'-.r',...
                                  tt,pOX2,'--r',tt,pYi,':r',tt,pY,'-.r');
                    else
                        p_un=plot(tt,pO,'--r', tt,pX,':r', tt,pX2,'-.r',...
                                  tt,pY,'-.r');
                    end
                    tmp = [' (' num2str(f_tr) ' / ' num2str(f_traj) ')'];
                    title(['Downscaling attempt: ', num2str(N_sc), tmp]);
                    pause(1e-9);
                end

                % AVERAGE KMC steps per time unit inside the DS window
                cost_DS(f_tr+1,N_sc) = counters_f_DS(11, f_tr, N_DS)/t_ahead_i;

                %% calculate Inter-arrival time error here to avoid saving all occurrence times
                % minimum # of events for reactions 1,2,7,8,9
                N_r(1,f_tr) = min([ counters_f_UN(1,N_sc) counters_f_DS(1,f_tr,N_DS) ]);
                N_r(2,f_tr) = min([ counters_f_UN(2,N_sc) counters_f_DS(2,f_tr,N_DS) ]);
                N_r(3,f_tr) = min([ counters_f_UN(7,N_sc) counters_f_DS(7,f_tr,N_DS) ]);
                N_r(4,f_tr) = min([ counters_f_UN(8,N_sc) counters_f_DS(8,f_tr,N_DS) ]);
                N_r(5,f_tr) = min([ counters_f_UN(9,N_sc) counters_f_DS(9,f_tr,N_DS) ]);
                IAT_err_Un(:,f_tr,N_sc) = calc_IAT_err(ocT_Un, ocT_DS, N_r(:,f_tr));

                %% calculate the Inter-Arrival Times of Downscaled trajectories
                R1_IAT_DS = diff( ocT_DS(1:N_r(1,f_tr), 1) );
                R2_IAT_DS = diff( ocT_DS(1:N_r(2,f_tr), 2) );
                R7_IAT_DS = diff( ocT_DS(1:N_r(3,f_tr), 3) );
                R8_IAT_DS = diff( ocT_DS(1:N_r(4,f_tr), 4) );
                R9_IAT_DS = diff( ocT_DS(1:N_r(5,f_tr), 5) );

                % arrays contain NaNs, pick the minimum # of IAT available
                % -1 because the IATs are 1 less than the # of firings
                IAT_dif{1,f_tr} = R1_IAT_DS - R1_IAT_Un(1:N_r(1,f_tr)-1);
                IAT_dif{2,f_tr} = R2_IAT_DS - R2_IAT_Un(1:N_r(2,f_tr)-1);
                IAT_dif{3,f_tr} = R7_IAT_DS - R7_IAT_Un(1:N_r(3,f_tr)-1);
                IAT_dif{4,f_tr} = R8_IAT_DS - R8_IAT_Un(1:N_r(4,f_tr)-1);
                IAT_dif{5,f_tr} = R9_IAT_DS - R9_IAT_Un(1:N_r(5,f_tr)-1);

            end % generation of DS trajectories
            pause(1e-9);

            %% error calculations
            for f_tr=1:f_traj
                set1 = f_unscaled(  :,[3 4 5 8]); % keep the 4 independent species
                set2 = f_DownScaled(:,[3 4 5 8],f_tr);
                [popul_error_Un(1:4,f_tr,N_sc), ...
                 popul_error_Un(6:9,f_tr,N_sc) ] = calc_popul_err(set1, set2);
            end

            % Plot Inter-Arrival Time (IAT) error norm
            if on_PC
                fg = figure('Position',p8); hold on; grid on; box on;
                mk = '-'; % marker format
                x = prod(df) .* d_factors; % take into account previous downscalings
%               plot(x, IAT_err_Un(1,:,N_sc),mk,'DisplayName',['R 1, Events: ' num2str(floor(mean(N_r(1,:))))]); % reaction 1
%               plot(x, IAT_err_Un(2,:,N_sc),mk,'DisplayName',['R 2, Events: ' num2str(floor(mean(N_r(2,:))))]); % reaction 2
%               plot(x, IAT_err_Un(3,:,N_sc),mk,'DisplayName',['R 7, Events: ' num2str(floor(mean(N_r(3,:))))]); % reaction 7
%               plot(x, IAT_err_Un(4,:,N_sc),mk,'DisplayName',['R 8, Events: ' num2str(floor(mean(N_r(4,:))))]); % reaction 8
%               plot(x, IAT_err_Un(5,:,N_sc),mk,'DisplayName',['R 9, Events: ' num2str(floor(mean(N_r(5,:))))]); % reaction 9
                plot(x, IAT_err_Un(1,:,N_sc),mk,'DisplayName','Reaction 1');
                plot(x, IAT_err_Un(2,:,N_sc),mk,'DisplayName','Reaction 2');
                plot(x, IAT_err_Un(3,:,N_sc),mk,'DisplayName','Reaction 7');
                plot(x, IAT_err_Un(4,:,N_sc),mk,'DisplayName','Reaction 8');
                plot(x, IAT_err_Un(5,:,N_sc),mk,'DisplayName','Reaction 9');
                xlabel('Downscale factor'); ylabel('Error norm, Normalised Cost');
                title(['Downscaling attempt: ', num2str(N_sc)]);
                lgd=legend('Location','SW','NumColumns',1);
                ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
%               yline(error_threshold,'--k','Threshold','Linewidth',2,...
%                   'LabelHorizontalAlignment','left');
                pause(1e-9);
            end

            % calculate and plot the Normalized cost
            cost_DS_normed = cost_DS(:,N_sc) ./ cost_DS(1,N_sc); % normalize cost
            % multiply by [df] to take into account previous downscalings
            x = prod(df) .* [1; d_factors]; % 1 accounts for the un-scaled cost
            y = cost_DS_normed ./ prod(df);
            if on_PC
                h1 = loglog(x,y,'-o','DisplayName','Cost');
                % Fill marker with the same color as the line
                set(h1,'MarkerFaceColor',get(h1,'Color'));
            end

            if on_PC
                % plot cost & 2nd derivative on the same Figure
                figure('Position',[10,260,560,420]);
                hold on; grid on; box on;
                h1 = plot(x, y, 'k-o');
                set(h1,'MarkerFaceColor',get(h1,'Color'));
                ax=gca; ax.XScale='log'; ax.YScale='log';
                xlabel('Downscale Factor');
                ylabel('Normalized Cost');
                title(['Downscaling attempt: ', num2str(N_sc)]);
            end

            % Calculate the 2nd derivative of the cost vs downscale_factor
            % The calculation is done on the Logarithm of cost & df
            % The ddy is 2 points shorter than the input data
            xx = log(x);
            yy = log(y);
            n = length(y);
            ddy = NaN(n,1);
            % i=1 corresponds to the unscaled trajectory
            % i=2 corersponds to the downscaled trajectory with df=1
            % i=3 is the cost for the first actual downscaled trajectory
            for i=3:n-1
                ddy(i) = ( yy(i+1)-2*yy(i)+yy(i-1) ) / log(down_factor)^2;
            end

            curv = ddy(2:end-1);
            if on_PC
                yyaxis right;
                % plot the curvature at the "middle" point
                plot(x,[NaN NaN curv(2:end)' NaN],'r-o','DisplayName','curvature');
                ylabel('2nd Derivative');
                axR = gca; axR.YColor = 'r';
                yline(0,'r','HandleVisibility','off');
                xlims = xlim();
                ct = 0.05; % curvature threshold
                patch('XData',[xlims(1)*[1 1] xlims(2)*[1 1]],...
                      'Ydata',[-ct ct ct -ct],'FaceColor','g','FaceAlpha',0.2,...
                      'EdgeColor','none','HandleVisibility','off')
                % send patch backwards wrt other entities in the plot
                set(gca,'children',flipud(get(gca,'children')))
            end

            % use curvature to find the point of interest
            % The value at index i is calculated via the formula
            % F{i} - 2*F{i+1} + F{i+2}
            curv_idx = find(abs(curv)>0.07, 1);
            % Once i has been found, the choices are:
            % point i  : the safest but also most conservative choice
            % point i+1: less risky choice
            % point i+2: linearity is already broken at this point
%           curv_max_df = x(curv_idx);   % F{i}
            curv_max_df = x(curv_idx+1); % F{i+1}
%           curv_max_df = x(curv_idx+2); % F{i+2}

            % convert df from absolute to relative
            curv_max_df = curv_max_df ./ prod(df);

            % Successive fit on Cost data-points. Data shifted/adjusted
            % in order to take into account previous downscalings
            x = log( prod(df) .* [1; d_factors] );
            y = log( cost_DS_normed ./ prod(df) );
            cfs = NaN(multi_DS-1, 8);
            for i=2:multi_DS % loop over the downscale factors
                % use at max 4 points in fitting: i-2 i-1 i i+1
                fp = max([2 i-2]);   % first point to fit
                lp = i+1; % last  point to fit
                [c, s] = fit(x(fp:lp), y(fp:lp), 'poly1');
                cfs(i-1,1) = fp;        % first point in fitting
                cfs(i-1,2) = lp;        % last  point in fitting
                cfs(i-1,3) = d_factors(i); % downscale factor of last point
                cfs(i-1,4) = c.p1;      % slope of this "partial" linear fit
                cfs(i-1,5) = c.p2;      % intercept of the linear fit
                cfs(i-1,6) = s.rsquare; % R^2 of the fit
                cfs(i-1,7) = s.sse;     % Sum of squares due to error
                cfs(i-1,8) = s.rmse;    % Root mean squared error (standard error)
            end

            % find point where linear decrease of cost breaks down
            for i=1:multi_DS-1
                % if slope goes > -0.9
                if cfs(i,4) > -0.90
                    % save the previous point as the max allowed downscale factor
                    grad_max_df = cfs(max([1 i-1]),3);
                    grad_max_df_idx = find(d_factors == grad_max_df);
                    break
                end
            end

            % The maximum allowed df is determined by 2 methods:
            % using gradients from the above fit & the 2nd derivative (curvature)
            max_allowed_df = max([curv_max_df grad_max_df]);
            max_allowed_df_idx = find(d_factors == max_allowed_df);

            % variables used for criteria
            cr1 = 0;
            if max_allowed_df_idx <= 2
                cr1 = 1; % too few points for a good decision
            end
            cr4 = 0;
            if max_allowed_df < 10
                cr4 = 1; % max_allowed_df is "small" so the risk is not worth it
            end

            % Linear fit on the IAT error of each slow reaction
            % use only those points that are inside the linear decrease of cost
            x   = log( prod(df) .* d_factors(1:max_allowed_df_idx) );
            yr1 = log( IAT_err_Un(1, 1:max_allowed_df_idx, N_sc))';
            yr2 = log( IAT_err_Un(2, 1:max_allowed_df_idx, N_sc))';
            yr7 = log( IAT_err_Un(3, 1:max_allowed_df_idx, N_sc))';
            yr8 = log( IAT_err_Un(4, 1:max_allowed_df_idx, N_sc))';
            yr9 = log( IAT_err_Un(5, 1:max_allowed_df_idx, N_sc))';
            [cf_R1, S_R1] = fit(x, yr1, 'poly1');
            [cf_R2, S_R2] = fit(x, yr2, 'poly1');
            [cf_R7, S_R7] = fit(x, yr7, 'poly1');
            [cf_R8, S_R8] = fit(x, yr8, 'poly1');
            [cf_R9, S_R9] = fit(x, yr9, 'poly1');

            % check if any of the error lines is flat
            gradients(1:5,N_sc) = [cf_R1.p1 cf_R2.p1 cf_R7.p1 cf_R8.p1 cf_R9.p1]';
            cr2 = 0;
            if (any(gradients(1:5,N_sc) < 0.2)) % flat error lines
                cr2 = 1; % reject this downscaling immediately
                disp('Reject because of flat error(s)');
            end

            % linear fit using ALL "accepted" IAT error points
            % the fixed number 5 corresponds to the 5 slow reactions
            x = prod(df) .* d_factors(1:max_allowed_df_idx)' .* ones(5,1);
            x = log(x);
            y = log(IAT_err_Un(:, 1:max_allowed_df_idx, N_sc));
            % the syntax (:) turns a row vec into a column vec
            % fit() needs columns
            [cf_avg, S_avg] = fit(x(:), y(:), 'poly1');

            % average individual gradients
            b1 = mean([cf_R1.p1 cf_R2.p1 cf_R7.p1 cf_R8.p1 cf_R9.p1]);
            b2 = mean([cf_R1.p2 cf_R2.p2 cf_R7.p2 cf_R8.p2 cf_R9.p2]);

            % Cross-check that gradient & intercept are correct
            res1 = abs(cf_avg.p1 - b1);
            res2 = abs(cf_avg.p2 - b2);
            if res1 > 1e-12 || res2 > 1e-12
                error('Inconsistent gradient or intercept of error fit')
            end

            % plot linear fits on cost and error
            x = prod(df) .* d_factors(1:max_allowed_df_idx);
            y = log( cost_DS_normed(2:max_allowed_df_idx+1) ./ prod(df) );
            z = log(x);
            [cf_cost, S_cost] = fit(z,y,'poly1');

            % Logarithmically spaced x for the objective function
            x_dense = logspace(0,log10(d_factors(max_allowed_df_idx)),1000)';
            x_dense = prod(df) .* x_dense;
            z_dense = log(x_dense);
            y_cost  = exp(cf_cost(z_dense));
            y_error = exp(cf_avg( z_dense));
            if on_PC
                figure(fg); % get the appropriate figure to plot on
                plot(x_dense, y_cost ,'r--','LineWidth',2,'HandleVisibility','off');
                plot(x_dense, y_error,'b--','LineWidth',2,'HandleVisibility','off');
            end
            % save the cost-error data
            CE_data{N_sc,1} = x_dense;
            CE_data{N_sc,2} = y_cost;
            CE_data{N_sc,3} = y_error;

            % calculate objective function
            obj{1,1} = 1; obj{1,2} =  0.5;
            obj{2,1} = 1; obj{2,2} =  1.0;
            obj{3,1} = 1; obj{3,2} =  2.0;
            obj{4,1} = 1; obj{4,2} = 10.0;

            pos = [.7 .12 .4 .2];
            for i=1:4
                obj{i,3} = obj{i,1} .* y_cost + obj{i,2} .* y_error;
                [mn, factor_idx] = min(obj{i,3});
                if on_PC
                    DN = ['\alpha = ' num2str(obj{i,1}) ', \beta = ' num2str(obj{i,2})];
                    pl = plot(x_dense,obj{i,3},'LineWidth',2,'DisplayName',DN);
                    cl = pl.Color;
                    xline(x_dense(factor_idx),'Color',cl,'HandleVisibility','off');
                    text(x_dense(factor_idx),pos(i),num2str(round(x_dense(factor_idx))),...
                        'HorizontalAlignment','center');
                end
                % Which Obj.Func to use to choose the df.
                % Here we use #3 (a=1, b=2), as defined by alpha_beta_pair
                if i == alpha_beta_pair
                    chosen_df_idx = factor_idx;
                    chosen_df = x_dense(chosen_df_idx);
                end
            end
            if on_PC
                % place legend outside and adjust figure size
                set(gcf,'units','pixels'); % set units
                figure_size =  get(gcf,'Position'); % get figure size
                set(lgd,'Location','NorthEastOutside')% place legend outside
                set(lgd,'units','pixels'); % set legend units
                legend_size = get(lgd,'Position'); % get legend size
                figure_size(3) = figure_size(3) + legend_size(3); % new width
                set(gcf, 'position', figure_size); % set new figure size
                xlims = xlim();
                xlim([1 xlims(2)])
            end
            clear pos DN pl mn cl % clear temporary variables

            %% calculate acceptable error with respect to error @ df=1
            % IAT error @ df=1. This is NOT equal to zero, even though
            % no downscaling is done
            err_df1 = y_error(1);
            % IAT error @ chosen df
            err_dfc = y_error(chosen_df_idx);
            % error difference as orders of magnitude
            err_dif_magn = abs( log10(err_df1) - log10(err_dfc) );

            % error increase rate (gradient in log-log plot)
            err_grad = cf_avg.p1;
            % the gradient of the average error is ~0.5
            % so error increases by 1 order of magn. with a df of ~100

            %% accept? reject ?
            % cr1: too few points for a good decision
            % cr2: at least one error line is flat
            % cr3: error of chosen df is below threshold
            % cr4: max_allowed_df is "small"
            % cr5: error increase wrt to df=1
            cr3 = err_dfc < error_threshold;
            cr5 = err_dif_magn < e_inc;

            % few points OR flat error(s) OR small max_allowed_df
            if (cr1 || cr2 || cr4) 
                accept = 0;
            elseif (cr3 && cr5) % overall fitted error is small enough
                accept = 1;
            else % High error. Enough points. Error lines OK.
                accept = 0;
            end

            %% is Down-scaling accepted or rejected ?
            if(accept == 1)
                %% Generate the trajectory for the chosen df
                % It is unlikely that the accepted df is among the few dfs
                % that trajectories exist fo. Generate it here
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % set the correct values to the downscaled rates
                % The "chosen_df" is the absolute (total) downscale factor
                % Rates are reduced gradually hence the relative df here
                chosen_df_rltv = chosen_df / prod(df);
                rates_ds(4) = rates(4) / chosen_df_rltv;
                rates_ds(5) = rates(5) / chosen_df_rltv;
                % set or reset the necessary simulation parameters
                ocT_DS = NaN(size(ocT_DS)); % erase previous data
                % RESTORE the state of the "slow" RNGs
                for i=1:size(GEN,1)
                    GEN{i}.State = GEN_states{i};
                end
                sampling_t_f = sampling_t; samples_f = 0;
                popul_fds=popul; t_f=t; % future downscaled population
                a_f=a; P_f=P; T_f=T; Dt_f=Dt; [D_f, m_f]=min(Dt_f);
                while(t_f < t + t_ahead_width(N_sc))
                    t_f = t_f + D_f;
                    % the 2nd condition prevents sampling 1 more point beyond the intented interval
                    while(t_f > sampling_t_f && samples_f+1 <= N_points)
                        samples_f = samples_f + 1;
                        f_DownScaled(samples_f,  1, multi_DS +1) = sampling_t_f;
                        f_DownScaled(samples_f,  2, multi_DS +1) = D_f;
                        f_DownScaled(samples_f,3:8, multi_DS +1) = popul_fds;
                        f_DownScaled(samples_f,9:10,multi_DS +1) = rates_ds(4:5);
                        sampling_t_f = sampling_t_f + sampling_dt;
                    end
                    % update population & counters
                    for i=1:6
                        popul_fds(i) = popul_fds(i) + r_mat(i,m_f);
                    end
                    counters_f_DS(m_f, multi_DS+1, N_DS) = counters_f_DS(m_f, multi_DS+1, N_DS) + 1;
                    % save occurrence time for slow reactions
                    if ~(m_f>2 && m_f<7) % m = 1, 2, 7, 8, 9
                        ocT_DS(counters_f_DS(m_f,multi_DS+1,N_DS), map(m_f)) = t_f;
                    end
                    % update Modified_Next_Reaction_Method quantities
                    [D_f, m_f, T_f, P_f, a_f, Dt_f] = ...
                        Mod_NRM(D_f, Dt_f, T_f, P_f, a_f, m_f, popul_fds,rates_ds,map,GEN);
                    jj=jj+1;
                end
                counters_f_DS(11, multi_DS+1, N_DS) = sum( counters_f_DS(1:9, multi_DS+1, N_DS) );
                % plot the final DownScaled trajectory
                if on_PC && plot_f_DS
                    % active figure changes randomly from time to time so...
                    figure(fig1); % retrieve figure I want to plot on
                    tt  = f_DownScaled(:,1,end);
                    pO  = f_DownScaled(:,3,end);
                    pX  = f_DownScaled(:,4,end);
                    pX2 = f_DownScaled(:,5,end);
                    pOX2= f_DownScaled(:,6,end); % dependent species
                    pYi = f_DownScaled(:,7,end); % dependent species
                    pY  = f_DownScaled(:,8,end);
                    p_un=plot(tt,pO,'k', tt,pX,'k', tt,pX2,'k',tt,pY,'k');
                    title('Downscaling attempt: Final Accepted Trajectory');
                    pause(1e-9);
                end
                f_traj = 1; % restore variable's value
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Accepted = Accepted + 1;
                Decision(Accepted+Rejected+Upscaled)=1;
                % update time, population and rates.
                t = t_f;
                popul = popul_fds;
                rates = rates_ds;
                sampling_t = sampling_t_f;
                % append the accepted Down-Scaled trajectory to
                % the array of the "accepted" data-points
                data(samples+1:samples+samples_f, :) = f_DownScaled(:,:,end);
                samples = samples + samples_f;
                % save occurrence time for slow reactions
                for ii=[1 2 7 8 9] % loop over the slow reactions: 1 2 7 8 9
                    % retrieve the Non-NaN occurrence times
                    xx = ocT_DS(~isnan(ocT_DS(:, map(ii))),map(ii));
                    xi = counters(ii) + 1;
                    xf = counters(ii) + size(xx,1);
                    ocT(xi:xf, map(ii)) = xx;
                end
                % Add the KMC step counters of the accepted DownScaled (DS) trajectory to the reference counters
                counters(:,1) = counters(:,1) + counters_f_DS(1:9,end,N_DS);
                disp([num2str(N_sc),' DOWN-SCALED Solution Accepted']);
                disp(' ');
                df(N_sc) = chosen_df_rltv; % store the downscale factor
                update_T_P
            else % Down_Scaling was rejected
                 disp(' ');
                 disp(num2str(N_sc));
                 disp(['Normalized Population / IAT error > ', num2str(error_threshold)]);
                 Reject
            end % accept DS
        end % early stop
    end % perform DS
end % main KMC loop

% disp(datetime('now'));
% dur = toc; % duration of run in seconds
% Simulation speed, including the downscaling attempts
% disp(['Speed: ' num2str(t_max/dur,'%.3f') ' KMC time units / wall sec']);

%% calculate concentrations from molecular populations
data_ppl = data; % copy population array
data_cnc = data; % create another copy
data_cnc(:,3:8) = data(:,3:8) ./ NV; % convert to concentration

% data = data_cnc; % plot concentrations

%% tidy-up things
ocT_tmp = cell(5,2);
for i=1:5
    ocT_tmp{i,1} = map_inv(i); % number of reaction
    ocT_tmp{i,2} = ocT(~isnan(ocT(:,i)),i); % occurrence times    
end
ocT = ocT_tmp; % overwrite to get rid of all NaNs
clear ocT_tmp ocT_Un ocT_DS % clear unnecessary variables

% Percentage of reaction execution
counters(:,2) = counters(:,1) / sum(counters(:,1)) * 100;

% create a line with zeros before the sums
counters(end+2,:) = sum(counters, 1); % dim=1, sum elements in columns

% sum all KMC steps, across all downscaling attempts.
counters_f_DS(end, end, N_DS+1) = sum( counters_f_DS(end, end, 1:N_DS) );
counters_f_UN(end,      N_UN+1) = sum( counters_f_UN(end,      1:N_UN) );

%% Clean variables
clear P_O P_X P_X2 P_Y O X X2 OX2 Yi Y
clear a12 a13 a14 a15 a16 a17 a18 a0 
clear a1 a2 a3 a4 a5 a6 a7 a8 a9
clear h1 h2 h3 h4 h5 h6 h7 h8 h9
clear j m r1 r2 t2 z tt avg sdv 
if ~DS
    clear counters_f_DS counters_f_UN
    clear popul_error_Un
end

%% Plot
if on_PC
    t=data(1:steps+1,1);

    fig01=figure('Position',[575 260 560 420]); hold on;
    %  3  4   5   6   7   8.     6 and 7 are dependent species
    %  O  X  X2  OX2  Yi  Y
    for i=[3 4 5 6 7 8]
        style='-'; lw=1;
        y = data(1:steps+1,i);
        if (i==6 || i==7)
            style='--'; lw=1;
        end
        plot(t,y,style,'LineWidth',lw);
    end
    hold off; grid on; box on;
    legend('O','X','X_2','OX_2','Y_i','Y','NumColumns',2,'AutoUpdate','off');
    % legend('O','X','X_2','Y','NumColumns',1,'AutoUpdate','off');
    xlabel('KMC time (min)');
    ylabel('Population');
    % ylabel('Concentration [nM]')
    xlim([0 t_max]);

    ylims = get(gca,'YLim');
    y = [ylims fliplr(ylims)];
    for i=1:N_sc
        ti_DS = t_DS_pos(i);

        t_from = ti_DS + sampling_dt;
        t_upto = ti_DS + t_ahead_width(i);
        x = [t_from t_from t_upto t_upto];
        if Decision(i)==1 % downscaling Accepted
            col = [0, 1, 0]; % bright green
        elseif Decision(i)==2 % downscaling Accepted after UP-scaling
            col = [0.13 0.55 0.13]; % Forest green (web)
        elseif Decision(i)==3 % downscaling Rejected
            col = [1, 0, 0]; % Red
        elseif Decision(i)==4 % UP-scaling
            col = [0, 0.28, 0.67];	% Cobalt Blue
        end
        patch('XData',x,'Ydata',y,'FaceColor',col,...
            'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off');
        text((t_upto+t_from)/2, ylims(2)-20-mod(i,2)*2, num2str(i),'Color','k',...
                         'FontSize',12,'HorizontalAlignment','center');
    end

    % Dt versus Time
    fig03=figure('Position',p5);
    semilogy(t,data(1:steps+1,2)); grid on;
    xlabel('KMC Time (min)');
    ylabel('KMC dt (min)');
    xlim([0 t_max]);

    % KMC steps vs KMC time
    fig=figure('Position',p7);
    % remove NaN from the steps_vs_t
    steps_vs_t = rmmissing( steps_vs_t );
    plot(steps_vs_t(:,1), steps_vs_t(:,2)); grid on;
    legend('Scaled','Location','SE');
    xlim([0 t_max]);
    title('Cumulative KMC steps vs KMC Time')
    xlabel('KMC Time (min)');
    ylabel('KMC Steps');
end

%% plot quasi-equilibrated quantities
if on_PC && plot_QE
    %  3  4   5   6   7   8   9   10
    %  O  X  X2  OX2  Yi  Y  chi  phi
    t  = data(:, 1);
    O  = data(:, 3);
    X  = data(:, 4);
    X2 = data(:, 5);
    OX2= data(:, 6);

    figure();
    plot(t, X .* X); hold on; grid on;
    plot(t, beta .* X2);
    legend('X * X','\beta * X_2');
    title('X*X = \beta*X_2');

    figure();
    plot(t, O .* X2); hold on; grid on;
    plot(t, alfa .* OX2);
    legend('O * X_2','\alpha * OX_2');
    title('O * X_2 = \alpha * OX_2');
end
