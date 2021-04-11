%--------------------------------------------------------------------------
% Comments refer to paper:

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO
%Wireless Powered Communication Networks," in IEEE Communications Letters,
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

% This is a simulation script for performance analysis. CVX should be
% installed in order to run this simulation. All functions are commented in
% detail.


%% Safety
clear
close all

%% HPC adjustment ,cvx and parpool

% This part is specific for our HPC server thus it is commented out.
% Reader might find it useful.


% setenv('TZ','Europe/Istanbul')
% % set_cvx
addpath('Required_codes');
% cls=parcluster;
% cls.NumWorkers=60;
% saveProfile(cls);
% parpool('local',60)

%% System constants
m_slot=5; % maximum number of slot
na=1;%number of how many user values is used (not active in this sim.)
max_it=1500; % maximum number of iteration
scaler=10^7; % For cvx safety
K=2^(4); % Number of antennas
RR=16; % Number of reciever antenna
N=(na)*10+50; % number of nodes
%% This part is not active for this simulation
k=1;% number of how many reciever antenna values is used
n=5;%number of how many reciever antenna values is used
%%
to=120;%number of different topologies
l=log2(K)+1; % number of how many rf chain values is used
s=2; % number of how many different rf chain values is used

rician_K_factor=0; % rician K factor - Taken as 0 to have rayleigh fading during simulations.
% Non-zero values do not correspond to Rician fading. User must alter the channel_creator.m function to 
% have Rician fading. This is a unused variable, remainder of a previous coding process.

fc=28*10^9; % Carrier frequency (Hz)
c=299792458; % speed of light(m/s)
dalgalamba= c/fc; % wavelength (m)
Pa=db2pow(40)/1000; % Total power of beams, Watt, 40 dBm
Ro=0.2; % Energy harvesting efficency
radiu=4; % Radius, meter

%% Matrix_decleration
simulation_step_am=zeros(n,l,s,to); % AM algorithm step size memory
simulation_objective=zeros(n,l,to); % SDR optimal value memory
simulation_objective_1=zeros(n,l,to); % SIMO case optimal value mem.
% SIMO stands for single input multiple output system.
simulation_objective_am=zeros(n,l,s,to); % AM algorithm optimal value mem.
simulation_error_am=zeros(n,l,s,to); % AM rank error memory
simulation_cvx_error_am=zeros(n,l,s,to); % CVX failure memory
%% Iterations

% This part uses parallel computing toolbox.
tic % time
parfor toa=1:120 % for every topology
    
    %% Energy requirement
    eng_scale=(10^(-6))*scaler; %Joule
    % Getting rid of numerical errors
    % I multiply both sides of the eqn.(2b) with scaler to avoid cvx
    % numerical errors.
    E=ones(N,1)*eng_scale; % Energy requirement
    %% Topology creation
    spacing=dalgalamba/4; % spacing between antennas
    d = topology_creator_3D(K,N,RR,radiu,spacing,1); % topology creation
    d_1 = topology_creator_3D(1,N,RR,radiu,spacing,1); % topology creation for SIMO CASE
    [Hi] = channel_creator(K, N ,RR,d,dalgalamba,rician_K_factor,scaler, Ro); % channel matrix creation
    [Hi_1] = channel_creator(1, N ,RR,d_1,dalgalamba,rician_K_factor,scaler, Ro);% channel matrix creation for SIMO CASE
    
    %% Calculations
    for la=1:5 % for 5 different rf chain values
        L=(2^(la-1)); % number of RF chain
        for sa=1:2  % for 2 different time slot values
            slota=1*(sa==1)+(sa==2)*(m_slot); % number of time slot
            
            if((la<=(log2(K)+1))) % this part is useful for different antenna trials
                % not active in this simulation.
                
                if(slota>1 && K>1) %if it is scheduled and not SIMO
                    
                    % solving the problem. Check the code.
                    [ am_optval, v_opt,remain_output,t_1,cvx_correct,step_m ] = am_algorithm_scheduled_latest( Hi,  K, N ,Pa ,L,E,slota,max_it);
                    
                    simulation_step_am(na,la,sa,toa)=step_m;% AM algorithm step size memory
                    simulation_objective_am(na,la,sa,toa)=am_optval;% AM algorithm optimal value mem.
                    simulation_error_am(na,la,sa,toa)=remain_output;% AM rank error memory
                    simulation_cvx_error_am(na,la,sa,toa)=cvx_correct;% CVX failure memory
                end
                
                if(slota==1) %if it is non-scheduled
                    
                    [ am_optval, v_opt,remain_output,t_1,cvx_correct,step_m ] = am_algorithm_nonscheduled_latest( Hi,  K, N ,Pa ,L,E,max_it );
                    
                    simulation_step_am(na,la,sa,toa)=step_m;% AM algorithm step size memory
                    simulation_objective_am(na,la,sa,toa)=am_optval;% AM algorithm optimal value mem.
                    simulation_error_am(na,la,sa,toa)=remain_output;% AM rank error memory
                    simulation_cvx_error_am(na,la,sa,toa)=cvx_correct;% CVX failure memory
                    
                    [ V,t_opt ,optimum_SDR,cx ] = Main_solver_schedule_am( K, N ,Pa ,L,Hi ,E,1); % solving the relaxed problem
                    simulation_objective(na,la,toa)=optimum_SDR; % SDR optimal value memory
                    V_n=max(E./Hi_1); % SIMO case
                    t_opt_n=V_n/Pa; % SIMO case
                    simulation_objective_1(na,la,toa)=sum(t_opt_n); % SIMO case optimal value mem.
                end
            else
                % This part is not active for this simulation. It can be
                % useful for simulations with different antenna values.
                simulation_error_am(na,la,sa,toa)=0;
                simulation_step_am(na,la,sa,toa)=0;
                simulation_objective_am(na,la,sa,toa)=0;
                simulation_objective(na,la,toa)=0;
                simulation_cvx_error_am(na,la,sa,toa)=0;
                
            end
        end
    end
end

t=toc; % simulation time
t=t/60; % simlation time in minutes
clear cls
%% Save workspace
formatSpec = 'simulation_user%d_to120.mat';
str = sprintf(formatSpec,N);
save(str);
