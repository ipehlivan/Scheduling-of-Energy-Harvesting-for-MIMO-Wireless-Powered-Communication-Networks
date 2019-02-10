% SCA ALGORITHM
function [ am_optval, v_opt,remain_output,t_1,cvx_correct,step_m ] = am_algorithm_nonscheduled_latest( Hi,  K, N ,Pa ,L,E,iter )
%--------------------------------------------------------------------------
% Comments refer to paper: 

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO 
%Wireless Powered Communication Networks," in IEEE Communications Letters, 
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function is an alternating minimization algorithm for non-scheduled 
%case that solves optimization problem at eqn.(3a).

%K: # antennas
%N: # users
%Pa: maximum transmit power
%L: # RF chain
%Hi: Channel matrix
%E: required energy
%S: number of RF chain
%iter: number of iteration


%% Initialization: Algorithm 1 (Line 1)
lambda=1; % Multiplier of the penalty function at eqn. (3a)
V_1_u=zeros(K,K); % solution at the (q-1), previous, iteration.
ei_V_1_u(1:3)=zeros(1,3); % storing first 3 eigenvalues for rank detection
mem=0; % memory of the optimal value
step_m=10; % memory of the step size
mem_t=0; % memory of the delay
ck=ones(1,iter); % error flag for every iteration
countt=1; % counting variable

%% Setting precision
% precision is being selected adaptively.

[ cvx_opt_prec, V_prec,t_prec,cvx_prec ] = am_algorithm_cvx_single( Hi,  K, N ,Pa ,L,E,1,V_1_u,lambda,mem_t );
prec=max(eigs(V_prec,1)*0.001,10^-7);
prect=max(t_prec*0.001);

%% Algorithm 1 starts
for ite=1:iter
    %% Solving problem : Algorithm 1 (Line 3)
    % Solving cvx problem. Check function itself for information.
    [ cvx_optval, V_1,t_1,cvx_status ] = am_algorithm_cvx_single( Hi,  K, N ,Pa ,L,E,1,V_1_u,lambda,mem_t ); % Solve the problem
    
    
    ei_V_1(1:3)=flip(sort(eigs(V_1(:,:),3)));% Obtaining first 3 eigenvalue
    % Important tip:
    % eigs returns ascending or descending eigenvalues depending on
    % platform dont forget to sort it.
    ei_rem=ei_V_1(2:3); % excluding biggest eigenvalue
    str = cvx_status; % status of the cvx solver
    
    %% If cvx fails precaution
    % Cvx sometimes fails, this part is for the case it fails.
    if(prod("Failed"==str(1:6)))
        lambda=lambda/2; % reduce the lambda
        % Solve it again. 
        [ cvx_optval, V_1,t_1,cvx_status2 ] = am_algorithm_cvx_single( Hi,  K, N ,Pa ,L,E,1,V_1_u,lambda,mem_t );
        str2=cvx_status2;
        if(prod("Failed"==str2(1:6)))
            ck(countt)=0; % If still failed, flag the error.
        end
    end
    %% Projection case: Algorithm 1 (Line 8)
    % In this part I project the result to rank 1 and continue searching
    % around it
    V_1_u=V_1; % Memory assignment
    step_m=real(mem-cvx_optval); % step assignment. Taking real for numerical issues.
    lo=((2*min(lambda)>10^6) && (step_m>=0)  && (step_m<=10^-9));
    proj_flag=0; % projected flag
    % If it is not rank 1 and at q-treshold. Project it.
    % Check the projection function
    if( ( (ite==floor(iter*0.5)) || lo ) && (max(ei_rem(:))>10^-10))     
        [V_1_u,t_1,lambda] = rank_1_proj(V_1,Hi,Pa,N,L,K,E,0);
        proj_flag=1;
    end
    %% Updating lambda: Algorithm 1 (Line 5)
    % Updating lambda if rank >1 and first eigenvalue is not increased
    % larger than the desired step.
    if(2*lambda<=10^6 &&(max(ei_V_1(2:3))>=10^-10)&& proj_flag==0 )
        lambda=2*lambda*((ei_V_1(1)<=(prec+ei_V_1_u(1)))&&max(ei_rem(:))>10^-10)+...
            not((ei_V_1(1)<=prec+ei_V_1_u(1))&&max(ei_rem(:))>10^-10)*lambda;
    end
    
    %% Memory variables update: Algorithm 1 (Line 3)
    % Updating memory variables.
    ei_V_1_u(1:3)=flip(sort(eigs(V_1_u(:,:),3)));
    remain=max(ei_V_1_u(2:3));
    step_m=real(mem-cvx_optval);
    mem=cvx_optval;
    mem_t=sum(t_1);
    
    %% Check convergence: Algorithm 1 (Line 10)
    % Checking convergence 
    if(step_m>=0 && step_m<=10^-6 && max(ei_V_1(2:3))<=10^-10)
        break;
    end
    %%
    % Counter update.
    countt=countt+1;
end

%% Last case projection: Algorithm 1 (Line 8)
% Projecting the solution to rank 1 if algorithm fails to find rank 1 sol.
% Last case precaution
if((remain>10^-10))
    [V_1,t_1,lambda] = rank_1_proj(V_1,Hi,Pa,N,L,K,E,1);
    ei_V_1(1:3)=flip(sort(eigs(V_1(:,:),3)));
    remain=max(ei_V_1(2:3));
    % update error variables
end

%% Output assignment: Algorithm 1 (Line 11)
% Assigning output
v_opt=V_1; % optimal matrix
am_optval=t_1; % optimal delay
remain_output=remain; % rank error
cvx_correct=prod(ck); % cvx fail flag
end