% SCA ALGORITHM
function [ am_optval, v_opt,remain_output,t_1,cvx_correct,step_m ] = am_algorithm_scheduled_latest( Hi,  K, N ,Pa ,L,E,S,iter )
%--------------------------------------------------------------------------
% Comments refer to paper: 

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO 
%Wireless Powered Communication Networks," in IEEE Communications Letters, 
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function is an alternating minimization algorithm for scheduled 
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
lambda=1*ones(1,S); % Multiplier of the penalty function at eqn. (3a)
lambda_mem=lambda; % lambda memory for cvx fail precaution
mem_t=0;  % memory of the delay
V_1_u=zeros(K,K,S);% solution at the (q-1), previous, iteration.
ei_V_1_u(:,1:3)=zeros(S,3);% storing first 3 eigenvalues for rank detection
mem=0;% memory of the optimal value
step_m=10;% memory of the step size
ck=ones(1,iter);% error flag for every iteration
countt=1;% counting variable

%% Setting precision
% precision is being selected adaptively.

[ cvx_opt_pr, V_pr,t_pr,cvx_pr ] = am_algorithm_cvx( Hi,  K, N ,Pa ,L,E,S,V_1_u,lambda,mem_t );
prec=zeros(1,S);
for sss=1:S
    prec(sss)=max(0.001*eigs(V_pr(:,:,sss),1),10^-7);
end
prec=ones(1,S)*min(prec);
prect=min(t_pr*0.001);

%% Algorithm 1 starts
for ite=1:iter

    %% Solving problem : Algorithm 1 (Line 3)
    % Solving cvx problem. Check function itself for information.
    [ cvx_optval, V_1,t_1,cvx_status ] = am_algorithm_cvx( Hi,  K, N ,Pa ,L,E,S,V_1_u,lambda,mem_t );
    ei_V_1=zeros(S,3);
    for sss=1:S % for every time slot
        % Obtaining first 3 eigenvalue
        ei_V_1(sss,1:3)=flip(sort(eigs(V_1(:,:,sss),3)));
        % Important tip:
        % eigs returns ascending or descending eigenvalues depending on
        % platform dont forget to sort it.
    end
    ei_rem(:,:)=ei_V_1(1:S,2:3);% excluding biggest eigenvalue
    str = cvx_status;% status of the cvx solver
    %% If cvx fails precaution
    % Cvx sometimes fails, this part is for the case it fails.
    if(prod("Failed"==str(1:6)))
        lambda=lambda_mem; % return to last solvable lambda
        % Solve it again. 
        [ cvx_optval, V_1,t_1,cvx_status2 ] = am_algorithm_cvx( Hi,  K, N ,Pa ,L,E,S,V_1_u,lambda,mem_t);
        for sss=1:S
            ei_V_1(sss,1:3)=flip(sort(eigs(V_1(:,:,sss),3)));
        end
        ei_rem(:,:)=ei_V_1(1:S,2:3);
        str2=cvx_status2;
        if(prod("Failed"==str2(1:6)))
            ck(countt)=0;% If still failed, flag the error.
        end
    end
    lambda_mem=lambda; % update lambda memory
    
    %% Updating lambda: Algorithm 1 (Line 5)
    % Updating lambda if rank >1 and first eigenvalue is not increased
    % larger than the desired step.
    for k=1:S % for every time slot
        if(2*lambda(k)<=10^6 &&(max(ei_V_1(k,2:3))>=10^-10) )
            lambda(k)=2*lambda(k)*((ei_V_1(k,1)<=(prec(k)+ei_V_1_u(k,1)))&&max(ei_rem(:))>10^-10)+...
                not((ei_V_1(k,1)<=prec(k)+ei_V_1_u(k,1))&&max(ei_rem(:))>10^-10)*lambda(k);
        end
    end
    
    %% Projection case: Algorithm 1 (Line 8)
    % In this part I project the result to rank 1 and continue searching
    % around it
    lo=( (2*min(lambda)>10^6) && (step_m>=0)  &&(step_m<=10^-9) ) ;
    V_1_u=V_1;% Memory assignment
    if( ( (ite==ceil(0.50*iter)) || lo) && (max(ei_rem(:))>10^-10) )
        % If it is not rank 1 and at q-treshold. Project it.
        % Check the projection function
        [V_1_u,t_1,lambda] = rank_1_proj_sch(V_1,t_1,prect,Hi,Pa,N,L,K,E,S,0);
    end
    
    %% Zero matrix cancelation
    
    % Some slots can appear redundant. In that case 
    % zero delay and zero beamforming matrix are assigned to these slots. 
    % To avoid oscillation, a selected non-zero matrix is also assigned to 
    % a redundant(!) slot and half of the selected non-zero matrice's time 
    % is assigned to both slots. Total delay remains the same. Two 
    % scenarios are equivalent. This assignment avoids oscillation and also 
    % preserves the degree of freedom. The algorithm might decide that slots 
    % are not redundant at the following iterations.
   
    for sss=1:S
        if(t_1(sss)<prect) % if time is smaller desired precision
            V_1_u(:,:,sss-1)=V_1_u(:,:,sss-1)/2;
            V_1_u(:,:,sss)=V_1_u(:,:,sss-1);
            t_1(sss-1)=t_1(sss-1)/2;
            t_1(sss)=t_1(sss-1);
            lambda(sss)=lambda(sss-1);
        end
        ei_V_1_u(sss,1:3)=flip(sort(eigs(V_1_u(:,:,sss),3)));
    end
    %% Memory variables update: Algorithm 1 (Line 3)
    % Updating memory variables.
    ei_rem_u(:,:)=ei_V_1_u(1:S,2:3);
    remain=max(ei_rem_u(:));
    step_m=real(mem-cvx_optval);
    mem=cvx_optval;
    mem_t=sum(t_1);
    %% Check convergence: Algorithm 1 (Line 10)
    % Checking convergence 
    if(step_m>=0 && step_m<=10^-6 && max(ei_rem(:))<=10^-10)
        break;
    end
    %%
    % Counter update.
    countt=countt+1;
end
%% Last case projection: Algorithm 1 (Line 8)
% Projecting the solution to rank 1 if algorithm fails to find rank 1 sol.
% Last case precaution
if( remain>10^-10)
    [V_1,t_1,lambda] = rank_1_proj_sch(V_1,t_1,prect,Hi,Pa,N,L,K,E,S,1);
    
    for sss=1:S
        if(t_1(sss)<prect)
            V_1(:,:,sss-1)=V_1(:,:,sss-1)/2;
            V_1(:,:,sss)=V_1(:,:,sss-1);
            t_1(sss-1)=t_1(sss-1)/2;
            t_1(sss)=t_1(sss-1);
        end
        ei_V_1(sss,1:3)=flip(sort(eigs(V_1(:,:,sss),3)));
    end

    ei_rem_e(:,:)=ei_V_1(1:S,2:3);
    remain=max(ei_rem_e(:));
end

%% Output assignment: Algorithm 1 (Line 11)
% Assigning output
v_opt=V_1; % optimal matrices
am_optval=sum(t_1); % optimal delay
remain_output=remain; % rank error
cvx_correct=prod(ck); % cvx fail flag


end