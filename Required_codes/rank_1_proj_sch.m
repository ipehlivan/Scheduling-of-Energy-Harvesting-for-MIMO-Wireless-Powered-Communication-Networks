function [V_1_u,t_1,lambda_out] = rank_1_proj_sch(V_1,t_1,prect,Hi,Pa,N,L,K,E,S,endo)
%--------------------------------------------------------------------------
% Comments refer to paper:

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO
%Wireless Powered Communication Networks," in IEEE Communications Letters,
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function projects given matrices to rank 1 matrices: Algorithm 1
% (Line 8)

%K: # antennas
%N: # users
%Pa: maximum transmit power
%L: # RF chain
%Hi: Channel matrix
%E: required energy
%S: number of RF chain
%V_1: Feasible matrices
%t_1: Feasible delays
%prect: desired precision
%endo: if endo is 0, the algorithm finds a suitable lambda.

%% Zero matrix cancelation
% Some slots can appear redundant. In that case
% zero delay and zero beamforming matrix are assigned to these slots.
% To avoid oscillation, a selected non-zero matrix is also assigned to
% a redundant(!) slot and half of the selected non-zero matrice's time
% is assigned to both slots. Total delay remains the same. Two
% scenarios are equivalent. This assignment avoids oscillation and also
% preserves the degree of freedom. The algorithm might decide that slots
% are not redundant at the following iterations.

prect=min(prect); % defining precision
for sss=1:S
    if(t_1(sss)<prect)
        V_1(:,:,sss-1)=V_1(:,:,sss-1)/2;
        V_1(:,:,sss)=V_1(:,:,sss-1);
        t_1(sss-1)=t_1(sss-1)/2;
        t_1(sss)=t_1(sss-1);
    end
end
%% Error calculation
for sss=1:S % for every time slot
    dum_mat(:,:)=V_1(:,:,sss);
    ei_V_1(sss,1:3)=flip(sort(eigs(dum_mat,3)));
end
ei_rem(:,:)=ei_V_1(1:S,2:3); % calculate rank violation error.
%% Initilization
V_1_u=V_1; % rank 1 matrix safey initialization
V_int=V_1; % scaled rank 1 matrix safey initialization
%% Projection
for sss=1:S % for every time slot
    if((max(ei_rem(sss,:))>=10^-10)) % rank error vs desired range
        [V,D] = eig(V_1(:,:,sss)); % eigenvalue decompostion
        Va=sum((V*abs(D)),2)/trace(abs(D)); % summation of eigenvectors
        Va=hybrid_beam_gen(Va,L); % reconstruct vector to satisfy eqn.(2d)
        V_1_u(:,:,sss)=Va*trace(abs(D))*Va'; % Reconstruct rank 1 matrix
        V_1_u(:,:,sss)=(V_1_u(:,:,sss)+V_1_u(:,:,sss)')/2; % Getting rid of numerical errors
    end
end
%% Scaling
% This part solves the problem in eqn.(2a) by scaling the generated matrices.
cvx_begin quiet
cvx_solver Sedumi
cvx_precision best
variable alfa(S) % scaling variable
expression V_int(K,K,S)
variable t_1(S)
minimize sum(t_1);

subject to

t_1==real(t_1);% eqn.(2g)
alfa==real(alfa);
for kkk=1:S % for every time slot
    if(kkk>=2)
        alfa(kkk)*trace(V_1_u(:,:,kkk))<=alfa(kkk-1)*trace(V_1_u(:,:,kkk-1));
        % ordering for implementation reasons. Check zero matrix
        % cancelation at am scheduling algorithm.
    end
    alfa(kkk)*trace(V_1_u(:,:,kkk))<=Pa*(t_1(kkk));% eqn.(2c)
    t_1(kkk)>=0;% eqn.(2g)
    alfa(kkk)>=0;
    V_int(:,:,kkk)=alfa(kkk)*V_1_u(:,:,kkk); % defining scaled matrices.
end

for j=1:N %for every user   
    Hi_e=Hi(1:K,1:K,j);
    real(sum_mat(Hi_e,V_int,S))>=E(j);% eqn.(2b)
end
cvx_end

V_1_u(:,:,:)=V_int(:,:,:); % assigning output.

%% Finding suitable lambda

  % This part finds suitable lambda for the algorithm to continue searching
  % around the projected matrix. Cvx may fail for high lambda values.This 
  % part  ensures cvx will not fail and lambda is the smallest lambda that
  % rank of the matrix is preserved.
  
lambda_mat(:,:)=ones(20*S,S); % lambda initialization
lambda_out=ones(1,S); % lambda output initialization
countt=1;% counting variable

while(min(lambda_mat(:))<10^6 && endo==0) % if endo is 0 and lambda is small enough.
    %% Solve the problem
    [ cvx_opt, V_mem,t_mm,cvx_s ] = am_algorithm_cvx( Hi,  K, N ,Pa ,L,E,S,V_1_u,lambda_mat(countt,:),t_1);
    str = cvx_s;
    
    if((prod("Failed"==str(1:6))))
        lambda_out=lambda_mat(countt-1,:);
        break;
        % if cvx fails. Return the last succesful lambda.
    else
        lambda_out=lambda_mat(countt,:);
    end
    %% Calculate error
    for sss=1:S
        ei_V_mem(sss,1:3)=flip(sort(eigs(V_mem(:,:,sss),3)));
    end
    ei_rem(:,:)=ei_V_mem(:,2:3);
    if(countt>1)
        if(max(ei_rem(:))<=(10^-6) || (prod(lambda_mat(countt,:)==lambda_mat(countt-1,:))))
            lambda_out=lambda_mat(countt,:);
            V_1_u=V_mem;
            t_1=t_mm;
            break;
            % if rank is preserved, break.
        end
    end
    %% Update lambda
    for sss=1:S
        if(max(ei_rem(sss,:))>(10^-6) && lambda_mat(countt,sss)<10^6)
            lambda_mat(countt+1,sss)=lambda_mat(countt,sss)*2;
        else
            lambda_mat(countt+1,sss)=lambda_mat(countt,sss);
        end
    end
    countt=countt+1;
end
end

