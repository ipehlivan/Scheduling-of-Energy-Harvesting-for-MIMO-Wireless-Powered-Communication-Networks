function [ V,t ,cvx_optval,cvx_correct ] = Main_solver_schedule_am( K, N ,Pa ,L,Hi ,E,S)
%--------------------------------------------------------------------------
% Comments refer to paper: 

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO 
%Wireless Powered Communication Networks," in IEEE Communications Letters, 
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------
%This function solves semidefinite relaxation of the problem (2).

%K: # antennas
%N: # users
%Pa: maximum power
%L: # RF chain
%Hi: Channel matrix
%E: required energy
%S: number of RF chain
%% Function starts
ant_pr_chain=K/L; % Antenna per chain

ck=1; % cvx failure flag

%% Scheduling case
% This part solves the problem with cvx. Constraints are indexed with 
% numbers at the article.

if(S>1) 
   cvx_begin quiet
    cvx_solver SeDuMi
    cvx_precision best
    variable V(K,K,S) hermitian
    variable t(S)
    minimize sum(t); % eqn.(2a)
    
    subject to
    
    t==real(t);% eqn.(2g)
    for kkk=1:S % for every slot
        
        V(:,:,kkk)==hermitian_semidefinite(K); % eqn.(2f)
        if(kkk>=2)
            % ordering for implementation reasons. Check zero matrix
            % cancelation at am scheduling algorithm.
            trace(V(:,:,kkk))<=trace(V(:,:,kkk-1));
        end
        for ga=1:L % for every RF chain
            for g=1:ant_pr_chain-1
                IA_e=zeros(K,K,S);
                IA=zeros(K);
                IA(ant_pr_chain*(ga-1)+1,ant_pr_chain*(ga-1)+1)=1;
                IA(ant_pr_chain*(ga-1)+1+(g),ant_pr_chain*(ga-1)+1+(g))=-1;
                trace(IA*V(:,:,kkk))==0; % eqn.(2d)
            end
        end
        
        trace(V(:,:,kkk))<=Pa*(t(kkk));% eqn.(2c)
        t(kkk)>=0;% eqn.(2g)
    end
    
    for j=1:N %for every user
        
        Hi_e=Hi(1:K,1:K,j);
        sum_mat(Hi_e,V,S)>=E(j);% eqn.(2b)
    end
    cvx_end
end
%% Non-scheduling case
if(S==1)
    cvx_begin quiet
    cvx_solver Sedumi
     cvx_precision best
    variable V(K,K) hermitian
    variable t
    minimize t; % eqn.(2a)
    subject to
    
    V==hermitian_semidefinite(K);% eqn.(2f)
    
    for j=1:N %for every user
        Hi_e=Hi(1:K,1:K,j);
        trace(Hi_e*V)>=E(j); % eqn.(2b)
    end
    
    
    for ga=1:L % for every RF chain
        for g=1:ant_pr_chain-1
            IA=zeros(K+1);
            IA(ant_pr_chain*(ga-1)+1,ant_pr_chain*(ga-1)+1)=1;
            IA(ant_pr_chain*(ga-1)+1+(g),ant_pr_chain*(ga-1)+1+(g))=-1;
            trace(IA(1:K,1:K)*V)==0; % eqn.(2d)
        end
    end
    trace(V)<=Pa*t; % eqn.(2c)
    t>=0; % eqn.(2g)
    cvx_end
end
%% Cvx solver precaution
% Collecting error information
   str = cvx_status;
    if(prod('Failed'==str(1:6)))
        ck=0;
    end
cvx_correct=ck; % cvx error flag
end

