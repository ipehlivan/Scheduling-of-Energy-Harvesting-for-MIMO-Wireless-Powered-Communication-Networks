% SCA ALGORITHM
function [ cvx_optval, V_1,t_1,cvx_status ] = am_algorithm_cvx( Hi,  K, N ,Pa ,L,E,S,V_1_u,lambda,t_mem )
%--------------------------------------------------------------------------
% Comments refer to paper: 

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO 
%Wireless Powered Communication Networks," in IEEE Communications Letters, 
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function solves the scheduled case of the optimization problem at
%eqn.(3a).

%K: # antennas
%N: # users
%Pa: maximum transmit power
%L: # RF chain
%Hi: Channel matrix
%E: required energy
%S: number of RF chain
%lambda: weights of the penalty function
%V_1_u: Memory optimal matrices
%t_mem: Memory optimal delays
%% Function starts

ant_pr_chain=K/L; % Antenna per chain
%% CVX starts
% This part solves the problem with cvx. Constraints are indexed with 
% numbers at the article.

    cvx_begin quiet
    cvx_solver Sedumi
    cvx_precision best
    variable V_1(K,K,S) hermitian
    variable t_1(S)
    minimize sum_trace(V_1_u,V_1,lambda)+sum(t_1)+sum(t_mem); % eqn.(3a)
    
    subject to
    
    t_1==real(t_1);% eqn.(2g)
    for kkk=1:S % for every slot
        
        V_1(:,:,kkk)==hermitian_semidefinite(K);% eqn.(2f)
        if(kkk>=2)
            % ordering for implementation reasons. Check zero matrix
            % cancelation at am scheduling algorithm.
            trace(V_1(:,:,kkk))<=trace(V_1(:,:,kkk-1));
        end
        for ga=1:L% for every RF chain
            for g=1:ant_pr_chain-1
                IA_e=zeros(K,K,S);
                IA=zeros(K);
                IA(ant_pr_chain*(ga-1)+1,ant_pr_chain*(ga-1)+1)=1;
                IA(ant_pr_chain*(ga-1)+1+(g),ant_pr_chain*(ga-1)+1+(g))=-1;
                trace(IA*V_1(:,:,kkk))==0;% eqn.(2d)
            end
        end
        
        trace(V_1(:,:,kkk))<=Pa*(t_1(kkk));% eqn.(2c)
        t_1(kkk)>=0;% eqn.(2g)
    end
    
    for j=1:N %for every user
        
        Hi_e=Hi(1:K,1:K,j);
        sum_mat(Hi_e,V_1,S)>=E(j);% eqn.(2b)
    end
    cvx_end
    

end