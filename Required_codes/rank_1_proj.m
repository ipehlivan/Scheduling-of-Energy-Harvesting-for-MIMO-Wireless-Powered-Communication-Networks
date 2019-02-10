function [V_1_u,t_1,lambda] = rank_1_proj(V_1,Hi,Pa,N,L,K,E,endo)
%--------------------------------------------------------------------------
% Comments refer to paper: 

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO 
%Wireless Powered Communication Networks," in IEEE Communications Letters, 
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function is projecting given matrix to rank 1 matrix: Algorithm 1 
% (Line 8)

%K: # antennas
%N: # users
%Pa: maximum transmit power
%L: # RF chain
%Hi: Channel matrix
%E: required energy
%V_1: Feasible matrix
%endo: if endo is 0, the algorithm finds suitable lambda.

%% Initialization 
[V,D] = eig(V_1); % eigenvalue decompostion
%% Projection
Va=sum((V*abs(D)),2)/trace(abs(D)); % summation of eigenvectors
Va=hybrid_beam_gen(Va,L); % reconstruct vector to satisfy eqn.(2d)
V_1_u=Va*trace(abs(D))*Va'; % Reconstruct rank 1 matrix
V_1_u=(V_1_u+V_1_u')/2; % Getting rid of numerical errors
%% Scaling
% This part solves the problem in eqn.(2a) by scaling the generated matrix.
  cvx_begin quiet
  cvx_solver Sedumi
  cvx_precision best 
  variable alfa   % scaling variable
  variable t_1
  minimize t_1
  
  subject to
  t_1==real(t_1);% eqn.(2g)
  alfa==real(alfa);
  alfa>=0;
  alfa*trace(V_1_u(:,:))<=Pa*(t_1);% eqn.(2c)
  t_1>=0;

  for j=1:N %for every user
      
      Hi_e=Hi(1:K,1:K,j);
      alfa*real(trace(Hi_e*V_1_u))>=E(j);% eqn.(2b)
  end
  cvx_end
  
    V_1_u=alfa*V_1_u; % scaling rank 1 matrix with alfa.assigning output.
  %% Finding suitable lambda
  
  % This part finds suitable lambda for the algorithm to continue searching
  % around the projected matrix. Cvx may fail for high lambda values. This 
  % part  ensures cvx will not fail and lambda is the smallest lambda that
  % rank of the matrix is preserved.
  

  lambda=1; % lambda initialization
  if(endo==0) % if endo is 0.
      for k=1:19
          %% Solve the problem
          [ cvx_opal, V_mem,tem,cvx_s ] = am_algorithm_cvx_single( Hi,  K, N ,Pa ,L,E,1,V_1_u,2^(k-1),t_1);
          ei_V_mem(1:3)=flip(sort(eigs(V_mem(:,:),3))); % Calculate error
          str = cvx_s; % cvx status
          lambda=2^(k-1);
          %Check rank preservation and cvx failure.
          if(max(ei_V_mem(2:3))<=(10^-6) || (prod("Failed"==str(1:6))))
              if(prod("Failed"==str(1:6)))
                  lambda=2^(k-2);
              else
                  lambda=2^(k-1); % Update output variables
                  V_1_u=V_mem;
                  t_1=tem;
              end
              break;
          end
      end
  end
end

