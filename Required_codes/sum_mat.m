function [y] = sum_mat(H,V,N)
%--------------------------------------------------------------------------
% Comments refer to paper:

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO
%Wireless Powered Communication Networks," in IEEE Communications Letters,
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function calculates sum of traces at % eqn.(2c).

% N: # time slots. Do not confuse with # users.
% H: The channel matrix
%% Sum of traces
y=0;
for k=1:N
    y=y+trace(H*V(:,:,k));
end
end

