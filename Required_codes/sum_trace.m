function [output] = sum_trace(V_1,V_2 ,ksi)
%--------------------------------------------------------------------------
% Comments refer to paper:

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO
%Wireless Powered Communication Networks," in IEEE Communications Letters,
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function calculates the penalty function at eqn.(3a).

%ksi: penalty function weights

%% Penalty function calculation
output=0;
for k=1:size(V_1,3)
output=output+ksi(k)*((trace(V_1(:,:,k))*trace(V_2(:,:,k)))-trace(V_1(:,:,k)*V_2(:,:,k)));
end
end

