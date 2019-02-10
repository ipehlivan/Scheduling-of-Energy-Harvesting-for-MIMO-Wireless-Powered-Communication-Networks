function [Y] = hybrid_beam_gen(X,num_chain)
%--------------------------------------------------------------------------
% Comments refer to paper:

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO
%Wireless Powered Communication Networks," in IEEE Communications Letters,
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

% This function generates a vector that satisfies eqn.(2d).

% num_chain: # RF chain
%%
X=X(:); % ensuring its a column vector.
X_mem=X; % memory
anpc=size(X)/num_chain; % antenna per chain
X=X./(abs(X)); % normalization
for s=1:num_chain % beamforming norms
    X((1+(s-1)*anpc):((s)*anpc))=X((1+(s-1)*anpc):((s)*anpc))*...
        mean(abs(X_mem((1+(s-1)*anpc):((s)*anpc))));
    % taking mean of weights of antennas that are assigned to same RF chain.
end
Y=X; % out
end

