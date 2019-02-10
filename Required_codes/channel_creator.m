function [Hi] = channel_creator(K, N ,RR,d,dalgalamba,rician_K_factor,scaler, Ro)
%--------------------------------------------------------------------------
% Comments refer to paper: 

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO 
%Wireless Powered Communication Networks," in IEEE Communications Letters, 
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function creates channel matrices.

%K: # antennas
%N: # users
%RR:# node(receiver) antennas
%d: antenna distance matrix
%dalgalamba: wavelength
%rician_K_factor: Rician K factor
%scaler: scaling factor of eqn.(2c). Eliminates cvx numerical errors.
%Ro: energy conversion efficiency

%% Channel creation

% Millimeter Wave Channel Modeling

power_gain=zeros(K,RR,N);
for nn=1:N
    power_gain(:,:,nn)=-((20*log10(4*pi*1/dalgalamba))+10*1.2*log10(d(:,:,nn))+normrnd(0,1.7,K,RR));
end
power_gain=db2pow(real(power_gain));
Hi_res=(randn(K,RR,N)+randn(K,RR,N)*1i).*(sqrt(power_gain/2))+rician_K_factor*power_gain; % Rayleigh channel coefficents

%% Channel matrix generation

Hi=zeros((K),(K),N); % For speed and safety
for jj=1:N %for every user
    Himat=zeros(K,K);
    Himat(1:K,1:K)=(Hi_res(:,:,jj)*Hi_res(:,:,jj)');
    Himat=scaler*(Himat'+Himat)/2; % Getting rid of numerical errors
    % We multiply both sides of the eqn.(2c) with scaler to avoid cvx
    % numerical errors.
    Hi(:,:,jj)=Ro*Himat;
end

end

