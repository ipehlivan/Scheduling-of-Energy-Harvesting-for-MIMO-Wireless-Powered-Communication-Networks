function [d] = topology_creator_3D(K,N,RR,radiu,spacing,w_case)
%--------------------------------------------------------------------------
% Comments refer to paper:

%I. Pehlivan and S. C. Ergen, "Scheduling of Energy Harvesting for MIMO
%Wireless Powered Communication Networks," in IEEE Communications Letters,
%vol. 23, no. 1, pp. 152-155, Jan. 2019.

%doi: 10.1109/LCOMM.2018.2881180

%--------------------------------------------------------------------------

%This function creates a random node distribution and calculates the antenna
%distance matrix.

%K: # antennas
%N: # users
%RR:# node(receiver) antennas
%radiu: Radius of the sphere
%spacing: spacing between antennas
%w_case: distributes on the sphere for 1 and in the sphere for 0.

%% Node distribution
absol=rand(1,N)*radiu; % generating radius of each user;
if(w_case==1)% On the sphere for 1.
    absol=ones(1,N)*radiu; % generating radius of each user;
end
azimuth=rand(1,N)*2*pi; % generating azimuth for each user
elevation=rand(1,N)*pi; % generating elavation for each user
if(K>1)% if base station has multiple antennas
    antendis_x=[([1:(K/2)]'-((K/2)+1)/2)*(spacing);zeros(1,K/2)']; % antenna distribution
    antendis_y=[zeros(1,K/2)';([1:(K/2)]'-((K/2)+1)/2)*(spacing)]; % antenna distribution
    antendis_z=zeros(K,1);
else
    antendis_x=0; % antenna distribution
    antendis_y=0; % antenna distribution
    antendis_z=0; % antenna distribution
end
[xu, yu,zu]=sph2cart(azimuth,elevation,absol); % polar to carthesian
userdis_x=zeros(RR,N); % user dist x
userdis_y=zeros(RR,N); % user dist y
userdis_z=zeros(RR,N); % user dist y
if(RR>1) % if nodes have multiple antennas
    for nn=1:N
        userdis_x(:,nn)=[([1:(RR/2)]'-((RR/2)+1)/2)*(spacing)+xu(nn);zeros(1,RR/2)'+xu(nn)];
        userdis_y(:,nn)=[zeros(1,RR/2)'+(yu(nn));([1:(RR/2)]'-((RR/2)+1)/2)*(spacing)+(yu(nn))];
        userdis_z(:,nn)=ones(RR,1)*zu(nn);
    end
else
    
    for nn=1:N
        userdis_x(:,nn)=ones(RR,1)*xu(nn);
        userdis_y(:,nn)=ones(RR,1)*yu(nn);
        userdis_z(:,nn)=ones(RR,1)*zu(nn);
    end
    
end
%% Scatter plot for checking
% You can observe your topology by uncommenting this line.

% for nn=1:N
%     scatter3(userdis_x(:,nn),userdis_y(:,nn),userdis_z(:,nn),'b');
%     hold on
%     scatter3(xu(nn),yu(nn),zu(nn),'r')
%     hold on
%     scatter3(antendis_x,antendis_y,antendis_z,'k')
% end
%% Distance matrix calculation

d=zeros(K,RR,N); % declaration for speed up
for nn=1:N
    dummy(:,:)=[userdis_x(:,nn),userdis_y(:,nn),userdis_z(:,nn)];
    for rr=1:RR
        dis_vec(:,1:3)=[ones(K,1)*dummy(rr,1)-antendis_x,ones(K,1)*dummy(rr,2)-antendis_y,ones(K,1)*dummy(rr,3)-antendis_z];
        d(:,rr,nn)=vec_norm(dis_vec,2);
        % check this function for comments.
    end
end

end

