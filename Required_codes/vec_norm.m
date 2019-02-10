function [out] = vec_norm(vector,dir)
% This is an implementation of Matlab's vecnorm function for lower versions 
% of  Matlab.
row_vec=size(vector,2);

coloum_vec=size(vector,1);
if(dir==1)
    out=zeros(1,size(vector,2));
    for k=1:row_vec
        out(k)=norm(vector(:,k));
    end
else
    out=zeros(size(vector,1),1);
    for k=1:coloum_vec
        out(k)=norm(vector(k,:));
    end
end
end

