function b = AB_multi(p,epsilon,Nr, Nt, dr, theta0)
e = ones(Nt,1);

p = reshape(p,Nr+1,Nr+1,Nt);
Rp_x = zeros(size(p));
Rp_y = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kt = 1:Nt
    p_temp = p(:,:,kt);
    
    Rp_pre_x = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
    Rp_pre_y = (p_temp(2:end,:) - p_temp(1:end-1,:))/dr;
    
    sinT = sin(theta0(kt));
    cosT = cos(theta0(kt));
    if sinT<=0
        Rp_y(1:end-1,:,kt) = sinT*Rp_pre_y;
        boundary_D(end,:,kt) = p_temp(end,:);
    else
        Rp_y(2:end,:,kt) = sinT*Rp_pre_y;
        boundary_D(1,:,kt) = p_temp(1,:);
    end
    if cosT>=0
        Rp_x(:,2:end,kt) = cosT*Rp_pre_x;
        boundary_D(:,1,kt) = p_temp(:,1);
    else
        Rp_x(:,1:end-1,kt) = cosT*Rp_pre_x;
        boundary_D(:,end,kt) = p_temp(:,end);
    end
end

for kx = 2:Nr
    for ky = 2:Nr
        p_temp = p(kx,ky,:); p_temp = p_temp(:);
        p_temp = p_temp - e*e'*p_temp/Nt;
        BCp(kx,ky,:) = p_temp;
    end
end

b = epsilon*Rp_x + epsilon* Rp_y + BCp + boundary_D;

b = b(:);

end