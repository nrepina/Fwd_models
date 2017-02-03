function [U_hat, U_k, alpha, beta] = Nestrov_forward_farfield(Vs,U_0,U_in,U_in_air,z_plane,iteration,tol)

global M N O fGreen E prop_wavenumber

F = @(x) fftn(x);
IF = @(x) ifftn(x);

U_k = zeros(M,N,O,iteration);
alpha = zeros(iteration,1);
beta = zeros(iteration,1);
t_k = 1;
k_count = 0;

error = zeros(iteration,1);
source = zeros(size(fGreen));
step_cal = 1; stepsize_pre = 0;

for Iter = 1:iteration
    
    source(:) = 0;
    if Iter == 1 
        source(1:N,1:M,1:O) = U_0.*Vs;
    else
        source(1:N,1:M,1:O) = U_k(:,:,:,k_count).*Vs;
    end
    
    source = IF(fGreen.*F(source));
    if Iter == 1
        res = U_0 - source(1:N,1:M,1:O) - U_in;
    else
        res = U_k(:,:,:,k_count) - source(1:N,1:M,1:O) - U_in;
    end
    error(Iter) = 0.5*norm(res(:))^2;
    
    if Iter > 1
        if error(Iter) > error(Iter-1)
            beta(k_count) = 0;
            U_k(:,:,:,k_count) = x_k;
            t_k = 1; 
            step_cal = 1;
%             disp('restart!');
            continue;
        elseif error(Iter) < tol*norm(U_in(:))^2
            break;
        end
    end
    
    k_count = k_count + 1;
    
    source(:) = 0;
    source(1:N,1:M,1:O) = res;
    source = IF(F(source).*conj(fGreen));
    grad_U = res - conj(Vs).*source(1:N,1:M,1:O);    

    if step_cal == 1
        source(:) = 0;
        source(1:N,1:M,1:O) = grad_U.*Vs;
        source = IF(fGreen.*F(source));
        A_grad_U = grad_U - source(1:N,1:M,1:O);
        alpha(k_count) = norm(grad_U(:))^2/norm(A_grad_U(:))^2;
        clear A_grad_U
        if abs(stepsize_pre-alpha(k_count))/alpha(k_count)<0.01
            step_cal = 0;
        end
    else
        alpha(k_count) = stepsize_pre;
    end
    
    % x_k update
    if Iter == 1
        x_k1 = U_0 - alpha(k_count)*grad_U;
    else
        x_k1 = U_k(:,:,:,k_count-1) - alpha(k_count)*grad_U;
    end
    stepsize_pre = alpha(k_count);
    
    % U_k update
    t_k1 = 0.5*(1+sqrt(1+4*t_k^2)); 
    beta(k_count) = (t_k-1)/t_k1;
    if beta(k_count) == 0
        U_k(:,:,:,k_count) = x_k1;
    else
        U_k(:,:,:,k_count) = x_k1 + beta(k_count)*(x_k1-x_k);
    end
    
    x_k = x_k1;
    t_k = t_k1;
    
%     hold on;
%     plot(1:Iter,log10(error(1:Iter)),['o','b']);axis square;drawnow;
%     disp(['forward iteration : ' num2str(Iter) '/' num2str(iteration)]);
end

U_k = U_k(:,:,:,1:k_count);
beta = beta(1:k_count);
alpha = alpha(1:k_count);

% farfield propagation
U_hat = zeros(M,N,length(z_plane));
fU_temp = zeros(2*M,2*N,O);
fU_temp(1:M,1:N,:) = U_k(:,:,:,end).*Vs;
fU_temp = sum(E.*fft2(fU_temp),3);

for zIdx = 1:length(z_plane)
     temp = ifft2(fU_temp.*exp(z_plane(zIdx)*prop_wavenumber));
     U_hat(:,:,zIdx) = temp(1:M,1:N);
end
U_hat = U_in_air + U_hat;
disp(['forward layer : ' num2str(k_count)]);

end