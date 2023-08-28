  It is Monte Carlo (MC) function with New equation instead Of  that we had in Petri's paper.
= 
In this Function the arccos equation of Chebyshev is used () for calculatiog base function.
 Next step was calculating of  and last step was study of Missmatch correction factor (F) from this equation 
function [F, delta] = MCEq2ChebyJuha_3(N, E_BD, u_c, V, S_rel)
wavelengthRegion = 360 : 830;
wavelengths = wavelengthRegion * 1e-9;
wavelengthNumber = length(wavelengthRegion);
Error Function:
lambda_1 = wavelengths(1);
lambda_2 = wavelengths(end);
f = ones(N + 1, wavelengthNumber);
C = sum(1 ./ (1 : N) .^ 2);

J = 2 * N;
T = ones(J, wavelengthNumber);

lambda = wavelengths;
x = (2 * lambda - lambda_1 - lambda_2) ./ (lambda_2 - lambda_1);

for j = 0 : J
    if j==0
        T(j+1,:)=1;
    elseif j==1
        T(j+1,:)= x;
    else 
        T(j+1,:)=2.*x.*T(j,:)-T(j-1,:);
    end
        
end
G = [ones(1, wavelengthNumber) ; T(2:J+1,:) ./ std(T(2:J+1,:), 0, 2)];
% G =  T(2:J+1,:) ./ std(T(2:J+1,:), 0, 2);
theta = zeros(1, N);
% R = ones(1, N + 1);
for i = 1 : N
    theta(i) = random('Uniform', 0, 2 * pi);
    theta_U(i) = random('Uniform', 0, 2 * pi);
    f(i + 1, :) =  (G(2 * i + 1, :) * cos(theta(i)) + G(2 * i, :) * sin(theta(i)));
    P = (cos(theta_U(i))./sqrt(cos(theta_U(i))^2+sin(theta_U(i))^2.*C));
    KJ = (sin(theta_U(i))./sqrt(cos(theta_U(i))^2+sin(theta_U(i))^2.*C));
    R(i+1) = P + (KJ / i);  
end

delta =  sum (R'.* f);
Modify Spectral Irradiance
 
Phi_e = (1 + delta .* u_c) .* E_BD; 
Spectral Mismatch Correction Factor

F = trapz(Phi_e .* V) ./ trapz(Phi_e .* S_rel);

end

