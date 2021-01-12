function X = MonteCarlo_Ising(J, h, T, L)
% Function to generate Monte Carlo samples of kinetic Ising model with
% possibly asymmetric interactions and synchronous updates

% Number of samples and number of burn-in steps:
% L = 1000000;
L_burn = 10000;

% Model parameters:
N = size(J,1);
% N = 100;
% T = 1;
% 
% J = normrnd(0, 1/sqrt(N), N, N);
% h = normrnd(0, 1, N, 1);

% Samples:
X = zeros(N, L);

% Initial configuration:
x = 2*double(rand(N,1) > .5) - 1; % Random

% Burn-in phase:
for i = 1:L_burn
    
    % Local fields:
    Hs = J*x + h;
    
    % Probabilities of each spin going "up":
    P_up = (1 + exp(-(2/T)*Hs)).^(-1);
    
    % New spin configuration:
    x = 2*double(rand(N,1) < P_up) - 1;
    
end

% Collection phase:
for i = 1:L
    
    % Local fields:
    Hs = J*x + h;
    
    % Probabilities of each spin going "up":
    P_up = (1 + exp(-(2/T)*Hs)).^(-1);
    
    % New spin configuration:
    x = 2*double(rand(N,1) < P_up) - 1;
    
    % Save:
    X(:,i) = x;

end

