clc; clear; close all;

%% Reaction Parameters
% Reaction: id = PGI, name = PGI	% Local Parameter:   id =  Keq, name = Keq
	reaction_PGI_Keq=0.36;
	% Local Parameter:   id =  KmF6P, name = KmF6P
	reaction_PGI_KmF6P=0.147;
	% Local Parameter:   id =  KmG6P, name = KmG6P
	reaction_PGI_KmG6P=0.28;
	% Local Parameter:   id =  KmPEP, name = KmPEP
	reaction_PGI_KmPEP=1.999;
	% Local Parameter:   id =  Vmax, name = Vmax
	reaction_PGI_Vmax_default = 2.32456;
	% Local Parameter:   id =  KmPGN, name = KmPGN
	reaction_PGI_KmPGN=0.515958;
    
    % Concentrations
    G6P = 0.861128;
    F6P = 0.261766;
    PEP = 0.997039;
    PGN = 0.131599;
    
%% Problem Definition

% costFunction = @(Vmax) ...
%     abs(MMTwoReactant(Vmax, reaction_PGI_Keq, G6P, reaction_PGI_KmG6P, F6P, reaction_PGI_KmF6P,...
%         PGN, reaction_PGI_KmPGN, 1, PEP, reaction_PGI_KmPEP,1) - 0.1683);     % Cost Function

costFunction = @(K) ...
    abs(MMTwoReactant(reaction_PGI_Vmax_default, reaction_PGI_Keq, G6P, K(1), F6P, K(2),...
        PGN, reaction_PGI_KmPGN, 1, PEP, reaction_PGI_KmPEP,1) - 0.1683);     % Cost Function


nVar =  2;                  % number of variables
VarSize = [ 1 nVar];     % matrix size of decision variables

VarMin = 0;             % Lower bound of decision variables
VarMax = 0.5;              % Ubber bound of decision variables


%% Parameters of PSO

% Constrictions Coeficients
kappa = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2*kappa /abs(2-phi-sqrt(phi^2-4*phi));

MaxIt = 1000;                % Maximum number of iterations
nPop  = 100;                 % Population size or Swarm size

w = chi;                        % The inertia Coeficient - Default
wDamp = 1;             % Damping Coeficient
c1 = chi*phi1;                       % Personal Acceleration Coeficient  
c2 = chi*phi2;                       % Global Acceleration Coeficient  

MaxVelocity = (VarMax - VarMin)*0.2;
MinVelocity = -MaxVelocity;
%% Initialization 

% the particle template.
empty_particle.Position = []; % Position of the particle in the search space
empty_particle.Velocity = [];  % Velocity of the particle
empty_particle.Cost = [];       % Cost function value
empty_particle.Best.Position =[]; % personal best
empty_particle.Best.Cost = [];

particle = repmat(empty_particle, nPop, 1);  % an array of empty particles.

% Initialize the Global Best
GlobalBest.Cost = Inf;


% Initialization of population members
for i=1:nPop
    % generate  a random solution
    % Update the solution
    particle(i).Position = unifrnd(VarMin , VarMax, VarSize);
    % Evaluation
    particle(i).Cost = costFunction(particle(i).Position);
    % Velocity
    particle(i).Velocity = zeros(VarSize);
    % Update the personal best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    % update global best
    if particle(i).Best.Cost < GlobalBest.Cost 
        GlobalBest = particle(i).Best;
    end
    
end

BestCosts = zeros(MaxIt,1); % hold best cost value in each iteration
BestPositions = zeros(MaxIt,nVar); % hold best cost value in each iteration

%% Main Loop

for It = 1:MaxIt
    
    for i=1:nPop
        % Update Velocity
        particle(i).Velocity = w* particle(i).Velocity + ...
            c1*rand(VarSize) .* (particle(i).Best.Position - particle(i).Position) + ...
            c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        % Apply lower bound Upper bound Limits on Velocity
        particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
        particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Apply lower bound Upper bound Limits on Position
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);
        
        % Evaluation
        particle(i).Cost = costFunction(particle(i).Position);
        
        % Update personal best
        if particle(i).Cost < particle(i).Best.Cost
            
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.Position = particle(i).Position;
        end
        
        % update global best
        if particle(i).Best.Cost < GlobalBest.Cost 
            GlobalBest = particle(i).Best;
        end
            
    end
    
    % Store the best cost value
    BestCosts(It) = GlobalBest.Cost;
    BestPositions(It,:) = GlobalBest.Position;
    % Damping Coeficients
    w = wDamp * w;
    % Display Information 
    
    disp(["Iteration: " num2str(It) "Cost Value= " num2str(BestCosts(It)) "Position= " GlobalBest.Position])
end

%% Plot 
figure;
grid on
subplot(2,1,1);
semilogy(BestCosts, 'lineWidth',2)
xlabel("Iteration")
ylabel("best cost")

subplot(2,1,2)
semilogy(BestPositions, 'lineWidth',2)
xlabel("Iteration")
ylabel("Best Position")














