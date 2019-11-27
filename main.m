%% APPM 2360 Project 2

%% 2.1 Questions

%2
Ps = 0.7;
Pe = 0.4;
Pi = 1;
Pr = 0.8;

transition_SEIR = [Ps, Pe, 0, 1-Pr;
                   1-Ps, 0, 0, 0;
                   0, 1/2*(1-Pe), (1-Pi), 0;
                   0, 1/2*(1-Pe), Pi, Pr];

fprintf('\n\n');
fprintf('The transistion matrix for the Markov Chain of the SEIR model:\n\n'); disp(transition_SEIR);

%3
exposed = [0; 1; 0; 0];

probability_e = transition_SEIR * exposed;

fprintf('\n\n');
fprintf('The probability an exposed individual is in each state after one day is:\n\n'); disp(probability_e);

%4
susceptible = [1; 0; 0; 0];

probability_s = transition_SEIR * susceptible;

fprintf('\n\n');
fprintf('The probability a susceptible individual is in each state after one day is:\n\n'); disp(probability_s);

%% 3.1 Questions
               
