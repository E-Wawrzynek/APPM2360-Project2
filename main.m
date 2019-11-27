%% APPM 2360 Project 2

%% 2.1 Questions

%% 2
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

%% 3
exposed = [0; 1; 0; 0];

probability_e = transition_SEIR * exposed;

fprintf('\n\n');
fprintf('The probability an exposed individual is in each state after one day is:\n\n'); disp(probability_e);

%% 4
susceptible = [1; 0; 0; 0];

probability_s = transition_SEIR * susceptible;

fprintf('\n\n');
fprintf('The probability a susceptible individual is in each state after one day is:\n\n'); disp(probability_s);

%% 3.1 Questions

%% 1

%(a)
days = 1:1:31;
prob_day = zeros(4, 1);
prob_s = zeros(1, 31);
prob_e = zeros(1, 31);
prob_i = zeros(1, 31);
prob_r = zeros(1, 31);


for n = 1:31
    prob_day = (transition_SEIR ^ n) * susceptible;
    prob_day = prob_day / sum(prob_day);
    
    prob_s(n) = prob_day(1);
    prob_e(n) = prob_day(2);
    prob_i(n) = prob_day(3);
    prob_r(n) = prob_day(4);
end

figure(1);
plot(days, 100*prob_s, '-g');
hold on
plot(days, 100*prob_e, '-r');
plot(days, 100*prob_i, '-m');
plot(days, 100*prob_r, '-b');
legend('susceptible', 'exposed', 'infected', 'recovered');
title('The Probability of being in a SEIR state on each day');
xlabel('Days (1-31)');
xlim([1, 31]);
ylabel('Probability');
    
% (b)
stat_dist = zeros(4, 1);
stat_dist(1) = prob_s(31);
stat_dist(2) = prob_e(31);
stat_dist(3) = prob_i(31);
stat_dist(4) = prob_r(31);

fprintf('\n\n');
fprintf('The stationary distribution is:\n\n'); disp(stat_dist);

%% 2

%(a)
susceptible_2 = [0.15; 0.85; 0; 0];
days_2 = 1:1:31;
prob_day_2 = zeros(4, 1);
prob_s_2 = zeros(1, 31);
prob_e_2 = zeros(1, 31);
prob_i_2 = zeros(1, 31);
prob_r_2 = zeros(1, 31);


for n = 1:31
    prob_day_2 = (transition_SEIR ^ n) * susceptible_2;
    prob_day_2 = prob_day_2 / sum(prob_day_2);
    
    prob_s_2(n) = prob_day_2(1);
    prob_e_2(n) = prob_day_2(2);
    prob_i_2(n) = prob_day_2(3);
    prob_r_2(n) = prob_day_2(4);
end

figure(2);
plot(days, 100*prob_s_2, '-g');
hold on
plot(days, 100*prob_e_2, '-r');
plot(days, 100*prob_i_2, '-m');
plot(days, 100*prob_r_2, '-b');
legend('susceptible', 'exposed', 'infected', 'recovered');
title('The Probability of being in a SEIR state on each day');
xlabel('Days (1-31)');
xlim([1, 31]);
ylabel('Probability');
    
% (b)
stat_dist_2 = zeros(4, 1);
stat_dist_2(1) = prob_s_2(31);
stat_dist_2(2) = prob_e_2(31);
stat_dist_2(3) = prob_i_2(31);
stat_dist_2(4) = prob_r_2(31);

fprintf('\n\n');
fprintf('The stationary distribution is:\n\n'); disp(stat_dist);

%% 4

