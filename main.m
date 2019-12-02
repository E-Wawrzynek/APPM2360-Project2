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

probability_s = transition_SEIR ^ 5 * susceptible;

fprintf('\n\n');
fprintf('The probability a susceptible individual is in state R after 5 days is:\n\n'); disp(probability_s(4));

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

% lots of math: x_inf = c1v1... working on finding v and c

[V, d] = eig(transition_SEIR);

lambdas = [d(1); d(6); d(11); d(16)];
eig_vec_i = V(1:4, 3);

x0_1 = [1; 0; 0; 0];
x0_2 = [0.15; 0.85; 0; 0];

c_1 = V \ x0_1;
c_2 = V \ x0_2;

c_i = c_1(3);
c_i_2 = c_2(3);

disp('Checking that the constant c is the same for initial conditions:');
disp('c_i:'); disp(c_i); disp('c_i_2'); disp(c_i_2);

%eig_vec_i = eig_vec_i / sum(eig_vec_i);

x_inf = c_i * eig_vec_i;
disp(x_inf);

abs_error = zeros(4, 31);
days_4 = 1:1:31;
prob_day_4 = zeros(4, 1);
abs_err_s = zeros(1, 31);
abs_err_e = zeros(1, 31);
abs_err_i = zeros(1, 31);
abs_err_r = zeros(1, 31);


for n = 1:31
    prob_day_4 = (transition_SEIR ^ n) * x0_2;
    prob_day_4 = prob_day_4 / sum(prob_day_4);
    
    abs_error(1:4, n) = abs(x_inf - prob_day_4);
    
    abs_err_s(n) = abs_error(1, n);
    abs_err_e(n) = abs_error(2, n);
    abs_err_i(n) = abs_error(3, n);
    abs_err_r(n) = abs_error(4, n);
end

figure(3);
semilogy(days_4, 100*abs_err_s);
hold on
semilogy(days_4, 100*abs_err_e);
semilogy(days_4, 100*abs_err_i);
semilogy(days_4, 100*abs_err_r);
%% 5

% (b)
Pim = 1;

transition_SEIR_Im = [Ps, Pe, 0, (1/2)*(1-Pr), 0;
                   1-Ps, 0, 0, 0, 0;
                   0, 1/2*(1-Pe), (1-Pi), 0, 0;
                   0, 1/2*(1-Pe), Pi, Pr, 0;
                   0, 0, 0, (1/2)*(1-Pr), Pim];
               
fprintf('\n\n');
fprintf('The transistion matrix for the Markov Chain of the SEIR-Im model:\n\n'); disp(transition_SEIR_Im); 

% (c)
%% 1

%(a)
susceptible_new = [1; 0; 0; 0; 0];
days_3 = 1:1:250;
prob_day_3 = zeros(5, 1);
prob_s_3 = zeros(1, 250);
prob_e_3 = zeros(1, 250);
prob_i_3 = zeros(1, 250);
prob_r_3 = zeros(1, 250);
prob_im = zeros(1, 250);


for n = 1:250
    prob_day_3 = (transition_SEIR_Im ^ n) * susceptible_new;
    prob_day_3 = prob_day_3 / sum(prob_day_3);
    
    prob_s_3(n) = prob_day_3(1);
    prob_e_3(n) = prob_day_3(2);
    prob_i_3(n) = prob_day_3(3);
    prob_r_3(n) = prob_day_3(4);
    prob_im(n) = prob_day_3(5);
end

figure(4);
plot(days_3, 100*prob_s_3, '-g');
hold on
plot(days_3, 100*prob_e_3, '-r');
plot(days_3, 100*prob_i_3, '-m');
plot(days_3, 100*prob_r_3, '-b');
plot(days_3, 100*prob_im, '-c');
legend('susceptible', 'exposed', 'infected', 'recovered', 'immune');
title('The Probability of being in a SEIR-Im state on each day');
xlabel('Days (1-250)');
%xlim([1, 250]);
ylabel('Probability');
    
% (b)
stat_dist_3 = zeros(5, 1);
stat_dist_3(1) = prob_s_3(250);
stat_dist_3(2) = prob_e_3(250);
stat_dist_3(3) = prob_i_3(250);
stat_dist_3(4) = prob_r_3(250);
stat_dist_3(5) = prob_im(250);

fprintf('\n\n');
fprintf('The stationary distribution is:\n\n'); disp(stat_dist_3);

% commenty comment