clear
clc
close all

l1 = 1.5;
l2 = .1;
sig = 10;

rat1 = [];
rat2 = [];
for t = 3 : 100
    num1 = l2^(2*t) + sig^2 * ((l2^(2*(t-1)) - 1)/ (l2^2 - 1));
    den1 = l1^(2*t) + sig^2 * ((l1^(2*(t-1)) - 1)/ (l1^2 - 1));
    rat1 = [rat1, num1 / den1];
    
    num2 = l2^(2*t);
    den2 = l1^(2*t);
    rat2 = [rat2, num2 / den2];
end

% l1 = 1.5;
% l2 = 1.1;
% sig = 1;
% t = 11;
% num = l2^(2*t) + sig^2 * ((l2^(2*(t-1)) - 1)/ (l2^2 - 1));
% den = l1^(2*t) + sig^2 * ((l1^(2*(t-1)) - 1)/ (l1^2 - 1));

figure;
semilogy(rat1, 'rs')
hold
semilogy(rat2, 'kd')