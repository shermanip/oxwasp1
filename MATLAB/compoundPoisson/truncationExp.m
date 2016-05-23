%SCRIPT FOR PLOTTING x VS E[Y|X=x]

%set parameters
psy = 5; %rate parameter
m = 3; %mean parameter
v = 0.1; %variance parameter

%declare array for plots
x_array = (-10):0.01:10;
yConditional_array = zeros(1,numel(x_array));

%for each x condition
for i = 1:numel(x_array)
    x = x_array(i); %get the x
    %work out the conditional expectation
    yConditional_array(i) = expectation_compoundPoisson(1000,x,m,v,psy);
end

%plot the conditional expectation
figure;
plot(x_array,yConditional_array);
%plot the discountinity at x=0
hold on;
scatter(0,0,'b','filled');
hold off;
%label the axis
xlabel('Conditional value');
ylabel('Conditional expectation');
