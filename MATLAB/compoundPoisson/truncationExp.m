%SCRIPT FOR PLOTTING x VS E[Y|X=x]

%set parameters
psy = 5;
m = 1;
v = 0.01;

%declare array for plots
x_array = (-10):0.01:10;
yConditional_array = zeros(1,numel(x_array));

for i = 1:numel(x_array)
    x = x_array(i);
    yConditional_array(i) = expectation_compoundPoisson(100,x,m,v,psy);
end

figure;
plot(x_array,yConditional_array);
hold on;
xlabel('Conditional value');
ylabel('Conditional expectation');
scatter(0,0,'b','filled');
scatter(0,yConditional_array(x_array==0),'b');