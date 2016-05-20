%SCRIPT FOR PLOTTING x VS E[Y|X=x]

%set parameters
psy = 5;
m = 1;
v = 1;
n_samples = 100;

%declare array for plots
x_array = 1:15;
yConditional_array = zeros(1,numel(x_array));
yConditionalError_array = zeros(1,numel(x_array));
rejection_array = zeros(n_samples,numel(x_array));

for i = 1:numel(x_array)
    x = x_array(i);
    [samples,n_reject] = rejection_compoundPoisson(n_samples,x,psy,m,v);
    yConditional_array(i) = mean(samples);
    yConditionalError_array(i) = std(samples);
    rejection_array(:,i) = n_reject;
end
yConditionalError_array = yConditionalError_array/sqrt(n_samples);

figure;
errorbar(x_array,yConditional_array,yConditionalError_array);
xlabel('Conditional value');
ylabel('Conditional expectation');
figure;
boxplot(rejection_array,x_array);
xlabel('Conditional value');
ylabel('Number of rejections per sample');