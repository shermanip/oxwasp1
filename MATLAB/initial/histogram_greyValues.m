clc;
clear all;
close all;

%load the data
X = load_stack('/data/tinamou/sip/block_images/orginial');
%reshape it in a column vector
X = reshape(X,numel(X),[]);

%plot the histogram
figure;
histogram(X,'Normalization','countdensity');
xlabel('Grey values (arb. unit)');
ylabel('Frequency density {(arb. unit^{-1})}');