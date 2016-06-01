clc;
clear all;
close all;

X = load_stack('/data/tinamou/sip/block_images/orginial');
X = reshape(X,numel(X),[]);

figure;
histogram(X,'Normalization','countdensity');
xlabel('Grey values (arb. unit)');
ylabel('Frequency density {(arb. unit^{-1})}');