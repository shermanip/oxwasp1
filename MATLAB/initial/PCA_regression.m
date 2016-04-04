%FUNCTION: PCA REGRESSION
    %Does PCA regression for 2 variables
%PARAMETERS:
    %x: vector of observations of one variable
    %y: vector of observations of another variable
%RETURN:
    %r_com: center of mass
    %u_1: gradient in vector form
    %y_error: root mean squared error in the y component
function [r_com,u_1,y_error] = PCA_regression(x,y)

    %get the mean and std and sample size of the data
    x_bar = mean(x);
    y_bar = mean(y);
    r_com = [x_bar;y_bar];
    std_x = std(x,1);
    std_y = std(y,1);
    n = numel(x);

    %normalize the data
    x = (x - x_bar)/std_x;
    y = (y - y_bar)/std_y;

    %estimate the covariance matrix
    X = [x,y];
    covariance = X'*X/n;
    %get the eigenvector with the largest eigenvalue
    [eigenvectors,eigenvalues] = eig(covariance);
    eigenvalues = diag(eigenvalues); %covert diagonal matrix into vector
    u_1 = eigenvectors(:,max(eigenvalues)==eigenvalues); %select eigenvector
    
    %if want y_error
    if nargout == 3
        u_2 = eigenvectors(:,min(eigenvalues)==eigenvalues); %select eigenvector
        %get perpendicular distance of each point to the line
        d = X*u_2;
        %get the root mean squared error of the perpendicular distance
        d = sqrt(mean(d.^2));
        %un-normalize the y component of the noise in the 2nd PC
        y_error = abs(d*u_2(2)*std_y);
    end

    %un-normalize the eigenvector
    u_1 = u_1.*[std_x;std_y];
   
end

