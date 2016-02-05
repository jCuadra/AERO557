clear all
xoi = [1 2 3 4 5 6 7 8 9 10]';
yoi = [1 2 2 3 4 4 5 7 9 9]';
H = [ones(length(xoi),1) xoi]
c = inv(H'*H)*H'*yoi;
alpha = c(1);
beta = c(2);
Pk = inv(H'*H);
fprintf('Using all the data points we find the equation of the line to be\n');
fprintf('y = %f + %fx with a confidence of %f %f\n', alpha, beta, Pk(1, 1), Pk(2, 2));
%% Editing Data
E = yoi-(alpha + beta*xoi);
twoSig = 2*rms(E);
for i = 1:length(xoi)
    if E(i) > twoSig
        ind(i) = 0;
    else
        ind(i) = 1;
    end
end
xoi = ind.*xoi
yoi = ind.*yoi
H = [ones(length(xoi),1) xoi]
c = inv(H'*H)*H'*yoi;
alpha = c(1);
beta = c(2);
Pk = inv(H'*H);
fprintf('Using all the data points we find the equation of the line to be\n');
fprintf('y = %f + %fx with a confidence of %f %f\n', alpha, beta, Pk(1, 1), Pk(2, 2));