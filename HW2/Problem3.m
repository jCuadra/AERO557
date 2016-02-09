%% Problem 3
clear all
%% Inputs
r = [5748.60114991; 2679.72874263; 3443.00728717];
v = [4.33046211; -1.92286233; -5.72656403];

atwa = [216.7346 327.1327 42.7977 35411.4705 55076.9922 4684.2902;
        327.1327 785.4208 87.5320 51191.8404 128003.7108 8758.9867;
        42.7977 87.5320 62.3664 5791.8322 11727.0485 9613.1211;
        35411.4705 51191.8404 5791.8322 5929294.3651 8844731.2848 597965.6811;
        55076.9922 128003.7108 11727.0485 8844731.2848 21446311.1547 1038893.8970;
        4684.2902 8758.9867 9613.1211 597965.6811 1038893.8970 1555196.0843];
atwb = [19.8969; -18.1124; -0.2737; 3815.5860; -3395.8956; 69.3738];
%% NOTE
%These numbers make the solution blow up becuase they are from the first
%iteration of Example 10-4 so I grabbed the fiinal P and delta xHat so that 
%I could back calculate the atwa and atwb but obviously there is rounding
%error becuase he only carries so many sig figs in what he records in his
%book.

% htwh = [136.7 240.2 73.7 9290.2 16859.2 3981.9;
%         240.2 1006.2 300.1 15580.9 56253.9 11843.9;
%         73.7 300.1 145.4 3906.1 12469.0 5665.1;
%         9290.2 15580.9 3906.1 767027.6 1322814.7 278780.4;
%         16859.2 56253.9 12469.0 1322814.7 4368107.9 764948.9;
%         3981.9 11843.9 5665.1 278780.4 764948.9 365885.3];
% htwb =[ 38634.2; 215383.4; 79253.9; 2260058.6; 10504371.5; 2995933.5];

%% Update State
%Taken from example 10-4
Pold = [.073098 .004509 -.008370 -8.6193e-3 -3.8277e-5 4.1945e-5;
        .004509 .040444 -.062260 -1.9932e-6 -4.4565e-4 6.6809e-4;
        -.00837 -.062260 .107335 5.2887e-5 6.9359e-4 -1.2223e-3;
        -.86193e-4 -1.9932e-6 5.2887e-5 1.3926e-5 -3.7863e-7 -9.578e-7;
        -3.8277e-5 -4.4565e-4 6.9359e-4 -3.7863e-7 5.9647e-6 -9.5101e-6;
        4.1945e-5 6.6809e-4 -1.223e-3 -9.578e-7 -9.5101e-6 1.9962e-5];
dxHatold = [-.000005; .000095; -.000005; 0; 0; 0];
atwbold = inv(Pold)*dxHatold;
xHat = [r;v];
dxHat = inv(atwa+inv(Pold))*(atwb+atwbold);
xHat = xHat +dxHat;
P = inv(atwa+inv(Pold));

rnew = r +dxHat(1:3);
vnew = v+ dxHat(4:6);

%% Results and Analysis
fprintf('The updated results after one iteration are: \n')
fprintf('rnew =\n')
confidence = diag(P);
ijk = 'ijk';
for i = 1:3
    fprintf('%f+/-%f %c km\n', rnew(i), confidence(i), ijk(i));
end
fprintf('vnew = \n')
for i = 1:3
    fprintf('%f+/-%f %c km/s \n',vnew(i), confidence(i+3), ijk(i));
end
fprintf('This is a %f%% improvement in the confidence of the state terms.\n',...
    abs((max(diag(P))-max(diag(Pold)))/max(diag(Pold)))*100);
fprintf('In order to get a better state value, we would need the initial \n')
fprintf('observation data so that we could iterate on delta xHat by continuing\n')
fprintf('to calculate new P, A and b matrices. As noted in the comments above,\n')
fprintf('I used the final P and delta xHat values from example 10-4 in the book \n')
fprintf('as the "kth" term before the new term I calculated. This resulted\n')
fprintf('in an answer more similar to the one in example 10-5. In order to\n')
fprintf('fully implement a sequential batch filter, you would need the observation\n')
fprintf('data and to continue to update your A, W and B matrices.\n')
