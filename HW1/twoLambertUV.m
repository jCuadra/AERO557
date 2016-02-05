function [ f1, g1] = lambertUV( r1, r2, r3, tau1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    mu = 398600;
    psi1 = 0;
    c2_1 = 1/2;
    c3_1 = 1/6;

    psi1_l = -4*pi;
    psi1_h = 4*pi^2;

    c_dv12 = dot(r1, r2)/norm(r1)/norm(r2);

    A12 = sqrt(norm(r1)*norm(r2)*(1+c_dv12));

    count = 0;
    err_dt = 10;
    while(err_dt >1e-6)
        y1 =getY( r1, r2, A12, psi1, c2_1, c3_1 );
        while(A12 >0 && y1<0)
            disp('adjusting')
            psi1_l = psi1_l +.01;
            psi1 = (psi1_l + psi1_h)/2;
            y1 =getY( r1, r2, A12, psi1, c2_1, c3_1 );
        end
       
        X1 = sqrt(y1/c2_1);
        dt1 = ((X1^3)*c3_1 + A12*sqrt(y1))/sqrt(mu);
        
        if dt1 <= tau1
            psi1_l = psi1;
        else
            psi1_h = psi1;
        end
        
        psi1 = (psi1_l + psi1_h)/2;
        
        [c2_1, c3_1] = lambertC(psi1);
        err_dt = abs(dt1-tau1);
        count = count +1;
        if count >50
            disp('quitting')
            break;
        end
    end
    f1 = 1-(y1/norm(r1));
    g1 = A12*sqrt(y1/mu);

end

