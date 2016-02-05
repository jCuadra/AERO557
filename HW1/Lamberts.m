function [ v1, v2 ] = Lamberts( r1_vect, r2_vect, inc,  t )
%Use lamberts solution to find v1 and v2 if you know r1, r2 (vectors)
%the inclination of the orbital plane in degrees and time in seconds 
%   Detailed explanation goes here
r1= norm(r1_vect);
r2 = norm(r2_vect);
u = 398600;
dTheta = acosd(dot(r1_vect, r2_vect)/(r1*r2));
if (inc<90)
    if(([0 0 1]' *cross(r1_vect, r2_vect))< 0)
        dTheta = 360-dTheta;
    end
else
     if(([0 0 1]' *cross(r1_vect, r2_vect))>= 0)  
           dTheta = 360 - dTheta;
     end
end
z = 0;
c =.5;
s = 1/6;
A = sin(dTheta)*((r1*r2)/(1-cos(dTheta)))^.5;

y = getY(r1, r2, A, z, s, c);
X = (y/c)^.5;
tol = 10^-9;
dt = (s*X^3)/(u^.5) +A*(y/u)^.5;
zl = -4*pi^2;
zu = 4*pi^2;
dtl = 0;
count = 0;
while(abs(dtl-tol)<t)
    if count >= 500
        disp('too manys')
        break;
    end
    
    c = C(z);
    s = S(z);
    y = Y(r1, r2, A, z, s, c);
    X = (y/c)^.5;
    dtl = dtl + (s*X^3)/(u^.5) +A*((y/u)^.5);
    if( dtl<=dt)
        zl = z;
    else
        zu = z;
    end
    z = (abs(zl+zu))/2;
    count = count +1;
end


f = 1- (y/r1);
g =A*((y/u)^.5);
gdot = 1-(y/r2);
fdot = (f*gdot - 1)/g;
v1 = (r2_vect-f*r1_vect)/g;
v2 = (fdot*r1_vect + gdot*v1);

end

