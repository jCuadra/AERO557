function [f, g] =  kepler  ( ro,vo, dtseco )
% Universal Varibale Function from Vallado modified By Morgan Yost
% -------------------------  implementation   -----------------
% set constants and intermediate printouts
constmath;
constastro;
numiter    =    50;

% --------------------  initialize values   -------------------
ktr   = 0;
xold  = 0.0;
znew  = 0.0;
errork = '      ok';
dtsec = dtseco;
mulrev = 0;
mu = 398600;

if ( abs( dtseco ) > small )
    magro = mag( ro );
    magvo = mag( vo );
    rdotv= dot( ro,vo );
    
    % -------------  find sme, alpha, and a  ------------------
    sme= ( (magvo^2)*0.5  ) - ( mu /magro );
    alpha= -sme*2.0/mu;
    
    if ( abs( sme ) > small )
        a= -mu / ( 2.0 *sme );
    else
        a= infinite;
    end
    if ( abs( alpha ) < small )   % parabola
        alpha= 0.0;
    end
    
    
    % ------------   setup initial guess for x  ---------------
    % -----------------  circle and ellipse -------------------
    
    period= twopi * sqrt( abs(a)^3.0/mu  );
    % ------- next if needed for 2body multi-rev ----------
    
end

ktr= 1;
dtnew = -10.0;

while ((abs(dtnew/sqrt(mu) - dtsec) >= small) & (ktr < numiter))
    xoldsqrd = xold*xold;
    znew     = xoldsqrd * alpha;
    
    % ------------- find c2 and c3 functions --------------
    [c2new,c3new] = findc2c3( znew );
    
    % ------- use a newton iteration for new values -------
    rval = xoldsqrd*c2new + rdotv/sqrt(mu) *xold*(1.0 -znew*c3new) + ...
        magro*( 1.0  - znew*c2new );
    dtnew= xoldsqrd*xold*c3new + rdotv/sqrt(mu)*xoldsqrd*c2new + ...
        magro*xold*( 1.0  - znew*c3new );
    
    % ------------- calculate new value for x -------------
    xnew = xold + ( dtsec*sqrt(mu) - dtnew ) / rval;
    
    % ------------------------------------------------------
    % check if the orbit is an ellipse and xnew > 2pi sqrt(a), the step
    % size must be changed.  this is accomplished by multiplying rval
    % by 10.0 .  note that 10.0  is arbitrary, but seems to produce good
    % results.  the idea is to keep xnew from increasing too rapidily.
    % ------------------------------------------------------
    %  including this doesn't work if you don't mod the dtsec
    %               if ( ( a > 0.0  ) and ( abs(xnew)>twopi*sqrt(a) ) and ( sme < 0.0  ) )
    %                   dx= ( dtsec-dtnew ) / rval  % *7.0   * 10.0
    %                   xnew = xold + dx / 7.0    % /(1.0  + dx)
    %                alternate method to test various values of change
    %                   xnew = xold + ( dtsec-dtnew ) / ( rval*10 chgamt  )
    %                 end
    
    ktr = ktr + 1;
    xold = xnew;
end
xnewsqrd = xnew*xnew;
f = 1.0  - ( xnewsqrd*c2new / magro );
g = dtsec - xnewsqrd*xnew*c3new/sqrt(mu);

end


