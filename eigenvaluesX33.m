function gama=eigenvaluescalculator(n,B1,B2)
    
if (B1*B2) < 1        
   
    if  ( ( ( B1 + B2 ) * ( B1 + B2 + ( B1 * B2 ) ) ) / ( ( 3 + B1 + B2 )^2 ) ) < 0.04 

        num1 = 3 * ( B1 + B2 + (B1*B2) ) ;
        den1 = (3 + B1 + B2);
        %gama(1) = sqrt(num1/den1);
        
        inicio1=sqrt(num1/den1);
        gama(1)=fzero(@func,inicio1);

    else
        
        num2 =  ((-45 -( 15 * ( B1 + B2 ) )) + sqrt( 255 * ( ( 3 + B1 + B2 )^2 ) + ( 180 * ( B1 + B2 ) ) * ( ( B1 * B2 ) + B1 + B2 ) )) ;
        den2 = 2 * ( B1 + B2 );
%       gama(1) = sqrt(num2/den2);
        
        inicio2=sqrt(num2/den2);
        gama(1)=fzero(@func,inicio2);

    end
        
else
    
    if ( (B1*B2) < (pi/2)^2 )
        
        num3 = ((pi/2) + sqrt( sqrt(pi/2) + 4 * ( B1 + B2 + 1 ) * ( B1 * B2 / ( B1 + B2 )^2 ) ));
        den3 = 2 * ( B1 + B2 + 1) / ( B1 + B2 );
%         gama(1) = num3/den3;

        inicio3=num3/den3;
        gama(1)=fzero(@func,inicio3);
    
    else
        
        if ( B1 + B2 + 3 ) < ( 16 * B1 * B2 / ( 3 * pi^2 ) )  
        
            num4 = B1*B2*pi;
            den4 = B1+B2+(B1*B2);
            %gama(1) = num4/den4;
            
            inicio4=num4/den4;
            gama(1)=fzero(@func,inicio4);

        else
        
            num5 = ((pi/2) + sqrt( sqrt(pi/2) + 4 * ( B1 + B2 + 1 ) * ( B1 * B2 / ( B1 + B2 )^2 ) ));
            den5 = 2 * ( B1 + B2 + 1) / ( B1 + B2 );
            %gama(1) = num5/den5;

            inicio5=num5/den5;
            gama(1)=fzero(@func,inicio5);
    
        end
        
    end
 
end
    
 for k=2:n
                 
            a=pi/2;
            b=3*pi/2;

            gama(k)=fzero(@func,[gama(k-1)+a, gama(k-1)+b]);
    
    end    
    
    function f=func(gama)
        
        f=((gama^2)-(B1*B2))*sin(gama) - gama*(B1+B2)*cos(gama);
        
    end

end
