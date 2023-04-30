function [U,Bi,V] = Bidiag_Francis_Step_Update_U_V( U,Bi,V )

 m = size(Bi,1);
    
    G = Givens_rotation( [ Bi(1,1) ^ 2 - (Bi(m,m) ^ 2 + Bi(m-1,m))
                                      Bi(1,2)* Bi(1,1)           ]);

    V(:,1:2) = V(:,1:2) * G;

    Bi(1:2,1:2) = Bi(1:2,1:2) * G;

    for index = 1:m-1

        F = Givens_rotation( [ Bi(index,index);
                               Bi(index+1,index)]);

        U(:,index:index+1) = U(:,index:index+1) * F;

        if index + 2 > m      
    
            Bi(index:index+1,index:index+1) = F' * Bi(index:index+1,index:index+1);
        else

            Bi(index:index+1,index:index+2) = F' * Bi(index:index+1,index:index+2);
        end

        if index + 2 <= m
    
            G_hat = Givens_rotation( [ Bi(index,index+1)
                                       Bi(index,index+2)] );
            V(:,index+1:index+2) = V(:,index+1:index+2) * G_hat;
    
            Bi(index:index+2,index+1:index+2) = Bi(index:index+2,index+1:index+2) * G_hat;
        end
    end

end