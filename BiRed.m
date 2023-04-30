function [ A_out, t_out, r_out ] = BiRed( A, t, r )
    [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );
    [ tT, ...
    tB ] = FLA_Part_2x1( t, ...
                         0, 'FLA_TOP' );
    [ rT, ...
    rB ] = FLA_Part_2x1( r, ...
                         0, 'FLA_TOP' );

    while ( size( ATL, 1 ) < size( A, 1 ))

    [ A00,  a01,     A02,  ...
      a10t, alpha11, a12t, ...
      A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                                    ABL, ABR, ...
                                                    1, 1, 'FLA_BR' );

    [ t0, ...
      tau1, ...
      t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
                                    tB, ...
                                    1, 'FLA_BOTTOM' );

    [ r0, ...
      rho1, ...
      r2 ] = FLA_Repart_2x1_to_3x1( rT, ...
                                    rB, ...
                                    1, 'FLA_BOTTOM' );

    %------------------------------------------------------------%
    
    [ alpha11, ...
      a21, tau1 ] = Housev( alpha11, ...
                            a21 );
    if size(a12t) >= 1

    w12t = (a12t + a21' * A22) / tau1;

    a12t = (a12t - w12t);

    A22 = (A22 - a21 * w12t);

    [u12,rho1] = Housev1(a12t');

    a12t = u12';

    u12(1) = 1;

    y12 = A22 * u12;

    A22 = A22 - (1/rho1) * y12 * u12';

    end

    
    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02, ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                            'FLA_TL' );

    [ tT, ...
      tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
                                       tau1, ...
                                       t2, ...
                                       'FLA_TOP' );

    [ rT, ...
      rB ] = FLA_Cont_with_3x1_to_2x1( r0, ...
                                       rho1, ...
                                       r2, ...
                                       'FLA_TOP' );

  end

  A_out = [ ATL, ATR
            ABL, ABR ];

  t_out = [ tT
            tB ];

  r_out = [ rT
            rB ];

return