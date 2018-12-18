function Param = getParam(Param,Problem_Name)

switch Problem_Name
    case 'Sphere'
        fitnessfun   = @Sphere;
        n_gene       = 100;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -5.12;
        ub(1:n_gene) =  5.12;
        
    case 'ScaledSphere'
        fitnessfun   = @ScaledSphere;
        n_gene       = 100;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -5.12;
        ub(1:n_gene) =  5.12;
        
    case 'Ellipsoid'
        fitnessfun   = @Ellipsoid;
        n_gene       = 100;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -5.12;
        ub(1:n_gene) =  5.12;
        
    case 'Cigar'
        fitnessfun   = @Cigar;
        n_gene       = 100;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -5.12;
        ub(1:n_gene) =  5.12;
        
    case 'k_tablet'
        fitnessfun   = @k_tablet;
        n_gene       = 100;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -5.12;
        ub(1:n_gene) =  5.12;
        
    case 'MMbenchmark'
        fitnessfun   = @MMbenchmark;
        n_gene       = 100;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -1;
        ub(1:n_gene) =  1;
        
    case 'Rosenbrock_star'
        fitnessfun   = @Rosenbrock_star;
        n_gene       = 50;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -2.048;
        ub(1:n_gene) =  2.048;
        
    case 'Rosenbrock_chain'
        fitnessfun   = @Rosenbrock_chain;
        n_gene       = 100;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -2.048;
        ub(1:n_gene) =  2.048;
        
    case 'Ackley'
        fitnessfun   = @Ackley;
        n_gene       = 50;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -32.768;
        ub(1:n_gene) =  32.768;
        
    case 'Bohachevsky'
        fitnessfun   = @Bohachevsky;
        n_gene       = 50;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -5.12;
        ub(1:n_gene) =  5.12;
        
    case 'Rastrigin'
        fitnessfun   = @Rastrigin;
        n_gene       = 20;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -5.12;
        ub(1:n_gene) =  5.12;
        
    case 'Schaffer'
        fitnessfun   = @Schaffer;
        n_gene       = 20;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -100;
        ub(1:n_gene) =  100;
        
    case 'Schwefel'
        fitnessfun   = @Schwefel;
        n_gene       = 10;
        n_constraint = 0;
        vtr          = 1e-6;
        lb(1:n_gene) = -512;
        ub(1:n_gene) =  512;
        
    case 'g01'
        fitnessfun   = @g01;
        n_gene       = 13;
        n_constraint = 9;
        vtr          = -15 * ( 1 - 1e-2 );
        lb( 1: 9)    = 0; 
        ub( 1: 9)    = 1;
        lb(10:12)    = 0;
        ub(10:12)    = 100;
        lb(13)       = 0;
        ub(13)       = 1;
        
    case 'g02'
        fitnessfun   = @g02;
        n_gene       = 20;
        n_constraint = 2;
        vtr          = -0.803619 * ( 1 - 1e-2 );
        lb(1:n_gene) = 0;
        ub(1:n_gene) = 10;
        
    case 'g03'
        fitnessfun   = @g03;
        n_gene       = 10;
        n_constraint = 1;
        vtr          = -1 * ( 1 - 1e-2 );
        lb(1:n_gene) = 0;
        ub(1:n_gene) = 1;
        
    case 'g04'
        fitnessfun   = @g04;
        n_gene       = 5;
        n_constraint = 6;
        vtr          = -30665.539 * ( 1 - 1e-2 );
        lb(1)        = 78; 
        ub(1)        = 102;
        lb(2)        = 33;
        ub(2)        = 45;
        lb(3:5)      = 27;
        ub(3:5)      = 45;
        
    case 'g05'
        fitnessfun   = @g05;
        n_gene       = 4;
        n_constraint = 5;
        vtr          = 5126.4981 * ( 1 + 1e-2 );
        lb(1:2)      = 0;
        ub(1:2)      = 1200;
        lb(3:4)      = -0.55;
        ub(3:4)      =  0.55;
        
    case 'g06'
        fitnessfun   = @g06;
        n_gene       = 2;
        n_constraint = 2;
        vtr          = -6961.81388 * ( 1 - 1e-2 );
        lb(1)        = 13;
        ub(1)        = 100;
        lb(2)        = 0;
        ub(2)        = 100;
        
    case 'g07'
        fitnessfun   = @g07;
        n_gene       = 10;
        n_constraint = 8;
        vtr          = 24.3062091 * ( 1 + 1e-2 );
        lb(1:n_gene) = -10;
        ub(1:n_gene) =  10;
        
    case 'g08'
        fitnessfun   = @g08;
        n_gene       = 2;
        n_constraint = 2;
        vtr          = -0.095825 * ( 1 - 1e-2 );
        lb(1:n_gene) =  0;
        ub(1:n_gene) = 10;
        
    case 'g09'
        fitnessfun   = @g09;
        n_gene       = 7;
        n_constraint = 4;
        vtr          = 680.6300573 * ( 1 + 1e-2 );
        lb(1:n_gene) = -10;
        ub(1:n_gene) =  10;
        
    case 'g10'
        fitnessfun   = @g10;
        n_gene       = 8;
        n_constraint = 6;
        vtr          = 7049.3307 * ( 1 + 1e-2 );
        lb(1)        = 2;
        ub(1)        = 4;
        lb(2:3)      = 3;
        ub(2:3)      = 4;
        lb(4:8)      = 1;
        ub(4:8)      = 3;
        
    case 'g11'
        fitnessfun   = @g11;
        n_gene       = 2;
        n_constraint = 1;
        vtr          = 0.75 * ( 1 + 1e-2 );
        lb(1:n_gene) = -1;
        ub(1:n_gene) =  1;
        
    case 'g12'
        fitnessfun   = @g12;
        n_gene       = 3;
        n_constraint = 1;
        vtr          = -1 * ( 1 - 1e-2 );
        lb(1:n_gene) =  0;
        ub(1:n_gene) = 10;
        
    case 'g13'
        fitnessfun   = @g13;
        n_gene       = 5;
        n_constraint = 3;
        vtr          = 0.0539498 * ( 1 + 1e-2 );
        lb(1:2)      = -2.3;
        ub(1:2)      =  2.3;
        lb(3:5)      = -3.2;
        ub(3:5)      =  3.2;
        
    otherwise
        error('Unexpected Problem_Name!');
end

Param.fitnessfun   = fitnessfun;
Param.n_gene       = n_gene;
Param.n_constraint = n_constraint;
Param.vtr          = vtr;
Param.lb           = lb;
Param.ub           = ub;
