function Param = getParam(Param,Problem_Name)

switch Problem_Name
    case 'Sphere'
        Param.fitnessfun   = @Sphere;
        Param.decodingfun  = @Sphere_decode;
        Param.n_gene       = 100;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'ScaledSphere'
        Param.fitnessfun   = @ScaledSphere;
        Param.decodingfun  = @ScaledSphere_decode;
        Param.n_gene       = 100;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Ellipsoid'
        Param.fitnessfun   = @Ellipsoid;
        Param.decodingfun  = @Ellipsoid_decode;
        Param.n_gene       = 100;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Cigar'
        Param.fitnessfun   = @Cigar;
        Param.decodingfun  = @Cigar_decode;
        Param.n_gene       = 100;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'k_tablet'
        Param.fitnessfun   = @k_tablet;
        Param.decodingfun  = @k_tablet_decode;
        Param.n_gene       = 100;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'MMbenchmark'
        Param.fitnessfun   = @MMbenchmark;
        Param.decodingfun  = @MMbenchmark_decode;
        Param.n_gene       = 100;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Rosenbrock_star'
        Param.fitnessfun   = @Rosenbrock_star;
        Param.decodingfun  = @Rosenbrock_star_decode;
        Param.n_gene       = 50;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Rosenbrock_chain'
        Param.fitnessfun   = @Rosenbrock_chain;
        Param.decodingfun  = @Rosenbrock_chain_decode;
        Param.n_gene       = 100;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Ackley'
        Param.fitnessfun   = @Ackley;
        Param.decodingfun  = @Ackley_decode;
        Param.n_gene       = 50;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Bohachevsky'
        Param.fitnessfun   = @Bohachevsky;
        Param.decodingfun  = @Bohachevsky_decode;
        Param.n_gene       = 50;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Rastrigin'
        Param.fitnessfun   = @Rastrigin;
        Param.decodingfun  = @Rastrigin_decode;
        Param.n_gene       = 20;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Schaffer'
        Param.fitnessfun   = @Schaffer;
        Param.decodingfun  = @Schaffer_decode;
        Param.n_gene       = 20;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'Schwefel'
        Param.fitnessfun   = @Schwefel;
        Param.decodingfun  = @Schwefel_decode;
        Param.n_gene       = 10;
        Param.n_constraint = 0;
        Param.vtr          = 1e-6;
        
    case 'g01'
        Param.fitnessfun   = @g01;
        Param.decodingfun  = @g01_decode;
        Param.n_gene       = 13;
        Param.n_constraint = 9;
        Param.vtr          = -15 * ( 1 - 1e-2 );
        
    case 'g02'
        Param.fitnessfun   = @g02;
        Param.decodingfun  = @g02_decode;
        Param.n_gene       = 20;
        Param.n_constraint = 2;
        Param.vtr          = -0.803619 * ( 1 - 1e-2 );
        
    case 'g03'
        Param.fitnessfun   = @g03;
        Param.decodingfun  = @g03_decode;
        Param.n_gene       = 10;
        Param.n_constraint = 1;
        Param.vtr          = -1 * ( 1 - 1e-2 );
        
    case 'g04'
        Param.fitnessfun   = @g04;
        Param.decodingfun  = @g04_decode;
        Param.n_gene       = 5;
        Param.n_constraint = 6;
        Param.vtr          = -30665.539 * ( 1 - 1e-2 );
        
    case 'g05'
        Param.fitnessfun   = @g05;
        Param.decodingfun  = @g05_decode;
        Param.n_gene       = 4;
        Param.n_constraint = 5;
        Param.vtr          = 5126.4981 * ( 1 + 1e-2 );
        
    case 'g06'
        Param.fitnessfun   = @g06;
        Param.decodingfun  = @g06_decode;
        Param.n_gene       = 2;
        Param.n_constraint = 2;
        Param.vtr          = -6961.81388 * ( 1 - 1e-2 );
        
    case 'g07'
        Param.fitnessfun   = @g07;
        Param.decodingfun  = @g07_decode;
        Param.n_gene       = 10;
        Param.n_constraint = 8;
        Param.vtr          = 24.3062091 * ( 1 + 1e-2 );
        
    case 'g08'
        Param.fitnessfun   = @g08;
        Param.decodingfun  = @g08_decode;
        Param.n_gene       = 2;
        Param.n_constraint = 2;
        Param.vtr          = -0.095825 * ( 1 - 1e-2 );
        
    case 'g09'
        Param.fitnessfun   = @g09;
        Param.decodingfun  = @g09_decode;
        Param.n_gene       = 7;
        Param.n_constraint = 4;
        Param.vtr          = 680.6300573 * ( 1 + 1e-2 );
        
    case 'g10'
        Param.fitnessfun   = @g10;
        Param.decodingfun  = @g10_decode;
        Param.n_gene       = 8;
        Param.n_constraint = 6;
        Param.vtr          = 7049.3307 * ( 1 + 1e-2 );
        
    case 'g11'
        Param.fitnessfun   = @g11;
        Param.decodingfun  = @g11_decode;
        Param.n_gene       = 2;
        Param.n_constraint = 1;
        Param.vtr          = 0.75 * ( 1 + 1e-2 );
        
    case 'g12'
        Param.fitnessfun   = @g12;
        Param.decodingfun  = @g12_decode;
        Param.n_gene       = 3;
        Param.n_constraint = 1;
        Param.vtr          = -1 * ( 1 - 1e-2 );
        
    case 'g13'
        Param.fitnessfun   = @g13;
        Param.decodingfun  = @g13_decode;
        Param.n_gene       = 5;
        Param.n_constraint = 3;
        Param.vtr          = 0.0539498 * ( 1 + 1e-2 );
        
    otherwise
        error('Unexpected Problem_Name!');
end
