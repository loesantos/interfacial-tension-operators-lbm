//=============================== E0 -- verificacao dos operadores de tensao interfacial ============================//
//
//  Confere, com as funcoes REAIS de Operadores_tensao.cpp e da biblioteca:
//
//   1) classe P: M0 = 0, M1 = 0, M2 = ( 2/9 ) A |g| ( n n - lambda I ) para normais aleatorias, e o
//      momento de 4a ordem contraido  M4_nnnn  contra o angulo no plano ( anisotropia de ordem superior );
//   2) colapso nz = 1: o operador P3( chi ) dobrado em z coincide com o membro chi + 2 da familia D2Q9;
//   3) fonte de Guo: soma S = 0 , soma S c = ( 1 - 1/2tau ) F ;
//   4) classe F: forca sobre um perfil tanh plano ( deve ser ~0 ) e sobre um disco ( integral radial
//      da forca = sigma / R ).
//
//  Compilar:  g++ -std=c++17 -O2 -fopenmp -I<SimBoltz_Functions> verifica_operadores.cpp -o verifica_operadores
//
//====================================================================================================================//

#define nvel 19
#define dim  3

#include "Definitions.cpp"
#include "LBM_functions.cpp"
#include "Collision_operators.cpp"
#include "Boundary_conditions.cpp"
#include "Initial_conditions.cpp"
#include "Other_functions.cpp"
#include "Operadores_tensao.cpp"

static int falhas = 0;

static void confere ( const char *o_que, double erro, double tol )
{
	const bool ok = fabs ( erro ) <= tol;
	if ( ! ok ) falhas++;
	printf ( "   %-58s  erro = %10.3e   %s\n", o_que, erro, ok ? "ok" : "FALHOU" );
}

int main ()
{
	LATTICE lattice;
	lattice.c_i = new double[ nvel * dim ];
	def_lattice_d3q19 ( lattice );

	const double A = 0.37, mg = 0.81;			// amplitude e |g| arbitrarios

	struct { OPERADOR op; double lambda; } casos[] = {
		{ { OP_P1, 0.5, 0, "P1" }, 0.0 },
		{ { OP_P3, -0.4, 0, "P3 chi=-2/5 (RP)" }, 1.0 },
		{ { OP_P3, 0.5,  0, "P3 chi=1/2 (LKR)" }, 1.0 },
		{ { OP_P3, 2.0,  0, "P3 chi=2 (LVK)" }, 1.0 },
	};

	mt19937 rng ( 7 );
	normal_distribution<double> N01;

	//------ 1) classe P -------------------------------------------------------------------------------------------//

	printf ( "\n1) Classe P: conservacao e segundo momento ( 200 normais aleatorias )\n" );

	for ( auto& cs : casos )
	{
		double e0 = 0, e1 = 0, e2 = 0;

		for ( int k = 0; k < 200; k++ )
		{
			double n[3] = { N01 ( rng ), N01 ( rng ), N01 ( rng ) };
			const double m = sqrt ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
			for ( auto& v : n ) v /= m;

			double g[3] = { mg * n[0], mg * n[1], mg * n[2] };
			double f[nvel] = { 0 };

			perturbacao ( f, A, g, cs.op, lattice );

			double M0 = 0, M1[3] = { 0 }, M2[3][3] = { { 0 } };

			for ( int i = 0; i < nvel; i++ )
			{
				const double *c = lattice.c_i + i * dim;
				M0 += f[i];
				for ( int a = 0; a < 3; a++ ) { M1[a] += f[i] * c[a]; for ( int b = 0; b < 3; b++ ) M2[a][b] += f[i] * c[a] * c[b]; }
			}

			e0 = max ( e0, fabs ( M0 ) );
			for ( int a = 0; a < 3; a++ )
			{
				e1 = max ( e1, fabs ( M1[a] ) );
				for ( int b = 0; b < 3; b++ )
				{
					const double alvo = 2.0 / 9.0 * A * mg * ( n[a] * n[b] - cs.lambda * ( a == b ) );
					e2 = max ( e2, fabs ( M2[a][b] - alvo ) );
				}
			}
		}

		printf ( "  %s  ( lambda = %.4f )\n", cs.op.nome.c_str (), cs.lambda );
		confere ( "M0", e0, 1e-14 );
		confere ( "M1", e1, 1e-14 );
		confere ( "M2 - (2/9) A |g| ( nn - lambda I )", e2, 1e-14 );
	}

	//------ anisotropia de 4a ordem no plano ------------------------------------------------------------------------//

	printf ( "\n   M4_nnnn / ( A |g| ) contra o angulo no plano ( 0, 15, 22.5, 30, 45 graus )\n" );

	for ( auto& cs : casos )
	{
		printf ( "   %-20s", cs.op.nome.c_str () );

		for ( double th : { 0.0, 15.0, 22.5, 30.0, 45.0 } )
		{
			const double t = th * M_PI / 180.0;
			double n[3] = { cos ( t ), sin ( t ), 0.0 };
			double g[3] = { mg * n[0], mg * n[1], 0.0 };
			double f[nvel] = { 0 };

			perturbacao ( f, A, g, cs.op, lattice );

			double M4 = 0;
			for ( int i = 0; i < nvel; i++ )
			{
				const double *c = lattice.c_i + i * dim;
				const double cn = c[0] * n[0] + c[1] * n[1];
				M4 += f[i] * cn * cn * cn * cn;
			}
			printf ( "  %9.6f", M4 / ( A * mg ) );
		}
		printf ( "\n" );
	}

	//------ 2) colapso nz = 1 -------------------------------------------------------------------------------------//

	printf ( "\n2) Colapso nz = 1: P3( chi ) dobrado em z  ==  familia D2Q9 com chi2 = chi + 2\n" );

	for ( double chi : { -0.4, 0.5, 2.0, 7.0 } )
	{
		OPERADOR op; op.tipo = OP_P3; op.chi = chi;

		double t = 0.3;
		double g[3] = { mg * cos ( t ), mg * sin ( t ), 0.0 };
		double f[nvel] = { 0 };

		perturbacao ( f, A, g, op, lattice );

		//  dobra em ( cx , cy )

		double f9[3][3] = { { 0 } };
		for ( int i = 0; i < nvel; i++ )
		{
			const double *c = lattice.c_i + i * dim;
			f9[ ( int ) c[0] + 1 ][ ( int ) c[1] + 1 ] += f[i];
		}

		//  D2Q9 direto: w = 4/9, 1/9, 1/36 ; B = familia com chi2

		const double c2x = chi + 2.0;
		const double B9[3] = { - c2x / ( 3.0 * ( c2x + 2.0 ) ), c2x / ( 6.0 * ( c2x + 2.0 ) ), 1.0 / ( 6.0 * ( c2x + 2.0 ) ) };

		double erro = 0;
		for ( int a = -1; a <= 1; a++ )
		for ( int b = -1; b <= 1; b++ )
		{
			const int s = abs ( a ) + abs ( b );
			const double w9 = ( s == 0 ) ? 4.0 / 9.0 : ( s == 1 ? 1.0 / 9.0 : 1.0 / 36.0 );
			const double cg = a * g[0] + b * g[1];
			const double om = A * mg * ( w9 * cg * cg / ( mg * mg ) - B9[s] );
			erro = max ( erro, fabs ( om - f9[a + 1][b + 1] ) );
		}

		char txt[80];
		snprintf ( txt, 80, "chi = %5.2f  ->  chi2 = %5.2f", chi, c2x );
		confere ( txt, erro, 1e-15 );
	}

	//------ 3) Guo --------------------------------------------------------------------------------------------------//

	printf ( "\n3) Fonte de Guo ( source da biblioteca )\n" );
	{
		const double tau = 0.8, rho = 1.03, F[3] = { 1.3e-3, -2.1e-3, 0.0 }, u[3] = { 0.02, -0.01, 0.0 };
		double S[nvel];
		source ( F[0], F[1], F[2], u[0], u[1], u[2], rho, tau, S, lattice );

		double M0 = 0, M1[3] = { 0 };
		for ( int i = 0; i < nvel; i++ ) { M0 += S[i]; for ( int a = 0; a < 3; a++ ) M1[a] += S[i] * lattice.c_i[ i * dim + a ]; }

		confere ( "soma S", M0, 1e-16 );
		confere ( "soma S c - ( 1 - 1/2tau ) F  ( x )", M1[0] - ( 1 - 0.5 / tau ) * F[0], 1e-16 );
		confere ( "soma S c - ( 1 - 1/2tau ) F  ( y )", M1[1] - ( 1 - 0.5 / tau ) * F[1], 1e-16 );
	}

	//------ 4) classe F sobre perfis analiticos ---------------------------------------------------------------------//

	printf ( "\n4) Classe F sobre perfis tanh ( sigma = 0.1, R = 60 )\n" );
	{
		const int nx = 256, ny = 256;
		const double sigma = 0.1, R = 60.0;

		GEOMETRY geo;
		geo.nx = nx; geo.ny = ny; geo.nz = 1; geo.fluid = nx * ny;
		geo.ini = new int[ nx * ny ];
		for ( int p = 0; p < nx * ny; p++ ) geo.ini[p] = p + 1;

		lattice.inif_R = new double[ ( size_t ) nx * ny * nvel ];
		lattice.inif_B = new double[ ( size_t ) nx * ny * nvel ];

		CAMPOS_F cf;
		aloca_campos_F ( cf, geo, lattice );
		double *F = new double[ 3 * nx * ny ];

		double xi = 0;

		auto semeia = [&] ( bool disco )
		{
			for ( int y = 0; y < ny; y++ )
			for ( int x = 0; x < nx; x++ )
			{
				const double s = disco ? R - hypot ( x - nx / 2.0, y - ny / 2.0 )
				                       : ( fabs ( x - nx / 2.0 ) < nx / 4.0 ? nx / 4.0 - fabs ( x - nx / 2.0 ) : - ( fabs ( x - nx / 2.0 ) - nx / 4.0 ) );
				const double ph = tanh ( s / xi );
				const int p = x + y * nx;
				for ( int i = 0; i < nvel; i++ )
				{
					lattice.inif_R[ ( size_t ) p * nvel + i ] = lattice.w[i] * 0.5 * ( 1 + ph );
					lattice.inif_B[ ( size_t ) p * nvel + i ] = lattice.w[i] * 0.5 * ( 1 - ph );
				}
			}
		};

		for ( double xi_ : { 1.569988, 3.0, 6.0 } )
		for ( int tipo : { OP_F1, OP_F2, OP_F3 } )
		{
			xi = xi_;
			OPERADOR op; op.tipo = tipo; op.xi = xi;
			op.nome = ( tipo == OP_F1 ) ? "F1 CSF" : ( tipo == OP_F2 ? "F2 div T" : "F3 mu grad phi" );

			//  plano: maior |F| relativo a sigma / xi ( escala natural da forca na interface )

			semeia ( false );
			calcula_forca ( cf, op, sigma, lattice, F );
			double fmax = 0;
			for ( int p = 0; p < nx * ny; p++ ) fmax = max ( fmax, fabs ( F[3 * p] ) + fabs ( F[3 * p + 1] ) );

			//  disco: integral radial ao longo de 16 raios,  int F . ( - r_hat ) dr  ~  sigma / R

			semeia ( true );
			calcula_forca ( cf, op, sigma, lattice, F );
			double soma = 0;
			for ( int y = 0; y < ny; y++ )
			for ( int x = 0; x < nx; x++ )
			{
				const double dx = x - nx / 2.0, dy = y - ny / 2.0, r = hypot ( dx, dy );
				if ( r < 1e-9 ) continue;
				const int p = x + y * nx;
				soma += - ( F[3 * p] * dx + F[3 * p + 1] * dy ) / r / ( 2.0 * M_PI * r );	// media azimutal * dr
			}

			//  F3 nao anula exatamente no plano: o laplaciano discreto erra O( 1 / xi^2 ) sobre o tanh.
			//  Isso e resultado ( H5 ), nao defeito; a tolerancia de F3 so exige a convergencia em xi.

			//  No disco sobra a correcao de espessura finita, O( ( xi / R )^2 ), em todas as forcas.

			const double tol_p = ( tipo == OP_F3 ) ? 0.2 / ( xi * xi ) : 1e-3;
			const double tol_d = 1.5 * ( xi / R ) * ( xi / R ) + 5e-4 + ( tipo == OP_F3 ? 0.25 / ( xi * xi ) : 0.0 );

			printf ( "  %s   ( xi = %.3f )\n", op.nome.c_str (), xi );
			confere ( "plano: max |F| / ( sigma / xi )", fmax / ( sigma / xi ), tol_p );
			confere ( "disco: ( int F_r dr ) / ( sigma / R ) - 1", soma / ( sigma / R ) - 1.0, tol_d );
		}
	}

	printf ( "\n%s ( %d falha%s )\n\n", falhas ? "HA FALHAS" : "TUDO CERTO", falhas, falhas == 1 ? "" : "s" );

	return falhas ? 1 : 0;
}
