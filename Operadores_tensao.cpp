//=============================== Operadores de tensao interfacial sob condicoes comuns ==============================//
//
//  Artigo "Interfacial Tension Paper".  Este arquivo reune, num lugar so, TODAS as formas de impor a tensao
//  interfacial que o artigo compara.  Tudo o mais e comum e vem da biblioteca SimBoltz_Functions:
//
//      rede          D3Q19 ( nz = 1 nos ensaios 2D )
//      colisao       BGK de tempo unico,  coll_BGK()  ( com forcamento de Guo quando ha forca )
//      recoloracao   Latva-Kokko & Rothman,  recolloring() ,  beta = recoll
//      gradiente     isotropico,  ( 1 / cs^2 ) soma_i w_i phi( x + c_i ) c_i  -- o mesmo dos mediadores
//      fase          rho^N = ( rho_R - rho_B ) / rho ,  em [ -1 , 1 ]
//
//  Nada aqui altera a biblioteca: os operadores novos ficam neste arquivo, ao lado do programa.
//
//  ----------------------------------------------------------------------------------------------------------------
//  CLASSE P -- perturbacao par, aplicada as populacoes pos-BGK
//
//      Omega_i = A | grad rho^N | [ w_i ( c_i . n )^2 - B_i ]
//
//  Conserva massa ( soma B_i = cs^2 ) e momento ( paridade ).  O segundo momento e
//
//      soma_i Omega_i c_i c_i = ( 2 / 9 ) A | grad rho^N | ( n n - lambda I )
//
//  e a tensao, pela integral mecanica, e  sigma = ( 4 / 9 ) A tau  para QUALQUER lambda e qualquer chi.
//
//      P1  Gunstensen/Rothman ponderado   B_i = cs^2 w_i                 lambda = 0
//      P3  familia de Liu-Valocchi-Kang   B_i( chi ), eq. (29)           lambda = 1
//
//  Em P3, chi = -2/5 , 1/2 e 2 sao Reis-Phillips ( D2Q9 ), Latva-Kokko/Leclaire e Liu-Valocchi-Kang:
//  com nz = 1 o membro chi do D3Q19 colapsa no membro chi + 2 da familia D2Q9.
//
//  CLASSE F -- forca de corpo F, aplicada pelo esquema de Guo dentro de coll_BGK
//
//      F1  CSF                     F = ( sigma / 2 ) kappa grad rho^N ,      kappa = - div n
//      F2  divergencia de tensao   F = div T ,   T = ( sigma / 2 ) ( |g| I - g g / |g| ) ,   g = grad rho^N
//      F3  potencial quimico       F = mu grad rho^N ,   mu = 4 b rho^N ( rho^N^2 - 1 ) - k lap rho^N
//
//  O fator 1/2 vem do salto 2 de rho^N.  Em F3,  k = 3 sigma xi / 4  e  b = 3 sigma / ( 8 xi ) , que
//  dao o perfil  rho^N = tanh( s / xi )  com tensao sigma.  O perfil de Latva-Kokko e exatamente esse
//  tanh no continuo, com  xi = 2 / K  ( K = 1.273895 em beta = 0.8, D3Q19  =>  xi = 1.569988 ).
//
//  Nas tres forcas sigma e ENTRADA ( parameters.A_fact ).  Nas perturbacoes A_fact e a amplitude A.
//
//====================================================================================================================//

#define OP_P1   1
#define OP_P3   3
#define OP_F1   11
#define OP_F2   12
#define OP_F3   13

struct OPERADOR
{
	int    tipo = OP_P3;				// um dos OP_* acima
	double chi  = 0.5;					// so para P3
	double xi   = 1.569988;				// so para F3: largura do tanh, 2 / K
	string nome = "P3";

	bool forca () const { return tipo >= OP_F1; }
};

//------------------ Leitura de operador.txt --------------------------------------------------------------------------//
//
//  Uma linha:   <P1|P3|F1|F2|F3>  [ chi ( P3 ) | xi ( F3 ) ]
//  Exemplos:    P3 0.5        P3 -0.4        F1        F3 1.569988

bool le_operador ( const string& arquivo, OPERADOR& op )
{
	ifstream in ( arquivo );

	if ( ! in ) return false;

	string s;

	if ( ! ( in >> s ) ) return false;

	op.nome = s;

	if      ( s == "P1" ) op.tipo = OP_P1;
	else if ( s == "P3" ) { op.tipo = OP_P3; double v; if ( in >> v ) op.chi = v; }
	else if ( s == "F1" ) op.tipo = OP_F1;
	else if ( s == "F2" ) op.tipo = OP_F2;
	else if ( s == "F3" ) { op.tipo = OP_F3; double v; if ( in >> v ) op.xi = v; }
	else return false;

	return true;
}

//------------------ Coeficientes B_i por camada de velocidade --------------------------------------------------------//
//
//  D3Q19: i = 0 repouso, 1..6 axiais ( |c|^2 = 1 ), 7..18 diagonais ( |c|^2 = 2 ).  A funcao usa |c|^2
//  lido da rede, nao o indice, para nao depender da ordem.


inline double coef_B ( const OPERADOR& op, double c2, double w )
{
	switch ( op.tipo )
	{
		case OP_P1:  return w / 3.0;										// cs^2 w_i
		default:
		{
			const double chi = op.chi;
			const double den = 6.0 * chi + 24.0;

			if ( c2 < 0.5 ) return - ( 2.0 + 2.0 * chi ) / ( 3.0 * chi + 12.0 );
			if ( c2 < 1.5 ) return chi / den;
			return 1.0 / den;
		}
	}
}

//------------------ Perturbacao da classe P ------------------------------------------------------------------------//
//
//  f   : populacao total pos-BGK ( entra e sai )
//  g   : gradiente de rho^N ( o sinal nao importa: so entram g.g e |g| )

void perturbacao ( double *f, double A, const double *g, const OPERADOR& op, LATTICE lattice )
{
	const double mod2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];

	if ( mod2 <= 1.0e-60 || A == 0.0 ) return;

	const double mod = sqrt ( mod2 );

	for ( int i = 0; i < nvel; i++ )
	{
		const double *c = lattice.c_i + i * dim;

		const double cg = c[0] * g[0] + c[1] * g[1] + c[2] * g[2];
		const double c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

		f[i] += A * mod * ( lattice.w[i] * cg * cg / mod2 - coef_B ( op, c2, lattice.w[i] ) );
	}
}

//------------------ Campos auxiliares da classe F ------------------------------------------------------------------//
//
//  Tudo por sitio fluido ( pto ), com a tabela de vizinhos viz[ pto * nvel + i ] = pto de x + c_i .
//  Com o dominio periodico todo fluido, estes estenceis sao EXATAMENTE os dos mediadores.

struct CAMPOS_F
{
	int     n    = 0;
	int    *viz  = nullptr;				// vizinho de x + c_i
	double *fase = nullptr;				// rho^N
	double *g    = nullptr;				// grad rho^N          ( 3 por sitio )
	double *aux  = nullptr;				// n ( F1, 3 ) | T ( F2, 6 ) | mu ( F3, 1 )
};

void aloca_campos_F ( CAMPOS_F& c, GEOMETRY geometry, LATTICE lattice )
{
	const int nx = geometry.nx, ny = geometry.ny, nz = geometry.nz;

	c.n    = geometry.fluid;
	c.viz  = new int   [ ( size_t ) c.n * nvel ];
	c.fase = new double[ c.n ];
	c.g    = new double[ ( size_t ) c.n * 3 ];
	c.aux  = new double[ ( size_t ) c.n * 6 ];

	for ( int z = 0; z < nz; z++ )
	for ( int y = 0; y < ny; y++ )
	for ( int x = 0; x < nx; x++ )
	{
		const int pos = x + y * nx + z * nx * ny;

		if ( ! geometry.ini[pos] ) continue;

		const int p = geometry.ini[pos] - 1;

		for ( int i = 0; i < nvel; i++ )
		{
			const double *ci = lattice.c_i + i * dim;

			const int xv = ( x + ( int ) ci[0] + nx ) % nx;
			const int yv = ( y + ( int ) ci[1] + ny ) % ny;
			const int zv = ( z + ( int ) ci[2] + nz ) % nz;

			const int pv = geometry.ini[ xv + yv * nx + zv * nx * ny ];

			c.viz[ ( size_t ) p * nvel + i ] = ( pv > 0 ) ? pv - 1 : p;	// solido: reflete no proprio sitio
		}
	}
}

//  Divergencia de um campo vetorial ( stride s, componentes a partir de off ) no sitio p:
//      div v = ( 1 / cs^2 ) soma_i w_i c_i . v( x + c_i )

inline double div_iso ( const CAMPOS_F& c, const double *v, int s, int p, LATTICE lattice )
{
	double d = 0.0;

	for ( int i = 1; i < nvel; i++ )
	{
		const double *ci = lattice.c_i + i * dim;
		const double *vv = v + ( size_t ) c.viz[ ( size_t ) p * nvel + i ] * s;

		d += lattice.w[i] * ( ci[0] * vv[0] + ci[1] * vv[1] + ci[2] * vv[2] );
	}

	return 3.0 * d;
}

void calcula_forca ( CAMPOS_F& c, const OPERADOR& op, double sigma, LATTICE lattice, double *F )
{
	const int n = c.n;

	//------ fase ----------------------------------------------------------------------------------------------------//

	#pragma omp parallel for
	for ( int p = 0; p < n; p++ )
	{
		const double rR = density ( lattice.inif_R + ( size_t ) p * nvel );
		const double rB = density ( lattice.inif_B + ( size_t ) p * nvel );

		c.fase[p] = ( rR - rB ) / ( rR + rB );
	}

	//------ gradiente -----------------------------------------------------------------------------------------------//

	#pragma omp parallel for
	for ( int p = 0; p < n; p++ )
	{
		double gx = 0.0, gy = 0.0, gz = 0.0;

		for ( int i = 1; i < nvel; i++ )
		{
			const double *ci = lattice.c_i + i * dim;
			const double  ph = lattice.w[i] * c.fase[ c.viz[ ( size_t ) p * nvel + i ] ];

			gx += ph * ci[0];  gy += ph * ci[1];  gz += ph * ci[2];
		}

		c.g[ 3 * p + 0 ] = 3.0 * gx;
		c.g[ 3 * p + 1 ] = 3.0 * gy;
		c.g[ 3 * p + 2 ] = 3.0 * gz;
	}

	//------ campo auxiliar ------------------------------------------------------------------------------------------//

	const double xi = op.xi;
	const double kf = 0.75 * sigma * xi;					// k da energia livre
	const double bf = 3.0 * sigma / ( 8.0 * xi );		// b da energia livre

	#pragma omp parallel for
	for ( int p = 0; p < n; p++ )
	{
		const double *g = c.g + 3 * p;
		double *a = c.aux + 6 * p;

		const double m = sqrt ( g[0] * g[0] + g[1] * g[1] + g[2] * g[2] );

		if ( op.tipo == OP_F1 )
		{
			//  normal unitaria; zero onde nao ha gradiente ( fora da interface )

			const double inv = ( m > 1.0e-12 ) ? 1.0 / m : 0.0;

			a[0] = g[0] * inv;  a[1] = g[1] * inv;  a[2] = g[2] * inv;
		}
		else if ( op.tipo == OP_F2 )
		{
			//  T = ( sigma / 2 ) ( |g| I - g g / |g| ) : xx yy zz xy xz yz.  Suave em |g| -> 0.

			const double h   = 0.5 * sigma;
			const double inv = ( m > 1.0e-30 ) ? 1.0 / m : 0.0;

			a[0] = h * ( m - g[0] * g[0] * inv );
			a[1] = h * ( m - g[1] * g[1] * inv );
			a[2] = h * ( m - g[2] * g[2] * inv );
			a[3] = - h * g[0] * g[1] * inv;
			a[4] = - h * g[0] * g[2] * inv;
			a[5] = - h * g[1] * g[2] * inv;
		}
		else
		{
			//  mu = 4 b phi ( phi^2 - 1 ) - k lap phi ,  lap phi = ( 2 / cs^2 ) soma_i w_i ( phi_i - phi )

			const double ph = c.fase[p];

			double lap = 0.0;

			for ( int i = 1; i < nvel; i++ )
				lap += lattice.w[i] * ( c.fase[ c.viz[ ( size_t ) p * nvel + i ] ] - ph );

			lap *= 6.0;

			a[0] = 4.0 * bf * ph * ( ph * ph - 1.0 ) - kf * lap;
		}
	}

	//------ forca ---------------------------------------------------------------------------------------------------//

	#pragma omp parallel for
	for ( int p = 0; p < n; p++ )
	{
		const double *g = c.g + 3 * p;
		double *Fp = F + 3 * p;

		if ( op.tipo == OP_F1 )
		{
			const double kappa = - div_iso ( c, c.aux, 6, p, lattice );

			Fp[0] = 0.5 * sigma * kappa * g[0];
			Fp[1] = 0.5 * sigma * kappa * g[1];
			Fp[2] = 0.5 * sigma * kappa * g[2];
		}
		else if ( op.tipo == OP_F2 )
		{
			//  F_a = d_b T_ab : tres divergencias das linhas do tensor

			double Fx = 0.0, Fy = 0.0, Fz = 0.0;

			for ( int i = 1; i < nvel; i++ )
			{
				const double *ci = lattice.c_i + i * dim;
				const double *T  = c.aux + ( size_t ) c.viz[ ( size_t ) p * nvel + i ] * 6;
				const double  w  = lattice.w[i];

				Fx += w * ( T[0] * ci[0] + T[3] * ci[1] + T[4] * ci[2] );
				Fy += w * ( T[3] * ci[0] + T[1] * ci[1] + T[5] * ci[2] );
				Fz += w * ( T[4] * ci[0] + T[5] * ci[1] + T[2] * ci[2] );
			}

			Fp[0] = 3.0 * Fx;  Fp[1] = 3.0 * Fy;  Fp[2] = 3.0 * Fz;
		}
		else
		{
			const double mu = c.aux[ 6 * p ];

			Fp[0] = mu * g[0];  Fp[1] = mu * g[1];  Fp[2] = mu * g[2];
		}
	}
}

//------------------ Colisao de um sitio, comum a todos os operadores -----------------------------------------------//
//
//  1) BGK da populacao total, com a forca F ( Guo ) quando o operador e da classe F;
//  2) perturbacao, quando o operador e da classe P;
//  3) recoloracao de Latva-Kokko com o gradiente dos mediadores.

void colisao_comum ( double *f_R, double *f_B, double *f_m, const double *F, const OPERADOR& op,
                     LATTICE lattice, PARAMETERS parameters )
{
	const double rho_R = density ( f_R );
	const double rho_B = density ( f_B );
	const double rho   = rho_R + rho_B;

	if ( rho <= 1.0e-30 ) return;

	const double conc_R = rho_R / rho;
	const double conc_B = 1.0 - conc_R;

	double gm[3];

	momentum ( f_m, gm[0], gm[1], gm[2], lattice );		// = - grad rho^N

	const double mod = sqrt ( gm[0] * gm[0] + gm[1] * gm[1] + gm[2] * gm[2] );

	double f[nvel];

	for ( int i = 0; i < nvel; i++ ) f[i] = f_R[i] + f_B[i];

	parameters.tau = parameters.tau_R;

	if ( op.forca () ) coll_BGK ( f, F[0] / rho, F[1] / rho, F[2] / rho, lattice, parameters );
	else
	{
		coll_BGK ( f, 0.0, 0.0, 0.0, lattice, parameters );

		perturbacao ( f, parameters.A_fact, gm, op, lattice );
	}

	recolloring ( f, f_R, f_B, gm, mod, conc_R, conc_B, rho, parameters.recoll, lattice );
}

//------------------ Modos de Fourier da interface ( isotropia da gota ) --------------------------------------------//
//
//  Raio da isolinha rho^N = 0 em N_ANG direcoes a partir do centro ( x0 , y0 ), por interpolacao bilinear
//  de rho^N ao longo do raio e bissecao.  Devolve R0 ( media ) e | a_k | / R0 para k = 2, 4, 8.

struct MODOS { double R0 = 0, a2 = 0, a4 = 0, a8 = 0, rmin = 0, rmax = 0; };

MODOS modos_interface ( GEOMETRY geometry, LATTICE lattice, double x0, double y0, double R_est )
{
	const int nx = geometry.nx, ny = geometry.ny;
	const int N_ANG = 256;

	auto fase_em = [&] ( double x, double y )
	{
		x = x - nx * floor ( x / nx );  y = y - ny * floor ( y / ny );

		const int i0 = ( int ) floor ( x ), j0 = ( int ) floor ( y );
		const double tx = x - i0, ty = y - j0;

		double v = 0.0;

		for ( int dj = 0; dj < 2; dj++ )
		for ( int di = 0; di < 2; di++ )
		{
			const int ii = ( i0 + di ) % nx, jj = ( j0 + dj ) % ny;
			const int p  = geometry.ini[ ii + jj * nx ] - 1;
			const double rR = density ( lattice.inif_R + ( size_t ) p * nvel );
			const double rB = density ( lattice.inif_B + ( size_t ) p * nvel );
			const double wgt = ( di ? tx : 1 - tx ) * ( dj ? ty : 1 - ty );

			v += wgt * ( rR - rB ) / ( rR + rB );
		}

		return v;
	};

	MODOS m;

	double c2 = 0, s2 = 0, c4 = 0, s4 = 0, c8 = 0, s8 = 0, soma = 0;

	m.rmin = 1e30;  m.rmax = 0;

	for ( int k = 0; k < N_ANG; k++ )
	{
		const double th = 2.0 * M_PI * k / N_ANG, ct = cos ( th ), st = sin ( th );

		double a = 0.3 * R_est, b = 1.7 * R_est;

		const double fa = fase_em ( x0 + a * ct, y0 + a * st );

		for ( int it = 0; it < 60; it++ )
		{
			const double mid = 0.5 * ( a + b );
			const double fm  = fase_em ( x0 + mid * ct, y0 + mid * st );

			if ( ( fm > 0 ) == ( fa > 0 ) ) a = mid; else b = mid;
		}

		const double r = 0.5 * ( a + b );

		soma += r;
		c2 += r * cos ( 2 * th );  s2 += r * sin ( 2 * th );
		c4 += r * cos ( 4 * th );  s4 += r * sin ( 4 * th );
		c8 += r * cos ( 8 * th );  s8 += r * sin ( 8 * th );

		m.rmin = min ( m.rmin, r );  m.rmax = max ( m.rmax, r );
	}

	m.R0 = soma / N_ANG;
	m.a2 = 2.0 * sqrt ( c2 * c2 + s2 * s2 ) / N_ANG / m.R0;
	m.a4 = 2.0 * sqrt ( c4 * c4 + s4 * s4 ) / N_ANG / m.R0;
	m.a8 = 2.0 * sqrt ( c8 * c8 + s8 * s8 ) / N_ANG / m.R0;

	return m;
}

//====================================================================================================================//
