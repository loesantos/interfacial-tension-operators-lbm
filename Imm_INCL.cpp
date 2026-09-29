//=============================== Lattice Boltzmann Method - Interface plana INCLINADA ( E3, isotropia ) =============//
//
//  Programa irmao de Imm_OPS.cpp, com o mesmo arcabouco comum ( D3Q19 com nz = 1, BGK de tempo unico,
//  recoloracao de Latva-Kokko, mediadores ) e os mesmos operadores de Operadores_tensao.cpp.  Muda
//  so a geometria e a medida: em vez da gota, uma faixa vermelha limitada por DUAS interfaces planas de
//  normal
//
//      n = ( p , q ) / sqrt( p^2 + q^2 ) ,      p , q inteiros ,      tan( theta ) = q / p
//
//  Com a caixa quadrada N x N x 1, o campo  c( x, y ) = f( ( p x + q y ) mod N )  e periodico nas duas
//  direcoes para quaisquer p, q inteiros.  O periodo ao longo da normal e  L = N / sqrt( p^2 + q^2 ),
//  e o vermelho ocupa metade dele:  0 <= ( p x + q y ) mod N < N / 2 .  As duas interfaces ficam a L/2
//  uma da outra e o comprimento total de interface na caixa e  2 N sqrt( p^2 + q^2 ) .
//
//  Numa interface plana e reta a fisica nao tem direcao preferida.  Tudo o que depender de theta vem da
//  rede -- do operador, do estencil do gradiente e da recoloracao, esta ultima comum a todos.  Medidas:
//
//    sigma_pop   integral mecanica  soma ( P_nn - P_T ) / comprimento de interface ,
//                P_T = ( P_tt + P_zz ) / 2 ,  P_ab = soma_i f_i c_ia c_ib  pos-colisao.
//                Mede a tensao que esta NAS POPULACOES: exata nas perturbacoes, ( 4/9 ) A tau em
//                theta = 0.  Nas forcas e ~0 por construcao ( sem curvatura, sem forca ).
//    int_grad    regra de soma  soma |grad rho^N| / comprimento  = 2  ( por interface )
//    K, residuo  ajuste logistico  ln[ c / ( 1 - c ) ] = K ( d - d0 )  contra a distancia d a interface
//                central, com c = rho_R / rho.  Largura 10-90 % = 2 ln 9 / K.  O residuo maximo mede o
//                quanto o perfil deixa de depender so de d ( degraus da rede ).
//    d0          deslocamento da interface ao longo da normal ( deriva ou fixacao na rede )
//    |u|max      corrente espuria; numa interface plana e reta deveria ser ZERO.  Tambem |u_n|max e
//                |u_t|max ( componentes normal e tangencial ).
//    soma F      forca liquida ( classe F )
//
//  Entradas:  data_in.txt, operador.txt ( como em Imm_OPS ), inclinacao.txt ( "p q" ), meio.vtk N x N x 1.
//  Saidas:    inclinada.dat ( serie temporal ), perfil.dat ( c contra d no fim ), resumo.dat ( uma linha ).
//
//====================================================================================================================//


//====================================================================================================================//
//   O P C O E S   D O   C A S O                                                                                      //
//====================================================================================================================//

const double ESPESSURA_INI = 3.45;		// largura 10-90 % do perfil de arranque
const int    PASSOS_DIAG   = 500;
const double U_LIMITE      = 0.4;

//====================================================================================================================//


#define nvel 	19
#define dim 	3

#include "Definitions.cpp"
#include "LBM_functions.cpp"
#include "Collision_operators.cpp"
#include "Boundary_conditions.cpp"
#include "Initial_conditions.cpp"
#include "Other_functions.cpp"
#include "Operadores_tensao.cpp"


//------------------ Medida da interface inclinada -------------------------------------------------------------------//

struct INCL
{
	double sigma_pop = 0, int_grad = 0, K = 0, res = 0, d0 = 0;
	double u_max = 0, u_rms = 0, un_max = 0, ut_max = 0;
	double massa_R = 0, massa_B = 0, qx = 0, qy = 0;
	double min_f = 1e300, SFx = 0, SFy = 0;
	long long neg_total = 0, neg_cor = 0;
	int n_ruim = 0, n_ajuste = 0;
};

//  distancia assinada a interface CENTRAL ( c = N/2 em unidades de p x + q y ), em [ -L/2 , L/2 )

inline double dist_central ( int x, int y, int p, int q, int N )
{
	const double m = sqrt ( ( double ) ( p * p + q * q ) );
	long long s = ( ( long long ) p * x + ( long long ) q * y ) % N;
	if ( s < 0 ) s += N;
	double d = ( ( double ) s - 0.5 * N );
	if ( d >= 0.5 * N ) d -= N;
	if ( d < -0.5 * N ) d += N;
	return - d / m;			// > 0 do lado vermelho ( s < N/2 )
}

bool mede_inclinada ( GEOMETRY geometry, LATTICE lattice, int p, int q, INCL& r, const double *forca,
                      bool gravar_perfil )
{
	const int N = geometry.nx;
	const double m = sqrt ( ( double ) ( p * p + q * q ) );
	const double nx_ = p / m, ny_ = q / m;			// normal
	const double tx_ = - ny_, ty_ = nx_;			// tangente no plano

	r = INCL ();

	double u2 = 0;
	double Sxx = 0, Sx = 0, Sxy = 0, Sy = 0;		// ajuste linear logit( c ) = a + K d
	int nfit = 0;

	vector<double> dd, ll;
	ofstream fp;
	if ( gravar_perfil ) { fp.open ( "perfil.dat" ); fp << "# d  c=rho_R/rho  rho  ux  uy\n" << setprecision ( 10 ); }

	for ( int y = 0; y < N; y++ )
	for ( int x = 0; x < N; x++ )
	{
		const int pto = geometry.ini[ x + y * N ] - 1;
		const double *fR = lattice.inif_R + ( size_t ) pto * nvel;
		const double *fB = lattice.inif_B + ( size_t ) pto * nvel;

		double rR = 0, rB = 0, jx = 0, jy = 0, Pxx = 0, Pyy = 0, Pxy = 0, Pzz = 0;

		for ( int i = 0; i < nvel; i++ )
		{
			const double *c = lattice.c_i + i * dim;
			const double f = fR[i] + fB[i];

			if ( ! std::isfinite ( fR[i] ) || ! std::isfinite ( fB[i] ) ) r.n_ruim++;
			if ( f < -1e-12 ) r.neg_total++;
			if ( fR[i] < -1e-12 ) r.neg_cor++;
			if ( fB[i] < -1e-12 ) r.neg_cor++;
			r.min_f = min ( r.min_f, f );

			rR += fR[i];  rB += fB[i];
			jx += f * c[0];  jy += f * c[1];
			Pxx += f * c[0] * c[0];  Pyy += f * c[1] * c[1];  Pxy += f * c[0] * c[1];  Pzz += f * c[2] * c[2];
		}

		const double rho = rR + rB;
		if ( ! ( rho > 0 ) ) { r.n_ruim++; continue; }

		if ( forca != nullptr )
		{
			jx -= 0.5 * forca[ 3 * pto ];  jy -= 0.5 * forca[ 3 * pto + 1 ];
			r.SFx += forca[ 3 * pto ];  r.SFy += forca[ 3 * pto + 1 ];
		}

		const double ux = jx / rho, uy = jy / rho;

		r.massa_R += rR;  r.massa_B += rB;  r.qx += jx;  r.qy += jy;

		const double u = hypot ( ux, uy ), un = fabs ( ux * nx_ + uy * ny_ ), ut = fabs ( ux * tx_ + uy * ty_ );
		r.u_max = max ( r.u_max, u );  r.un_max = max ( r.un_max, un );  r.ut_max = max ( r.ut_max, ut );
		u2 += u * u;

		//  anisotropia mecanica ( sem o termo convectivo, desprezivel e simetrico )

		const double Pnn = nx_ * nx_ * Pxx + 2 * nx_ * ny_ * Pxy + ny_ * ny_ * Pyy;
		const double Ptt = tx_ * tx_ * Pxx + 2 * tx_ * ty_ * Pxy + ty_ * ty_ * Pyy;
		r.sigma_pop += Pnn - 0.5 * ( Ptt + Pzz );

		//  regra de soma com o gradiente dos mediadores

		double gx, gy, gz;
		momentum ( lattice.inif_m + ( size_t ) pto * nvel, gx, gy, gz, lattice );
		r.int_grad += sqrt ( gx * gx + gy * gy + gz * gz );

		//  perfil ao redor da interface central

		const double d = dist_central ( x, y, p, q, N );
		const double cR = rR / rho;

		if ( gravar_perfil && fabs ( d ) < 12.0 ) fp << d << " " << cR << " " << rho << " " << ux << " " << uy << "\n";

		if ( cR > 0.02 && cR < 0.98 && fabs ( d ) < 8.0 )
		{
			const double l = log ( cR / ( 1.0 - cR ) );
			dd.push_back ( d );  ll.push_back ( l );
			Sx += d;  Sy += l;  Sxx += d * d;  Sxy += d * l;  nfit++;
		}
	}

	const double comprimento = 2.0 * N * m;		// duas interfaces

	r.sigma_pop /= comprimento;
	r.int_grad  /= comprimento;
	r.u_rms = sqrt ( u2 / ( N * N ) );
	r.n_ajuste = nfit;

	if ( nfit > 3 )
	{
		const double K = ( nfit * Sxy - Sx * Sy ) / ( nfit * Sxx - Sx * Sx );
		const double a = ( Sy - K * Sx ) / nfit;
		r.K  = K;
		r.d0 = - a / K;
		for ( size_t k = 0; k < dd.size (); k++ ) r.res = max ( r.res, fabs ( ll[k] - ( a + K * dd[k] ) ) );
	}

	return r.n_ruim == 0 && std::isfinite ( r.u_max );
}


int main ()
{
	//=========================== Le o arquivo de entrada ============================================================//

	GEOMETRY geometry;
	PARAMETERS parameters;

	read_data ( geometry, parameters );

	OPERADOR op;
	if ( ! le_operador ( "operador.txt", op ) )
	{
		cerr << "\noperador.txt ausente ou invalido.  Uma linha: P1 | P3 chi | F1 | F2 | F3 xi\n";
		return 1;
	}

	int p = 1, q = 0;
	{
		ifstream in ( "inclinacao.txt" );
		if ( ! ( in >> p >> q ) || ( p == 0 && q == 0 ) ) { cerr << "inclinacao.txt: \"p q\" inteiros, nao ambos nulos.\n"; return 1; }
	}

	if ( fabs ( parameters.tau_R - parameters.tau_B ) > 1e-14 ) { cerr << "tau_R e tau_B tem de ser iguais.\n"; return 1; }

	//=========================== Initializing variables =============================================================//

	const int n_steps = parameters.n_steps;

	int set_threads = parameters.n_threads;
	if ( set_threads == 0 ) set_threads = omp_get_max_threads ();
	omp_set_num_threads ( set_threads );

	//=========================== Read the geometry file (.vtk) ======================================================//

	{
		ifstream file_geo ( geometry.file );
		if ( ! file_geo ) { cout << "\nNao foi possivel abrir " << geometry.file << endl; return 1; }
		string line, dump; stringstream data;
		for ( int i = 0; i < 4; i++ ) getline ( file_geo, dump );
		getline ( file_geo, line );
		data << line;
		data >> dump >> geometry.nx >> geometry.ny >> geometry.nz;
	}

	const int N = geometry.nx;

	geometry.ini   = new int[ geometry.nx * geometry.ny * geometry.nz ];
	geometry.fluid = read_geo ( geometry.file, geometry.ini, 0, 0 );

	const int n_fluid = geometry.fluid;

	if ( geometry.nz != 1 || geometry.ny != N || n_fluid != N * N )
	{
		cerr << "\nEste programa requer caixa quadrada N x N x 1, toda fluida.\n";
		return 1;
	}

	//------------------ Geometria da faixa ---------------------------------------------------------------------------//

	const double m     = sqrt ( ( double ) ( p * p + q * q ) );
	const double L     = N / m;					// periodo ao longo da normal
	const double theta = atan2 ( ( double ) q, ( double ) p ) * 180.0 / M_PI;
	const double tau   = parameters.tau_R;
	const double A_in  = parameters.A_fact;
	const double sigma_nominal = op.forca () ? A_in : 4.0 / 9.0 * A_in * tau;

	cout << "\nInterface plana inclinada -- arcabouco comum ( D3Q19, BGK, Latva-Kokko )" << endl;
	cout << "   operador          = " << op.nome << endl;
	cout << "   caixa             = " << N << " x " << N << endl;
	cout << "   ( p , q )         = ( " << p << " , " << q << " )   theta = " << theta << " graus" << endl;
	cout << "   distancia entre interfaces = " << 0.5 * L << endl;
	cout << "   sigma nominal     = " << sigma_nominal << endl;

	if ( 0.5 * L < 8.0 * ESPESSURA_INI )
		cout << "   ATENCAO: interfaces a menos de 8 espessuras uma da outra; aumente N." << endl;

	//====================== Allocates the memory used by the other arrays ===========================================//

	LATTICE lattice;

	const int size_to_alloc = n_fluid * nvel;

	lattice.inif_R     = new double[ size_to_alloc ];
	lattice.inif_R_new = new double[ size_to_alloc ];
	lattice.inif_B     = new double[ size_to_alloc ];
	lattice.inif_B_new = new double[ size_to_alloc ];
	lattice.inif_m     = new double[ size_to_alloc ];
	lattice.inif_m_new = new double[ size_to_alloc ];
	lattice.ini_stream = new int [ size_to_alloc ];
	lattice.ini_solid  = new bool[ size_to_alloc ];

	double *forca = nullptr;
	if ( op.forca () )
	{
		forca = new double[ 3 * n_fluid ]();
		lattice.ini_force = forca;
	}

	//====================== Define the lattice ======================================================================//

	lattice.c_i = new double[ nvel * dim ];
	def_lattice_d3q19 ( lattice );
	lattice.Q_i = new double[ nvel * dim * dim ];
	calc_Q ( lattice );

	//====================== Positions to propagate ==================================================================//

	def_dir_prop ( geometry, lattice );

	CAMPOS_F campos;
	if ( op.forca () ) aloca_campos_F ( campos, geometry, lattice );

	//====================== Geometria da faixa: condicao inicial ====================================================//
	//
	//  Equilibrio em repouso, rho = 1, fracao vermelha logistica na distancia a interface mais proxima.
	//  Os mediadores sao emitidos no inicio de cada passo, a partir das densidades: nao precisam de valor.

	{
		const double k = 2.0 * log ( 9.0 ) / ESPESSURA_INI;

		for ( int y = 0; y < N; y++ )
		for ( int x = 0; x < N; x++ )
		{
			long long s = ( ( long long ) p * x + ( long long ) q * y ) % N;
			if ( s < 0 ) s += N;

			//  interfaces em s = 0 ( = N ) e s = N/2;  distancia assinada a mais proxima, > 0 no vermelho

			const double dred = ( 2 * s < N )
			                  ? min ( ( double ) s, 0.5 * N - s ) / m
			                  : - min ( s - 0.5 * N, ( double ) ( N - s ) ) / m;

			const double cR  = 1.0 / ( 1.0 + exp ( - k * dred ) );
			const int pto = geometry.ini[ x + y * N ] - 1;

			for ( int i = 0; i < nvel; i++ )
			{
				lattice.inif_R[ ( size_t ) pto * nvel + i ] = lattice.w[i] * cR;
				lattice.inif_B[ ( size_t ) pto * nvel + i ] = lattice.w[i] * ( 1.0 - cR );
				lattice.inif_m[ ( size_t ) pto * nvel + i ] = 0.0;
			}
		}
	}

	//====================== Arquivos de saida =======================================================================//

	ofstream f_inc ( "inclinada.dat" );
	f_inc << setprecision ( 10 );
	f_inc << "# 1:passo  2:sigma_pop  3:int_grad  4:K  5:espessura_10_90  6:residuo  7:d0  8:|u|max  9:|u|rms"
	         "  10:|u_n|max  11:|u_t|max  12:massa_R  13:massa_B  14:qx  15:qy  16:soma_Fx  17:soma_Fy"
	         "  18:min_f  19:neg_total  20:neg_cor\n";

	INCL r;
	double s_sig = 0, s_um = 0, s_un = 0, s_ut = 0, s_K = 0, s_res = 0, s_ig = 0;
	int n_med = 0;
	long long neg_t = 0, neg_c = 0;
	double min_f = 1e300, d0_ini = NAN, d0_fim = NAN;

	string estado = "estavel";
	int passo_falha = -1;

	const auto t0 = chrono::steady_clock::now ();
	int step = 0;

	//============================== Looping principal ===============================================================//

	for ( step = 0; step < n_steps; step++ )
	{
		if ( step % 100 == 0 ) cout << "\rStep : " << step << flush;

		//------------------ Emissao e propagacao dos mediadores -----------------------------------------------------//

		#pragma omp parallel for
		for ( int pto = 0; pto < n_fluid; pto++ )
		{
			const double rho_R = density ( lattice.inif_R + pto * nvel );
			const double rho_B = density ( lattice.inif_B + pto * nvel );

			emite_mediadores ( lattice.inif_m + pto * nvel, ( rho_R - rho_B ) / ( rho_R + rho_B ), lattice );

			propag_site_med ( lattice, parameters, pto );
		}

		//------------------ Troca dos mediadores --------------------------------------------------------------------//

		{ double *t = lattice.inif_m; lattice.inif_m = lattice.inif_m_new; lattice.inif_m_new = t; }

		//------------------ Forca de tensao interfacial ( classe F ) ------------------------------------------------//

		if ( op.forca () ) calcula_forca ( campos, op, A_in, lattice, forca );

		//------------------ Colisao e propagacao --------------------------------------------------------------------//

		#pragma omp parallel for
		for ( int pto = 0; pto < n_fluid; pto++ )
		{
			static const double zero[3] = { 0.0, 0.0, 0.0 };
			const double *F = op.forca () ? forca + 3 * pto : zero;

			colisao_comum ( lattice.inif_R + pto * nvel, lattice.inif_B + pto * nvel, lattice.inif_m + pto * nvel,
			                F, op, lattice, parameters );

			double mx = 0., my = 0., mz = 0.;
			propag_site ( lattice, pto, mx, my, mz );
		}

		//------------------ Diagnostico -----------------------------------------------------------------------------//

		const bool ultimo = ( step == n_steps - 1 );

		if ( step % PASSOS_DIAG == 0 || ultimo )
		{
			const bool ok = mede_inclinada ( geometry, lattice, p, q, r, forca, ultimo );

			f_inc << step << " " << r.sigma_pop << " " << r.int_grad << " " << r.K << " "
			      << ( r.K > 0 ? 2.0 * log ( 9.0 ) / r.K : 0.0 ) << " " << r.res << " " << r.d0 << " "
			      << r.u_max << " " << r.u_rms << " " << r.un_max << " " << r.ut_max << " "
			      << r.massa_R << " " << r.massa_B << " " << r.qx << " " << r.qy << " "
			      << r.SFx << " " << r.SFy << " " << r.min_f << " " << r.neg_total << " " << r.neg_cor << "\n";
			f_inc.flush ();

			if ( step == 0 ) d0_ini = r.d0;
			d0_fim = r.d0;

			neg_t = max ( neg_t, r.neg_total );  neg_c = max ( neg_c, r.neg_cor );  min_f = min ( min_f, r.min_f );

			if ( step >= n_steps / 2 )
			{
				s_sig += r.sigma_pop;  s_um += r.u_max;  s_un += r.un_max;  s_ut += r.ut_max;
				s_K += r.K;  s_res += r.res;  s_ig += r.int_grad;
				n_med++;
			}

			cout << "\rstep = " << step << "   sigma_pop = " << r.sigma_pop << "   K = " << r.K
			     << "   |u|max = " << r.u_max << "   neg = " << r.neg_total << "        " << endl;

			if ( ! ok || r.u_max > U_LIMITE )
			{
				estado = "instavel";  passo_falha = step;
				cout << "\nSimulacao instavel no passo " << step << endl;
				break;
			}
		}

		//------------------ Atualiza ( fnovo => f ) -----------------------------------------------------------------//

		{
			double *t = lattice.inif_R; lattice.inif_R = lattice.inif_R_new; lattice.inif_R_new = t;
			t = lattice.inif_B;         lattice.inif_B = lattice.inif_B_new; lattice.inif_B_new = t;
		}
	}

	//====================== Resumo ==================================================================================//

	const double elapsed = chrono::duration<double> ( chrono::steady_clock::now () - t0 ).count ();
	const double mlups   = ( elapsed > 0 ) ? ( double ) n_fluid * step / ( elapsed * 1e6 ) : 0.0;
	auto med = [&] ( double s ) { return n_med ? s / n_med : NAN; };

	ofstream f_res ( "resumo.dat" );
	f_res << setprecision ( 10 );
	f_res << "# operador chi xi tau beta fat_R_B N p q theta sigma_nominal estado passo_falha passos "
	         "sigma_pop_med int_grad_med K_med residuo_med umax_med un_max_med ut_max_med umax/sigma_nominal "
	         "d0_ini d0_fim neg_total_max neg_cor_max min_f MLUPS\n";
	f_res << op.nome << " " << op.chi << " " << op.xi << " " << tau << " " << parameters.recoll << " " << A_in << " "
	      << N << " " << p << " " << q << " " << theta << " " << sigma_nominal << " "
	      << estado << " " << passo_falha << " " << step << " "
	      << med ( s_sig ) << " " << med ( s_ig ) << " " << med ( s_K ) << " " << med ( s_res ) << " "
	      << med ( s_um ) << " " << med ( s_un ) << " " << med ( s_ut ) << " " << med ( s_um ) / sigma_nominal << " "
	      << d0_ini << " " << d0_fim << " " << neg_t << " " << neg_c << " " << min_f << " " << mlups << "\n";

	cout << "\n\n================ Resultado ================" << endl;
	cout << "   estado            = " << estado << endl;
	cout << "   theta             = " << theta << " graus" << endl;
	cout << "   sigma_pop         = " << med ( s_sig ) << "   ( nominal " << sigma_nominal << ";  ~0 nas forcas )" << endl;
	cout << "   int |grad|        = " << med ( s_ig ) << "   ( deve ser 2 )" << endl;
	cout << "   K / espessura     = " << med ( s_K ) << " / " << 2.0 * log ( 9.0 ) / med ( s_K ) << endl;
	cout << "   |u|max            = " << med ( s_um ) << "   ( |u_n| " << med ( s_un ) << ", |u_t| " << med ( s_ut ) << " )" << endl;
	cout << "   d0 inicial / final= " << d0_ini << " / " << d0_fim << endl;
	cout << "   MLUPS             = " << mlups << endl;

	return 0;
}
