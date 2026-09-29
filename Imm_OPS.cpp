//=============================== Lattice Boltzmann Method - Operadores de tensao interfacial ( gota ) ===============//
//
//  Programa unico do artigo "Interfacial Tension Paper".  Uma gota cilindrica de fluido vermelho em fluido
//  azul, caixa periodica 2D ( nz = 1 ), parada ( v0 = 0 ) ou em translacao uniforme ( v0 > 0 ).  O que
//  varia entre execucoes e SO a forma de impor a tensao interfacial, lida de operador.txt:
//
//      P1          perturbacao, Gunstensen/Rothman ponderada          ( lambda = 0 )
//      P3 chi      perturbacao, familia de Liu-Valocchi-Kang          ( lambda = 1 )
//      F1          forca CSF                    ( Guo )
//      F2          forca div T  ( SHC 2010 )    ( Guo )
//      F3 xi       forca mu grad rho^N          ( Guo )
//
//  Todo o resto e fixo e comum: D3Q19, BGK de tempo unico ( tau = tau_R, exige tau_R = tau_B ), equilibrio
//  classico, recoloracao de Latva-Kokko com beta = recoll, gradiente isotropico dos mediadores.  Ver
//  Operadores_tensao.cpp para as formulas.
//
//  Entradas:
//      data_in.txt     formato de sempre.  fat_R_B = A ( classe P ) ou sigma ( classe F ).
//      operador.txt    uma linha, ex.: "P3 0.5" , "F1" , "F3 1.569988"
//      velocity.txt    v0 ( opcional; ausente => 0 )
//      meio.vtk        caixa nx x ny x 1 toda fluida
//
//  Saidas:
//      laplace.dat     as mesmas 22 colunas dos programas de With Mediators
//      galilean.dat    diagnostico no referencial da gota, incluindo contagem de populacoes negativas
//      iso.dat         modos de Fourier do contorno rho^N = 0:  R0, |a2|/R0, |a4|/R0, |a8|/R0, rmin, rmax
//      resumo.dat      uma linha: operador, parametros, estado final ( estavel / instavel ) e medias
//
//  Compilacao ( ver Makefile ):
//      g++ -std=c++17 -O3 -march=native -fopenmp -I<SimBoltz_Functions> Imm_OPS.cpp -o Imm_OPS
//
//====================================================================================================================//


//====================================================================================================================//
//   O P C O E S   D O   C A S O                                                                                      //
//====================================================================================================================//

const double RAIO_GOTA     = 0.0;		// 0 => R = nx / 4
const double ESPESSURA_INI = 3.45;		// largura 10-90 % do perfil de arranque ( a de equilibrio em beta = 0.8 )
const double FRAC_INTERNA  = 0.5;
const double FRAC_EXTERNA  = 2.0;

const int    PASSOS_DIAG   = 500;		// intervalo do diagnostico
const double U_LIMITE      = 0.4;		// acima disto: instavel

#define GRAVA_CAMPOS_FINAIS 1			// rho_R, rho_B e velocidade no ultimo passo

//====================================================================================================================//


#define nvel 	19
#define dim 	3

#include "Definitions.cpp"
#include "LBM_functions.cpp"
#include "Collision_operators.cpp"
#include "Boundary_conditions.cpp"
#include "Initial_conditions.cpp"
#include "Other_functions.cpp"
#include "Galilean_diagnostics.cpp"
#include "Operadores_tensao.cpp"


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

	double v0 = 0.0;
	{
		ifstream vin ( "velocity.txt" );
		if ( vin && ! ( vin >> v0 ) ) v0 = 0.0;
		if ( ! std::isfinite ( v0 ) || v0 < 0.0 ) { cerr << "velocity.txt invalido.\n"; return 1; }
	}

	if ( fabs ( parameters.tau_R - parameters.tau_B ) > 1e-14 )
	{
		cerr << "\nEste programa usa BGK de tempo UNICO: tau_R e tau_B tem de ser iguais.\n";
		return 1;
	}

	//=========================== Initializing variables =============================================================//

	const int n_steps     = parameters.n_steps;
	const int n_files     = parameters.n_files;
	const int passos_diag = PASSOS_DIAG;

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

	const int nx = geometry.nx, ny = geometry.ny, nz = geometry.nz;

	geometry.ini   = new int[ nx * ny * nz ];
	geometry.fluid = read_geo ( geometry.file, geometry.ini, 0, 0 );
	geometry.phi   = ( double ) geometry.fluid / ( double ) ( nx * ny * nz );

	const int n_fluid = geometry.fluid;

	if ( nz != 1 || n_fluid != nx * ny )
	{
		cerr << "\nEste programa requer caixa 2D periodica toda fluida ( nz = 1 ).\n";
		return 1;
	}

	//====================== Forma da gota e raio inicial ============================================================//

	const bool   cilindro = true;
	const double raio_ini = ( RAIO_GOTA > 0.0 ) ? RAIO_GOTA : 0.25 * ( double ) nx;
	const int    x0_ini = nx / 2, y0_ini = ny / 2, z0_ini = 0;

	//------------------ Resumo do caso ------------------------------------------------------------------------------//

	const double tau  = parameters.tau_R;
	const double A_in = parameters.A_fact;
	const double sigma_nominal = op.forca () ? A_in : 4.0 / 9.0 * A_in * tau;

	cout << "\nOperadores de tensao interfacial -- arcabouco comum ( D3Q19, BGK, Latva-Kokko )" << endl;
	cout << "   operador          = " << op.nome;
	if ( op.tipo == OP_P3 ) cout << "   chi = " << op.chi;
	if ( op.tipo == OP_F3 ) cout << "   xi = " << op.xi;
	cout << endl;
	cout << "   caixa             = " << nx << " x " << ny << " x " << nz << endl;
	cout << "   raio inicial      = " << raio_ini << endl;
	cout << "   tau               = " << tau << "   ( nu = " << ( tau - 0.5 ) / 3.0 << " )" << endl;
	cout << "   beta              = " << parameters.recoll << endl;
	cout << "   fat_R_B           = " << A_in << ( op.forca () ? "   ( = sigma )" : "   ( = A )" ) << endl;
	cout << "   sigma nominal     = " << sigma_nominal << endl;
	cout << "   v0                = " << v0 << endl;

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

	//  A forca por sitio existe so na classe F; nas perturbacoes ini_force = nullptr desliga a correcao da
	//  meia-forca nos diagnosticos da biblioteca.

	double *forca = nullptr;

	if ( op.forca () )
	{
		forca = new double[ 3 * n_fluid ];
		for ( int k = 0; k < 3 * n_fluid; k++ ) forca[k] = 0.0;
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

	//====================== Initializing the distribution function ==================================================//

	const double k_perfil = 2.0 * log ( 9.0 ) / ESPESSURA_INI;

	initial_conditions_drop ( raio_ini, k_perfil, 0.0, x0_ini, y0_ini, z0_ini,
	                          v0, 0., 0., cilindro, geometry, lattice, parameters, true, EQ_CLASSICO );

	//====================== Arquivos de saida =======================================================================//

	ofstream f_lap ( "laplace.dat" );
	f_lap << setprecision ( 10 );
	f_lap << "# 1:passo  2:raio_area  3:raio_grad  4:raio_linha  5:p_in  6:p_out  7:delta_p"
	         "  8:sigma  9:|u-v0|max  10:|u-v0|rms  11:Ca  12:massa_R  13:massa_B  14:|q|"
	         "  15:sd_p_in  16:sd_p_out  17:sigma_2pontos"
	         "  18:K_perfil  19:r_perfil  20:espessura_10_90  21:residuo_ajuste  22:n_camadas\n";

	ofstream f_iso ( "iso.dat" );
	f_iso << setprecision ( 10 );
	f_iso << "# 1:passo  2:R0  3:|a2|/R0  4:|a4|/R0  5:|a8|/R0  6:rmin  7:rmax\n";

	ofstream f_gi ( "galilean.dat" );
	gi_header ( f_gi );

	const double nu = lattice.c_s2 * ( tau - 0.5 );

	DROP drop;
	drop.cilindro  = cilindro;
	drop.frac_in   = FRAC_INTERNA;
	drop.frac_out  = FRAC_EXTERNA;
	drop.espessura = ESPESSURA_INI;

	GI_DATA gi;

	//  Medias da metade final da serie ( o que os estudos de With Mediators reportam )

	double s_sigma = 0, s_umax = 0, s_urms = 0, s_a4 = 0, s_a8 = 0, s_K = 0;
	int    n_med = 0;
	long long neg_total_max = 0, neg_cor_max = 0;
	double min_f_global = 1e30;

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

		{
			double *t = lattice.inif_m; lattice.inif_m = lattice.inif_m_new; lattice.inif_m_new = t;
		}

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
		//
		//  Populacoes pos-colisao ( a troca vem depois ); a meia-forca de Guo e descontada pelos diagnosticos.

		const bool medir = ( step % passos_diag == 0 ) || ( step == n_steps - 1 );

		if ( medir )
		{
			//--------------- A medida, na biblioteca ---------------------------------------------------------------//

			const bool ok = calc_galilean ( geometry, lattice, drop, gi, v0, step, x0_ini, y0_ini );

			if ( ok ) calc_perfil_logistico ( geometry, lattice, drop );

			MODOS md;
			if ( ok ) md = modos_interface ( geometry, lattice, drop.x0, drop.y0, drop.raio );

			//--------------- Saida ---------------------------------------------------------------------------------//

			if ( ok )
			{
				const double Ca = ( drop.sigma > 0.0 ) ? nu * drop.u_max / drop.sigma : 0.0;
				const double q  = sqrt ( drop.qx * drop.qx + drop.qy * drop.qy + drop.qz * drop.qz );

				f_lap << step << " " << drop.raio << " " << drop.raio_grad << " " << drop.raio_linha << " "
				      << drop.p_in << " " << drop.p_out << " " << drop.delta_p << " " << drop.sigma << " "
				      << drop.u_max << " " << drop.u_rms << " " << Ca << " "
				      << drop.massa_R << " " << drop.massa_B << " " << q << " "
				      << drop.sd_in << " " << drop.sd_out << " " << drop.sigma_linha << " "
				      << drop.K_perfil << " " << drop.r_perfil << " " << drop.esp_perfil << " "
				      << drop.res_perfil << " " << drop.n_perfil << "\n";
				f_lap.flush ();

				gi_write ( f_gi, step, v0, lattice.c_s2, drop, gi );

				f_iso << step << " " << md.R0 << " " << md.a2 << " " << md.a4 << " " << md.a8 << " "
				      << md.rmin << " " << md.rmax << "\n";
				f_iso.flush ();

				neg_total_max = max ( neg_total_max, gi.negative_total );
				neg_cor_max   = max ( neg_cor_max,   gi.negative_color );
				min_f_global  = min ( min_f_global,  gi.min_f );

				if ( step >= n_steps / 2 )
				{
					s_sigma += drop.sigma;  s_umax += drop.u_max;  s_urms += drop.u_rms;
					s_a4 += md.a4;  s_a8 += md.a8;  s_K += drop.K_perfil;
					n_med++;
				}

				cout << "\rstep = " << step << "   R = " << drop.raio << "   sigma = " << drop.sigma
				     << "   |u|max = " << drop.u_max << "   a4 = " << md.a4
				     << "   neg = " << gi.negative_total << "        " << endl;
			}

			//--------------- Avisos e paradas ----------------------------------------------------------------------//

			if ( ! ok || drop.n_ruim > 0 || drop.u_max > U_LIMITE )
			{
				estado = "instavel";
				passo_falha = step;
				cout << "\nSimulacao instavel no passo " << step << endl;
				break;
			}
		}

		//------------------ Grava os campos -------------------------------------------------------------------------//

#if GRAVA_CAMPOS_FINAIS
		if ( step == n_steps - 1 || ( n_files > 0 && step % max ( 1, n_steps / n_files ) == 0 && step > 0 ) )
		{
			rec_density ( "rho_R", geometry, lattice.inif_R, ( unsigned int ) step );
			rec_density ( "rho_B", geometry, lattice.inif_B, ( unsigned int ) step );
			rec_velocity ( geometry, lattice, ( unsigned int ) step );
		}
#endif

		//------------------ Atualiza ( fnovo => f ) -----------------------------------------------------------------//

		{
			double *t = lattice.inif_R; lattice.inif_R = lattice.inif_R_new; lattice.inif_R_new = t;
			t = lattice.inif_B;         lattice.inif_B = lattice.inif_B_new; lattice.inif_B_new = t;
		}
	}

	//====================== Resumo ==================================================================================//

	const double elapsed = chrono::duration<double> ( chrono::steady_clock::now () - t0 ).count ();
	const double mlups   = ( elapsed > 0 ) ? ( double ) n_fluid * step / ( elapsed * 1e6 ) : 0.0;

	const double m_sigma = n_med ? s_sigma / n_med : NAN;
	const double m_umax  = n_med ? s_umax  / n_med : NAN;

	ofstream f_res ( "resumo.dat" );
	f_res << setprecision ( 10 );
	f_res << "# operador chi xi tau beta fat_R_B v0 nx R_ini sigma_nominal estado passo_falha passos "
	         "sigma_med umax_med urms_med umax/sigma Ca_s=nu*umax/sigma a4_med a8_med K_med "
	         "neg_total_max neg_cor_max min_f MLUPS\n";
	f_res << op.nome << " " << op.chi << " " << op.xi << " " << tau << " " << parameters.recoll << " "
	      << A_in << " " << v0 << " " << nx << " " << raio_ini << " " << sigma_nominal << " "
	      << estado << " " << passo_falha << " " << step << " "
	      << m_sigma << " " << m_umax << " " << ( n_med ? s_urms / n_med : NAN ) << " "
	      << m_umax / m_sigma << " " << nu * m_umax / m_sigma << " "
	      << ( n_med ? s_a4 / n_med : NAN ) << " " << ( n_med ? s_a8 / n_med : NAN ) << " "
	      << ( n_med ? s_K / n_med : NAN ) << " "
	      << neg_total_max << " " << neg_cor_max << " " << min_f_global << " " << mlups << "\n";

	cout << "\n\n================ Resultado ================" << endl;
	cout << "   estado            = " << estado << endl;
	cout << "   sigma ( media )   = " << m_sigma << "   ( nominal " << sigma_nominal << " )" << endl;
	cout << "   |u|max ( media )  = " << m_umax << "   |u|max / sigma = " << m_umax / m_sigma << endl;
	cout << "   populacoes < 0    = " << neg_total_max << " ( total )  " << neg_cor_max << " ( cor )" << endl;
	cout << "   MLUPS             = " << mlups << endl;

	return 0;
}
