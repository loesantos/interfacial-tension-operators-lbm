//------ Paralelização: OpenACC (GPU) e/ou OpenMP (CPU) ------------------------------------------//

#ifdef _OPENACC
#include <openacc.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#else
//	Permite compilar sem OpenMP (execução serial ou apenas OpenACC)
static int  omp_get_max_threads ( void )  { return 1; }
static void omp_set_num_threads ( int n ) { ( void ) n; }
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <array>
#include <vector>
#include <chrono>

using namespace std;

//=============================== Define Struct Geometry =============================================================//

struct GEOMETRY 
{	
	string file;		// Name of the geometry file
	
	double ftesc;   	// Dimensions of the pixels

	int nx, ny, nz;		// Dimensions (pixels) of the geometry 
	
	double phi;			// Porosity	

	int fluid;			// Number of fluid points

	int* ini;			// First adress of the geometry (pointer) 
};


//=============================== Define Struct Lattice ==============================================================//

struct LATTICE 
{		
	double *inif = nullptr;				// Distribution functions (initial adress)
	double *inif_R = nullptr;			// Distribution functions (initial adress - Red fluid)
	double *inif_B = nullptr;			// Distribution functions (initial adress - Blue fluid)
	double *inif_m = nullptr;			// Distribution functions (initial adress - mediators)
	
	double *inif_new = nullptr;			// Distribution functions (initial adress) t + delta_t
	double *inif_R_new = nullptr;		// Distribution functions (initial adress - Red fluid) t + delta_t
	double *inif_B_new = nullptr;		// Distribution functions (initial adress - Blue fluid) t + delta_t
	double *inif_m_new = nullptr;		// Distribution functions (initial adress - mediators) t + delta_t
	
	double *ini_psi = nullptr;			// Points to a scalar (phase) field
	
	//------ Campos auxiliares do modelo Montessori-Hegele-Lauricella (um valor por sitio fluido) --//
	//
	//  ini_phi  : campo de fase phi = soma de g_i
	//  ini_grad : gradiente de phi, dim componentes por sitio
	//
	//  Sao escalares por sitio (4 doubles), e nao populacoes: o gradiente e a curvatura sao obtidos
	//  lendo os vizinhos pelo mapa de propagacao ini_stream, o que custa muito menos memoria do que
	//  carregar mais um par de populacoes.
	
	//------ Campo de tensao interfacial do modelo de Spencer-Halliday-Care -----------------------//
	//
	//  Seis componentes independentes de Dsigma'_ab por sitio fluido, na ordem
	//  xx, yy, zz, xy, xz, yz.  So e alocado por quem usa coll_SHC_2010_guo, que precisa da
	//  DIVERGENCIA do tensor -- e divergencia nao e local, tem de ler os vizinhos.

	double *ini_stress = nullptr;		// Dsigma'_ab (Spencer-Halliday-Care 2010, eq. 5)

	//------ Forca de corpo por sitio, para a correcao da meia-forca de Guo ------------------------//
	//
	//  Tres componentes por sitio fluido.  So precisa existir em programas com FORCAMENTO DE GUO,
	//  e serve a um proposito: dizer as medidas ( calc_drop_radius, calc_laplace_tension ) qual
	//  forca subtrair do momento cru para chegar a velocidade fisica.  nullptr = sem correcao.

	double *ini_force = nullptr;		// F_a por sitio ( esquema de Guo )

	double *ini_phi = nullptr;			// Phase field (MHL model)
	
	double *ini_grad = nullptr;			// Gradient of the phase field (MHL model)
	
	int *ini_stream = nullptr;					// Points to the propagation site
	
	bool *ini_solid	= nullptr;			// Indicates propagation to a solid site
	
	double *ini_mom_x = nullptr;		// Momentum x lost during the propagation
	double *ini_mom_y = nullptr;		// Momentum y lost during the propagation
	double *ini_mom_z = nullptr;		// Momentum z lost during the propagation
		
	double *c_i = nullptr;						// velocity vectors c_i; 
	
	double w[nvel];						// Weights of the lattice
	
	double c_s2;						// square of the sound velocity
	
	double one_over_c_s2;				// Inverse of the square of the sound velocity
	
	double *Q_i = nullptr;						// Tensor Q_i_alpha_beta = ci_alpha * ci_bet - delta_alpha_beta
	
	double *K2_i = nullptr;						// Tensor K2_i_alpha_beta = w_i/(2cs^4)Q_i
};


//=============================== Define Struct Collision ============================================================//

struct PARAMETERS 
{		
	double tau;							// Relaxation time
	
	double visc;						// kinematic viscosity
	
	double tau_R;						// Relaxation time (Red fluid)
	
	double visc_R;						// kinematic viscosity (Red fluid)
	
	double tau_B;						// Relaxation time (Blue fluid)

	double visc_B;						// kinematic viscosity (Blue fluid)
	
	double tau_m;						// Relaxation time (mixed fluid)
	
	double rho_ini;						// Initial density
	
	double rho_ini_R;					// Initial density (Red fluid)
	
	double rho_ini_B;					// Initial density (Blue fluid)
	
	double tau_nh;						// Relaxation time (non-hidrodynamic)
	
	double A_fact;						// Interaction factor
	
	double A_fact_R;					// Interaction factor (Red fluid)
	
	double A_fact_B;					// Interaction factor (Blue fluid)
	
	double recoll;						// Recolloring factor
	
	double wett_R;						// Wettability (Red fluid)

	//------ Modelo color-gradient de Wen, Li, Yu & Luo (2019) ------------------------------------//
	//
	//  alfa_R e alfa_B fixam as velocidades do som de cada fluido: cs_k^2 = ( 1 - alpha_k ) / 2.
	//  O equilibrio de pressao entre as fases exige rho_R^in/rho_B^in = (1-alpha_B)/(1-alpha_R),
	//  de modo que so alpha_R e livre.  delta_tau e a largura da emenda parabolica de tau_v
	//  atraves da interface ( eq. 59 do artigo; o valor usual e 0.98 ).

	double alfa_R = 1. / 3.;			// alpha_R  ( 1/3 => cs^2 = 1/3, caso de densidades iguais )
	double alfa_B = 1. / 3.;			// alpha_B  ( derivado de alfa_R e da razao de densidades )
	double delta_tau = 0.98;			// delta da suavizacao de tau_v atraves da interface
	double lyl_melhorado = 1.0;			// 1 = modelo melhorado ( termo de ordem alta + correcao )
										// 0 = modelo original ( equilibrio da eq. 6, sem correcao )

	//------ Modelo de Liu, Valocchi & Kang (2012) ------------------------------------------------//
	//
	//  chi e o parametro livre da familia de B_i da eq. (29).  A tensao interfacial NAO depende
	//  dele -- todos os chi satisfazem as mesmas tres condicoes da eq. (28) --, so a parte da
	//  perturbacao que nao contribui para o tensor S.  O artigo usa chi = 2; chi = 1/2 reproduz
	//  exatamente o interf_tension_RPL() desta biblioteca.  alfa_R e alfa_B sao os mesmos de cima.

	double chi_LVK = 2.0;				// chi da eq. (29) de Liu-Valocchi-Kang

	//------ Modelos de Spencer, Halliday & Care --------------------------------------------------//
	//
	//  QUAL RECOLORACAO.  O artigo de 2010 nao escreve a regra na secao de dois fluidos -- diz so
	//  "we therefore employ identical segregation rules to those analyzed previously [14]", e [14] e
	//  Halliday, Hollis & Care, Phys. Rev. E 76, 026708 (2007).  Ela so aparece explicita na eq. (69),
	//  a versao de N componentes, e la e  t_i ( rho_n rho_m / rho^2 ) c_i . n  --  SEM dividir por
	//  | c_i | .  Latva-Kokko e D'Ortona aparecem no artigo apenas como linhagem, nao como a regra
	//  usada.  Esta arvore sempre usou a de Latva-Kokko, que divide.
	//
	//  A diferenca nao e cosmetica: no mesmo beta o fluxo de segregacao de Latva-Kokko vale
	//  0.90237 ( D2Q9 ), 0.80474 ( D3Q19 ) ou 0.82286 ( D3Q27 ) do de Halliday, e a interface fica
	//  proporcionalmente mais grossa.  Alem disso a forma de Halliday degenera exatamente de nz = 1
	//  para o D2Q9 e a de Latva-Kokko nao.
	//
	//  0 ( padrao ) preserva tudo o que ja foi medido nesta arvore.

	int recoll_halliday = 0;			// 0 = Latva-Kokko ( recolloring )   1 = Halliday, eq. (69)

	//------ Modelo de Saito et al. (2023) --------------------------------------------------------//
	//
	//  Aqui a velocidade do som de cada fluido entra DIRETO, e nao por um alpha: a pressao e
	//  p = rho_R cs2_R + rho_B cs2_B ( eq. 12 ), e a razao de densidades de equilibrio e
	//  gamma = cs2_B / cs2_R ( eq. 13 ).  Com cs2_R = cs2_B = 1/3 recai no caso classico.
	//
	//  Na convencao alpha do artigo, cs_k^2 = xi ( 1 - alpha_k ) com xi = 9/19 no D3Q27 --
	//  alpha_k = 8/27 da cs_k^2 = 1/3.  Guardar cs^2 evita ter de carregar o xi de cada rede.

	double cs2_R = 1. / 3.;				// ( cs^R )^2
	double cs2_B = 1. / 3.;				// ( cs^B )^2

	//  LIMITADOR DE POSITIVIDADE DA RECOLORACAO.  A recoloracao de Halliday ( eqs. 18-19 ) soma a
	//  cada cor  d_i = beta ( rho_r rho_b / rho^2 ) p w_i ( c_i . n ) / cs^2 .  Para f_b_i >= 0 e
	//  preciso  f_i / ( rho w_i )  >=  beta ( rho_r / rho ) ( p / rho cs^2 ) ( c_i . n ) .  Em
	//  repouso o lado esquerdo vale 1 e a condicao e folgada; numa gota que se move com velocidade
	//  u ele cai como  1 - 3 u + 3 u^2  nos elos contra o movimento, e a partir de u ~ 0.1 a
	//  recoloracao passa a produzir populacoes de cor NEGATIVAS.  Dai vem phi fora de [ -1, 1 ],
	//  | grad phi | de ruido no seio do fluido, uma normal n aleatoria em todo o dominio e, por
	//  fim, uma forca CSF espuria que freia e deforma a gota.
	//
	//  recoll_positivo = 1 multiplica d_i por um UNICO lambda <= 1 por sitio -- o maior que mantem
	//  as duas cores nao negativas.  Como e um fator global do sitio, soma_i d_i continua zero e a
	//  massa de cada cor segue exatamente conservada.  Em repouso lambda = 1 em todo o dominio e o
	//  resultado e BIT A BIT o mesmo de recoll_positivo = 0 -- ver SAITO.md, secao 8.
	//
	//  recoll_positivo = 0 ( padrao ) reproduz a recoloracao crua do artigo.

	//  PADRAO 0: fiel ao artigo.  Cada programa escolhe explicitamente ( RECOLL_POSITIVO ), de modo
	//  que nenhum resultado mude sem que alguem tenha pedido.

	int recoll_positivo = 0;			// 1 = limita a recoloracao para nao gerar cor negativa

	//  QUAL RECOLORACAO no modelo de Saito.  0 = Halliday, eqs. (18)-(19), a do artigo.
	//  1 = Latva-Kokko, que divide por | c_i |.  O limitador acima vale para as duas.
	//  No mesmo beta a de Latva-Kokko segrega 0.8229 do que a de Halliday no D3Q27 ( 0.8047 no
	//  D3Q19 ), entao os beta NAO sao comparaveis entre as duas -- ver PLANO_DO_ARTIGO.md.

	int recoll_saito = 0;				// 0 = Halliday ( artigo )   1 = Latva-Kokko
	
	
	int n_steps;						// Number of time steps
	
	int n_files;						// Number of files 
	
	int n_threads;						// Number of threads
};


//=============================== Define Struct Drop =================================================================//
//
//  Medidas de uma gota isolada ( calc_drop_radius e calc_laplace_tension, em Other_functions.cpp ).
//
//  Os quatro primeiros campos sao ENTRADAS, ajustadas pelo programa antes da primeira chamada; o
//  resto e saida.  Os valores por omissao servem a uma gota 2D com interface fina.
//
//====================================================================================================================//

struct DROP
{
	//------ Entradas -----------------------------------------------------------------------------//

	bool   cilindro  = true;			// true  : gota cilindrica de eixo z ( 2D )  =>  sigma = R dp
										// false : gota esferica ( 3D )              =>  sigma = R dp / 2

	double frac_in   = 0.5;				// p_in  : media da pressao sobre  r <  frac_in  * raio
	double frac_out  = 2.0;				// p_out : media da pressao sobre  r >= frac_out * raio

	double espessura = 0.0;				// largura 10%-90% da interface; afasta a regiao externa dela

	//  Velocidade do som de cada fluido, para modelos em que a pressao NAO e c_s^2 rho.  No de
	//  Liu-Valocchi-Kang ( e no de Wen-Li-Yu-Luo ) cada fluido tem a sua,  cs_k^2 = ( 1 - alpha_k )/2,
	//  e a pressao do sitio e  p = rho_R cs_R^2 + rho_B cs_B^2 .  Negativo ( o padrao ) manda usar
	//  lattice.c_s2 para os dois, que e o que todos os outros modelos da arvore querem.

	double cs2_R = -1.0;				// cs^2 do fluido vermelho ( < 0 => lattice.c_s2 )
	double cs2_B = -1.0;				// cs^2 do fluido azul

	//------ Curva interfacial de Latva-Kokko & Rothman ( 2005 ), eqs. (15)-(17) -------------------//
	//
	//  A recoloracao daqueles autores leva o perfil de cor a uma logistica,
	//
	//      d phi / ds = K phi ( 1 - phi )       =>       phi( s ) = 1 / ( 1 + exp[ K ( s - s_0 ) ] )
	//
	//  com s a distancia radial e K uma constante de comprimento INVERSO, proporcional a beta.
	//  calc_perfil_logistico() ajusta  ln[ phi / ( 1 - phi ) ]  contra r por minimos quadrados: a
	//  inclinacao e -K e o coeficiente linear da a posicao da interface.

	//------ Interface PLANA ( calc_flat_tension ) ------------------------------------------------//
	//
	//  A tensao interfacial pela definicao mecanica,  sigma = integral ( P_N - P_T ) dx , com a
	//  interface normal a x.  E a medida complementar a lei de Laplace: sem curvatura, e portanto
	//  sem correcao de Tolman e sem extrapolacao -- o valor sai direto.
	//
	//  plana = true faz calc_perfil_logistico() ajustar o perfil ao longo de x em vez de radialmente.

	bool   plana       = false;			// geometria da medida do perfil
	double sigma_pop   = 0.0;			// sigma pela anisotropia das POPULACOES
	double integral_grad = 0.0;			// integral de | grad( fase ) | atraves da interface
	double x_interface = 0.0;			// posicao da interface, pelo ajuste logistico
	double min_f_R = 0.0, min_f_B = 0.0;	// menores populacoes de cada cor no diagnostico plano
	double min_f_total = 0.0;			// menor f_R + f_B; limite conservador de admissibilidade

	double K_perfil  = 0.0;				// K medido ( inverso da largura da interface )
	double r_perfil  = 0.0;				// raio onde phi = 1/2, pelo mesmo ajuste
	double res_perfil = 0.0;			// maior residuo de ln[phi/(1-phi)] no ajuste
	double esp_perfil = 0.0;			// largura 10%-90% implicada:  2 ln(9) / K
	int    n_perfil  = 0;				// numero de camadas radiais usadas

	//------ Centro e tamanho ---------------------------------------------------------------------//

	double x0 = 0.0, y0 = 0.0, z0 = 0.0;	// centro, pelo centroide circular de phi = rho_R / rho

	double area = 0.0;					// soma_sitios phi  ( area * nz no caso 2D, volume no 3D )

	double raio = 0.0;					// raio equimolar, a partir de 'area'
	double raio_grad = 0.0;				// primeiro momento radial de rho_R rho_B ( estimador independente )
	double raio_linha = 0.0;			// estimador de uma linha so ( compatibilidade; e fragil )

	//------ Conservacao --------------------------------------------------------------------------//

	double massa_R = 0.0, massa_B = 0.0;

	double qx = 0.0, qy = 0.0, qz = 0.0;

	//------ Correntes espurias -------------------------------------------------------------------//

	double u_max = 0.0, u_rms = 0.0;

	//------ Pressao e tensao interfacial ---------------------------------------------------------//

	double p_in = 0.0, p_out = 0.0;		// medias sobre as duas regioes
	double sd_in = 0.0, sd_out = 0.0;	// desvios padrao: a barra de erro da medida

	double delta_p = 0.0;
	double sigma = 0.0;

	double p_centro = 0.0, p_canto = 0.0;	// medida de dois pontos, do jeito antigo
	double sigma_linha = 0.0;				// sigma que sairia dela, com raio_linha

	int n_in = 0, n_out = 0;			// sitios em cada regiao de media

	//------ Diagnostico --------------------------------------------------------------------------//

	int n_fluid = 0;					// sitios de fluido varridos
	int n_ruim = 0;						// sitios com densidade ou velocidade nao finita
};


//====================================================================================================================//
//
//   Medida de uma ONDA CAPILAR ( calc_onda_capilar, em Other_functions.cpp ).
//
//   Geometria do ensaio: a interface e normal a x e a onda corre em y.  O vermelho ocupa x < x_int,
//   o azul x > x_int, a caixa e periodica em y e tem derivada nula nas duas faces x.  A interface
//   fica em
//
//       x_int( y, t ) = x_med( t ) + a( t ) cos( k y + fase )          k = 2 pi modo / ny
//
//   A amplitude modal a( t ) sai por projecao de Fourier da posicao medida da interface, e nao de
//   ler a altura numa coluna so: a projecao usa as ny colunas e o ruido cai com sqrt( ny ).
//
//   Teoria de referencia, dois fluidos semi-infinitos de mesma densidade e mesma viscosidade
//   cinematica, no limite de amortecimento fraco:
//
//       w_0^2 = sigma k^3 / ( rho_R + rho_B )              ( Lamb, onda capilar sem viscosidade )
//       Gamma = 2 nu k^2                                   ( taxa de decaimento da amplitude )
//       a( t ) = a_0 exp( - Gamma t ) cos( w t )           com  w^2 = w_0^2 - Gamma^2
//
//   Vale enquanto Gamma << w_0.  Fora disso a referencia correta e a solucao de Prosperetti (1981),
//   que nao e uma exponencial vezes cosseno.
//
//====================================================================================================================//

struct WAVE
{
	//------ Entradas -----------------------------------------------------------------------------//

	int    modo = 1;					// comprimentos de onda na caixa;  k = 2 pi modo / ny

	double cs2_R = -1.0;				// idem DROP: < 0 manda usar lattice.c_s2 para os dois
	double cs2_B = -1.0;

	//------ Posicao da interface -----------------------------------------------------------------//

	double k = 0.0;						// numero de onda usado na projecao
	double x_med = 0.0;					// posicao media da interface ( deve ficar parada )

	double a_cos = 0.0, a_sin = 0.0;	// componentes de Fourier do deslocamento, no modo excitado
	double amp = 0.0;					// amplitude   sqrt( a_cos^2 + a_sin^2 )
	double fase = 0.0;					// fase        atan2( a_sin, a_cos ), em radianos

	double amp_res = 0.0;				// o que sobra fora do modo excitado: harmonicos e ruido.
										// Se chegar perto de 'amp', a onda deixou de ser linear.

	int    n_linhas = 0;				// colunas y em que a interface foi localizada
	int    n_faltou = 0;				// colunas em que NAO houve cruzamento ( interface perdida )

	//------ Conservacao e estabilidade -----------------------------------------------------------//

	double massa_R = 0.0, massa_B = 0.0;
	double qx = 0.0, qy = 0.0, qz = 0.0;

	double u_max = 0.0, u_rms = 0.0;

	int    n_fluid = 0;
	int    n_ruim = 0;
};

//====================================================================================================================//


//=============================== Declarations (Functions_LBM) =======================================================//

//	Read the inicialization file
void read_data ( GEOMETRY&, PARAMETERS& );

//	Read the geometry file
int read_geo ( string, int*, int, int );

// 	Define a D3Q19 lattice
void def_lattice_d3q19 ( LATTICE& );

//	Define the directions used in the propagation step
void def_dir_prop ( GEOMETRY geometry, LATTICE& lattice );

//	Return the dot product
#pragma acc routine seq
double dot_product ( double*, double* );

//	Return the equilibrium distribution function 
#pragma acc routine seq
void dist_eq ( double*, double, double, double, double, LATTICE );

//	Return the equilibrium distribution function  (sixth order)
#pragma acc routine seq
void dist_eq_sixth ( double*, double, double, double, double, LATTICE );

//	Calculate density and velocities
#pragma acc routine seq
void calcula ( double*, double&, double&, double&, double&, LATTICE );

//  Compute momentum
#pragma acc routine seq
void momentum ( double*, double&, double&, double&, LATTICE );

//	Propagation step for one site
#pragma acc routine seq
void propag_site ( LATTICE, int, double&, double&, double& );

//	Emission of the mediators for one site  ( M_i = w_i fase / cs^2 )
#pragma acc routine seq
void emite_mediadores ( double*, double, LATTICE );

//	Propagation step for one site for mediators
#pragma acc routine seq
void propag_site_med ( LATTICE, PARAMETERS, int );

//	Record the velocity field (monophasic) - time is integer
void rec_velocity ( GEOMETRY, LATTICE, unsigned int );

//Record the density field (one of two fluids)
void rec_density ( string, GEOMETRY, double*, unsigned int );

//	Returns the density of a site 
#pragma acc routine seq
double density ( double* );

//	Calculates the force term (source)
#pragma acc routine seq
void source ( double, double, double, double, double, double, double, double, double*, LATTICE );

// Calculates the tensor Q_i = ci_alpha * ci_bet - delta_alpha_beta
void calc_Q ( LATTICE lattice );

//	Recolloring step (Latva-Koko)
#pragma acc routine seq
void recolloring ( double*, double*, double*, double*, double, double, double, double, double, LATTICE );


//============================== Declarations (Functions_Collision) ===================================================//

//	Collision step using the BGK model
#pragma acc routine seq
void coll_BGK ( double*, double, double, double, LATTICE, PARAMETERS );

//	Inverse central-moment transform, k_abc -> f
#pragma acc routine seq
void cm_saito_inv ( const double*, double, double, double, double* );

//	Generalized equilibrium of eq. (56), written through its velocity-independent CMs (57)-(63)
#pragma acc routine seq
void dist_eq_saito ( double*, double, double, double, double, double, LATTICE );

//	Modelo color-gradient 3D melhorado - Wen, Li, Yu & Luo, Phys. Rev. E 100, 023301 (2019)

//	Runtime check of the LYL model: M^-1 M = I, M f^eq = m^eq, and eqs. (20) and (23).  Host only.
#if nvel == 19
#endif

#pragma acc routine seq
void dist_eq_LYL ( double*, double, double, double, double, double, double, LATTICE );

//	Initial Condition of a Drop with a diffuse interface (immiscible fluids)
//	( radius, k_perfil, delta_rho, x0, y0, z0, vx, vy, vz, cilindro, ... )
//	Modo de semeadura do equilibrio em initial_conditions_drop.  EQ_ALFA vale 1, entao um
//	'true' passado por programas antigos continua significando o equilibrio de alpha variavel.

#define EQ_CLASSICO  0
#define EQ_ALFA      1
#define EQ_SAITO     2

void initial_conditions_drop ( double, double, double, int, int, int, double, double, double, bool,
							GEOMETRY, LATTICE, PARAMETERS, bool fase_simetrica = false,
							int eq_modo = EQ_CLASSICO );


//============================== Declarations (Functions_Diverses) ===================================================//

//	Round a number
int round_number ( double );

//	Returns the radius of an isolated drop (fills DROP: center, radii, masses, spurious currents)
double calc_drop_radius ( GEOMETRY, LATTICE, DROP& );

//	Returns the interfacial tension of an isolated drop by Laplace's law (fills DROP)
double calc_laplace_tension ( GEOMETRY, LATTICE, DROP& );

//	Fits the interfacial colour profile to the logistic of Latva-Kokko & Rothman, eq. (17)
double calc_perfil_logistico ( GEOMETRY, LATTICE, DROP& );
