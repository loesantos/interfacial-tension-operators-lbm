


//================================ Etapa de colisão usando BGK =======================================================//
//
//      Input: distribution function, lattice vectors, acceleration in the x,y,z directions
//      Output: pos-collisional distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void coll_BGK ( double *f, double acc_x, double acc_y, double acc_z, LATTICE lattice, PARAMETERS parameters )
{
	double tau = parameters.tau;

	double f_eq[nvel];

	double op_col[nvel];
	
	double S[nvel];

	double vx, vy, vz, rho;

	double one_over_tau = 1.0 / tau;

	calcula ( f, vx, vy, vz, rho, lattice );

	double vx_alt = vx + 0.5 * acc_x;
	double vy_alt = vy + 0.5 * acc_y;
	double vz_alt = vz + 0.5 * acc_z;
	
	double Fx = acc_x * rho;    
	double Fy = acc_y * rho;
	double Fz = acc_z * rho;
	
	source ( Fx, Fy, Fz, vx_alt, vy_alt, vz_alt, rho, tau, S, lattice );
	
	if ( nvel == 77 ) dist_eq_sixth ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
	
	else dist_eq ( f_eq, vx_alt, vy_alt, vz_alt, rho, lattice );
	
	for ( int i = 0; i < nvel; i++ )
	{
		op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau;

		f[i] = f[i] + op_col[i] + S[i];
	}
	
}

//====================================================================================================================//


//=============================== Matrizes de transformacao do modelo MRT ============================================//
//
//  As matrizes sao constantes: ficam em escopo de arquivo e sao enviadas uma unica vez ao acelerador
//  por 'acc declare copyin'. Dentro de uma 'acc routine seq' o gcc rejeita dados estaticos declarados
//  na propria funcao ('requires a declare directive for use in a routine function'); alem disso, a
//  versao local reinicializava 722 doubles (D3Q19) a cada chamada da colisao na CPU.
//
//  Ordem dos momentos (D3Q19): rho, e, eps, jx, qx, jy, qy, jz, qz, 3pxx, 3pi_xx, p_ww, pi_ww,
//                              p_xy, p_yz, p_xz, mx, my, mz
//  Fonte: d'Humieres-Lallemand-Luo (2002), reordenada para a numeracao de def_lattice_d3q19.
//
//====================================================================================================================//

static const double MRT_D3Q19_M[19][19] = { 
{ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. },
{ -30., -11., -11., -11., -11., -11., -11., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8. },
{ 12., -4., -4., -4., -4., -4., -4., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. },
{ 0., 1., -1., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0. },
{ 0., -4., 4., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., 1., -1., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., -1., 1., -1., 1. },
{ 0., 0., 0., -4., 4., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., -1., 1., -1., 1. },
{ 0., 0., 0., 0., 0., 1., -1., 0., 0., 0., 0., 1., -1., -1., 1., -1., 1., 1., -1. },
{ 0., 0., 0., 0., 0., -4., 4., 0., 0., 0., 0., 1., -1., -1., 1., -1., 1., 1., -1. },
{ 0., 2., 2., -1., -1., -1., -1., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2. },
{ 0., -4., -4., 2., 2., 2., 2., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2. },
{ 0., 0., 0., 1., 1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., -2., -2., 2., 2., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1. },
{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., 1., -1., 1., -1., -1., 1., -1., 1., 0., 0., 0., 0. },
{ 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., -1., 0., 0., 0., 0., -1., 1., -1., 1. },
{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 1., -1., -1., 1. }
};

static const double MRT_D3Q19_M_INV[19][19] = {
{ 1./19.,-5./399., 1./21., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 1./10.,-1./10., 0., 0., 0., 0., 1./18.,-1./18., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63.,-1./10., 1./10., 0., 0., 0., 0., 1./18.,-1./18., 0., 0., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0., 1./10.,-1./10., 0., 0.,-1./36., 1./36., 1./12.,-1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0.,-1./10., 1./10., 0., 0.,-1./36., 1./36., 1./12.,-1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0., 0., 0., 1./10.,-1./10.,-1./36., 1./36.,-1./12., 1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19.,-11./2394.,-1./63., 0., 0., 0., 0.,-1./10., 1./10.,-1./36., 1./36.,-1./12., 1./12., 0., 0., 0., 0., 0., 0. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40., 1./10., 1./40., 0., 0., 1./36., 1./72., 1./12., 1./24., 1./4., 0., 0., 1./8.,-1./8., 0. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40.,-1./10.,-1./40., 0., 0., 1./36., 1./72., 1./12., 1./24., 1./4., 0., 0.,-1./8., 1./8., 0. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40.,-1./10.,-1./40., 0., 0., 1./36., 1./72., 1./12., 1./24.,-1./4., 0., 0., 1./8., 1./8., 0. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40., 1./10., 1./40., 0., 0., 1./36., 1./72., 1./12., 1./24.,-1./4., 0., 0.,-1./8.,-1./8., 0. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40., 0., 0., 1./10., 1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0., 1./4.,-1./8., 0., 1./8. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40., 0., 0.,-1./10.,-1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0., 1./4., 1./8., 0.,-1./8. },
{ 1./19., 4./1197., 1./252., 1./10., 1./40., 0., 0.,-1./10.,-1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0.,-1./4.,-1./8., 0.,-1./8. },
{ 1./19., 4./1197., 1./252.,-1./10.,-1./40., 0., 0., 1./10., 1./40., 1./36., 1./72.,-1./12.,-1./24., 0., 0.,-1./4., 1./8., 0., 1./8. },
{ 1./19., 4./1197., 1./252., 0., 0.,-1./10.,-1./40.,-1./10.,-1./40.,-1./18.,-1./36., 0., 0., 0., 1./4., 0., 0.,-1./8., 1./8. },
{ 1./19., 4./1197., 1./252., 0., 0., 1./10., 1./40., 1./10., 1./40.,-1./18.,-1./36., 0., 0., 0., 1./4., 0., 0., 1./8.,-1./8. },
{ 1./19., 4./1197., 1./252., 0., 0.,-1./10.,-1./40., 1./10., 1./40.,-1./18.,-1./36., 0., 0., 0.,-1./4., 0., 0.,-1./8.,-1./8. },
{ 1./19., 4./1197., 1./252., 0., 0., 1./10., 1./40.,-1./10.,-1./40.,-1./18.,-1./36., 0., 0., 0.,-1./4., 0., 0., 1./8., 1./8. }
};

static const double MRT_D2Q9_M[9][9] = {
{ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 },	// m0 = rho
{ -4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0 },	// m1 = e
{  4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0 },	// m2 = epsilon
{  0.0,  1.0, -1.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0 },	// m3 = jx
{  0.0, -2.0,  2.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0 },	// m4 = qx
{  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0 },	// m5 = jy
{  0.0,  0.0,  0.0, -2.0,  2.0,  1.0, -1.0, -1.0,  1.0 },	// m6 = qy
{  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0 },	// m7 = pxx
{  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0 }	// m8 = pxy
};

static const double MRT_D2Q9_M_INV[9][9] = {

{  1./9.,   -1./9.,   1./9.,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 },				// f0
{  1./9.,  -1./36., -1.0/18.0,  1./6.,  -1./6.,   0.0,   0.0,   1./4.,   0.0 },			// f1  
{  1./9.,  -1./36., -1.0/18.0, -1./6.,   1./6.,   0.0,   0.0,   1./4.,   0.0 },			// f2  
{  1./9.,  -1./36., -1.0/18.0,  0.0,   0.0,   1./6.,  -1./6.,  -1./4.,   0.0 },			// f3  
{  1./9.,  -1./36., -1.0/18.0,  0.0,   0.0,  -1./6.,   1./6.,  -1./4.,   0.0 },			// f4  
{  1./9.,   1./18.,   1./36.,   1./6.,   1./12.,   1./6.,   1./12.,  0.0,   1./4. },	// f5  
{  1./9.,   1./18.,   1./36.,  -1./6.,  -1./12.,  -1./6.,  -1./12.,  0.0,   1./4. },	// f6  
{  1./9.,   1./18.,   1./36.,   1./6.,   1./12.,  -1./6.,  -1./12.,  0.0,  -1./4. },	// f7  
{  1./9.,   1./18.,   1./36.,  -1./6.,  -1./12.,   1./6.,   1./12.,  0.0,  -1./4. }		// f8  
};

//====================================================================================================================//


//=============================== Modelo color-gradient 3D melhorado - Wen, Li, Yu & Luo (2019) ======================//
//
//  Z. X. Wen, Q. Li, Y. Yu e Kai H. Luo, "Improved three-dimensional color-gradient lattice Boltzmann
//  model for immiscible two-phase flows", Phys. Rev. E 100, 023301 (2019).
//
//  Modelo de gradiente de cor com tres operadores, colisao MRT:
//
//      Omega_i = ( Omega_i )^(3) [ ( Omega_i )^(1) + ( Omega_i )^(2) ]
//
//      (1) colisao monofasica MRT com termo de correcao   -- eq. (44)
//      (2) operador de perturbacao ( tensao interfacial ) -- eqs. (14) e (60)
//      (3) recoloracao de Latva-Kokko                     -- eqs. (9) e (10)
//
//  O QUE E "MELHORADO":  nos modelos anteriores a pressao entra no equilibrio como p_k = rho_k cs_k^2
//  no SEGUNDO momento, mas o TERCEIRO momento continua carregando rho_k c^2/3.  Essa inconsistencia
//  gera termos de erro na equacao de momento -- eq. (17) do artigo -- que quebram a invariancia
//  galileana quando as densidades sao diferentes.  O modelo corrige isso de duas maneiras:
//
//    a) um termo de ordem alta no equilibrio, eq. (22), que poe p_k tambem nos elementos FORA da
//       diagonal do terceiro momento;
//
//    b) um termo de correcao G_i na colisao, eqs. (24) e (52), que trata os elementos DA diagonal,
//       que a simetria da D3Q19 nao permite corrigir pelo equilibrio.
//
//  Equilibrio ( eq. 22, com c = 1 ):
//
//      f^eq_i = rho_k [ phi^k_i + w_i ( 3 (e.u) + 9/2 (e.u)^2 - 3/2 |u|^2
//                                       + 3 (e.u) ( 3 cs_k^2 - 1 ) ( 3 |e_i|^2 - 5 ) ) ]
//
//      phi^k_0 = alpha_k ,  phi^k_{1..6} = (1-alpha_k)/12 ,  phi^k_{7..18} = (1-alpha_k)/24
//      cs_k^2  = ( 1 - alpha_k ) / 2       p_k = rho_k cs_k^2
//
//  NOTA sobre o coeficiente do termo de ordem alta:  o artigo o imprime como 3(e.u)/(2c^2)(...)(...),
//  mas esse valor nao satisfaz a eq. (23), que e o proprio objetivo do termo.  Impondo a eq. (23) --
//  terceiro momento fora da diagonal igual a p_k, e nao a rho_k c^2/3 -- o coeficiente sai unico e
//  vale 3, nao 3/2.  Com 3 as quinze condicoes das eqs. (20) e (23) sao satisfeitas exatamente
//  ( verificado em algebra simbolica ); com 3/2, seis delas falham.  Ver LI_YU_LUO.md.
//
//  Correcao G_i, em espaco de momentos ( eq. 52 ):  so os momentos 4, 5 e 6 sao nao nulos,
//
//      C_4 = Qx + Qy + Qz ,   C_5 = 2 Qx - Qy - Qz ,   C_6 = Qy - Qz
//
//      Q_alpha = d_alpha [ rho_k u_alpha ( c^2 - 3 cs_k^2 ) ]        eqs. (39)-(41)
//
//  Q_alpha e calculado fora desta funcao ( precisa dos vizinhos ) e chega pelos vetores Q_R e Q_B.
//  Repare que  c^2 - 3 cs_k^2 = 0  quando alpha_k = 1/3: com densidades iguais a correcao some, que
//  e exatamente o que o artigo diz na secao II C.
//
//  Perturbacao ( eq. 14 ), com n = grad(rho^N) / |grad(rho^N)| :
//
//      ( Omega^k_i )^(2) = ( A_k / 2 ) |grad rho^N| [ w_i (e_i.n)^2 - B_i ]
//      B_0 = -1/3 ,  B_{1..6} = 1/18 ,  B_{7..18} = 1/36
//
//  Em espaco de momentos esse vetor tem forma fechada -- nao e preciso multiplicar por M:
//
//      mP_4 = -4/9                mP_7 = 2 nx ny / 9        mP_16 = - nz^2 / 9
//      mP_5 = 2 ( 3 nx^2 - 1 )/9  mP_8 = 2 nx nz / 9        mP_17 = - ny^2 / 9
//      mP_6 = 2 ( ny^2 - nz^2 )/9 mP_9 = 2 ny nz / 9        mP_18 = - nx^2 / 9
//
//  e todos os outros sao zero -- em particular mP_0 = mP_1 = mP_2 = mP_3 = 0, ou seja a perturbacao
//  conserva massa e quantidade de movimento exatamente.  Seguindo a eq. (60), ela e multiplicada
//  pela matriz de relaxacao, o que torna a tensao interfacial independente de tau:
//
//      sigma = 2 ( A_R + A_B ) c^4 dt / 9
//
//  Recoloracao de Latva-Kokko ( eqs. 9 e 10 ):
//
//      f^{R,+}_i = (rho_R/rho) f*_i + beta (rho_R rho_B / rho^2) cos(phi_i) soma_k rho_k phi^k_i
//      f^{B,+}_i = (rho_B/rho) f*_i - beta (rho_R rho_B / rho^2) cos(phi_i) soma_k rho_k phi^k_i
//
//      cos(phi_i) = ( e_i . grad rho^N ) / ( |e_i| |grad rho^N| )
//
//  Parametros, lidos de PARAMETERS:
//
//      tau_R , tau_B    tau_v de cada fluido puro     ( nu_k = cs_k^2 ( tau_v^k - 1/2 ) )
//      tau_m            tau_e = tau_q = tau_pi        ( o artigo usa 1.0 )
//      A_fact_R/_B      A_R e A_B da tensao interfacial
//      recoll           beta da recoloracao
//      rho_ini_R/_B     rho_R^in e rho_B^in           ( definicao de rho^N, eq. 12 )
//      alfa_R           alpha_R;  alpha_B sai da condicao de equilibrio de pressao
//
//====================================================================================================================//

//  As matrizes MRT e a colisao de Li-Yu-Luo sao especificas do D3Q19.  Num programa D3Q27
//  ( Saito ) elas nao existem -- e indexar 19x19 com nvel = 27 seria erro de memoria.

#if nvel == 19

static const double MRT_LYL_M[19][19] =
{
	{        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1. },
	{        0.,        1.,       -1.,        0.,        0.,        0.,        0.,        1.,       -1.,        1.,       -1.,        1.,       -1.,        1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        1.,       -1.,        0.,        0.,        1.,       -1.,       -1.,        1.,        0.,        0.,        0.,        0.,       -1.,        1.,       -1.,        1. },
	{        0.,        0.,        0.,        0.,        0.,        1.,       -1.,        0.,        0.,        0.,        0.,        1.,       -1.,       -1.,        1.,       -1.,        1.,        1.,       -1. },
	{        0.,        1.,        1.,        1.,        1.,        1.,        1.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2.,        2. },
	{        0.,        2.,        2.,       -1.,       -1.,       -1.,       -1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,        1.,       -2.,       -2.,       -2.,       -2. },
	{        0.,        0.,        0.,        1.,        1.,       -1.,       -1.,        1.,        1.,        1.,        1.,       -1.,       -1.,       -1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,       -1.,       -1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,       -1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,       -1.,       -1. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,       -1.,        1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,        1.,       -1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,       -1.,        1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,       -1.,        1.,       -1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,       -1.,        1.,        1.,       -1. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,       -1.,        1.,       -1.,        1. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1.,        1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1.,        1.,        0.,        0.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1.,        1. },
};

static const double MRT_LYL_M_INV[19][19] =
{
	{        1.,        0.,        0.,        0.,       -1.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        1.,        1.,        1. },
	{        0.,     1./2.,        0.,        0.,     1./6.,     1./6.,        0.,        0.,        0.,        0.,        0.,    -1./2.,        0.,    -1./2.,        0.,        0.,    -1./2.,    -1./2.,        0. },
	{        0.,    -1./2.,        0.,        0.,     1./6.,     1./6.,        0.,        0.,        0.,        0.,        0.,     1./2.,        0.,     1./2.,        0.,        0.,    -1./2.,    -1./2.,        0. },
	{        0.,        0.,     1./2.,        0.,     1./6.,   -1./12.,     1./4.,        0.,        0.,        0.,    -1./2.,        0.,        0.,        0.,        0.,    -1./2.,    -1./2.,        0.,    -1./2. },
	{        0.,        0.,    -1./2.,        0.,     1./6.,   -1./12.,     1./4.,        0.,        0.,        0.,     1./2.,        0.,        0.,        0.,        0.,     1./2.,    -1./2.,        0.,    -1./2. },
	{        0.,        0.,        0.,     1./2.,     1./6.,   -1./12.,    -1./4.,        0.,        0.,        0.,        0.,        0.,    -1./2.,        0.,    -1./2.,        0.,        0.,    -1./2.,    -1./2. },
	{        0.,        0.,        0.,    -1./2.,     1./6.,   -1./12.,    -1./4.,        0.,        0.,        0.,        0.,        0.,     1./2.,        0.,     1./2.,        0.,        0.,    -1./2.,    -1./2. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,     1./4.,     1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,    -1./4.,    -1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,    -1./4.,     1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,     1./4.,    -1./4.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,     1./4.,     1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,    -1./4.,    -1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,    -1./4.,     1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,     1./4.,    -1./4.,        0.,        0.,        0.,     1./4.,        0. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,        0.,    -1./4.,    -1./4.,        0.,        0.,     1./4. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,     1./4.,        0.,        0.,        0.,        0.,     1./4.,     1./4.,        0.,        0.,     1./4. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,        0.,     1./4.,    -1./4.,        0.,        0.,     1./4. },
	{        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,        0.,    -1./4.,        0.,        0.,        0.,        0.,    -1./4.,     1./4.,        0.,        0.,     1./4. },
};

#pragma acc declare copyin( MRT_LYL_M, MRT_LYL_M_INV )

#endif   // nvel == 19

//------------------ Equilibrio do modelo, eq. (22) -----------------------------------------------------------------//

#pragma acc routine seq
void dist_eq_LYL ( double *f_eq, double vx, double vy, double vz, double rho_k, double alfa_k,
                   double hi_order, LATTICE lattice )
{
	const double cs2_k = 0.5 * ( 1.0 - alfa_k );

	//  hi_order = 1 liga o termo de ordem alta da eq. (22);  0 recai no equilibrio da eq. (6),
	//  que e o "modelo original" com que o artigo compara.

	const double fator = hi_order * 3.0 * ( 3.0 * cs2_k - 1.0 );

	const double u2 = vx * vx + vy * vy + vz * vz;

	for ( int i = 0; i < nvel; i++ )
	{
		const double cx = lattice.c_i[ i * dim + 0 ];
		const double cy = lattice.c_i[ i * dim + 1 ];
		const double cz = ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0;

		const double eu = cx * vx + cy * vy + cz * vz;
		const double e2 = cx * cx + cy * cy + cz * cz;

		//  phi^k_i : alpha_k no repouso, (1-alpha_k)/12 nos eixos, (1-alpha_k)/24 nas diagonais

		double phi_i;

		if      ( i == 0 ) phi_i = alfa_k;
		else if ( i <  7 ) phi_i = ( 1.0 - alfa_k ) / 12.0;
		else               phi_i = ( 1.0 - alfa_k ) / 24.0;

		f_eq[i] = rho_k * ( phi_i + lattice.w[i] * ( 3.0 * eu + 4.5 * eu * eu - 1.5 * u2
		                                             + fator * eu * ( 3.0 * e2 - 5.0 ) ) );
	}
}

//------------------ Colisao completa de um sitio -------------------------------------------------------------------//
//
//      grad_N   grad( rho^N ), calculado fora ( precisa dos vizinhos )
//      Q_R,Q_B  Q_alpha de cada fluido, eqs. (39)-(41), tambem calculados fora
//      gx,gy,gz aceleracao de corpo ( entra pelo esquema de Guo em espaco de momentos )

#if nvel == 19

#endif   // nvel == 19

//====================================================================================================================//


//=============================== Modelo de Saito et al. ( 2023 ), D3Q27 =============================================//
//
//      Saito, Takada, Baba, Someya & Ito, Phys. Rev. E 108, 065305 (2023),
//      "Generalized equilibria for color-gradient lattice Boltzmann model based on higher-order
//       Hermite polynomials: A simplified implementation with central moments".
//
//      Este bloco implementa o modelo inteiro.  Ele difere dos outros color-gradient desta arvore
//      em tres pontos, e vale separar o que e novidade do artigo do que e so heranca:
//
//      1) A TENSAO INTERFACIAL NAO VEM DE OPERADOR DE PERTURBACAO.  Vem de uma forca de corpo do
//         tipo continuum surface force, eq. (15):
//
//             F_s = ( 1/2 ) sigma kappa grad( phi ) ,      kappa = - div( n ) ,   n = grad(phi)/|grad(phi)|
//
//         Isto e: sigma e ENTRADA, nao resultado de calibracao.  Nao ha "sigma = C A" para medir --
//         o teste da gota estatica serve para CONFERIR o modelo, nao para calibra-lo.  A curvatura
//         nao e local ( precisa dos vizinhos de n ), entao quem a calcula e o programa.
//
//      2) O EQUILIBRIO GENERALIZADO, eq. (56), que e a contribuicao do artigo.  Em espaco de
//         populacoes ele e um polinomio de Hermite ate SEXTA ordem em u, com dezenas de termos.  Em
//         espaco de MOMENTOS CENTRAIS ele e isto, e nada mais ( eqs. 57-63 ):
//
//             k000 = rho ;   k200 = k020 = k002 = p ;
//             k220 = k202 = k022 = p cs^2 ;   k222 = p cs^4 ;   todos os outros = 0
//
//         INDEPENDENTE DA VELOCIDADE.  E dai que sai a invariancia galileana: nos equilibrios
//         anteriores ( modelos A, B e C do artigo ) sobram termos proporcionais a ( p - rho cs^2 )
//         vezes potencias de u nos momentos de terceira ordem e acima, e sao eles que deformam uma
//         gota que se move.
//
//         Por isso esta implementacao NAO escreve a eq. (56): escreve os sete numeros acima e
//         volta ao espaco de populacoes pela transformada inversa.  E o mesmo equilibrio, exato, e
//         cabe em cinco linhas em vez de trinta termos.
//
//      3) A RECOLORACAO e a de Halliday, eqs. (18)-(19):
//
//             f_i^{r,*} = ( rho_r / rho ) f_i^*  +  ( w_i / cs^2 ) c_i . R ,    R = beta ( rho_r rho_b / rho^2 ) p n
//
//         Repare que e  c_i . R , e nao  |c_i| cos(phi_i) : em relacao a recoloracao de
//         Latva-Kokko usada nos outros modelos desta arvore, sobra um fator | c_i | nas direcoes
//         diagonais.  E a amplitude e p, a pressao, e nao rho w_i.
//
//      A pressao e  p = soma_k rho_k ( cs^k )^2 , eq. (12), com uma velocidade do som por fluido --
//      e dela que vem a razao de densidades, eq. (13):  gamma = ( cs^b )^2 / ( cs^r )^2 .
//
//      A rede e D3Q27 ( eq. 1 ).  Numa caixa com nz = 1 ela degenera EXATAMENTE no D2Q9:
//      axial 2/27 + 2 x 1/54 = 1/9 , diagonal 1/54 + 2 x 1/216 = 1/36 , repouso 8/27 + 2 x 2/27 = 4/9.
//      O teste da gota estatica do artigo ( secao IV A, Tabela II ) e justamente bidimensional.
//
//====================================================================================================================//

//------------------ Indice da rede a partir do trio ( a, b, c ) ----------------------------------------------------//
//
//  a, b, c em { 0, 1, 2 }  correspondem a  c_alpha em { 0, +1, -1 }.  As 27 velocidades do D3Q27 sao
//  o produto tensorial dos tres eixos, e e isso que faz a transformada de momentos centrais ser
//  SEPARAVEL: em vez de uma matriz 27 x 27, tres passagens de tamanho 3, uma por eixo.
//
//  SAITO_MAP[ a * 9 + b * 3 + c ] = indice i da velocidade ( c_x, c_y, c_z ) nesta biblioteca.

#define SAITO_NM 27

static const int SAITO_MAP[27] =
{
	 0,  5,  6,
	 3, 16, 18,
	 4, 17, 15,
	 1, 11, 13,
	 7, 19, 24,
	 9, 25, 22,
	 2, 14, 12,
	10, 21, 26,
	 8, 23, 20
};

//------------------ Transformada inversa: momentos centrais -> populacoes ------------------------------------------//
//
//  Num eixo:   m0 = k0 ,  m1 = k1 + u k0 ,  m2 = k2 + 2 u k1 + u^2 k0
//              f0 = m0 - m2 ,  f+ = ( m2 + m1 ) / 2 ,  f- = ( m2 - m1 ) / 2

#pragma acc routine seq
void cm_saito_inv ( const double *k, double ux, double uy, double uz, double *f )
{
	double t[27], s[27];

	//------ eixo z ---------------------------------------------------------------------------------//

	for ( int a = 0; a < 3; a++ )
	{
		for ( int b = 0; b < 3; b++ )
		{
			const double k0 = k[ a * 9 + b * 3 + 0 ];
			const double k1 = k[ a * 9 + b * 3 + 1 ];
			const double k2 = k[ a * 9 + b * 3 + 2 ];

			const double m0 = k0;
			const double m1 = k1 + uz * k0;
			const double m2 = k2 + 2.0 * uz * k1 + uz * uz * k0;

			t[ a * 9 + b * 3 + 0 ] = m0 - m2;
			t[ a * 9 + b * 3 + 1 ] = 0.5 * ( m2 + m1 );
			t[ a * 9 + b * 3 + 2 ] = 0.5 * ( m2 - m1 );
		}
	}

	//------ eixo y ---------------------------------------------------------------------------------//

	for ( int a = 0; a < 3; a++ )
	{
		for ( int c = 0; c < 3; c++ )
		{
			const double k0 = t[ a * 9 + 0 * 3 + c ];
			const double k1 = t[ a * 9 + 1 * 3 + c ];
			const double k2 = t[ a * 9 + 2 * 3 + c ];

			const double m0 = k0;
			const double m1 = k1 + uy * k0;
			const double m2 = k2 + 2.0 * uy * k1 + uy * uy * k0;

			s[ a * 9 + 0 * 3 + c ] = m0 - m2;
			s[ a * 9 + 1 * 3 + c ] = 0.5 * ( m2 + m1 );
			s[ a * 9 + 2 * 3 + c ] = 0.5 * ( m2 - m1 );
		}
	}

	//------ eixo x ---------------------------------------------------------------------------------//

	for ( int b = 0; b < 3; b++ )
	{
		for ( int c = 0; c < 3; c++ )
		{
			const double k0 = s[ 0 * 9 + b * 3 + c ];
			const double k1 = s[ 1 * 9 + b * 3 + c ];
			const double k2 = s[ 2 * 9 + b * 3 + c ];

			const double m0 = k0;
			const double m1 = k1 + ux * k0;
			const double m2 = k2 + 2.0 * ux * k1 + ux * ux * k0;

			f[ SAITO_MAP[ 0 * 9 + b * 3 + c ] ] = m0 - m2;
			f[ SAITO_MAP[ 1 * 9 + b * 3 + c ] ] = 0.5 * ( m2 + m1 );
			f[ SAITO_MAP[ 2 * 9 + b * 3 + c ] ] = 0.5 * ( m2 - m1 );
		}
	}
}

//------------------ Equilibrio generalizado, eq. (56) --------------------------------------------------------------//
//
//  Escrito pelos momentos centrais das eqs. (57)-(63), que sao independentes da velocidade, e
//  trazido de volta ao espaco de populacoes.  E EXATAMENTE a eq. (56) -- a transformada e bijetiva.
//
//      k000 = rho ;  k200 = k020 = k002 = p ;  k220 = k202 = k022 = p cs^2 ;  k222 = p cs^4
//
//  Conferido em runtime por check_saito(): soma f^eq = rho, soma f^eq c = rho u, e os momentos de
//  terceira ordem da eq. (64).
//
//      Input: velocidade, densidade total, pressao, lattice
//      Output: f_eq[27]

#pragma acc routine seq
void dist_eq_saito ( double *f_eq, double ux, double uy, double uz, double rho, double p, LATTICE lattice )
{
	const double cs2 = lattice.c_s2;

	double k[27];

	for ( int n = 0; n < 27; n++ ) k[n] = 0.0;

	k[ 0 * 9 + 0 * 3 + 0 ] = rho;						// k000

	k[ 2 * 9 + 0 * 3 + 0 ] = p;							// k200
	k[ 0 * 9 + 2 * 3 + 0 ] = p;							// k020
	k[ 0 * 9 + 0 * 3 + 2 ] = p;							// k002

	k[ 2 * 9 + 2 * 3 + 0 ] = p * cs2;					// k220
	k[ 2 * 9 + 0 * 3 + 2 ] = p * cs2;					// k202
	k[ 0 * 9 + 2 * 3 + 2 ] = p * cs2;					// k022

	k[ 2 * 9 + 2 * 3 + 2 ] = p * cs2 * cs2;				// k222

	cm_saito_inv ( k, ux, uy, uz, f_eq );
}

//====================================================================================================================//
