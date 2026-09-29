
//=============================== round a number =====================================================================//


int round_number ( double  num )
{
    int num_int;

    if ( num < 0 ) num_int = ( int ) ( num - 0.5 );
    else num_int = ( int ) ( num + 0.5 );

    return num_int;
}

//=============================== Medidas de uma gota ================================================================//
//
//      calc_drop_radius     : centro, raio e correntes espurias de uma gota isolada
//      calc_laplace_tension : salto de pressao e tensao interfacial pela lei de Laplace
//
//      As duas trabalham sobre a struct DROP ( Definitions.cpp ), que leva as entradas
//      ( cilindro, frac_in, frac_out, espessura ) e recebe todas as saidas.  Rodam no HOSPEDEIRO:
//      num programa OpenACC, faca 'acc update self' de inif_R e inif_B antes de chamar.  O custo e
//      desprezivel perto dos milhares de passos entre dois diagnosticos.
//
//====================================================================================================================//

//=============================== Centro, raio e correntes espurias de uma gota ======================================//
//
//      Preenche, em DROP:
//
//          x0, y0, z0    centro da gota, pelo CENTROIDE CIRCULAR de phi = rho_R / ( rho_R + rho_B ).
//                        Circular porque a caixa e periodica: a media aritmetica de x erraria assim
//                        que a gota encostasse na borda.
//
//          area          soma_sitios phi.  Como soma_sitios rho_R e conservada exatamente, esta e a
//                        medida estavel do tamanho da gota.
//
//          raio          raio equimolar:  raiz( area / ( pi nz ) )  ( cilindro )
//                                         ( 3 area / 4 pi )^(1/3)   ( esfera )
//
//          raio_grad     primeiro momento radial de rho_R rho_B sobre o dominio inteiro: o raio
//                        medio da superficie de tensao.  Estimador independente do de cima; os dois
//                        concordam dentro de ( espessura / raio )^2.  Se divergirem, a gota nao esta
//                        redonda ou nao convergiu.
//
//          raio_linha    o estimador de uma linha so ( primeiro momento de rho_R rho_B ao longo da
//                        linha que passa pelo centro ), guardado por compatibilidade com programas
//                        antigos.  E FRAGIL: usa poucos sitios e soma a linha inteira, inclusive o
//                        campo distante, onde rho_R rho_B e minusculo mas o braco de alavanca e
//                        grande.  Nao use para calibrar nada.
//
//          massa_R, massa_B, qx, qy, qz, u_max, u_rms, n_ruim
//
//      Entradas usadas: drop.cilindro.
//
//      Retorna drop.raio.
//
//====================================================================================================================//

double calc_drop_radius ( GEOMETRY geometry, LATTICE lattice, DROP& drop )
{
	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	int n_voxels = nx * ny * nz;

	const double pi = 3.14159265358979323846;

	const double tx = 2.0 * pi / ( double ) nx;
	const double ty = 2.0 * pi / ( double ) ny;
	const double tz = 2.0 * pi / ( double ) nz;

	//--------------- Passo 1: massas, area, centroide e correntes espurias ---------------------------//

	double massa_R = 0.0, massa_B = 0.0;
	double qx = 0.0, qy = 0.0, qz = 0.0;
	double area = 0.0, u2_sum = 0.0, u_max = 0.0;
	double cx_c = 0.0, cx_s = 0.0, cy_c = 0.0, cy_s = 0.0, cz_c = 0.0, cz_s = 0.0;

	int n_ruim = 0, n_fluid = 0;

	#pragma omp parallel for \
		reduction(+:massa_R,massa_B,qx,qy,qz,area,u2_sum,cx_c,cx_s,cy_c,cy_s,cz_c,cz_s,n_ruim,n_fluid) \
		reduction(max:u_max)

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		int pto = geo[pos] - 1;

		double vxR, vyR, vzR, rho_R;
		double vxB, vyB, vzB, rho_B;

		calcula ( lattice.inif_R + pto * nvel, vxR, vyR, vzR, rho_R, lattice );
		calcula ( lattice.inif_B + pto * nvel, vxB, vyB, vzB, rho_B, lattice );

		n_fluid = n_fluid + 1;

		//  Teste de finitude escrito sem isnan() para poder ser reaproveitado em codigo de
		//  dispositivo compilado com matematica rapida.

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) || ! ( vxR == vxR ) || ! ( vxB == vxB ) )
		{
			n_ruim = n_ruim + 1;

			continue;
		}

		massa_R = massa_R + rho_R;
		massa_B = massa_B + rho_B;

		qx = qx + rho_R * vxR + rho_B * vxB;
		qy = qy + rho_R * vyR + rho_B * vyB;
		qz = qz + rho_R * vzR + rho_B * vzB;

		double rho = rho_R + rho_B;

		if ( rho <= 0.0 ) continue;

		double ux = ( rho_R * vxR + rho_B * vxB ) / rho;
		double uy = ( rho_R * vyR + rho_B * vyB ) / rho;
		double uz = ( rho_R * vzR + rho_B * vzB ) / rho;

		//  CORRECAO DA MEIA-FORCA.  Num esquema forcado a Guo a velocidade fisica NAO e o momento
		//  cru das populacoes:
		//
		//      antes da colisao    rho u = soma_i f_i c_i  +  F / 2
		//      depois da colisao   rho u = soma_i f_i c_i  -  F / 2
		//
		//  Esta medida roda DEPOIS da colisao, entao o sinal e o de baixo.  Sem isto, uma gota em
		//  equilibrio perfeito ( u = 0 em toda parte ) aparece com uma "corrente espuria" igual a
		//  | F | / 2 rho -- radial, axissimetrica, concentrada na interface e proporcional a
		//  tensao interfacial.  Nao e corrente nenhuma: e a metade da forca.
		//
		//  Medido na lei de Laplace, SHC com forcamento de Guo, R = 40, beta = 0.8, sigma = 5.36e-03:
		//  |u|max caiu de 3.54e-04 ( = |F|max / 2 rho a 0.7% ) para 1.97e-05 quando esta linha
		//  passou a existir.  Um fator 18.
		//
		//  lattice.ini_force = nullptr ( o padrao ) desliga a correcao: os modelos que nao usam
		//  forca de corpo -- Gunstensen, Santos, a versao H -- nao precisam dela.

		if ( lattice.ini_force != nullptr )
		{
			const double *F = lattice.ini_force + ( long long ) pto * 3;

			ux = ux - 0.5 * F[0] / rho;
			uy = uy - 0.5 * F[1] / rho;
			uz = uz - 0.5 * F[2] / rho;

			qx = qx - 0.5 * F[0];
			qy = qy - 0.5 * F[1];
			qz = qz - 0.5 * F[2];
		}

		double u2 = ux * ux + uy * uy + uz * uz;

		u2_sum = u2_sum + u2;

		double u = sqrt ( u2 );

		if ( u > u_max ) u_max = u;

		double phi = rho_R / rho;

		area = area + phi;

		int z = pos / ( nx * ny );
		int y = ( pos - z * nx * ny ) / nx;
		int x = pos - z * nx * ny - y * nx;

		cx_c = cx_c + phi * cos ( tx * x );   cx_s = cx_s + phi * sin ( tx * x );
		cy_c = cy_c + phi * cos ( ty * y );   cy_s = cy_s + phi * sin ( ty * y );
		cz_c = cz_c + phi * cos ( tz * z );   cz_s = cz_s + phi * sin ( tz * z );
	}

	double x0 = atan2 ( cx_s, cx_c ) / tx;   if ( x0 < 0.0 ) x0 += nx;
	double y0 = atan2 ( cy_s, cy_c ) / ty;   if ( y0 < 0.0 ) y0 += ny;
	double z0 = atan2 ( cz_s, cz_c ) / tz;   if ( z0 < 0.0 ) z0 += nz;

	if ( nz == 1 ) z0 = 0.0;

	double raio;

	if ( drop.cilindro ) raio = ( area > 0.0 ) ? sqrt ( area / ( pi * ( double ) nz ) ) : 0.0;

	else                 raio = ( area > 0.0 ) ? cbrt ( 3.0 * area / ( 4.0 * pi ) )     : 0.0;

	//--------------- Passo 2: raio da interface, pelo momento de rho_R rho_B -------------------------//

	double som_w = 0.0, som_wr = 0.0;
	double sx_e = 0.0, sw_e = 0.0, sx_d = 0.0, sw_d = 0.0;

	//  A LINHA DO ESTIMADOR raio_linha.
	//
	//  Num CILINDRO o campo e uniforme em z e o centroide circular em z e indefinido: cz_s e cz_c
	//  sao exatamente zero no papel e ruido de arredondamento ( ~1e-13 ) no computador, e o atan2
	//  de dois ruidos devolve um angulo qualquer em [ 0, nz ).  Com z0 > nz - 0.5 o arredondamento
	//  dava  z_lin = nz , que nao e plano nenhum: nenhum sitio casava, sw_e = sw_d = 0 e tanto
	//  raio_linha quanto sigma_2pontos saiam ZERO -- de forma intermitente, so em alguns passos.
	//  Num cilindro qualquer plano serve, entao fixamos z = 0.  Na esfera, o clamp evita o mesmo
	//  arredondamento na borda.

	int y_lin = ( int ) ( y0 + 0.5 );   if ( y_lin >= ny ) y_lin -= ny;

	int z_lin = ( nz == 1 || drop.cilindro ) ? 0 : ( int ) ( z0 + 0.5 );   if ( z_lin >= nz ) z_lin -= nz;

	#pragma omp parallel for reduction(+:som_w,som_wr,sx_e,sw_e,sx_d,sw_d)

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		int pto = geo[pos] - 1;

		double rho_R = density ( lattice.inif_R + pto * nvel );
		double rho_B = density ( lattice.inif_B + pto * nvel );

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) ) continue;

		int z = pos / ( nx * ny );
		int y = ( pos - z * nx * ny ) / nx;
		int x = pos - z * nx * ny - y * nx;

		//  Imagem minima: a caixa e periodica.

		double dx = ( double ) x - x0;   if ( dx >  0.5 * nx ) dx -= nx;   if ( dx < -0.5 * nx ) dx += nx;
		double dy = ( double ) y - y0;   if ( dy >  0.5 * ny ) dy -= ny;   if ( dy < -0.5 * ny ) dy += ny;
		double dz = ( double ) z - z0;   if ( dz >  0.5 * nz ) dz -= nz;   if ( dz < -0.5 * nz ) dz += nz;

		if ( drop.cilindro ) dz = 0.0;

		double r = sqrt ( dx * dx + dy * dy + dz * dz );

		double w = rho_R * rho_B;

		som_w  = som_w  + w;
		som_wr = som_wr + w * r;

		if ( y == y_lin && z == z_lin )
		{
			if ( dx < 0.0 ) { sx_e = sx_e + w * ( double ) x;   sw_e = sw_e + w; }

			else            { sx_d = sx_d + w * ( double ) x;   sw_d = sw_d + w; }
		}
	}

	//--------------- Devolve --------------------------------------------------------------------------//

	drop.x0 = x0;   drop.y0 = y0;   drop.z0 = z0;

	drop.area = area;

	drop.raio = raio;

	drop.raio_grad = ( som_w > 0.0 ) ? som_wr / som_w : 0.0;

	drop.raio_linha = ( sw_e > 0.0 && sw_d > 0.0 ) ? 0.5 * ( sx_d / sw_d - sx_e / sw_e ) : 0.0;

	drop.massa_R = massa_R;   drop.massa_B = massa_B;

	drop.qx = qx;   drop.qy = qy;   drop.qz = qz;

	drop.u_max = u_max;

	drop.u_rms = ( n_fluid > 0 ) ? sqrt ( u2_sum / ( double ) n_fluid ) : 0.0;

	drop.n_fluid = n_fluid;

	drop.n_ruim = n_ruim;

	return raio;
}

//====================================================================================================================//


//=============================== Curva interfacial de Latva-Kokko & Rothman ( 2005 ) ================================//
//
//      Latva-Kokko & Rothman, Phys. Rev. E 71, 056702 (2005), eqs. (15)-(17).  A recoloracao
//      daqueles autores nao so separa as cores: ela leva o perfil a uma forma FECHADA.  O balanco
//      entre a difusao de cor e o termo de segregacao da
//
//          d phi / ds = K phi ( 1 - phi )                                        eq. (15)
//
//      cuja solucao e uma logistica
//
//          phi( s ) = 1 / ( 1 + exp[ K ( s - s_0 ) ] )                           eq. (17)
//
//      com  K  uma constante de comprimento INVERSO, proporcional a beta.  O artigo mostra isso no
//      modelo D1Q2 ( onde K = 4 beta / N ) e afirma que a FORMA vale para qualquer rede que use a
//      recoloracao das eqs. (9)-(10) -- so o prefator K( beta ) muda.  E o que esta funcao mede.
//
//      Aqui  s  e a distancia radial ao centro da gota, e o perfil sai da media azimutal de
//      phi = rho_R / rho  sobre cascas de meio sitio de largura.  O ajuste e linear:
//
//          ln[ phi / ( 1 - phi ) ] = - K ( r - r_0 )
//
//      entao a inclinacao da K e o zero da o raio da interface.  So entram as camadas com
//      0.02 < phi < 0.98 : fora dessa faixa phi e dominado pela cauda exponencial e pelo ruido de
//      arredondamento, e o logito estoura.
//
//      O que sai em DROP:
//
//          K_perfil    K medido
//          r_perfil    raio onde phi = 1/2 ( terceiro estimador independente de raio )
//          res_perfil  maior residuo de ln[phi/(1-phi)] -- e ele que diz se a logistica cabe
//          esp_perfil  largura 10%-90% implicada, 2 ln(9) / K
//          n_perfil    quantas camadas entraram no ajuste
//
//      Devolve K_perfil.  Roda no hospedeiro, com OpenMP, como as outras medidas.
//
//====================================================================================================================//

double calc_perfil_logistico ( GEOMETRY geometry, LATTICE lattice, DROP& drop )
{
	constexpr double eps = 1.0e-30;

	const int nx = geometry.nx;
	const int ny = geometry.ny;
	const int nz = geometry.nz;

	const int n_voxels = nx * ny * nz;

	const int *geo = geometry.ini;

	//  Cascas de meio sitio.  O alcance vai ate a diagonal da caixa, o que basta e sobra.

	const double d_casca = 0.5;

	const int n_casca = ( int ) ( sqrt ( ( double ) ( nx * nx + ny * ny + nz * nz ) ) / d_casca ) + 2;

	double *soma_phi = new double[ n_casca ]();
	double *soma_n   = new double[ n_casca ]();

	const double x0 = drop.x0, y0 = drop.y0, z0 = drop.z0;

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		const int pto = geo[pos] - 1;

		const double rho_R = density ( lattice.inif_R + pto * nvel );
		const double rho_B = density ( lattice.inif_B + pto * nvel );

		const double rho = rho_R + rho_B;

		if ( rho <= eps ) continue;

		const int z = pos / ( nx * ny );
		const int y = ( pos - z * nx * ny ) / nx;
		const int x = pos - z * nx * ny - y * nx;

		//  Mesma convencao periodica das outras medidas.

		double dx = ( double ) x - x0;   if ( dx >  0.5 * nx ) dx -= nx;   if ( dx < -0.5 * nx ) dx += nx;
		double dy = ( double ) y - y0;   if ( dy >  0.5 * ny ) dy -= ny;   if ( dy < -0.5 * ny ) dy += ny;
		double dz = ( double ) z - z0;   if ( dz >  0.5 * nz ) dz -= nz;   if ( dz < -0.5 * nz ) dz += nz;

		if ( drop.cilindro ) dz = 0.0;

		//  drop.plana = true : a interface e um plano normal a x, e a coordenada do perfil e a
		//  distancia |x - x0| ao centro da laje -- o analogo unidimensional do raio.

		const double r = drop.plana ? fabs ( dx ) : sqrt ( dx * dx + dy * dy + dz * dz );

		const int k = ( int ) ( r / d_casca );

		if ( k < 0 || k >= n_casca ) continue;

		soma_phi[k] += rho_R / rho;
		soma_n  [k] += 1.0;
	}

	//  Ajuste de minimos quadrados de  ln[ phi / ( 1 - phi ) ]  contra r.

	double S = 0.0, Sr = 0.0, Sl = 0.0, Srr = 0.0, Srl = 0.0;

	int n_usadas = 0;

	for ( int k = 0; k < n_casca; k++ )
	{
		if ( soma_n[k] < 1.0 ) continue;

		const double phi = soma_phi[k] / soma_n[k];

		if ( phi <= 0.02 || phi >= 0.98 ) continue;

		const double r = ( ( double ) k + 0.5 ) * d_casca;

		const double l = log ( phi / ( 1.0 - phi ) );

		S += 1.0;  Sr += r;  Sl += l;  Srr += r * r;  Srl += r * l;

		n_usadas++;
	}

	double K = 0.0, r_meio = 0.0, res_max = 0.0;

	if ( n_usadas >= 3 )
	{
		const double den = S * Srr - Sr * Sr;

		if ( fabs ( den ) > eps )
		{
			const double a = ( S * Srl - Sr * Sl ) / den;		// inclinacao
			const double b = ( Sl - a * Sr ) / S;				// intercepto

			K = -a;

			if ( fabs ( a ) > eps ) r_meio = -b / a;			// onde ln[...] = 0, isto e phi = 1/2

			//  Residuo: quanto a logistica deixa de explicar.

			for ( int k = 0; k < n_casca; k++ )
			{
				if ( soma_n[k] < 1.0 ) continue;

				const double phi = soma_phi[k] / soma_n[k];

				if ( phi <= 0.02 || phi >= 0.98 ) continue;

				const double r = ( ( double ) k + 0.5 ) * d_casca;

				const double d = log ( phi / ( 1.0 - phi ) ) - ( a * r + b );

				if ( fabs ( d ) > res_max ) res_max = fabs ( d );
			}
		}
	}

	delete[] soma_phi;
	delete[] soma_n;

	drop.K_perfil   = K;
	drop.r_perfil   = r_meio;

	if ( drop.plana ) drop.x_interface = drop.x0 + r_meio;
	drop.res_perfil = res_max;
	drop.esp_perfil = ( K > eps ) ? 2.0 * log ( 9.0 ) / K : 0.0;
	drop.n_perfil   = n_usadas;

	return K;
}

//====================================================================================================================//


double calc_laplace_tension ( GEOMETRY geometry, LATTICE lattice, DROP& drop )
{
	calc_drop_radius ( geometry, lattice, drop );

	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	int n_voxels = nx * ny * nz;

	double x0 = drop.x0, y0 = drop.y0, z0 = drop.z0;

	double raio = drop.raio;

	//--------------- Limites das duas regioes -------------------------------------------------------//

	double r_in_max = drop.frac_in * raio;

	double r_out_min = drop.frac_out * raio;

	if ( r_out_min < raio + 4.0 * drop.espessura ) r_out_min = raio + 4.0 * drop.espessura;

	double r_caixa = 0.45 * ( double ) ( ( nx < ny ) ? nx : ny );

	if ( r_out_min > r_caixa ) r_out_min = r_caixa;

	//--------------- Varredura ----------------------------------------------------------------------//

	double p_in_s = 0.0, p_out_s = 0.0, p_in_s2 = 0.0, p_out_s2 = 0.0;
	double p_ctr = 0.0, p_cnt = 0.0;

	int n_in = 0, n_out = 0, n_ctr = 0;

	#pragma omp parallel for reduction(+:p_in_s,p_out_s,p_in_s2,p_out_s2,n_in,n_out,p_ctr,p_cnt,n_ctr)

	for ( int pos = 0; pos < n_voxels; pos++ )
	{
		if ( geo[pos] == 0 ) continue;

		int pto = geo[pos] - 1;

		double rho_R = density ( lattice.inif_R + pto * nvel );
		double rho_B = density ( lattice.inif_B + pto * nvel );

		if ( ! ( rho_R == rho_R ) || ! ( rho_B == rho_B ) ) continue;

		int z = pos / ( nx * ny );
		int y = ( pos - z * nx * ny ) / nx;
		int x = pos - z * nx * ny - y * nx;

		double dx = ( double ) x - x0;   if ( dx >  0.5 * nx ) dx -= nx;   if ( dx < -0.5 * nx ) dx += nx;
		double dy = ( double ) y - y0;   if ( dy >  0.5 * ny ) dy -= ny;   if ( dy < -0.5 * ny ) dy += ny;
		double dz = ( double ) z - z0;   if ( dz >  0.5 * nz ) dz -= nz;   if ( dz < -0.5 * nz ) dz += nz;

		if ( drop.cilindro ) dz = 0.0;

		double r = sqrt ( dx * dx + dy * dy + dz * dz );

		//  Pressao do sitio.  Ver DROP::cs2_R / cs2_B: modelos de alpha variavel dao a cada fluido
		//  a sua velocidade do som, e ai a pressao nao e c_s^2 rho.

		double p = ( drop.cs2_R >= 0.0 ) ? ( rho_R * drop.cs2_R + rho_B * drop.cs2_B )
		                                 : ( lattice.c_s2 * ( rho_R + rho_B ) );

		if      ( r <  r_in_max  ) { p_in_s  = p_in_s  + p;   p_in_s2  = p_in_s2  + p * p;   n_in  = n_in  + 1; }

		else if ( r >= r_out_min ) { p_out_s = p_out_s + p;   p_out_s2 = p_out_s2 + p * p;   n_out = n_out + 1; }

		if ( r < 1.5 )   { p_ctr = p_ctr + p;   n_ctr = n_ctr + 1; }

		if ( pos == 0 )    p_cnt = p_cnt + p;
	}

	//--------------- Devolve --------------------------------------------------------------------------//

	double p_in  = ( n_in  > 0 ) ? p_in_s  / ( double ) n_in  : 0.0;
	double p_out = ( n_out > 0 ) ? p_out_s / ( double ) n_out : 0.0;

	double var_in  = ( n_in  > 0 ) ? p_in_s2  / ( double ) n_in  - p_in  * p_in  : 0.0;
	double var_out = ( n_out > 0 ) ? p_out_s2 / ( double ) n_out - p_out * p_out : 0.0;

	double fator_forma = drop.cilindro ? 1.0 : 0.5;

	drop.p_in  = p_in;    drop.sd_in  = ( var_in  > 0.0 ) ? sqrt ( var_in  ) : 0.0;
	drop.p_out = p_out;   drop.sd_out = ( var_out > 0.0 ) ? sqrt ( var_out ) : 0.0;

	drop.n_in = n_in;   drop.n_out = n_out;

	drop.delta_p = p_in - p_out;

	drop.sigma = fator_forma * raio * drop.delta_p;

	drop.p_centro = ( n_ctr > 0 ) ? p_ctr / ( double ) n_ctr : 0.0;

	drop.p_canto = p_cnt;

	drop.sigma_linha = fator_forma * drop.raio_linha * ( drop.p_centro - drop.p_canto );

	return drop.sigma;
}

//====================================================================================================================//
