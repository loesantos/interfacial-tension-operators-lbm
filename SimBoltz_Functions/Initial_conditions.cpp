
//====================================================================================================================//


//=============================== Initial Condition of a Drop with a diffuse interface ===============================//
//
//      Gota de fluido vermelho ( raio 'radius', centro x0,y0,z0 ) num banho de fluido azul, ja com o
//      perfil de interface do modelo, em vez de um degrau:
//
//          phi(r) = 1 / ( 1 + exp( k_perfil ( r - radius ) ) )         phi = rho_R / ( rho_R + rho_B )
//
//      Para os modelos com mediadores ( Santos imiscivel ), o perfil de equilibrio e a logistica de
//      inclinacao  k = lambda_ef / cs^2 = 3 lambda_ef ( D3Q19 ), com
//
//          lambda_ef = lambda / ( 1 - 1 / ( 2 tau_m ) )        lambda = parameters.A_fact
//
//      Partir dai poupa o transiente acustico do degrau, que numa caixa periodica reverbera por
//      milhares de passos.  k_perfil <= 0 volta a dar o degrau ( equivale a initial_conditions_sphere ).
//
//      'cilindro' = true trata a gota como um cilindro de eixo z ( o caso 2D, nz = 1 ): a distancia
//      e medida so em x e y.  Com false a gota e uma esfera.
//
//      'delta_rho' e um acrescimo de densidade DENTRO da gota, para ja entrar com o salto de Laplace
//      ( delta_rho = delta_p / cs^2 ).  ATENCAO: densidade a mais e massa a mais, e a massa se
//      conserva -- a gota de equilibrio fica MAIOR que 'radius'.  delta_rho = 0 e o caso normal.
//
//      'fase_simetrica' escolhe a convencao dos mediadores: false ( padrao ) semeia
//      phi = rho_R / rho, em [0,1], que e a do modelo de Santos; true semeia
//      rho^N = ( rho_R - rho_B ) / rho, em [-1,1], que e a de Spencer-Halliday-Care.
//
//      Input : radius, k_perfil, delta_rho, x0, y0, z0, vx, vy, vz, cilindro, geometry, lattice,
//              parameters, fase_simetrica
//      Output: distribution functions ( inif_R, inif_B ), mediators ( inif_m ), scalar field ( ini_psi )
//
//====================================================================================================================//

void initial_conditions_drop ( double radius, double k_perfil, double delta_rho, int x0, int y0, int z0,
							double vx, double vy, double vz, bool cilindro, GEOMETRY geometry, LATTICE lattice,
							PARAMETERS parameters, bool fase_simetrica, int eq_modo )
{
	int *geo = geometry.ini;

	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;

	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{
				int pos = x + y * nx + z * nx * ny;

				//  No solido nao ha populacao, mas o campo escalar existe em todos os voxels e e
				//  lido pelos vizinhos em gradient(): recebe a molhabilidade, como fazem as outras
				//  condicoes iniciais deste arquivo.

				if ( geo[pos] == 0 )
				{
					if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = parameters.wett_R;

					continue;
				}

				double dx = ( double ) ( x - x0 );
				double dy = ( double ) ( y - y0 );
				double dz = cilindro ? 0.0 : ( double ) ( z - z0 );

				double r = sqrt ( dx * dx + dy * dy + dz * dz );

				//  Perfil.  O expoente e limitado para nao estourar exp() longe da interface.

				double phi;

				if ( k_perfil > 0.0 )
				{
					double arg = k_perfil * ( r - radius );

					if ( arg >  60.0 ) arg =  60.0;
					if ( arg < -60.0 ) arg = -60.0;

					phi = 1.0 / ( 1.0 + exp ( arg ) );
				}

				else phi = ( r < radius ) ? 1.0 : 0.0;

				//  A densidade total interpola entre os dois fluidos puros; delta_rho e o salto de
				//  Laplace, que so existe dentro da gota.

				double rho_tot = parameters.rho_ini_B + phi * ( parameters.rho_ini_R - parameters.rho_ini_B )
							   + phi * delta_rho;

				int pto = geo[pos] - 1;

				//  eq_modo escolhe COM QUAL EQUILIBRIO as populacoes sao semeadas.  Importa quando
				//  o modelo nao usa o equilibrio classico: com w_i a pressao inicial sairia rho/3
				//  em vez da pressao do modelo, e a caixa inteira arrancaria com uma onda acustica
				//  que reverbera por milhares de passos numa caixa periodica.
				//
				//      EQ_CLASSICO : dist_eq , o de sempre;
				//      EQ_ALFA     : equilibrio de alpha variavel ( eq. (8) de Liu-Valocchi-Kang,
				//                    que e dist_eq_LYL com hi_order = 0 ) -- so difere se alpha != 1/3;
				//      EQ_SAITO    : equilibrio generalizado de Saito ( eq. 56 ), pelos momentos
				//                    centrais.  Precisa da pressao, p = rho_R cs2_R + rho_B cs2_B.

				double *fR = lattice.inif_R + pto * nvel;
				double *fB = lattice.inif_B + pto * nvel;

				const double rho_R = phi * rho_tot;
				const double rho_B = ( 1.0 - phi ) * rho_tot;

				if ( eq_modo == EQ_ALFA )
				{
					dist_eq_LYL ( fR, vx, vy, vz, rho_R, parameters.alfa_R, 0.0, lattice );
					dist_eq_LYL ( fB, vx, vy, vz, rho_B, parameters.alfa_B, 0.0, lattice );
				}

				else if ( eq_modo == EQ_SAITO )
				{
					//  Em Saito o equilibrio e da MISTURA: um so f^eq, com rho e p totais, repartido
					//  entre as duas cores pela concentracao.  E o que a recoloracao pressupoe.

					const double rho_t = rho_R + rho_B;

					const double p_t = rho_R * parameters.cs2_R + rho_B * parameters.cs2_B;

					double f_eq[nvel];

					dist_eq_saito ( f_eq, vx, vy, vz, rho_t, p_t, lattice );

					const double cR = ( rho_t > 0.0 ) ? rho_R / rho_t : 0.0;

					for ( int i = 0; i < nvel; i++ )
					{
						fR[i] = cR * f_eq[i];
						fB[i] = ( 1.0 - cR ) * f_eq[i];
					}
				}

				else
				{
					dist_eq ( fR, vx, vy, vz, rho_R, lattice );
					dist_eq ( fB, vx, vy, vz, rho_B, lattice );
				}

				//  Convencoes do resto do arquivo: ini_psi guarda ( rho_R - rho_B ) / rho, em [-1,1].
				//  Ja os mediadores mudam de modelo para modelo:
				//
				//      fase_simetrica = false  ( Santos ) :  phi = rho_R / rho , em [0,1]
				//      fase_simetrica = true   ( SHC )    :  rho^N = ( rho_R - rho_B ) / rho , em [-1,1]
				//
				//  Nos dois casos os mediadores sao reemitidos do zero a cada passo, entao esta
				//  semeadura so importa para um diagnostico feito antes da primeira colisao.

				double fase_med = fase_simetrica ? ( 2.0 * phi - 1.0 ) : phi;

				if ( lattice.inif_m != nullptr )
				{
					double *f_m = lattice.inif_m + pto * nvel;

					emite_mediadores ( f_m, fase_med, lattice );
				}

				if ( lattice.ini_psi != nullptr ) lattice.ini_psi[pos] = 2.0 * phi - 1.0;
			}
		}
	}
}

//====================================================================================================================//
