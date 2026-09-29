
//=============================== Read the inicialization file =======================================================//
//
//      Input:
//      Output: geometry file name, pixel dimension, steps, files, relaxation time,
//              initial density
//
//====================================================================================================================//

void read_data ( GEOMETRY& geometry, PARAMETERS& parameters )
{
    //----------------------------------------------------------------------------------------------------------------//

    string name_in = "data_in.txt";

    ifstream f_in( name_in );

	//----------------------------------------------------------------------------------------------------------------//

    string name_out = "dat_out.txt";

    cout << "\nName of the output file: " << name_out << endl;

    //----------------------------------------------------------------------------------------------------------------//

    ofstream fdat( name_out );

    fdat << "Name of the output file: " << name_out << endl;

    //----------------------------------------------------------------------------------------------------------------//

    f_in >> geometry.file;

    cout << "\nName of the geometry file: " << geometry.file << endl;
    fdat << "\nName of the geometry file: " << geometry.file << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_ftesc;
    
    f_in >> st_ftesc;

    f_in >> geometry.ftesc;  // read the pixel's dimension ( m )

    cout << "\nPixel dimension = " << geometry.ftesc << " m" << endl; 
    fdat << "\nPixel dimension = " << geometry.ftesc << " m" << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_steps;
    
    f_in >> st_steps;

    f_in >> parameters.n_steps;

    cout << "\nNumber of steps: " << parameters.n_steps << endl;
    fdat << "\nNumber of steps: " << parameters.n_steps << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_files;
    
    f_in >> st_files;

    f_in >> parameters.n_files;

    cout << "\nNumber of files: " << parameters.n_files << endl;
    fdat << "\nNumber of files: " << parameters.n_files << endl;
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_threads;
    
    f_in >> st_threads;

    f_in >> parameters.n_threads;

    cout << "\nNumber of threads: " << parameters.n_threads << endl;
    fdat << "\nNumber of threads: " << parameters.n_threads << endl;

    //----------------------------------------------------------------------------------------------------------------//

    string st_tau;
    
    f_in >> st_tau;

    f_in >> parameters.tau;
    
    parameters.visc = (1./3.) * ( parameters.tau - 0.5 );

    cout << "\nRelaxation time: " << parameters.tau << endl;
    fdat << "\nRelaxation time: " << parameters.tau << endl;
    
    cout << "\nKinematic viscosity = " << parameters.visc << endl;
    fdat << "\nKinematic viscosity = " << parameters.visc << endl;
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_nh;
    
    f_in >> st_tau_nh;

    f_in >> parameters.tau_nh;

    cout << "\nRelaxation time (non-hidrodynamic): " << parameters.tau_nh << endl;
    fdat << "\nRelaxation time (non-hidrodynamic): " << parameters.tau_nh << endl;

	//----------------------------------------------------------------------------------------------------------------//
	
	string st_rho;
    
    f_in >> st_rho;

    f_in >> parameters.rho_ini;

    cout << "\nDensity: " << parameters.rho_ini << endl;
    fdat << "\nDensity: " << parameters.rho_ini << endl;
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_R;
    
    f_in >> st_tau_R;

    f_in >> parameters.tau_R;

    parameters.visc_R = ( 1.0 / 3.0 ) * ( parameters.tau_R - 0.5 );

    if ( parameters.tau_R )
    {
		cout << "\nRelaxation time (Red): " << parameters.tau_R << " => viscosity R = " << parameters.visc_R << endl;		
		fdat << "\nRelaxation time (Red): " << parameters.tau_R << " => viscosity R = " << parameters.visc_R << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_B;
    
    f_in >> st_tau_B;

    f_in >> parameters.tau_B;

    parameters.visc_B = ( 1.0 / 3.0 ) * ( parameters.tau_B - 0.5 );

    if ( parameters.tau_B )
    {
		cout << "\nRelaxation time (Blue): " << parameters.tau_B << " => viscosity B = " << parameters.visc_B << endl;		
		fdat << "\nRelaxation time (Blue): " << parameters.tau_B << " => viscosity B = " << parameters.visc_B << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_tau_m;
    
    f_in >> st_tau_m;

    f_in >> parameters.tau_m;
    
    if ( parameters.tau_R )
    {
		cout << "\nRelaxation time (mixture): " << parameters.tau_m << endl;		
		fdat << "\nRelaxation time (mixture): " << parameters.tau_m << endl;
	}
	
    //----------------------------------------------------------------------------------------------------------------//

    string st_rho_ini_R;
    
    f_in >> st_rho_ini_R;

    f_in >> parameters.rho_ini_R;
    
    if ( parameters.tau_R )
    {
		cout << "\nInitial density (Red): " << parameters.rho_ini_R << endl;		
		fdat << "\nInitial density (Red): " << parameters.rho_ini_R << endl;
	}
	
	//----------------------------------------------------------------------------------------------------------------//

    string st_rho_ini_B;
    
    f_in >> st_rho_ini_B;

    f_in >> parameters.rho_ini_B;
    
    if ( parameters.tau_B )
    {
		cout << "\nInitial density (Blue): " << parameters.rho_ini_B << endl;		
		fdat << "\nInitial density (Blue): " << parameters.rho_ini_B << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_fat_int;
    
    f_in >> st_fat_int;

    f_in >> parameters.A_fact;
    
    if ( parameters.tau_R )
    {
		cout << "\nInterfacial factor[0;0.4]: " << parameters.A_fact << endl;		
		fdat << "\nInterfacial factor[0;0.4]: " << parameters.A_fact << endl;
	}
	
	//----------------------------------------------------------------------------------------------------------------//

    string st_fat_int_R;
    
    f_in >> st_fat_int_R;

    f_in >> parameters.A_fact_R;
    
    if ( parameters.tau_R )
    {
		cout << "\nInteraction factor (Red) [0;0.4]: " << parameters.A_fact_R << endl;		
		fdat << "\nInteraction factor (Red) [0;0.4]: " << parameters.A_fact_R << endl;
	}
	
	//----------------------------------------------------------------------------------------------------------------//

    string st_fat_int_B;
    
    f_in >> st_fat_int_B;

    f_in >> parameters.A_fact_B;
    
    if ( parameters.tau_B )
    {
		cout << "\nInteraction factor (Blue) [0;0.4]: " << parameters.A_fact_B << endl;		
		fdat << "\nInteraction factor (Blue) [0;0.4]: " << parameters.A_fact_B << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

	string st_recoll;

    f_in >> st_recoll;

    f_in >> parameters.recoll;
    
    if ( parameters.tau_R )
    {
		cout << "\nRecolloring factor[0;1]: " << parameters.recoll << endl;		
		fdat << "\nRecolloring factor[0;1]: " << parameters.recoll << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    string st_wett_R;
    
    f_in >> st_wett_R;
    
    f_in >> parameters.wett_R;
    
    if ( parameters.tau_R )
    {
		cout << "\nRed_wall interaction [0.0 ; 1.0] : " << parameters.wett_R << endl;		
		fdat << "\nRed_wall interaction [0.0 ; 1.0] : " << parameters.wett_R << endl;
	}
    
    //----------------------------------------------------------------------------------------------------------------//

    f_in.close();

    fdat.close();
}

//====================================================================================================================//


//=============================== Read the geometry file =============================================================//
//
//      Input: geometry file name, pointer to adress of begining of the geometry
//      Output: number of fluid points
//
//====================================================================================================================//

int read_geo ( string nome_geo, int *meio, int pts_in, int pts_out )
{
    int nx, ny, nz;

    ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;

    cout << "\nTamanho:  x = " << nx << ";  y = "   << ny << ";  z = "   << nz << endl;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    //--------------- Acrescenta layers ----------------------------------------------------------//

    nx = nx + pts_in + pts_out;

    //--------------------------------------------------------------------------------------------//

    int poros = 0;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;

                if ( x >= pts_in && x < nx - pts_out )
                {
                    fmatriz >> meio[pos];

                    if ( meio[pos] )
                    {
                        poros++;

                        meio[pos] = poros;
                    }
                }
                else if ( x == 0 )
                {
                    meio[pos] = 0;
                }
                
                else if ( x == nx - 1 )
                {
                    meio[pos] = 0;
                }
                
                else
                {
                    poros++;

                    meio[pos] = poros;
                }
            }
        }
    }

    fmatriz.close();

    return poros;
}

//====================================================================================================================//


//=============================== Defines a D3Q19 lattice ============================================================//
//
//      Input: pointer to the vectors
//      Output:
//
//====================================================================================================================//

void def_lattice_d3q19 ( LATTICE& lattice )
{
    //------------- |ci| = 0 ----------------------//
    int i = 0;
	
    lattice.c_i[ i * dim + 0 ] = 0; 
    lattice.c_i[ i * dim + 1 ] = 0;
    lattice.c_i[ i * dim + 2 ] = 0;

    //------------- |ci| = 1 ----------------------//
	i = 1;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 2;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 3;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 4;

    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 5;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 6;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

    //------------ |c| = sqrt(2) -----------------//

	i = 7;

    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 8;

    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 9;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 10;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  0;

	i = 11;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 12;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 13;
	
    lattice.c_i[ i * dim + 0 ] =  1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] = -1;

	i = 14;
	
    lattice.c_i[ i * dim + 0 ] = -1;
    lattice.c_i[ i * dim + 1 ] =  0;
    lattice.c_i[ i * dim + 2 ] =  1;
    
    i = 15;
    
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] = -1;
    
    i = 16;
    	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] =  1;
    
	i = 17;    
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] = -1;
    lattice.c_i[ i * dim + 2 ] =  1;

	i = 18;
	
    lattice.c_i[ i * dim + 0 ] =  0;
    lattice.c_i[ i * dim + 1 ] =  1;
    lattice.c_i[ i * dim + 2 ] = -1;

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    lattice.w[0] = 1. / 3.;

    for ( int i = 1; i < 7; i++ ) lattice.w[i] =  1. / 18.;

    for ( int i = 7; i < nvel; i++ ) lattice.w[i] = 1. / 36.; 
    
    lattice.c_s2 = 1./3.; // sound velocity
    
    lattice.one_over_c_s2 = 3.; // inverse of the sound velocity
}

//====================================================================================================================//


//=============================== Define sites used in the streaming process =========================================//
//
//      Input: geometry, lattice
//      Output: addresses, *ini_dir
//
//====================================================================================================================//

void def_dir_prop ( GEOMETRY geometry, LATTICE& lattice )
{
	int *ini_meio = geometry.ini; 
	int *ini_dir = lattice.ini_stream;
	
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
	bool* solid = nullptr;
	
	double *qlost_x = nullptr;
	double *qlost_y = nullptr;
	double *qlost_z = nullptr;
	
	double *ini_qlost_x = nullptr;
	double *ini_qlost_y = nullptr;
	double *ini_qlost_z = nullptr;
	
	if ( lattice.ini_mom_x != nullptr ) ini_qlost_x = lattice.ini_mom_x;
	if ( lattice.ini_mom_y != nullptr ) ini_qlost_y = lattice.ini_mom_y;
	if ( lattice.ini_mom_z != nullptr ) ini_qlost_z = lattice.ini_mom_z;

	if ( dim == 2 )
	{
		//---------------- Define os passos para propagação ------------------------------------------//

		double  step_x[nvel];
		double  step_y[nvel];
		
		int steps[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			double cx = lattice.c_i[ i * dim + 0 ];
			double cy = lattice.c_i[ i * dim + 1 ];

			int cx_int =  round_number ( cx );
			int cy_int =  round_number ( cy );

			steps[i] = abs ( cx_int );

			if ( abs ( cy_int ) > steps[i] )
			{
				steps[i] = abs ( cy_int );
			}

			step_x[i] = cx / steps[i];
			step_y[i] = cy / steps[i];
		}

		//-------------- Encontra as direções contrárias ---------------------------------------------//

		int i_op[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			int cx_i = round_number ( lattice.c_i[ i * dim + 0 ] );
			int cy_i = round_number ( lattice.c_i[ i * dim + 1 ] );

			for ( int j = 1; j < nvel; j++ )
			{
				int cx_j = round_number ( lattice.c_i[ j * dim + 0 ] );
				int cy_j = round_number ( lattice.c_i[ j * dim + 1 ] );

				if ( cx_j == -cx_i  && cy_j == -cy_i )
				{
					i_op[i] = j;
				}
			}
		}

		//--------------------------------------------------------------------------------------------//

		for ( int y = 0; y < ny; y++ )
		{
			for ( int x = 0; x < nx; x++ )
			{
				int* meio_local = ini_meio + x + y * nx;

				if ( *meio_local )
				{
					if ( lattice.ini_mom_x != nullptr )
					{
						qlost_x = ini_qlost_x + ( *meio_local - 1 ) * nvel;
						qlost_x[0] = 0;
					}
					if ( lattice.ini_mom_y != nullptr )
					{
						qlost_y = ini_qlost_y + ( *meio_local - 1 ) * nvel;
						qlost_y[0] = 0;
					}
					
					if ( lattice.ini_solid != nullptr )
					{
						solid = lattice.ini_solid + ( *meio_local - 1 ) * nvel;
						solid[0] = 0;
					}					 
										
					int* dir = ini_dir + ( *meio_local - 1 ) * nvel;

					int* meio_prop = meio_local;
					
					dir[0] = ( *meio_prop - 1 ) * nvel;

					//--------------------------------------------------------------------------------//

					for ( int i = 1; i < nvel; i++ )
					{
						if ( lattice.ini_solid != nullptr ) solid[i] = 0; // fluid
						
						int inv = 0; // número de inversões

						double  stpx = step_x[i];
						double  stpy = step_y[i];

						int x_prop = x;
						int y_prop = y;

						double  x_prop_f = ( double ) x;
						double  y_prop_f = ( double ) y;

						for ( int stp = 1; stp < steps[i] + 1; stp++ )
						{
							x_prop_f = x_prop_f + stpx;
							y_prop_f = y_prop_f + stpy;

							x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
							y_prop = ( round_number ( y_prop_f ) + ny ) % ny;

							meio_prop = ini_meio + x_prop + y_prop * nx;

							if ( *meio_prop == 0 )
							{
								stpx = -stpx;
								stpy = -stpy;

								x_prop_f = x_prop_f + stpx;
								y_prop_f = y_prop_f + stpy;

								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;

								meio_prop = ini_meio + x_prop + y_prop * nx;

								inv++;
								
								if ( lattice.ini_solid != nullptr ) solid[i] = 1; // solid
							}
						}

						if ( inv % 2 == 0 )
						{
							dir[i] = i + ( *meio_prop - 1 ) * nvel;
						}
						else
						{
							dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
							
							if ( lattice.ini_mom_x != nullptr )	qlost_x[i] = - 2. * stpx;
							if ( lattice.ini_mom_y != nullptr )	qlost_y[i] = - 2. * stpy;
						}
					}

					//----------------------------------------------------------------------------//
				}
			}
		}
	}
	
	else
	{
		//---------------- Define os passos para propagação ------------------------------------------//

		double  step_x[nvel];
		double  step_y[nvel];
		double  step_z[nvel];

		int steps[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			double cx = lattice.c_i[ i * dim + 0 ];
			double cy = lattice.c_i[ i * dim + 1 ];
			double cz = lattice.c_i[ i * dim + 2 ];

			int cx_int =  round_number ( cx );
			int cy_int =  round_number ( cy );
			int cz_int =  round_number ( cz );

			steps[i] = abs ( cx_int );

			if ( abs ( cy_int ) > steps[i] )
			{
				steps[i] = abs ( cy_int );
			}

			if ( abs ( cz_int ) > steps[i] )
			{
				steps[i] = abs ( cz_int );
			}

			step_x[i] = cx / steps[i];
			step_y[i] = cy / steps[i];
			step_z[i] = cz / steps[i];
		}

		//-------------- Encontra as direções contrárias ---------------------------------------------//

		int i_op[nvel];

		for ( int i = 1; i < nvel; i++ )
		{
			int cx_i = round_number ( lattice.c_i[ i * dim + 0 ] );
			int cy_i = round_number ( lattice.c_i[ i * dim + 1 ] );
			int cz_i = round_number ( lattice.c_i[ i * dim + 2 ] );

			for ( int j = 1; j < nvel; j++ )
			{
				int cx_j = round_number ( lattice.c_i[ j * dim + 0 ] );
				int cy_j = round_number ( lattice.c_i[ j * dim + 1 ] );
				int cz_j = round_number ( lattice.c_i[ j * dim + 2 ] );

				if ( cx_j == -cx_i  && cy_j == -cy_i && cz_j == -cz_i )
				{
					i_op[i] = j;
				}
			}
		}

		//--------------------------------------------------------------------------------------------//

		for ( int z = 0; z < nz; z++ )
		{
			for ( int y = 0; y < ny; y++ )
			{
				for ( int x = 0; x < nx; x++ )
				{
					int* meio_local = ini_meio + x + y * nx + z * nx * ny;

					if ( *meio_local )
					{
						if ( lattice.ini_mom_x != nullptr )
						{
							qlost_x = ini_qlost_x + ( *meio_local - 1 ) * nvel;
							qlost_x[0] = 0;
						}
						if ( lattice.ini_mom_y != nullptr )
						{
							qlost_y = ini_qlost_y + ( *meio_local - 1 ) * nvel;
							qlost_y[0] = 0;
						}
						if ( lattice.ini_mom_z != nullptr )
						{
							qlost_z = ini_qlost_z + ( *meio_local - 1 ) * nvel;
							qlost_z[0] = 0;
						}
						
						if ( lattice.ini_solid != nullptr )
						{
							solid = lattice.ini_solid + ( *meio_local - 1 ) * nvel;
							solid[0] = 0;
						}
					
						int* dir = ini_dir + ( *meio_local - 1 ) * nvel;

						int* meio_prop = meio_local;
					
						dir[0] = ( *meio_prop - 1 ) * nvel;

						//----------------------------------------------------------------------------//

						for ( int i = 1; i < nvel; i++ )
						{
							if ( lattice.ini_solid != nullptr ) solid[i] = 0; // fluid
							
							int inv = 0; // número de inversões

							double  stpx = step_x[i];
							double  stpy = step_y[i];
							double  stpz = step_z[i];

							int x_prop = x;
							int y_prop = y;
							int z_prop = z;

							double  x_prop_f = ( double ) x;
							double  y_prop_f = ( double ) y;
							double  z_prop_f = ( double ) z;

							for ( int stp = 1; stp < steps[i] + 1; stp++ )
							{
								x_prop_f = x_prop_f + stpx;
								y_prop_f = y_prop_f + stpy;
								z_prop_f = z_prop_f + stpz;

								x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
								y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
								z_prop = ( round_number ( z_prop_f ) + nz ) % nz;

								meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

								if ( *meio_prop == 0 )
								{
									stpx = -stpx;
									stpy = -stpy;
									stpz = -stpz;

									x_prop_f = x_prop_f + stpx;
									y_prop_f = y_prop_f + stpy;
									z_prop_f = z_prop_f + stpz;

									x_prop = ( round_number ( x_prop_f ) + nx ) % nx;
									y_prop = ( round_number ( y_prop_f ) + ny ) % ny;
									z_prop = ( round_number ( z_prop_f ) + nz ) % nz;

									meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

									inv++;
					
									if ( lattice.ini_solid != nullptr ) solid[i] = 1;
								}
							}

							if ( inv % 2 == 0 )
							{
								dir[i] = i + ( *meio_prop - 1 ) * nvel;
							}
							else
							{
								dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
																
								if ( lattice.ini_mom_x != nullptr )	qlost_x[i] = - 2. * stpx;
								if ( lattice.ini_mom_y != nullptr )	qlost_y[i] = - 2. * stpy;
								if ( lattice.ini_mom_z != nullptr )	qlost_z[i] = - 2. * stpz;
							}
						}

						//----------------------------------------------------------------------------//
					}
				}
			}
		}
	}
}

//====================================================================================================================//


//=============================== Return the dot product =============================================================//
//
//      Input: two vetors
//      Output: dot product
//
//====================================================================================================================//

#pragma acc routine seq
double  dot_product ( double  vet_1[dim], double  vet_2[dim] )
{

    double  result = 0.0;

    for ( int i = 0; i < dim; i++ )
    {
        result = result + vet_1[i] * vet_2[i];
    }

    return result;
}

//====================================================================================================================//


//=============================== Return the equilibrium distribution function =======================================//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

#pragma acc routine seq
void dist_eq ( double* feq, double vx, double vy, double vz, double rho, LATTICE lattice )
{
	if ( dim == 2 )
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;

		double vquad = ( vx * vx + vy * vy );
		
		double c_i[dim];
		
		double one_over_c_s2 = lattice.one_over_c_s2;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2
									+ 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2
									- 0.5 * vquad * one_over_c_s2 );
		}
	}
	
	else
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;

		double vquad = ( vx * vx + vy * vy + vz * vz );
		
		double c_i[dim];
		
		double one_over_c_s2 = lattice.one_over_c_s2;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];        
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2
									+ 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2
									- 0.5 * vquad * one_over_c_s2 );
		}
	}
}

//====================================================================================================================//


//=============================== Return the equilibrium distribution function =======================================//
//
//      Input: velocities, density, lattice
//      Output: equilibrium distribution
//
//====================================================================================================================//

#pragma acc routine seq
void dist_eq_sixth ( double* feq, double vx, double vy, double vz, double rho, LATTICE lattice )
{
	if ( dim == 2 )
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		
		double c_i[dim];
		
		double u_u = ( vx * vx + vy * vy );
		
		double one_over_c_s2 = lattice.one_over_c_s2;    
		
		double one_over_c_s4 = one_over_c_s2 * one_over_c_s2;
		
		double one_over_c_s6 = one_over_c_s4 * one_over_c_s2;
		
		double one_over_c_s8 = one_over_c_s4 * one_over_c_s4;

		double one_over_c_s10 = one_over_c_s8 * one_over_c_s2;
		
		double one_over_c_s12 = one_over_c_s8 * one_over_c_s4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2 + 0.5 * cv*cv * one_over_c_s4
							- 0.5 * u_u * one_over_c_s2 + (1./6.)*cv*cv*cv*one_over_c_s6
							- 0.5 * cv * u_u * one_over_c_s4 + (1./24.)*cv*cv*cv*cv * one_over_c_s8 
							- 0.25 * cv*cv*u_u * one_over_c_s6 + (1./8.)*u_u*u_u * one_over_c_s4 
							+ (1./120.)*cv*cv*cv*cv*cv * one_over_c_s10- (1./12.) * u_u * cv*cv*cv * one_over_c_s8 
							+ (1./8.)*u_u*u_u *cv * one_over_c_s6	+ (1./720.) * cv*cv*cv*cv*cv*cv * one_over_c_s12 
							- (1./48.) * u_u * cv*cv*cv*cv * one_over_c_s10 + (1./16.) * u_u*u_u * cv*cv * one_over_c_s8 
							- (1./48.) * u_u*u_u*u_u * one_over_c_s6 );
		}
	}
	
	else
	{
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;
		
		double c_i[dim];
		
		double u_u = ( vx * vx + vy * vy + vz * vz );
		
		double one_over_c_s2 = lattice.one_over_c_s2;    
		
		double one_over_c_s4 = one_over_c_s2 * one_over_c_s2;
		
		double one_over_c_s6 = one_over_c_s4 * one_over_c_s2;
		
		double one_over_c_s8 = one_over_c_s4 * one_over_c_s4;
		
		double one_over_c_s10 = one_over_c_s8 * one_over_c_s2;
		
		double one_over_c_s12 = one_over_c_s8 * one_over_c_s4;

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];        
			
			double cv = dot_product ( c_i, v );

			feq[i] = lattice.w[i] * rho * ( 1.0 + cv * one_over_c_s2 + 0.5 * cv*cv * one_over_c_s4
							- 0.5 * u_u * one_over_c_s2 + (1./6.)*cv*cv*cv*one_over_c_s6
							- 0.5 * cv * u_u * one_over_c_s4 + (1./24.)*cv*cv*cv*cv * one_over_c_s8 
							- 0.25 * cv*cv*u_u * one_over_c_s6 + (1./8.)*u_u*u_u * one_over_c_s4 
							+ (1./120.)*cv*cv*cv*cv*cv * one_over_c_s10- (1./12.) * u_u * cv*cv*cv * one_over_c_s8 
							+ (1./8.)*u_u*u_u *cv * one_over_c_s6	+ (1./720.) * cv*cv*cv*cv*cv*cv * one_over_c_s12 
							- (1./48.) * u_u * cv*cv*cv*cv * one_over_c_s10 + (1./16.) * u_u*u_u * cv*cv * one_over_c_s8 
							- (1./48.) * u_u*u_u*u_u * one_over_c_s6 );							  
		}
	}
	
}

//===================================================================================================================//


//=============================== Compute momentum ===================================================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//====================================================================================================================//

#pragma acc routine seq
void momentum ( double *f, double& mx, double& my, double& mz, LATTICE lattice )
{
	if ( dim == 2 )
	{
		mz = 0.0;
		
		mx = 0.0;
		my = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{
			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
		}
	}
	 
	else
	{
		mx = 0.0;
		my = 0.0;
		mz = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{
			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
			mz = mz + lattice.c_i[ i * dim + 2 ] * f[i];
		}
	}

}

//====================================================================================================================//


//=============================== Calculate density and velocities ===================================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//====================================================================================================================//

#pragma acc routine seq
void calcula ( double *f, double& vx, double& vy, double& vz, double& rho, LATTICE lattice )
{
	if ( dim == 2 )
	{
		vz = 0.0;
		
		double mx = 0.0;
		double my = 0.0;

		double one_over_rho;

		rho = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{

			rho = rho + f[i];

			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
		}

		if ( rho )
		{
			one_over_rho = 1.0 / ( rho );

			vx = mx * one_over_rho;
			vy = my * one_over_rho;
		}
		else
		{
			vx = 0.;
			vy = 0.;
		}
	}
	 
	else
	{
		double mx = 0.0;
		double my = 0.0;
		double mz = 0.0;

		double one_over_rho;

		rho = 0.0;

		for ( int i = 0 ; i < nvel; i++ )
		{

			rho = rho + f[i];

			mx = mx + lattice.c_i[ i * dim + 0 ] * f[i];
			my = my + lattice.c_i[ i * dim + 1 ] * f[i];
			mz = mz + lattice.c_i[ i * dim + 2 ] * f[i];
		}

		if ( rho )
		{
			one_over_rho = 1.0 / ( rho );

			vx = mx * one_over_rho;
			vy = my * one_over_rho;
			vz = mz * one_over_rho;
		}
		else
		{
			vx = 0.;
			vy = 0.;
			vz = 0.;
		}
	}
}

//====================================================================================================================//


//=============================== Propagation step for one site ======================================================//
//
//      Input: lattice, site
//      Output: 
//
//====================================================================================================================//

#pragma acc routine seq
void propag_site ( LATTICE lattice, int pto, double &mx, double &my, double &mz )
{

	double *qlost_x = lattice.ini_mom_x;
	double *qlost_y = lattice.ini_mom_y;
	double *qlost_z = lattice.ini_mom_z;
	
	double sum_qlost_x = 0.0;
	double sum_qlost_y = 0.0;
	double sum_qlost_z = 0.0;

	if ( lattice.inif != nullptr ) // Monophasic
	{
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f = lattice.inif + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;
		
		if ( lattice.ini_mom_x != nullptr ) qlost_x = lattice.ini_mom_x + ( pto ) * nvel;
		if ( lattice.ini_mom_y != nullptr ) qlost_y = lattice.ini_mom_y + ( pto ) * nvel;
		if ( lattice.ini_mom_z != nullptr ) qlost_z = lattice.ini_mom_z + ( pto ) * nvel;

		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_new = lattice.inif_new + dir[i];

			*f_new = f[i];
			
			if ( lattice.ini_mom_x != nullptr ) sum_qlost_x = sum_qlost_x + f[i] * qlost_x[i];
			if ( lattice.ini_mom_y != nullptr ) sum_qlost_y = sum_qlost_y + f[i] * qlost_y[i];
			if ( lattice.ini_mom_z != nullptr ) sum_qlost_z = sum_qlost_z + f[i] * qlost_z[i];
		}
	}
    else 
    {
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f_R = lattice.inif_R + ( pto ) * nvel;
		double *f_B = lattice.inif_B + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;
		
		if ( lattice.ini_mom_x != nullptr ) qlost_x = lattice.ini_mom_x + ( pto ) * nvel;
		if ( lattice.ini_mom_y != nullptr ) qlost_y = lattice.ini_mom_y + ( pto ) * nvel;
		if ( lattice.ini_mom_z != nullptr ) qlost_z = lattice.ini_mom_z + ( pto ) * nvel;

		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_R_new = lattice.inif_R_new + dir[i];
			double* f_B_new = lattice.inif_B_new + dir[i];

			*f_R_new = f_R[i];
			*f_B_new = f_B[i];
			
			if ( lattice.ini_mom_x != nullptr ) sum_qlost_x = sum_qlost_x + ( f_R[i] + f_B[i] ) * qlost_x[i];
			if ( lattice.ini_mom_y != nullptr ) sum_qlost_y = sum_qlost_y + ( f_R[i] + f_B[i] ) * qlost_y[i];
			if ( lattice.ini_mom_z != nullptr ) sum_qlost_z = sum_qlost_z + ( f_R[i] + f_B[i] ) * qlost_z[i];
		}
	}
    
    mx = sum_qlost_x;
    my = sum_qlost_y;
    mz = sum_qlost_z;
	
}

//====================================================================================================================//


//=============================== Emission of the mediators for one site ==============================================//
//
//      Input: mediator populations of the site, the phase carried by them, lattice
//      Output: mediator populations, ready to be propagated
//
//      Emite  M_i = w_i fase / cs^2 .  Depois da propagacao, M_i no sitio x vale  w_i fase( x - c_i ) ,
//      e entao
//
//          soma_i M_i c_i  =  ( 1 / cs^2 ) soma_i w_i fase( x - c_i ) c_i  =  - grad( fase )
//
//      porque  soma_i w_i c_i_alfa c_i_beta = cs^2 delta_alfa_beta .  O 1/cs^2 e o que faz o momento
//      dos mediadores ser o GRADIENTE, e nao cs^2 vezes ele.
//
//      ATENCAO -- correcao de 2026-08-26.  Ate esta data a emissao era  M_i = w_i fase , sem o
//      1/cs^2, e o momento dos mediadores saia cs^2 = 1/3 vezes menor que o gradiente.  Como as
//      sobrecargas dos operadores que recebem ( x, y, z, geometry ) usam gradient(), que ja devolve
//      o gradiente verdadeiro, as duas rotas do MESMO operador davam tensoes interfaciais na razao
//      de 3 para 1 com o mesmo fat_R_B.  Medido num operador de perturbacao: sigma = 0.014944 pelos
//      mediadores contra 0.044832 por gradient(), razao 3.0002.
//
//      Consequencia: fat_R_B passou a valer 3 vezes mais nos programas com mediadores.  Para
//      reproduzir uma tensao interfacial obtida antes desta data, divida o fat_R_B por 3.
//
//      O fator sai de lattice.one_over_c_s2, entao vale para qualquer rede -- D2Q9, D3Q19, D3Q27 e
//      as multivelocidade, onde cs^2 nao e 1/3.
//
//====================================================================================================================//

#pragma acc routine seq
void emite_mediadores ( double *f_m, double fase, LATTICE lattice )
{
	double fat = fase * lattice.one_over_c_s2;

	for ( int i = 0; i < nvel; i++ )   f_m[i] = lattice.w[i] * fat;
}

//====================================================================================================================//


//=============================== Propagation step for one site (mediators) ==========================================//
//
//      Input: lattice, site
//      Output: 
//
//      Nas direcoes que dao em solido, injeta a molhabilidade no lugar do vizinho.  O valor injetado
//      segue a MESMA normalizacao da emissao ( ver emite_mediadores ): a parede se comporta como um
//      fluido de fase wett_R.  Se um dos dois nao tiver o 1/cs^2, wett_R muda de significado e o
//      angulo de contato vai junto.
//
//====================================================================================================================//

#pragma acc routine seq
void propag_site_med ( LATTICE lattice, PARAMETERS parameters, int pto )
{

	if ( lattice.inif_m != nullptr ) // Mediators
	{
		//--------------- Aponta os ponteiros --------------------------------------------------------//

		double *f_m = lattice.inif_m + ( pto ) * nvel;

		int *dir = lattice.ini_stream + ( pto ) * nvel;
		
		bool *solid = lattice.ini_solid + ( pto ) * nvel;
		
		double *W = lattice.w;
		
		//--------------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double* f_m_new = lattice.inif_m_new + dir[i];

			if ( solid[i] )	*f_m_new = W[i] * parameters.wett_R * lattice.one_over_c_s2;  

			else *f_m_new = f_m[i];
			
		}
	}
}

//====================================================================================================================//


//=============================== Record the velocity field (monophasic) =============================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_velocity ( GEOMETRY geometry, LATTICE lattice, unsigned int passo )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    string nomevel = "vel_" + to_string( passo ) + ".vtk";

    ofstream fvel ( nomevel );

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx * ny * nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
				//int x_desloc = ( x + nx / 2 ) % nx;
				//int y_desloc = ( y + 2*ny / 3 ) % ny;
				
				int x_desloc = x;
				int y_desloc = y;
				
                int *meio = geometry.ini + x_desloc + y_desloc * nx + z * ny * nx;

                if ( *meio )
                {
                    if ( lattice.inif != nullptr )
					{
						double *f = lattice.inif + ( *meio - 1 ) * nvel;

						double vx, vy, vz, rho;

						calcula ( f, vx, vy, vz, rho, lattice );

						fvel << vx  << " " << vy << " " << vz << " ";
					}
					else
					{
						double *f_R = lattice.inif_R + ( *meio - 1 ) * nvel;
						double *f_B = lattice.inif_B + ( *meio - 1 ) * nvel;

						double vx_R, vy_R, vz_R, rho_R;

						calcula ( f_R, vx_R, vy_R, vz_R, rho_R, lattice );
						
						double vx_B, vy_B, vz_B, rho_B;

						calcula ( f_B, vx_B, vy_B, vz_B, rho_B, lattice );
						
						double concR = rho_R / ( rho_R + rho_B );
						double concB = 1. - concR;
						
						double vx = concR * vx_R + concB * vx_B;
						double vy = concR * vy_R + concB * vy_B;
						double vz = concR * vz_R + concB * vz_B;

						//  Correcao da meia-forca de Guo, se o programa tiver gravado a forca:
						//  pos-colisao,  rho u = soma_i f_i c_i - F / 2 .  Sem isto o campo
						//  gravado e dominado por |F|/2rho na interface, que nao e velocidade.

						if ( lattice.ini_force != nullptr )
						{
							double rho = rho_R + rho_B;

							const double *F = lattice.ini_force + ( long long ) ( *meio - 1 ) * 3;

							if ( rho > 0.0 )
							{
								vx = vx - 0.5 * F[0] / rho;
								vy = vy - 0.5 * F[1] / rho;
								vz = vz - 0.5 * F[2] / rho;
							}
						}

						fvel << vx  << " " << vy << " " << vz << " ";
					}
                }
                else
                {
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvel << endl;
        }

    }
    fvel.close();
}
//====================================================================================================================//


//=============================== Record the density field (one of two fluids) =======================================//
//
//      Input: geometry, lattice, step
//      Output:
//
//====================================================================================================================//

void rec_density ( string name, GEOMETRY geometry, double* ini_f, unsigned int passo )
{
	int nx = geometry.nx;
	int ny = geometry.ny;
	int nz = geometry.nz;
	
    string name_rho = name + to_string( passo ) + ".vtk";

    ofstream frho ( name_rho );

    frho << "# vtk DataFile Version 2.0" << endl;
    frho << "Densidade" << endl;
    frho << "ASCII" << endl;
    frho << "DATASET STRUCTURED_POINTS" << endl;
    frho << "DIMENSIONS " << geometry.nx << " " << geometry.ny << " " << geometry.nz << endl;
    frho << "ASPECT_RATIO 1 1 1" << endl;
    frho << "ORIGIN 0 0 0" << endl;
    frho << "POINT_DATA " << geometry.nx * geometry.ny * geometry.nz << endl;
    frho << "SCALARS densidade double" << endl;
    frho << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
				//int x_desloc = ( x + nx / 2 ) % nx;
				//int y_desloc = ( y + 2*ny / 3 ) % ny;
				
				int x_desloc = x;
				int y_desloc = y;
				
                int *meio = geometry.ini + x_desloc + y_desloc * nx + z * ny * nx;

                if ( *meio )
                {
					double *f = ini_f + ( *meio - 1 ) * nvel;

					double rho = density( f );

					frho << rho << " ";					
                }
                else
                {
                    frho << 0.0 << " ";
                }
            }

            frho << endl;
        }

    }
    frho.close();
}
//====================================================================================================================//


//====================================================================================================================//
//                                                                                                                    //
//   G R A V A C A O   E M   V T K   B I N A R I O                                                                    //
//                                                                                                                    //
//   As funcoes rec_*_bin abaixo escrevem exatamente os mesmos campos das rec_* em ASCII, no formato                   //
//   VTK legacy binario.  Regras do formato (VTK File Formats, legacy):                                                //
//                                                                                                                    //
//      * o cabecalho continua sendo texto, uma linha por item, terminada em '\n';                                     //
//      * a terceira linha e BINARY em vez de ASCII;                                                                   //
//      * logo depois de "LOOKUP_TABLE default" (escalares) ou de "VECTORS ..." (vetores) vem o bloco                  //
//        de dados cru, sem separadores, seguido de um '\n';                                                           //
//      * os dados sao gravados em BIG-ENDIAN, independentemente da maquina;                                           //
//      * o tipo declarado (float / double / int) tem de bater com o numero de bytes gravados.                         //
//                                                                                                                    //
//   Precisao: por omissao grava-se float (4 bytes).  Passando double_prec = true grava-se double                      //
//   (8 bytes).  Para visualizar, float e mais do que suficiente -- o campo tem 7 digitos significativos,              //
//   enquanto o ASCII gravava 6.  Ou seja, o arquivo binario em float e ao mesmo tempo menor e um pouco                //
//   mais preciso do que o ASCII que substitui.                                                                        //
//                                                                                                                    //
//   Tamanho, por voxel:                                                                                               //
//                                                                                                                    //
//      escalar ASCII   ~ 9 a 13 bytes        escalar binario float   4 bytes                                          //
//      vetor  ASCII    ~ 27 a 39 bytes       vetor  binario float   12 bytes                                          //
//                                                                                                                    //
//   Os arquivos sao lidos pelo ParaView, VisIt, meshio e pelo vtk do Python sem nenhuma mudanca --                     //
//   o leitor reconhece o formato pela linha BINARY.                                                                    //
//                                                                                                                    //
//====================================================================================================================//


//------------------ A maquina e little-endian? ---------------------------------------------------------------------//

inline bool vtk_host_little_endian ()
{
	const unsigned int um = 1u;

	return *( const unsigned char * ) &um == 1u;
}

//------------------ Empilha um valor no buffer, em big-endian ------------------------------------------------------//

template < class T >
inline void vtk_push ( vector<char> &buf, T valor )
{
	const char *p = ( const char * ) &valor;

	if ( vtk_host_little_endian() )
	{
		for ( int i = ( int ) sizeof ( T ) - 1; i >= 0; i-- ) buf.push_back ( p[ i ] );
	}
	else
	{
		for ( size_t i = 0; i < sizeof ( T ); i++ ) buf.push_back ( p[ i ] );
	}
}

//====================================================================================================================//


//=================================== Returns the density of a site ==================================================//
//
//      Input: distribution function
//      Output: density
//
//====================================================================================================================//

#pragma acc routine seq
double density ( double *f )
{

    double rho = 0.0;

    for ( int i = 0 ; i < nvel; i++ ) rho = rho + f[i];

    return rho;
}

//====================================================================================================================//


//================================ Calculates the force term =========================================================//
//
//      Input: velocities, Forces, velocities, density, relaxation time, lattice
//      Output: S (source)
//
//====================================================================================================================//

#pragma acc routine seq
void source ( double Fx, double Fy, double Fz, double vx, double vy, double vz, double rho, double tau, double S[nvel], 
				LATTICE lattice )
{
	if ( dim == 2 )
	{
		double one_ov_cs_sqd = lattice.one_over_c_s2;
		
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		
		double F[dim];
		
		F[0] = Fx;
		F[1] = Fy;
		
		double factor = ( 1. - 1. / ( 2. * tau ) );
			
		double vF = dot_product ( F, v );
		
		double c_i[dim];

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			  
			double cF = dot_product ( c_i, F );
			
			double cv = dot_product ( c_i, v );

			S[i] = factor * lattice.w[i] *  one_ov_cs_sqd * ( cF + one_ov_cs_sqd * cF * cv - vF );
		}
	}
	else
	{
		double one_ov_cs_sqd = lattice.one_over_c_s2;
		
		double v[dim];

		v[0] = vx;
		v[1] = vy;
		v[2] = vz;
		
		double F[dim];
		
		F[0] = Fx;
		F[1] = Fy;
		F[2] = Fz;
		
		double factor = ( 1. - 1. / ( 2. * tau ) );
			
		double vF = dot_product ( F, v );
		
		double c_i[dim];

		for ( int i = 0; i < nvel; i++ )
		{
			c_i[0] = lattice.c_i[ i * dim + 0 ];
			c_i[1] = lattice.c_i[ i * dim + 1 ];
			c_i[2] = lattice.c_i[ i * dim + 2 ];           
			  
			double cF = dot_product ( c_i, F );
			
			double cv = dot_product ( c_i, v );

			S[i] = factor * lattice.w[i] *  one_ov_cs_sqd * ( cF + one_ov_cs_sqd * cF * cv - vF );
		}
	}
}

//====================================================================================================================//


//============================= Calculates the tensor Q_i = c_alpha c_beta - delta_alpha_beta ========================//
//
//      Input: lattice vectors
//      Output: tensor Q_i
//
//====================================================================================================================//

void calc_Q ( LATTICE lattice )
{
    for ( int i = 0 ; i < nvel; i++ )
    {
        for ( int alpha = 0; alpha < dim; alpha++ )
        {
            for ( int beta = 0; beta < dim; beta++ )
            {

                lattice.Q_i[ alpha + beta * dim + i * dim * dim ] 
                = lattice.c_i[ i*dim + alpha ] * lattice.c_i[ i*dim + beta ] - ( lattice.c_s2 ) * ( alpha == beta );
            }
        }
    }
}

//====================================================================================================================//


//===================== Etapa de recoloração (Latva-Koko) ============================================================//
//
//      Input: distribution functions, lattice vectors, gradient, magnitude of gradient,
//              mass fractions, density, recolloring factor
//      Output: distribution function
//
//====================================================================================================================//

#pragma acc routine seq
void recolloring ( double* f, double* f_R, double* f_B, double* grad, double mod_mM, double conc_R, double conc_B, 
					double rho, double beta, LATTICE lattice )
{
			
    double fator = beta * conc_R * conc_B;

	double c_i[3];
	
	double cos_phi = 0.;
	
    //----------------------------------------------------------------------------------------------------------------//
	
    f_R[0] = conc_R * f[0];
    f_B[0] = conc_B * f[0];

    for ( int i = 1; i < nvel; i++ )
    {

        c_i[0] = lattice.c_i[ i * dim + 0 ];
		c_i[1] = lattice.c_i[ i * dim + 1 ];
		c_i[2] = ( dim == 3 ) ? lattice.c_i[ i * dim + 2 ] : 0.0;           
			  
		double c_2 = dot_product( c_i, c_i );

        double prod_vm_ci = dot_product( c_i, grad );

        if ( mod_mM ) cos_phi = prod_vm_ci / ( mod_mM * sqrt ( c_2 ) );

        else cos_phi = 0.0;

        f_R[i] = conc_R * f[i] - fator * lattice.w[i] * rho * cos_phi;
        f_B[i] = conc_B * f[i] + fator * lattice.w[i] * rho * cos_phi;
    }
}

//====================================================================================================================//


