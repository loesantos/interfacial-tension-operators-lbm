


//====================================================================================================================//
//
//   A GEOMETRIA NO DISPOSITIVO
//
//   As fronteiras deste arquivo sao as unicas rotinas da biblioteca que saem do indice de fluido
//   e voltam ao indice de voxel: elas leem  geometry.ini  para achar, em cada coluna, o sitio e o
//   seu vizinho.  Os kernels do laco principal nunca precisam disso, entao um programa tipico
//   NAO mapeia  geometry.ini  -- e ai uma clausula  present  sobre ele falha em execucao, com
//
//       FATAL ERROR: data in PRESENT clause was not found on device
//
//   A struct GEOMETRY inteira nao pode ir para o dispositivo, porque contem um std::string; so o
//   ponteiro do mapa de indices pode.  Como aqui  geometry.ini  e apenas LIDO, a solucao segura e
//   declara-lo com  copyin :  do OpenACC 2.5 em diante  copyin  e "presente-ou-copia", isto e,
//   reaproveita a copia do dispositivo quando ela existe e faz uma temporaria quando nao existe.
//
//   As populacoes continuam em  present , e isso e proposital: elas sao ESCRITAS, e uma copia
//   temporaria receberia as escritas e as jogaria fora ao fim da regiao.
//
//   Copiar o mapa a cada chamada funciona mas e desperdicio.  O aviso abaixo sai uma unica vez e
//   diz como evitar: mapear a geometria uma vez, no programa.
//
//====================================================================================================================//

#ifdef _OPENACC

static void confere_geometria_no_dispositivo ( const int *geometry_ini, int geometry_size )
{
	static bool avisado = false;

	if ( avisado ) return;

	if ( acc_is_present ( ( void * ) geometry_ini,
	                      ( size_t ) geometry_size * sizeof ( int ) ) ) return;

	avisado = true;

	cerr << "\n   AVISO: geometry.ini nao esta no dispositivo.  As condicoes de fronteira vao"
	        "\n          copia-lo a cada chamada; o resultado esta certo, mas custa uma"
	        "\n          transferencia por chamada.  Para evitar, acrescente ao programa, junto"
	        "\n          dos outros  enter data :"
	        "\n"
	        "\n              int *geo_ini = geometry.ini;"
	        "\n              #pragma acc enter data copyin( geo_ini[ 0 : ( long long ) nx * ny * nz ] )"
	        "\n"
	        "\n          e o  exit data delete  correspondente ao fim.\n" << endl;
}

#define CONFERE_GEOMETRIA( p, n )   confere_geometria_no_dispositivo ( ( p ), ( n ) )

#else

#define CONFERE_GEOMETRIA( p, n )   ( ( void ) 0 )

#endif

//====================================================================================================================//
