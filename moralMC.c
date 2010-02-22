//******************************************************************************
//moralMC.cpp
//Este codigo implementa uma simulação MC para uma 
//sociedade de perceptrons lineares com uma hamiltoneana de 
//influencia V(h,b)=-hb[1-(1-\delta)\Theta(hb)]
//Linhas de transiccao ordem-desordem sao calculadas no plano
// \delta vs P  com parametro de ordem dado por (1/N) \sum_k |J_k.Z|/(|J_k||Z|)
//Esta versao constroi curvas de histerese
//19.05.09 
//R. Vicente
//******************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <time.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_erf.h>
#include <igraph/igraph.h>

#define SIZE N*N	   
#define DIMENSION 2        //lattice dimension
#define TORUS 0            //periodic boundary conditions (=1)  
#define NBINS 200 	   //number of histogram bins 
#define EPS 5
#define EPS2 0.00
#define STEPTEMP 0.1
#define NTEMP 41  

void help();

main (int argc,char *argv[])
{
	int ia,ib,ic,id,it,inow,ineigh,icont,ix,iy,itemp;
	int in,ia2,ia3,irun,icurrent,iflag;
	int P,L,N,TMAX,NRUNS;
	double u,field1,field2,S,q,aux1,aux2;
	double erro_mean, erro_sigma;
	double alfa,sigma,aux,Q1,Q2,QZ,RZQ,rho,R;
	double pm,D,wmax,mQ,wx,wy,h_sigma,h_mean;	
	double h_overlap_mean,h_overlap_sigma;
	double DELTA,RP,BETA,TEMP,TEMPMIN,TEMPMAX;	
	double H,Delta,E,DeltaE,E_new,E_old,TEMPSTEP;
	long seed;

//***********************************
// Verificaccao de sintaxe
//***********************************
	if (argc<8){
		help();
		return(1);
	}
	else{
		DELTA = atof(argv[1]);
		P = atoi(argv[2]);
		RP = atof(argv[3]); 
		L = atof(argv[4]); 
		N = atof(argv[5]);
		NRUNS =atoi(argv[6]); 
		TMAX = atoi(argv[7]);
	}  	  

//*********************
// Aloca matrizes                                              
//*********************
	gsl_matrix * sociedade = gsl_matrix_alloc(SIZE,L);
	gsl_matrix * issue = gsl_matrix_alloc(P,L);
	gsl_vector * current_issue = gsl_vector_alloc(L);
	gsl_vector * perceptron1 = gsl_vector_alloc(L);
	gsl_vector * perceptron2 = gsl_vector_alloc(L);
	gsl_vector * perceptron_trial = gsl_vector_alloc(L);
	gsl_vector * rnd_vector = gsl_vector_alloc(L);
	gsl_vector * Z = gsl_vector_alloc(L);  
	gsl_vector * opinions = gsl_vector_alloc(SIZE);
	gsl_vector * Tschedule = gsl_vector_alloc(NTEMP);
	gsl_vector * Dmean = gsl_vector_alloc(NTEMP);
	gsl_vector * Dsigma = gsl_vector_alloc(NTEMP); 
	gsl_vector * Dmean2 = gsl_vector_alloc(NTEMP);
	gsl_vector * Dsigma2 = gsl_vector_alloc(NTEMP); 

//*********************
// Inicialization                                               
//*********************
	const gsl_rng_type * T;
	gsl_rng * r; 

	gsl_rng_env_setup();  
	T = gsl_rng_default;
	r=gsl_rng_alloc (T);

	seed = time (NULL) * getpid();
	//seed = 13188839657852;
	gsl_rng_set(r,seed);

	igraph_t graph;
	igraph_vector_t neighbors;
	igraph_vector_t result;
	igraph_vector_t dim_vector;
	igraph_real_t res;
	igraph_bool_t C;

	igraph_vector_init(&neighbors,1000);
	igraph_vector_init(&result,0);
	igraph_vector_init(&dim_vector,DIMENSION);
	igraph_vector_fill(&dim_vector,N); 

	gsl_histogram * h_overlap = gsl_histogram_alloc (NBINS);
	gsl_histogram_set_ranges_uniform (h_overlap,-1.01,1.01);

	gsl_histogram * h_dim = gsl_histogram_alloc (NBINS);
	gsl_histogram_set_ranges_uniform (h_dim,0,L+1);

//*******************
// Social Graph                                             
//*******************
//Regular lattice with periodic boundary conditions
	igraph_lattice(&graph,&dim_vector,0,0,0,TORUS);

/*  for (inow=0;inow<SIZE;inow++){
	igraph_neighbors(&graph,&neighbors,inow,IGRAPH_OUT);
	printf("%d ",inow);
	for(ic=0;ic<igraph_vector_size(&neighbors);ic++)
	{	
		ineigh=(int)VECTOR(neighbors)[ic];
		printf("%d ",ineigh);
	}
	printf("\n");	
	}*/

	//Random rewiring with probability RP
		igraph_rewire_edges(&graph,RP);

	//*************************************************
	//Temperature schedule
	//*************************************************
	TEMP=STEPTEMP;
	gsl_vector_set(Tschedule,0,TEMP);	
	gsl_vector_set(Dmean,0,0);
	gsl_vector_set(Dsigma,0,0);

	TEMP+=STEPTEMP;
	for(ia=1;ia<NTEMP;ia++)
	{
		gsl_vector_set(Tschedule,ia,TEMP);  	
		TEMP+=STEPTEMP;	
		gsl_vector_set(Dmean,0,0);
		gsl_vector_set(Dsigma,0,0);    	
	}	

	for (in=0;in<NRUNS;in++) 
	{	
		//**************************************************
		//Quenched issues set and Zeitgeist
		//**************************************************	
		gsl_vector_set_zero(Z);             
		for (ib=0;ib<P;ib++)
		{
			for(ic=0;ic<L;ic++)
			{
				u= gsl_ran_ugaussian(r);
				gsl_matrix_set(issue,ib,ic,u);
			}   
			gsl_matrix_get_row(perceptron1,issue,ib);
			gsl_vector_add(Z,perceptron1);
		}
		gsl_blas_ddot(Z,Z,&QZ);	

	//**********************************************
		// Initial ordered configuration                                            
		//**********************************************
		for(ia=0;ia<SIZE;ia++)
		{	
			gsl_vector_memcpy(perceptron1,Z);
			gsl_blas_ddot(perceptron1,perceptron1,&Q1);
			aux=1/sqrt(Q1);
			gsl_vector_scale(perceptron1,aux);

			gsl_vector_memcpy(perceptron_trial,perceptron1);	
			for(ic=0;ic<L;ic++)
			{
				u= gsl_ran_ugaussian(r);
			//printf("u=%f\n",u);
				gsl_vector_set(rnd_vector,ic,u);
				gsl_blas_ddot(rnd_vector,rnd_vector,&Q1);
				aux=EPS2/sqrt(Q1);
				gsl_vector_scale(rnd_vector,aux);
			}   

			gsl_vector_add(perceptron_trial,rnd_vector);
			gsl_blas_ddot(perceptron_trial,perceptron_trial,&Q1);
			aux=1/sqrt(Q1);
			gsl_vector_scale(perceptron_trial,aux);	
			gsl_matrix_set_row(sociedade,ia,perceptron_trial);						
		}	 		

		//Temperature annealing
		for (itemp=0;itemp<NTEMP;itemp++)	
		{ 	 
			TEMP=gsl_vector_get(Tschedule,itemp);
			BETA=1/TEMP;

		//**************************************************************************	 		      	
			//MC steps
		//**************************************************************************
			E=0;	
			it=0;
			iflag=0;
			gsl_histogram_reset(h_dim);				
			do
			{								
				for(inow=0;inow<SIZE;inow++)
				{
			//*********************************
			//current perceptron and neighbors
			//*********************************
					gsl_matrix_get_row(perceptron1,sociedade,inow);
					igraph_neighbors(&graph,&neighbors,inow,IGRAPH_OUT);

			//********************************
			//current energy
			//********************************
					E_old=0;
					for (icurrent=0;icurrent<P;icurrent++)
					{
				//current issue
						gsl_matrix_get_row(current_issue,issue,icurrent);

				//influenced field
						gsl_matrix_get_row(perceptron1,sociedade,inow);
						gsl_blas_ddot(perceptron1,current_issue,&field1);	

						H=0;	
						for(ic=0; ic<igraph_vector_size(&neighbors); ic++)
						{

							ineigh=(int)VECTOR(neighbors)[ic];
					//influencer opinion			
							gsl_matrix_get_row(perceptron2,sociedade,ineigh);
							gsl_blas_ddot(perceptron2,current_issue,&field2);

							if(field2*field1>0) alfa = DELTA/(double)P;
							else alfa =1/(double)P;
							H+=alfa*field2;	
						} 

						E_old+=-field1*H;
					}									

			//********************************		
			//trial
			//********************************
					gsl_vector_memcpy(perceptron_trial,perceptron1);	
					for(ic=0;ic<L;ic++)
					{
						u= gsl_ran_ugaussian(r);
				//printf("u=%f\n",u);
						gsl_vector_set(rnd_vector,ic,u);
						gsl_blas_ddot(rnd_vector,rnd_vector,&Q1);
						aux=EPS/sqrt(Q1);
						gsl_vector_scale(rnd_vector,aux);
					}   

					gsl_vector_add(perceptron_trial,rnd_vector);
					gsl_blas_ddot(perceptron_trial,perceptron_trial,&Q1);
					aux=1/sqrt(Q1);
					gsl_vector_scale(perceptron_trial,aux);

			//********************************
			//new energy
			//********************************
					E_new=0;
					for (icurrent=0;icurrent<P;icurrent++)
					{
				//current issue
						gsl_matrix_get_row(current_issue,issue,icurrent);

				//influenced field
						gsl_blas_ddot(perceptron_trial,current_issue,&field1);	

						H=0;	
						for(ic=0; ic<igraph_vector_size(&neighbors); ic++)
						{				
							ineigh=(int)VECTOR(neighbors)[ic];
					//influencer opinion			
							gsl_matrix_get_row(perceptron2,sociedade,ineigh);
							gsl_blas_ddot(perceptron2,current_issue,&field2);

							if(field2*field1>0) alfa = DELTA/(double)P;
							else alfa =1/(double)P;
							H+=alfa*field2;	
						} 
						E_new+=-field1*H;
					}

			//**********************************
			//accept new configuration
			//**********************************

					DeltaE=E_new-E_old;			
					u=gsl_rng_uniform(r);		

					if( DeltaE<0) {
						gsl_matrix_set_row(sociedade,inow,perceptron_trial);
						E+=DeltaE;		
					}
					else if ( u<exp(-BETA*DeltaE) ) { 
				//printf("u=%f %f p=%f\n",u,BETA,exp(-BETA*DeltaE));      	
						gsl_matrix_set_row(sociedade,inow,perceptron_trial);
						E+=DeltaE;			
					}	
				//printf("%f\n",E);			

				}//inow

				it++;
				if(it>TMAX){

				//*******************************	
					//Gather statistics
					//*******************************	 
					for(ia=0;ia<SIZE;ia++)
					{
						gsl_matrix_get_row(perceptron1,sociedade,ia);
						gsl_blas_ddot(perceptron1,perceptron1,&Q1);

						gsl_blas_ddot(perceptron1,Z,&RZQ);

						D = fabs(RZQ)/sqrt(Q1*QZ);
						gsl_histogram_increment(h_dim,D);
					}	
				} 

				if(it>TMAX+500) iflag=1;


			}while(iflag==0);   

			h_sigma=gsl_histogram_sigma(h_dim);
			h_mean=gsl_histogram_mean(h_dim);  	

			aux1=gsl_vector_get(Dmean,itemp);
			aux2=gsl_vector_get(Dsigma,itemp);
			aux1+=h_mean/(double)NRUNS;
			aux2+=h_sigma/(double)NRUNS;
			gsl_vector_set(Dmean,itemp,aux1);
			gsl_vector_set(Dsigma,itemp,aux2);

			aux1=gsl_vector_get(Dmean2,itemp);
			aux2=gsl_vector_get(Dsigma2,itemp);
			aux1+=h_mean*h_mean/(double)NRUNS;
			aux2+=h_sigma*h_sigma/(double)NRUNS;
			gsl_vector_set(Dmean2,itemp,aux1);
			gsl_vector_set(Dsigma2,itemp,aux2);

		//printf("%d %f\n",in,TEMP);
		}//TEMP 	 
	}//NRUNS

	//*********************
	//  Imprime resultados                                                  
	//*********************	

	for (itemp=0;itemp<NTEMP;itemp++)	
	{ 	 
		TEMP=gsl_vector_get(Tschedule,itemp); 	
		BETA=1/TEMP;
		h_mean=gsl_vector_get(Dmean,itemp);
		h_sigma=gsl_vector_get(Dsigma,itemp);
		erro_mean=sqrt(fabs(gsl_vector_get(Dmean2,itemp)-h_mean*h_mean));
		erro_sigma=sqrt(fabs(gsl_vector_get(Dsigma2,itemp)-h_sigma*h_sigma)); 	   
		printf("%f %f %f %f %f %f\n",TEMP,BETA,h_mean,erro_mean,h_sigma,erro_sigma);
	}		

	//****************
	// Finalizacoes                                                  
	//****************
	igraph_destroy(&graph);
	igraph_vector_destroy(&neighbors);
	igraph_vector_destroy(&result);
	gsl_matrix_free(issue);
	gsl_vector_free(current_issue);
	gsl_vector_free(opinions);
	gsl_vector_free(perceptron1);
	gsl_vector_free(perceptron2);
	gsl_vector_free(perceptron_trial);
	gsl_matrix_free(sociedade);  
	gsl_vector_free(Tschedule);  
	gsl_vector_free(Dmean);
	gsl_vector_free(Dsigma); 
	gsl_vector_free(Dmean2);
	gsl_vector_free(Dsigma2);    	   	
	gsl_rng_free (r);

	return(0);
}

void help()
{
	printf( "\n");
	printf( "USAGE: moralMC  DELTA  ISSUES RP DIM SIZE NRUNS TMAX \n\n");	
	printf( "DELTA: confirmation bias\n");
	printf( "ISSUES: number of simultaneously debated issues\n");
	printf( "RP: rewiring probability\n");
	printf( "DIM: moral matrix dimension\n");
	printf( "SIZE: lattice length (N_AGENTS=SIZE*SIZE) \n");
	printf( "NRUNS\n");
	printf( "TMAX:  number of lattice steps \n\n");
	printf( "OUTPUT  TEMP,BETA,h_mean,erro_mean,h_sigma,erro_sigma \n\n\n");
	exit(1);
}




