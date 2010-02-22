//******************************************************************************
//moralWL.cpp
//Este codigo implementa uma simulação MC Wang-Landau para uma 
//sociedade de perceptrons lineares com uma hamiltoneana de 
//influencia V(h,b)=-hb[1-(1-\delta)\Theta(hb)]
//29.06.09 
//R. Vicente
//******************************************************************************

#include "includes.h"
#include "moral.h" 

#define SIZE N*N	   
#define DIMENSION 2        //lattice dimension
#define TORUS 1            //periodic boundary conditions (=1)  
#define NBINS 400 	   //number of histogram bins 
#define CHECK 10000
#define DE 3.0
#define NSWEEPS 100000
#define MAXSWEEPS 100000000
#define NDE 20
#define WINDOW 10.0

void help();

main (int argc,char *argv[])
{
    int ia,ib,ic,id,it,inow,ineigh,icont;
    int in,ia2,ia3,irun,icurrent,ORTOGONALFLAG;
    int P,L,N,NRUNS,next,sweep,SHOWFLAG;
    double u,field1,field2,field0,q,aux1,aux2;
    double alfa,aux,Q1,Q2,QZ,RZQ,rho,R;
    double pm,D,wmax,mQ,wx,wy,h_sigma,h_mean;	
    double TOL,MINLOGF,E;
    double DELTA,RP;	
    double E_new,Ex,DeltaE,ER;
    double EW,meanhist,hvalue,wE,aratio;
    double logG_old,logG_new,lf;	
    size_t  i_old,i_new;	
    long seed;
    double lGvR,lGv,DlG;
    size_t iL,iR,i1,i2;
    int I_endpoint[NBINS];
    double lower,upper;
    size_t i0;		
	
    FILE * wlsrange; 
    FILE * dos;
    FILE * thermodynamics;
    FILE * canonical; 
    FILE * logfile;
    	     	
//***********************************
// Verificaccao de sintaxe
//*********************************** 
   if (argc<15){
	   help();
	   return(1);
   }
   else{
    	DELTA = atof(argv[1]);
   	P = atoi(argv[2]);
    	RP = atof(argv[3]); 
	L = atoi(argv[4]); 
	N = atoi(argv[5]);
	TOL = atof(argv[6]);
	MINLOGF = atof(argv[7]);
   }  	  
   wlsrange=fopen(argv[8],"w");	
   dos=fopen(argv[9],"w");  
   thermodynamics=fopen(argv[10],"w");
   canonical=fopen(argv[11],"w");
   logfile=fopen(argv[12],"w");	
   SHOWFLAG = atoi(argv[13]);
   ORTOGONALFLAG = atoi(argv[14]);

   if ((ORTOGONALFLAG==1) && (P>L)) P=L; //maximum number of orthogonal issues 

   if (SHOWFLAG==1){
  	printf("# parameters are DELTA=%1.2f P=%d ",DELTA,P);
        printf("D=%d L=%d RP=%1.2f TOL=%1.2f MINLOGF=%g \n",L,N,RP,TOL,MINLOGF); 
   }

  fprintf(logfile,"# parameters are DELTA=%1.2f P=%d D=%d",DELTA,P,D);
  fprintf(logfile,"L=%d RP=%1.2f TOL=%1.2f MINLOGF=%g\n",L,N,RP,TOL,MINLOGF); 	

//******************************************************************************
// Aloca matrizes                                              
//******************************************************************************
    gsl_matrix * sociedade = gsl_matrix_alloc(SIZE,L);
    gsl_matrix * issue = gsl_matrix_alloc(P,L);
    gsl_vector * current_issue = gsl_vector_alloc(L);
    gsl_vector * v0 = gsl_vector_alloc(L);
    gsl_vector * v1 = gsl_vector_alloc(L);
    gsl_vector * Z = gsl_vector_alloc(L);  
    gsl_vector * E_borda = gsl_vector_alloc(NBINS);  	
     
//******************************************************************************
// Inicialization                                                
//******************************************************************************
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
    for(ic=0;ic<DIMENSION;ic++) VECTOR(dim_vector)[ic]=N;

    gsl_histogram * HE = gsl_histogram_alloc (NBINS);
    gsl_histogram * logG = gsl_histogram_alloc (NBINS);
    gsl_histogram * LG = gsl_histogram_alloc (NBINS);

//******************************************************************************
// Social Graph                                             
//******************************************************************************
   //Regular lattice with periodic boundary conditions
   igraph_lattice(&graph,&dim_vector,0,0,0,TORUS);
    
    //Random rewiring with probability RP
    int NR=(int)(RP*(double)SIZE);
    igraph_rewire(&graph,NR,IGRAPH_REWIRING_SIMPLE);
   

//******************************************************************************
//Quenched issues set and Zeitgeist
//******************************************************************************	 
    gsl_vector_set_zero(Z);  
    gera_config(Z,issue,P,L,r,1.0);
     
    if (ORTOGONALFLAG==1) gsl_matrix_set_identity(issue);
   	 
    for (ib=0;ib<P;ib++)
    {
	gsl_matrix_get_row(current_issue,issue,ib);
	gsl_blas_ddot(current_issue,current_issue,&Q1);
	gsl_vector_scale(current_issue,1/sqrt(Q1));	
	gsl_vector_add(Z,current_issue);
     }
     gsl_blas_ddot(Z,Z,&QZ);
     gsl_vector_scale(Z,1/sqrt(QZ));			  
	
//******************************************************************************
//Determina a energia do ground state
//******************************************************************************
     double E0; 	
     gera_config(Z,sociedade,SIZE,L,r,0);  									
     E0 = hamiltoneana(sociedade,issue,SIZE,L,P,DELTA,graph,neighbors);
    
     double EMIN=E0;	
     double EMAX=-E0; 	
     double E_old=E0;
     		
     gsl_histogram_set_ranges_uniform (HE,EMIN,EMAX);
     gsl_histogram_set_ranges_uniform (logG,EMIN,EMAX); 
       
     if (SHOWFLAG==1) printf("# ground state: %3.0f\n",E0);
     fprintf(logfile,"# ground state: %3.0f\n",E0); 

//******************************************************************************		      	
// Determina intervalo de amostragem
//******************************************************************************    
     //printf("#finding the sampling interval\n");    

     lf=1;
     sweep=0;
     icont=0;	
     int iflag=0;
     int TMAX=NSWEEPS;	

     while(sweep<=TMAX){
	if (icont==10000) {
			//printf("%d sweeps\n",sweep);
			icont=0;
	}	
	for(it=0;it<SIZE;it++){
			 //escolhe sítio aleatoriamente
			  do{
              		 	inow=gsl_rng_uniform_int(r,SIZE);
			 }while((inow<0)||(inow>=SIZE)); 
	      		 gsl_matrix_get_row(v1,sociedade,inow);
	      		 igraph_neighbors(&graph,&neighbors,inow,IGRAPH_OUT); 
		
	      		 //gera vetor aleatorio v1
	      		 gsl_vector_memcpy(v0,v1);	
			 gera_vetor(v1,L,r);

			 //calcula variaccao de energia na trasicao v0->v1
			 //no sitio inow
			 DeltaE=variacaoE(v0,v1,inow,sociedade,
					  issue,N,L,P,DELTA,graph,neighbors);	      		
			 E_new=E_old+DeltaE;
			
			 //WL: aceita dentro da janela [EMIN,EMAX]
			 if ((E_new>EMIN) && (E_new<EMAX))
	      		 {
		   		gsl_histogram_find(logG,E_old,&i_old);
		   		logG_old=gsl_histogram_get(logG,i_old);
		   		gsl_histogram_find(logG,E_new,&i_new);
		   		logG_new=gsl_histogram_get(logG,i_new); 	
		  		wE = GSL_MIN(exp(logG_old-logG_new),1);		
		   		if (gsl_rng_uniform(r)<wE){
					E_old=E_new;
					gsl_matrix_set_row(sociedade,inow,v1);
		   		}
	      		 }
			 //WL: atualiza histogramas
			 gsl_histogram_increment(HE,E_old); 
	     		 gsl_histogram_accumulate(logG,E_old,lf); 
		 }	
		sweep++;
		icont++;	
     }		 
     	
     gsl_histogram_fprintf(wlsrange,HE,"%g","%g");
    
     double maxH=gsl_histogram_max_val(HE);
     	
     //printf("ok\n");
     Ex=0;
     hvalue=maxH;		
     while((hvalue>TOL*maxH)&&(Ex>EMIN)){
	gsl_histogram_find(HE,Ex,&i0);
	hvalue=gsl_histogram_get(HE,i0);
	Ex-=1;
	if(Ex<=EMAX)TMAX+=10000;
     }		
     EMIN=Ex;
	
     Ex=0;	
     hvalue=maxH;	
     while((hvalue>TOL*maxH)&&(Ex<EMAX)) {
	gsl_histogram_find(HE,Ex,&i0);
	hvalue=gsl_histogram_get(HE,i0);
	Ex+=1;
	if(Ex>=EMAX)TMAX+=10000;
     }		
     EMAX=Ex;	   
     EMAX=GSL_MIN(10.0,Ex);
     if (SHOWFLAG==1) 
       printf("# the sampling interval is [%3.0f,%3.0f] found in %d sweeps \n\n"
	                                                      ,EMIN,EMAX,sweep);

     fprintf(logfile,
	"# the sampling interval is [%3.0f,%3.0f] found in %d sweeps \n\n"
	                                                      ,EMIN,EMAX,sweep);
     
     gsl_histogram_set_ranges_uniform (HE,EMIN-1,EMAX+1);
     gsl_histogram_set_ranges_uniform (logG,EMIN-1,EMAX+1); 
     gsl_histogram_set_ranges_uniform (LG,EMIN-1,EMAX+1); 		
     
//******************************************************************************		      	
// WLS
//******************************************************************************		
     int iE,itera=0;
     double endpoints[NBINS];		
     double w = WINDOW; //(EMAX-EMIN)/10.0;	
     //printf("W=%f\n",w);	
     lf=1;

//RESOLUTION ---->                                        <------RESOLUTION*****            
     do{
	int iverify=0,iborda=0,flat=0; 
	sweep=0;
	Ex=EMAX; 
	EW=EMAX;
	E_old=EMAX+1;
	iE=0;
	endpoints[iE]=EMAX;
	iE++;
	gsl_histogram_reset(LG);		

//WINDOWS -->                                                  <--WINDOWS*******           
	while((Ex>EMIN)&&(sweep<MAXSWEEPS)){	 
	   //configuracao inicial
	   gera_config(Z,sociedade,SIZE,L,r,0);
	   E_old = hamiltoneana(sociedade,issue,SIZE,L,P,DELTA,graph,neighbors);
	   while( (E_old<EMIN+1)||(E_old>Ex) ){
		//printf("%d %3.1f\n",E_old);
		do{
              		inow=gsl_rng_uniform_int(r,SIZE);
		  }while((inow<0)||(inow>=SIZE)); 	
	      	gsl_matrix_get_row(v1,sociedade,inow);
	   	gera_vetor(v1,L,r); 	
		gsl_matrix_set_row(sociedade,inow,v1);					
                E_old = hamiltoneana(sociedade,issue,SIZE,L,
				     P,DELTA,graph,neighbors);
	   }
	   
	   if (SHOWFLAG==1){
		printf("# sampling [%f,%f]\n",EMIN,Ex);
	   	printf("# walking from E=%3.0f\n",E_old);
	   }
	   
	   fprintf(logfile,"# sampling [%f,%f]\n",EMIN,Ex);
	   fprintf(logfile,"# walking from E=%3.0f\n",E_old);

	   do{	//FLAT WINDOW------>                     <------FLAT WINDOW*****
//MC sweep ---->                                         <------MC sweep********	
		
	    	for(it=0;it<SIZE;it++){
			 //escolhe sítio aleatoriamente
			 do{
              		 	inow=gsl_rng_uniform_int(r,SIZE);
			 }while((inow<0)||(inow>=SIZE)); 
	      		 gsl_matrix_get_row(v1,sociedade,inow);
	      		 igraph_neighbors(&graph,&neighbors,inow,IGRAPH_OUT); 
		
	      		 //gera vetor aleatorio v1
	      		 gsl_vector_memcpy(v0,v1);	
			 gera_vetor(v1,L,r);

			 //calcula variaccao de energia na trasicao 
			 //v0->v1 no sitio inow
			 DeltaE=variacaoE(v0,v1,inow,sociedade,issue,
					  N,L,P,DELTA,graph,neighbors);	      		
			 E_new=E_old+DeltaE;
			
			 //WL: aceita dentro da janela [EMIN,Ex]
			 if ((E_new>EMIN) && (E_new<Ex))
	      		 {
		   		gsl_histogram_find(logG,E_old,&i_old);
		   		logG_old=gsl_histogram_get(logG,i_old);
		    		gsl_histogram_find(logG,E_new,&i_new);
		   		logG_new=gsl_histogram_get(logG,i_new); 	
		  		wE = GSL_MIN(exp(logG_old-logG_new),1);		
		   		if (gsl_rng_uniform(r)<wE){
					E_old=E_new;
					gsl_matrix_set_row(sociedade,inow,v1);
		   		}
	      		 }
			 //WL: atualiza histogramas
			 gsl_histogram_increment(HE,E_old); 
	     		 gsl_histogram_accumulate(logG,E_old,lf); 
			 itera++;
		 }
//MC sweep ---->                                           <--------MC sweep**** 
		sweep++; iverify++;   

		if( (EMAX-EMIN)<NDE*DE ) {
			EW=EMIN;
		}else{	    
	   		EW=GSL_MAX(Ex-w,EMIN);
		}	

	   	if (iverify==CHECK){//Verify flatness 
			if (SHOWFLAG==1)  
			    printf(" #verificando flatness em [%f,%f]\n",EW,Ex);
	
			fprintf(logfile," #verificando flatness em [%f,%f]\n"
				                                        ,EW,Ex);
	   		iverify=0;
			flat=flatness(HE,EW,Ex,TOL,itera,meanhist,hvalue);
			if (SHOWFLAG==1) 
			    printf("#minH= %8.0f\t k<H>=%8.0f\t %d sweeps\t ",
				                hvalue,TOL*meanhist,sweep,flat);

			fprintf(logfile,
				"#minH= %8.0f\t k<H>=%8.0f\t %d sweeps\t ",
			                        hvalue,TOL*meanhist,sweep,flat);
	   	}

	    }while(flat==0);//                          <------FLAT WINDOW******	 	  
            flat=0;
  
	    //Find ER
	    if( (EMAX-EMIN)<NDE*DE ) {
		Ex=EMIN;
		endpoints[iE]=EMIN;
	    } 
	    else {		
	    	if (EW>EMIN){
			 ER=flatwindow(HE,EW,TOL,meanhist);
			 if (SHOWFLAG==1)  
			      printf("# extendendo flatness a [%f,%f]\n",ER,Ex);

			 fprintf(logfile,
				"# extendendo flatness a [%f,%f]\n",ER,Ex);

			 if((ER-EMIN)<1){
				ER=EMIN;
				Ex=EMIN;
				endpoints[iE]=EMIN;
		 	}else{
		 		endpoints[iE]=GSL_MIN(ER+DE,EMAX);
				Ex=GSL_MIN(ER+2*DE,EMAX);
			}
	     	}
	     	else{
			 endpoints[iE]=EMIN;
			 Ex=EMIN;	
			 ER=EMIN;	   			
	    	} 	 	
	    }	    
	   
	    if (SHOWFLAG==1) 
		 printf("# window %d [%3.0f,%3.0f] is flat after %d sweeps \n",
					iE,endpoints[iE],endpoints[iE-1],sweep);

	  fprintf(logfile,"# window %d [%3.0f,%3.0f] is flat after %d sweeps\n",
			iE,endpoints[iE],endpoints[iE-1],sweep);	
	   		
	   	     
	    //registra histograma
	    if (iE==1){
		gsl_histogram_find(logG,endpoints[iE],&i1);
		gsl_histogram_find(logG,endpoints[iE-1],&i2);
		for(i0=i1;i0<=i2;i0++){
			lGv=gsl_histogram_get(logG,i0);
			gsl_histogram_get_range(logG,i0,&lower,&upper);
			E=0.5*(upper+lower);
			gsl_histogram_accumulate(LG,E,lGv);
		}				
	    }else{
		gsl_histogram_find(logG,endpoints[iE],&i1);
		gsl_histogram_find(logG,endpoints[iE-1],&i2);
		lGv=gsl_histogram_get(logG,i2);
		lGvR=gsl_histogram_get(LG,i2);
		DlG=lGvR-lGv;
	     //printf("i1=%d i2=%d lGv=%f lGvR=%f DlG=%f\n",i1,i2,lGv,lGvR,DlG);
		for(i0=i1;i0<i2;i0++){
			lGv=gsl_histogram_get(logG,i0);
			lGv=lGv+DlG;
			gsl_histogram_get_range(logG,i0,&lower,&upper);
			E=(upper+lower)*0.5;
			//printf("i0=%d E=%f lGv=%f\n",i0,E,lGv);
			gsl_histogram_accumulate(LG,E,lGv);
		}		
	    }	
				 
	    //printf("#########################################\n");		
	    //gsl_histogram_fprintf(stdout,LG,"%g","%g");		
	    //printf("#########################################\n");			

	    iE++;
	    if(Ex>EMIN) {
		if (SHOWFLAG==1) 
		     printf("# random walk is now restricted to [%3.0f,%3.0f]\n"
			    					      ,EMIN,Ex);
	    fprintf(logfile,"# random walk is now restricted to [%3.0f,%3.0f]\n"
			                                              ,EMIN,Ex);
	    }
	    gsl_histogram_reset(HE);
	   		     	 		
      }  
//WINDOWS -->

     if(sweep<MAXSWEEPS){	
     	if (SHOWFLAG==1) 
		  printf("# log(f)=%f converged within %d sweeps\n\n",lf,sweep);	
	fprintf(logfile,"# log(f)=%f converged within %d sweeps\n\n",lf,sweep);		 	
     	lf=lf/2.0;	 
     	gsl_histogram_reset(HE);
     	gsl_histogram_memcpy(logG,LG);
     }else {
	if (SHOWFLAG==1) 
		printf("# FAILED: no convergence has been attained.");
	fprintf(logfile,
           "# FAILED: no convergence has been attained. Simulation ABANDONED.");
	return(1);
     }	 			
	     
    }while(lf>MINLOGF); 
//RESOLUTION -->                                            <-----RESOLUTION****

     //************************************************************************* 	
     //Density of states	     	
     //*************************************************************************
     double minlogG=gsl_histogram_min_val(logG);	
     gsl_histogram_shift(logG,-minlogG);	
     gsl_histogram_fprintf(dos,logG,"%g","%g");	

     //************************************************************************* 
     //Thermodynamics    	
     //************************************************************************* 
     double beta,A,wT,Zmin_beta;
     double lGvalue,maxA,betaC,CTMAX=0;	
     double Z_beta,U,U2,CT,F,S;
     
     for (beta=0.01;beta<=15;beta+=0.01)
     {
	//********************************************************************** 	
	//Energy, free-energy, entropy, specific heat and Tc
 	//********************************************************************** 
	maxA=0;
	for (ia2=0; ia2<NBINS;ia2++)
	{
		lGvalue=gsl_histogram_get(logG,ia2);
		gsl_histogram_get_range(logG,ia2,&lower,&upper);
		E=(lower+upper)/2;
		A=lGvalue-beta*E;
		if (A>maxA) maxA=A;
	}       

        gsl_histogram_find(logG,EMIN,&i0); 

	Z_beta=0;U=0;U2=0;
	for (ia2=0; ia2<NBINS;ia2++)
	{
			lGvalue=gsl_histogram_get(logG,ia2);
			gsl_histogram_get_range(logG,ia2,&lower,&upper);
			E=(lower+upper)/2;
			A=lGvalue-beta*E-maxA;
			Z_beta+=exp(A);  
			U+=E*exp(A);
			U2+=E*E*exp(A);  
			if(ia2==i0) Zmin_beta=exp(A);
	}	
	wT=Zmin_beta/Z_beta;

	F=-log(Z_beta)/beta - maxA/beta; 
	U=U/Z_beta;
	S= (U-F)*beta;
	U2=U2/Z_beta;
	CT=(U2-U*U)*beta*beta;	
			
	if(CT>CTMAX){
		CTMAX=CT;
		betaC=beta;
	}

	fprintf(thermodynamics,"%f %f %f %f %f %f %f \n"
		,beta,1/beta,F/(double)(SIZE),S/(double)(SIZE),
		 U/(double)(SIZE),CT/(double)(SIZE),wT);
     }
     
     if (SHOWFLAG==1) printf("# BETAc: %f  Tc:%f \n",betaC,1/betaC); 
     fprintf(logfile,"# BETAc: %f  Tc:%f \n",betaC,1/betaC); 	
		

    //**************************************************************************	
    //canonical distribuition at Tc
    //**************************************************************************	
     beta=betaC; 
     double distr_canonica;	
     
      maxA=0;
      for (ia2=0; ia2<NBINS;ia2++)
      {
		lGvalue=gsl_histogram_get(logG,ia2);
		gsl_histogram_get_range(logG,ia2,&lower,&upper);
		E=(lower+upper)/2;
		A=lGvalue-beta*E;
		if (A>maxA) maxA=A;
      }       

      for (ia2=0; ia2<NBINS;ia2++)
      {
	  lGvalue=gsl_histogram_get(logG,ia2);
	  gsl_histogram_get_range(logG,ia2,&lower,&upper);
	  E=(lower+upper)/2;
	  A=lGvalue-beta*E-maxA;
	  distr_canonica=exp(A);
	  fprintf(canonical,"%f %f %f\n",
		  E/(double)(SIZE),distr_canonica,A);  
      }			

     //*************************************************************************
     // Finalizacoes                                                     
     //*************************************************************************
     igraph_destroy(&graph);
     igraph_vector_destroy(&neighbors);
     igraph_vector_destroy(&result);  
     gsl_matrix_free(issue);
     gsl_vector_free(current_issue);
     gsl_vector_free(v1);
     gsl_vector_free(v0);
     gsl_matrix_free(sociedade);     	   	
     gsl_rng_free(r);
	
     fclose(wlsrange);
     fclose(dos);
     fclose(thermodynamics);
     fclose(canonical);   
     fclose(logfile);   	
   
     return(0);
}


//******************************************************************************
// help
//******************************************************************************
void help()
{
  printf( "\n");
  printf("USAGE: moralWLwindows  DELTA ISSUES RP DIM SIZE FLAT MINLOGF ");
  printf("FILE_WLSrange FILE_dos FILE_thermo FILE_canonical FILE_log SHOWFLAG");
  printf(" ORTOGONALFLAG\n\n");
  printf("-----------------------------------------------------\n");		
  printf("DELTA: confirmation bias\n");
  printf("ISSUES: number of simultaneously debated issues\n");
  printf("RP: rewiring probability\n");
  printf("DIM: moral matrix dimension\n");
  printf("SIZE: lattice length (N_AGENTS=SIZE*SIZE) \n");
  printf("FLAT: flat histogram criterium \n");	
  printf("MINLOGF: minimum log(f)\n");	
  printf("SHOWFLAG=1 log on the screen, ");
  printf("ORTHOGONALFLAG=1 P<D orthogonal issues\n");	
  printf("-----------------------------------------------------\n");	
  printf("OUTPUT FILES: \n");
  printf("thermo.dat:  beta T f s u c R \n");
  printf("WLSrange.dat:  E E+DE H(E) for at least 100000 sweeps\n");
  printf("dos.dat:  E E+DE g(E) \n");
  printf("canonical.dat:  E  g(E)exp(-E/T) log(g)-E/T\n");
  printf("logfile.dat \n");	
  printf("-----------------------------------------------------\n\n");			
  exit(1);
}

