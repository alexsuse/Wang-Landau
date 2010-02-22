//******************************************************************************
// Gera SIZE vetores aleatórios com dimensao L como 
//((1-eps)v0 + eps v)/|(1-eps)v0 + eps v| -> v 
// eps=0 -> v=v0  eps=1 -> v aleatorio
//******************************************************************************
void gera_config (gsl_vector * v0, gsl_matrix * M, 
		  int N, int L, gsl_rng * r, double eps)    	
{
	int ia,ib;
	double u,Q;
	gsl_vector * v = gsl_vector_alloc(L);
	
	gsl_vector_scale(v0,1-eps);
	for(ia=0;ia<N;ia++){

		for(ib=0;ib<L;ib++){
  			u=gsl_ran_ugaussian(r);
        		gsl_vector_set(v,ib,u);
       		}
		gsl_blas_ddot(v,v,&Q);
		gsl_vector_scale(v,eps/sqrt(Q));
		gsl_vector_add(v,v0);

		gsl_blas_ddot(v,v,&Q);
		gsl_vector_scale(v,1/sqrt(Q));
		gsl_matrix_set_row(M,ia,v); 
    	} 
	gsl_vector_free(v);
}

//******************************************************************************
// Hamiltoneana
//******************************************************************************
double hamiltoneana( gsl_matrix * sociedade, gsl_matrix * issue, 
		     int N, int L, int P, double DELTA, 
		     igraph_t graph,  igraph_vector_t neighbors)
{
    double E=0,H;
    double field1,field2,alfa; 	
    int inow,icurrent,ineigh,ic;
    gsl_vector * v1 = gsl_vector_alloc(L);	
    gsl_vector * v2 = gsl_vector_alloc(L);	
    gsl_vector * v_issue = gsl_vector_alloc(L);
	
    for(inow=0;inow<N;inow++)
    {
	gsl_matrix_get_row(v1,sociedade,inow);
	igraph_neighbors(&graph,&neighbors,inow,IGRAPH_OUT);

	for (icurrent=0;icurrent<P;icurrent++)
	{
		gsl_matrix_get_row(v_issue,issue,icurrent);	
		gsl_matrix_get_row(v1,sociedade,inow);
	   	gsl_blas_ddot(v1,v_issue,&field1);	

	    	H=0;	
	   	for (ic=0; ic<igraph_vector_size(&neighbors); ic++)
	    	{
			ineigh=(int)VECTOR(neighbors)[ic];			
		   	gsl_matrix_get_row(v2,sociedade,ineigh);
			gsl_blas_ddot(v2,v_issue,&field2);
		
			alfa=1;		
			if(field2*field1>0) alfa = DELTA;
			H+=alfa*field2/(double)P;	
		 } 
		 E+=-field1*H;
	}									
     }
     E=E/2;  
     gsl_vector_free(v1);
     gsl_vector_free(v2);	
     gsl_vector_free(v_issue);
     return(E);	
}

//******************************************************************************
//gera vetor aleatorio normalizado
//******************************************************************************
void gera_vetor(gsl_vector * v, int L,gsl_rng * r)
{
	double u,Q;
	int ic;
	for(ic=0;ic<L;ic++)
	{
  		u=gsl_ran_ugaussian(r);
		gsl_vector_set(v,ic,u);
	}   
	gsl_blas_ddot(v,v,&Q);
	gsl_vector_scale(v,1/sqrt(Q));
}

//******************************************************************************
//Calcula variação de energia ao substituirmos v0 por v1 no sitio inow
//******************************************************************************
double variacaoE (gsl_vector * v0, gsl_vector * v1, int inow, 
		  gsl_matrix * sociedade, gsl_matrix * issue, 
	          int N, int L, int P, double DELTA, 
		  igraph_t graph, igraph_vector_t neighbors)

{
	      double DeltaE;
	      double E0=0,E1=0,H0,H1;
    	      double field0,field1,field2,alfa; 	
    	      int    icurrent,ineigh,ic;	

              gsl_vector * v2 = gsl_vector_alloc(L);	
              gsl_vector * v_issue = gsl_vector_alloc(L); 
 					
	      for (icurrent=0;icurrent<P;icurrent++)
	      {
		     gsl_matrix_get_row(v_issue,issue,icurrent);
		     gsl_blas_ddot(v0,v_issue,&field0);	
		     gsl_blas_ddot(v1,v_issue,&field1);
		     H0=H1=0;	
		     for(ic=0; ic<igraph_vector_size(&neighbors); ic++)
		     {
			  ineigh=(int)VECTOR(neighbors)[ic];			
	   		  gsl_matrix_get_row(v2,sociedade,ineigh);
			  gsl_blas_ddot(v2,v_issue,&field2);

			  alfa=1;
			  if(field2*field0>0) alfa = DELTA;
			  H0+=alfa*field2/(double)P;
				
			  alfa=1;
			  if(field2*field1>0) alfa = DELTA;
			  H1+=alfa*field2/(double)P;
		     } 
		     E0+=-field0*H0;
		     E1+=-field1*H1;
              }
	      DeltaE=E1-E0;
	      gsl_vector_free(v2);	
	      gsl_vector_free(v_issue);			
											
	      return(DeltaE);
}

//******************************************************************************
//  Verifica se histograma é flat no intervalo [E1,E2]
//******************************************************************************
int flatness(gsl_histogram * H, double E1, double E2, double TOL, 
	     int itera, double &meanhist, double &hvalue)
{
	double nbins=0,auxhist;
	size_t i1,i2,ih;
	int flat;		

	meanhist=0;
	hvalue=gsl_histogram_max_val(H);
	gsl_histogram_find(H,E1,&i1);
	gsl_histogram_find(H,E2-0.1,&i2);

	ih=i1;
	//printf("ih=%d i2=%d\n",ih,i2);
	while(ih<=i2)
	{
		auxhist=gsl_histogram_get(H,ih);
		meanhist+=auxhist;
		if(auxhist<hvalue) hvalue=auxhist;
		//printf("H[%d]=%f minH=%f \n",ih,auxhist,hvalue);
		nbins++;
		ih++;
	}
	meanhist=meanhist/(double)nbins;
	
	if ( (hvalue>TOL*meanhist) && (hvalue>0) )  flat=1;
	else flat=0;	

	return(flat);	
}

//******************************************************************************
//Encontra menor valor Ex tal que [Ex,E2] 
//******************************************************************************
double flatwindow(gsl_histogram * H, double EW, double TOL,double &meanhist)
{
	double aux;
	double lower,upper,E;
	size_t iW,ih;
	int iflag;		
	gsl_histogram_find(H,EW,&iW);
	ih=iW;
	do{
		aux=gsl_histogram_get(H,ih);
		//printf("H[%d]=%f k<H>=%f \n",ih,aux,meanhist*TOL);
		if((aux>TOL*meanhist)&&(ih>0)){
			 iflag=1;
			 ih--;
		}	
		else iflag=0;
	}while(iflag==1);
	ih++;
	gsl_histogram_get_range(H,ih,&lower,&upper);
	E=lower;
	return(E);
}






