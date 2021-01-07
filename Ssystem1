			      // FILE 2
// This File is to implement the S-System and Runge Kutta 4 Method

double power(double a, double b) // This is a modified power function
{
	 double result;

	 if (a < 0)
	       result = - pow(-a,b);
	 else if (fabs(a) < EPS && b<0)
	       result =  1111111111.11;//MAXFLOAT;
	 else result = (double)pow(a,b);

	 if(result>BIG || result <EPS)
	   result = 0.0;

	 return (result);

}

double CalculateRK(double *PrevData)  // Solving the S-System Equation
{
  double Temp1, Temp2;
  int j;

	 	Temp1=New.Alpha;
		Temp2=New.Beta;
		for(j=1;j<=N;j++)
		{
			Temp1 = Temp1* (double)power(PrevData[j], New.G[j]);
			Temp2 = Temp2* (double)power(PrevData[j], New.H[j]);			
		}
  return Temp1-Temp2;

}


double RKMethod(double *PrevSample, int DSet, int Tstamp) // Calculating New Sample Values
{
   int i;
   double Temp[N+1], D, halfstep= h*0.5;
   double k1,k2,k3,k4;

   for(i=1;i<=N;i++)
   {
	Temp[i] = PrevSample[i];
   }

   k1 = CalculateRK(Temp);
   
   for(i=1;i<=N;i++)
   {
	  Temp[i] = T1[DSet][i][Tstamp];
   }

   Temp[GenePos] = PrevSample[GenePos] +halfstep*k1;
   
   k2= CalculateRK(Temp);

   for(i=1;i<=N;i++)
   {
	  Temp[i]= T2[DSet][i][Tstamp];
   }
   Temp[GenePos] = PrevSample[GenePos] + halfstep*k2;

   k3 = CalculateRK(Temp);

   for(i=1;i<=N;i++)
   {
	   Temp[i] =	T3[DSet][i][Tstamp];
   }
   
   Temp[GenePos] = PrevSample[GenePos] + h*k3;

   k4 = CalculateRK(Temp);


   D = PrevSample[GenePos] + h/6.0*(k1 + 2 *(k2 + k3) + k4);
   

   return D;
}

 // From this function, RK will be called
// Once this function will be called, one set of data values (for N Genes, For T Timestamps)
// will be generated
void S_System(int Set, double &F,int InRefinement=0)
{
    int i,j;
    double PrevData[N+1], D;

	
    for(i=0;i<=T-1;i++)	
    {
		for(j=1;j<=N;j++)
		{
			PrevData[j]=DataSet1[Set][j][i];
		}

		D = RKMethod(PrevData, Set, i);
		if(D>=0 && D<=BIG)
		{
			DataSet2[Set][GenePos][i+1] = D ;
		}
	    else if(D>BIG)
	    {
			DataSet2[Set][GenePos][i+1] = BIG;
	    }
	    else
	    {
		    DataSet2[Set][GenePos][i+1] = 0.0;
	    }
	    
		if(DataSet2[Set][GenePos][i+1]<EPS)		// this might not be necessary
		    DataSet2[Set][GenePos][i+1] = EPS;		
if(InRefinement)
	F += fabs(DataSet2[Set][GenePos][i+1] - DataSet1[Set][GenePos][i+1]);
else
	 F += pow((DataSet2[Set][GenePos][i+1] - DataSet1[Set][GenePos][i+1])/DataSet1[Set][GenePos][i+1],2.00);
	//	
	}   
	
	/*for(i=1;i<=N;i++)
	{
		if(New.G[i]>0 && New.H[i]<0)
		{
			F=F*10;
			break;
		}
		else if(New.G[i]<0 && New.H[i]>0)
		{
			F=F*10;
			break;
		}
	}*/

}

int sort_function( const void *a, const void *b)
{
   double a1 = *(double *)a;
   double b1 = *(double *)b;
   if ( a1-b1==0)
     return 0;
   else if(a1-b1<0)
    return -1;
   else
    return 1;
}

double NomanIba()
{
	double ForSort[2*N], Sum=0;
	int j,k;

	// for copying the G and H values	
		    
		    for(j=1;j<=N;j++)
			{
				ForSort[j-1]=fabs(New.G[j]);  // should be started from 0, but G is started from 1, that is why j-1 is used, but 
				ForSort[N+j-1]=fabs(New.H[j]);				
			}
			qsort((void *)ForSort, 2*N, sizeof(ForSort[0]), sort_function);

			for(k=0;k<2*N-I;k++)
			{
				Sum+= ForSort[k];	   
			}		
   return Sum*c;
}


void CopyMicroArray()
{
 int  j,k;
  for(k=1;k<=M;k++)
      for(j=1;j<=N;j++)	 		  
			  DataSet2[k][j][0]=DataSet1[k][j][0];

}

int NewAhsan()
{
	int i, ZCount=0;
	double DD,S;
	for(i=1;i<= N;i++)
	{
		if(New.G[i]==0)
			ZCount++;
		if(New.H[i]==0)
			ZCount++;
	}
	if(ZCount==0)
		return N;
	else if(ZCount<2*N-I)
	{
		DD= (2*N)/(double)ZCount;
		return DD*N;// extra penalty
	}
	else if(ZCount>=2*N-I && (2*N-ZCount)>=J)
	{
		DD= (2*N)/(double)ZCount;
		//return DD;// no penalty
		S = NomanIba();
		return S*DD;
	}
	else
	{
		DD= (2*N)/(double)ZCount;
		return DD*N;// Huge penalty
	}
}


void CreateDataSet(double &F,int InRefinement=0)
{
 int k;
 double S;
 double SS;
 
 F = 0.0;
 CopyMicroArray();
 S=SS=0;
 
 for(k=1;k<=M;k++)
 {    
	     S_System(k, F,InRefinement);	 		
 } 
 
 if(!InRefinement)
 {
	
	SS = NewAhsan();
	F=N*N*N*F+SS;
 }
 /*else 
 {
	 F=F+NomanIba();
 }*/
// F=F+S+SS;
//printf("%.3lf",F);
 //F=F*(1-(SS-I)/(N));
 //printf("->%.2lf\n",SS);
 
}
