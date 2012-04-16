# include <R.h>
# include <math.h>
# include <Rmath.h>
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include <R_ext/Utils.h>
     
 void R_CheckUserInterrupt(void);

 
/* these functions are used to implement the R version to make the fitting procedure faster*/

/* Pay attention to the arrays of result, it should be initiated with appropriate length*/


/*this function calculates sum , and corrected sum of squares of a matrix respect to a given vector N*/
/*this function has 2 outputs,  meanXN and cssXN as a result, ttention should be paied to initialize the result vectos with appropriate length,
it is nN times ncol(X) */


void covmaker(double sigma2, double tau2eta, double tau2theta, double r[], int nr, double *Cov)
{
int i, j, k,rsum=0,counter=0,rdone,rremain,ddone,dremain;
double offone=tau2theta,offtwo=tau2theta+tau2eta,diag=tau2theta+tau2eta+sigma2;
for (i=0;i<nr;i++){rsum+=r[i];}
rdone=0;
rremain=rsum-r[0];
for (i=0;i<nr;i++)
	{
	ddone=0;
	dremain=r[i];
	for (j=0;j<r[i];j++)
		{
		if (rdone>0)
			{
			for (k=0;k<(rdone);k++)
				{
				Cov[counter]=offone;
				counter++;
				}
			}	
		if (ddone>0) 
			{
			for (k=0;k<(ddone);k++)
				{
				Cov[counter]=offtwo;
				counter++;
				}
			}
		Cov[counter]=diag;
		counter++;
		dremain--;
		ddone++;
		if (dremain>0)
			{
			for (k=0;k<(dremain);k++)
				{
				Cov[counter]=offtwo;
				counter++;
				}
			}
		if(rremain>0)
			{
			for (k=0;k<(rremain);k++)
				{
				Cov[counter]=offone;
				counter++;
				}
			}
		}
rdone=rdone+r[i];
if (i<nr-1) {rremain=rremain-r[i+1];}
	}
}



void Rcovmaker(double *sigma2, double *tau2eta, double *tau2theta, double *r, int *nr, double *Cov)
{
	covmaker(*sigma2, *tau2eta, *tau2theta, r, *nr, Cov);
}

void productmaker (double b[], int nb,double d[], double CholA[], double *bAinvb, double *bAinvd, double *dAinvd)
{
	char myL='L',myN='N';
	int info=0, KD=nb, N=nb, LD=nb+1, LDB=nb, LDA=nb, one=1, i=0;
	double bstar[nb],dstar[nb],bb=0,dd=0,bd=0;

for (i=0;i<nb;i++)
	{
	bstar[i]=b[i];
	dstar[i]=d[i];
	}

/*F77_CALL(dpbtrf)(&myL, &N, &KD , Chol, &LD,&info);/* Now Chol is Cholesky factor of A */
F77_CALL(dtrtrs)(&myL, &myN, &myN, &N, &one, CholA, &LDA,  bstar,&LDB, &info); /* Now bstar is inverse of Chol times b */
F77_CALL(dtrtrs)(&myL, &myN, &myN, &N, &one, CholA, &LDA,  dstar,&LDB, &info); /* Now dstar is inverse of Chol times b */

for (i=0;i<nb;i++)
	{
	bb+=bstar[i]*bstar[i];
	dd+=dstar[i]*dstar[i];
	bd+=bstar[i]*dstar[i];
	}
*bAinvb=bb;
*dAinvd=dd;
*bAinvd=bd;
}


void Rproductmaker (double *b, int *nb,double *d, double *A, double *bAinvb, double *bAinvd, double *dAinvd)
{
	productmaker (b, *nb,d, A, bAinvb, bAinvd, dAinvd);
}

void dmvnorm(double X[], int nX, double Mu[], double Sigma[],double *logdensity)
{
char myL='L',myN='N';
int info=0,N=nX,KD=nX,LD=KD+1,one=1,LDA=nX,LDB=nX,i=0;
double Chol[(nX*nX)];	/*used as Cholesky factor of Sigma*/
double Centered[nX]; /*used to centered density values*/
double d2=0;/*gives Mahalanobis distance*/
double sumlogeigen=0; /*gives sum of log eigen values of Cholesky factor, useful as to be traslated as log determinant of Sigma*/
for (i=0;i<(nX*nX);i++)
	{
	Chol[i]=Sigma[i];
	}
for (i=0;i<nX;i++)
	{
	Centered[i]=X[i]-Mu[i];
	}
	

/*SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )*/
F77_CALL(dpbtrf)(&myL, &N, &KD , Chol, &LD,&info);
/*Now Sigma is B, a lower triangular matrix that Sigma =B'B */

/*DTRTRS(          UPLO, TRANS, DIAG,   N, NRHS, A,      LDA,      B, LDB, INFO )*/
F77_CALL(dtrtrs)(&myL, &myN, &myN, &N, &one, Chol, &LDA,  Centered ,&LDB, &info);
/* Now Centred is inverse of Chol times (X-Mu)*/
for (i=0;i<nX;i++)
	{
	d2 += Centered[i]*Centered[i];
	sumlogeigen += log(Chol[i*nX+i]);
	}

*logdensity=-0.5*(nX*log(2*M_PI)+d2)-sumlogeigen;
}






void Rdmvnorm(double *X, int *nX, double *Mu, double *Sigma,double *logdensity)
{
dmvnorm(X, *nX, Mu, Sigma, logdensity);
}




void cssofmatrix (double X[],int nX, double N[], int nN, double *meanXN, double *cssXN)
	{
	int i=0,j,k,l=0;
	double sumx=0, sumx2=0;
	while (i<nX)
		{
		for (j=0;j<nN;j++)
			{
			for (k=0;k<N[j];k++)
				{
				sumx   += X[i];
				sumx2 += pow(X[i],2);
				i++;
				}
			meanXN[l] = sumx / N[j];
			cssXN[l]  = sumx2 - pow(sumx,2) / N[j];
			l++;
			sumx  = 0;
			sumx2= 0;
			}
		}
	}

/* matsum is like cssofmatrix, but with the result in integer instead of double and does not give the corrected sum of squares*/
void matsum (double X[],int nX, double N[], int nN, double *result)
	{
	int i=0,j,k,l=0;
	double sumx=0;
	while (i<nX)
		{
		for (j=0;j<nN;j++)
			{
			for (k=0;k<N[j];k++)
				{
				sumx   += X[i];
				i++;
				}
			result[l] = sumx;
			l++;
			sumx  = 0;
			}
		}
	}


/*
void sumab(double a[], double b[], int n, double *result)
      {
      *result=a[n]+b[n];
      }

void Rsumab(double *a, double *b,int *n, double *result)
      {
      sumab(a, b, *n, result);
      }

*/

void logmarg0tmG (double sigma2, double tau2eta, double mu, double mean[], double css[], int n,double r[], int nr, double  *result) /* here r works like repno */
{
	int i=0,j;
	while (i<n)
		{
		for (j=0; j<nr;j++)
			{
			result[i] = - r[j]/2 * log(2 * M_PI) - (r[j]-1)/2 * log (sigma2)- log(r[j]*tau2eta+sigma2)/2  - 
			css[i]/(2*sigma2)-pow(mean[i]-mu,2)/(2*(tau2eta+sigma2/r[j]));
			i++;
			}
		}
}

void logmarg0tmAL (double sigma2, double tau2eta, double mu, double mean[], double css[], int n,double r[], int nr, double  *result) /*
it is the same as Gaussian effect model */
{
	int i=0,j;
	while (i<n)
		{
		for (j=0; j<nr;j++)
			{
			result[i] = - r[j]/2 * log(2 * M_PI) - (r[j]-1)/2 * log (sigma2)- log(r[j]*tau2eta+sigma2)/2  - 
			css[i]/(2*sigma2)-pow(mean[i]-mu,2)/(2*(tau2eta+sigma2/r[j]));
			i++;
			}
		}
}

void Rlogmarg0tmAL (double *sigma2, double *tau2eta, double *mu, double *mean, double *css, int *n,double *r, int *nr, double  *result) 
{
	logmarg0tmAL (*sigma2, *tau2eta, *mu, mean, css, *n,r, *nr, result) ;
}

void logmarg1tmG (double sigma2, double tau2eta, double tau2theta, double mu, double mean[], double css[], int n, double r[], int nr, double  *result)
{
logmarg0tmG(sigma2, tau2eta+tau2theta, mu, mean, css, n, r, nr, result);
}


void logmarg1tmAL (double sigma2, double tau2eta, double tau2thetaL, double tau2thetaR, double mu, double mean[], double css[], int n, double r[], int nr, double  *result)
{
double loginside1[n],loginside2[n];

	int i=0,j;
	while (i<n)
		{
		for (j=0; j<nr;j++)
			{
/*			inside1[i]=1/(2*sqrt(tau2thetaL))*
					exp(tau2eta/(2*tau2thetaL)+sigma2/(2*r[j]*tau2thetaL)+(mean[i]-mu)/sqrt(tau2thetaL))*
					pnorm((mu-mean[i]-sigma2/(r[j]*sqrt(tau2thetaL))-tau2eta/sqrt(tau2thetaL))/sqrt(tau2eta+sigma2/r[j])
				,0,1,1,0);*/


			loginside1[i]=-log(2)-0.5*log(tau2thetaL)+tau2eta/(2*tau2thetaL)+sigma2/(2*r[j]*tau2thetaL)+(mean[i]-mu)/sqrt(tau2thetaL)+
			pnorm((mu-mean[i]-sigma2/(r[j]*sqrt(tau2thetaL))-tau2eta/sqrt(tau2thetaL))/sqrt(tau2eta+sigma2/r[j])
				,0,1,1,1);				
				
/*			inside2[i]=1/(2*sqrt(tau2thetaR))*
					 exp(tau2eta/(2*tau2thetaR)+sigma2/(2*r[j]*tau2thetaR)+(mu-mean[i])/sqrt(tau2thetaR))*
					(1-pnorm( (mu-mean[i]+sigma2/(r[j]*sqrt(tau2thetaR))+tau2eta/sqrt(tau2thetaR))/sqrt(tau2eta+sigma2/r[j])
				,0,1,1,0));*/

				
			loginside2[i]=-log(2)-0.5*log(tau2thetaR)+tau2eta/(2*tau2thetaR)+sigma2/(2*r[j]*tau2thetaR)+(mu-mean[i])/sqrt(tau2thetaR)+
					pnorm( (mu-mean[i]+sigma2/(r[j]*sqrt(tau2thetaR))+tau2eta/sqrt(tau2thetaR))/sqrt(tau2eta+sigma2/r[j])
				,0,1,0,1);
				



			if (loginside2[i]>loginside1[i])
			result[i]=-r[j]/2 *log(2*M_PI*sigma2)+0.5*log(2*M_PI*sigma2/r[j])-css[i]/(2*sigma2)+
			loginside2[i]+log(1+exp(loginside1[i]-loginside2[i]));
			else
			result[i]=-r[j]/2 *log(2*M_PI*sigma2)+0.5*log(2*M_PI*sigma2/r[j])-css[i]/(2*sigma2)+
			loginside1[i]+log(1+exp(loginside2[i]-loginside1[i]));
			i++;
			}
		}
}


void Rlogmarg1tmAL(double *sigma2, double *tau2eta, double *tau2thetaL,double *tau2thetaR, double *mu, double *mean, 
		double *css, int *n, double *r, int *nr, double  *result)
	{
		logmarg1tmAL(*sigma2, *tau2eta, *tau2thetaL,*tau2thetaR, *mu, mean, 
		css, *n, r, *nr, result);
	}

	
void logmargtmG (double theta[], double mean[], double css[], int n, double r[], int nr, double *result)
{
double sigma2,tau2eta,tau2theta,mu,p,l0[n], l1[n],suml=0;
int i;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2theta=exp(theta[2]);
mu=theta[3];
p=1/(1+exp(-theta[4]));
logmarg0tmG(sigma2, tau2eta, mu, mean, css, n, r, nr, l0);
logmarg1tmG(sigma2, tau2eta,tau2theta, mu, mean, css, n, r, nr, l1);
	for (i=0;i<n;i++)
	{
	if ( (l1[i]-l0[i])>0) 
	suml+=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	suml+=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
*result=-suml;
}



void logmargtmvsG(double theta[], double mean[], double css[], int nrow,int  ncol , double r[], int nr, double *result)
{
int i, n=nrow*ncol;
double sigma2,tau2eta,tau2theta,mu,p,q,l0[n], l1[n],lvs0[ncol], lvs1[ncol],suml=0,l[n],nrowpointer[1];
nrowpointer[0]=(double) nrow;	
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2theta=exp(theta[2]);
mu=theta[3];
p=1/(1+exp(-theta[4]));
q=1/(1+exp(-theta[5]));

logmarg0tmG(sigma2, tau2eta,mu, mean, css, n, r, nr, l0);
logmarg1tmG(sigma2, tau2eta,tau2theta, mu, mean, css, n, r, nr, l1);

	for (i=0;i<n;i++)
	{
	if ( (l1[i]-l0[i])>0) 
	l[i]=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	l[i]=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}

matsum(l ,n, nrowpointer, 1, lvs1); 
matsum(l0 ,n, nrowpointer, 1, lvs0);
suml=0;

	for (i=0;i<ncol;i++)
	{
	if ( (lvs1[i]-lvs0[i])>0) 
	suml+=lvs1[i] + log(q + (1-q) * exp(lvs0[i]-lvs1[i]));
	else 
	suml+=lvs0[i] + log((1-q) + q * exp(lvs1[i]-lvs0[i])); 
	}
*result=-suml;	
}

void logmargtmvsAL(double theta[], double mean[], double css[], int nrow,int  ncol , double r[], int nr, double *result)
{
int i, n=nrow*ncol;
double sigma2,tau2eta,tau2thetaL,tau2thetaR,mu,p,q,l0[n], l1[n],lvs0[ncol], lvs1[ncol],suml=0,l[n],nrowpointer[1];
nrowpointer[0]=(double) nrow;	
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2thetaL=exp(theta[2]);
tau2thetaR=exp(theta[3]);
mu=theta[4];
p=1/(1+exp(-theta[5]));
q=1/(1+exp(-theta[6]));
logmarg0tmAL(sigma2, tau2eta, mu, mean, css, n, r, nr, l0);
logmarg1tmAL(sigma2, tau2eta,tau2thetaL,tau2thetaR, mu, mean, css, n, r, nr, l1);
	for (i=0;i<n;i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	l[i]=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	l[i]=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
matsum(l ,n, nrowpointer, 1, lvs1); /* sum of l over types*/
matsum(l0 ,n, nrowpointer, 1, lvs0); /*sum of l0 over types*/

	for (i=0;i<ncol;i++)
	{
	if ( (lvs1[i]-lvs0[i])>0) 
	suml+=lvs1[i] + log(q + (1-q) * exp(lvs0[i]-lvs1[i]));
	else 
	suml+=lvs0[i] + log((1-q) + q * exp(lvs1[i]-lvs0[i])); 
	}

*result=-suml;	
}


void RlogmargtmvsG(double *theta, double *mean, double *css, int *nrow,int  *ncol , double *r, int *nr, double *result)
{
	logmargtmvsG(theta, mean, css, *nrow,*ncol , r, *nr, result);
}
void RlogmargtmvsAL(double *theta, double *mean, double *css, int *nrow,int  *ncol , double *r, int *nr, double *result)
{
	logmargtmvsAL(theta, mean, css, *nrow,*ncol , r, *nr, result);
}




void logmargtmAL (double theta[], double mean[], double css[], int n, double r[], int nr, double *result)
{
double sigma2,tau2eta,tau2thetaL,tau2thetaR,mu,p,l0[n], l1[n],suml=0;
int i;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2thetaL=exp(theta[2]);
tau2thetaR=exp(theta[3]);
mu=theta[4];
p=1/(1+exp(-theta[5]));
logmarg0tmAL(sigma2, tau2eta, mu, mean, css, n, r, nr, l0);
logmarg1tmAL(sigma2, tau2eta,tau2thetaL,tau2thetaR, mu, mean, css, n, r, nr, l1);
	for (i=0;i<n;i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	suml+=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	suml+=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
*result=-suml;
}


void RlogmargtmAL (double *theta, double *mean, double *css, int *n, double *r, int *nr, double *result)
{
logmargtmAL (theta, mean, css, *n, r, *nr,result);
}

void RlogmargtmG (double *theta, double *mean, double *css, int *n, double *r, int *nr, double *result)
{
logmargtmG (theta, mean, css, *n, r, *nr,result);
}


void typenoreorder(double vec[],int nvec, int noi,int noj, double *result)
{
int l=0 , j;
result[l]=vec[noi-1];l++;
result[l]=vec[noj-1];l++;
for (j=0; j<nvec; j++)
	{
	if ( (l != (noi-1)) || (l != (noj-1)) ) 
	{result[l] = vec[l];}
	l++;
	}
}


int asint(double d) 
{
int i;
i=d;
return i;
}

int asdouble(int i) 
{
double d;
d=i;
return d;
}

void pickuprowsofmat(double y[],  int nrowy, int ncoly, int irowstart,int nrowpick, int startfill, double *result)
{
int i,j,l=0;
	for (j=1; j<=nrowpick;j++)
	{
	for (i=0 ; i<ncoly; i++ )
		{
		result[l+startfill]=y[i*nrowy+j+irowstart-2];
		l++;
		}
	}
}

void eliminaterowsofmat(double y[],  int nrowy, int ncoly, int irow1,int nrowpick1, 
int irow2, int nrowpick2,int startfill, double *result)
{
int i,j,l=0;
	for (i=1; i<=nrowy;i++)
	{
	for (j=1 ; j<=ncoly; j++ )
		{
		if( ((i>=irow1) && (i< (irow1+nrowpick1))) || ((i>=irow2) && (i< (irow2+nrowpick2))) ) 
			{break;}
			else {
				result[l+startfill]=y[((i-1)+nrowy*(j-1) )];
				l++;
			}
		}
	}
}

void transpose(double y[],int nrowy,int ncoly,double *result)
{
int i,j,l=0;
	for (i=0;i<nrowy;i++)
	{
		for (j=0;j<ncoly;j++)
		{
		result[i+j*nrowy]=y[l];
		l++;
		}
	}
}/* It is re-indexing of a vector in an appropriate way works like transpose re-indexing*/

void yreorder(double y[], int nrowy, int ncoly, double repno[], int nrepno, double typeno[], int ntypeno, int noi, int  noj, double *result)
{
	/* for noi=1 or noj=1 it won't work*/
double repnosum[ntypeno],tresult[nrowy*ncoly], irowstart1=1,irowstart2=1;
int i,j, startfill=0;
matsum (repno, nrepno, typeno, ntypeno, repnosum);

for (i=0;i<(noi-1);i++)
	{
	irowstart1 += repnosum[i];
	}
pickuprowsofmat(y,  nrowy, ncoly, asint(irowstart1), asint(repnosum[(noi-1)] ) , startfill, tresult); 
/*it picks up the first part related to noi  and put it in the first line of "result" matrix*/

for (i=0;i<(noj-1);i++)
	{
	irowstart2 += repnosum[i];
	}
startfill=asint(repnosum[(noi-1)] ) * ncoly;
/*it picks up the second part related to noj  and put it in "result" immidiately after first part taken in the previous step*/
pickuprowsofmat(y,  nrowy, ncoly, asint(irowstart2), asint(repnosum[(noj-1)] ) , startfill, tresult);

startfill=asint(repnosum[(noi-1)]+repnosum[(noj-1)] ) * ncoly;

/*it eiminates the values taken in the previous 2 steps and out it at last*/
	eliminaterowsofmat(y,nrowy,ncoly, irowstart1, asint(repnosum[(noi-1)] ), 
irowstart2, asint(repnosum[(noj-1)] ) ,startfill, tresult);

/*we have arranged by rows but we read the matrix by columns so we need to transpose the result matrix*/
transpose(tresult, nrowy, ncoly, result);
}

void Ryreorder(double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, double *typeno, int *ntypeno, int *noi, int  *noj, double *result)
{
yreorder(y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, *noi, *noj, result);
}

void repnoreorder(double repno[], int nrepno, double typeno[], int ntypeno, int noi, int noj,double  *result)
{
double irowstart1=1, irowstart2=1;
int i;
	
for (i=0;i<(noi-1);i++)
	{
	irowstart1 += typeno[i];
	}
pickuprowsofmat(repno,  nrepno, 1, asint(irowstart1) , asint(typeno[(noi-1)] ) , 0, result); 

for (i=0;i<(noj-1);i++)
	{
	irowstart2 += typeno[i];
	}

pickuprowsofmat(repno,  nrepno, 1, asint(irowstart2) , asint(typeno[(noj-1)] ) , asint(typeno[(noi-1)] ) , result); 

eliminaterowsofmat(repno,nrepno,1,asint(irowstart1) , asint(typeno[(noi-1)]), asint(irowstart2), asint(typeno[(noj-1)]) ,asint(typeno[(noi-1)]+typeno[(noj-1)] )  , result);
}


void Rrepnoreorder(double *repno, int *nrepno, double *typeno, int *ntypeno, int *noi, int *noj,double  *result)
{
	repnoreorder(repno, *nrepno, typeno, *ntypeno, *noi, *noj,result);
}



void logmarg0dataG (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, mu, ybar[nrepno*ncoly],css[nrepno*ncoly], l0[nrepno*ncoly];
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
mu=theta[3];
/*tau2theta=exp(theta[2]); 
p=1/(1+exp(-theta[4])); they are useless*/
cssofmatrix (y, (nrowy*ncoly), repno, nrepno, ybar, css);
logmarg0tmG(sigma2, tau2eta, mu, ybar, css, (nrepno*ncoly ) , repno , nrepno, l0);
matsum (l0, (nrepno*ncoly) , typeno, ntypeno, result);
}


void logmarg0datavsG (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, mu, ybar[nrepno*ncoly],css[nrepno*ncoly], l0[nrepno*ncoly];
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
mu=theta[3];
/*tau2theta=exp(theta[2]); 
p=1/(1+exp(-theta[4])); they are useless*/
cssofmatrix (y, (nrowy*ncoly), repno, nrepno, ybar, css);
logmarg0tmG(sigma2, tau2eta, mu, ybar, css, (nrepno*ncoly ) , repno , nrepno, l0);
matsum (l0, (nrepno*ncoly) , typeno, ntypeno, result);
}

void logmarg0dataAL (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, mu, ybar[nrepno*ncoly],css[nrepno*ncoly], l0[nrepno*ncoly];
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
mu=theta[4];
/*tau2theta=exp(theta[2]); 
p=1/(1+exp(-theta[4])); they are useless*/
cssofmatrix (y, (nrowy*ncoly), repno, nrepno, ybar, css);
logmarg0tmAL(sigma2, tau2eta, mu, ybar, css, (nrepno*ncoly ) , repno , nrepno, l0);
matsum (l0, (nrepno*ncoly) , typeno, ntypeno, result);
}

void logmarg0datavsAL (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, mu, ybar[nrepno*ncoly],css[nrepno*ncoly], l0[nrepno*ncoly];
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
mu=theta[4];
/*tau2theta=exp(theta[2]); 
p=1/(1+exp(-theta[4])); they are useless*/
cssofmatrix (y, (nrowy*ncoly), repno, nrepno, ybar, css);
logmarg0tmAL(sigma2, tau2eta, mu, ybar, css, (nrepno*ncoly ) , repno , nrepno, l0);
matsum (l0, (nrepno*ncoly) , typeno, ntypeno, result);
}
void Rlogmarg0datavsAL (double *y, int *nrowy, int *ncoly, double *repno , int *nrepno, double *typeno, int *ntypeno, double *theta,double  *result)
{
logmarg0datavsAL (y, *nrowy, *ncoly, repno ,*nrepno, typeno, *ntypeno, theta,result);
}

/*
void logmarg1dataG (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, tau2theta, mu, ybar[ntypeno*ncoly],css[ntypeno*ncoly], newrepno[ntypeno];
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2theta=exp(theta[2]); 
mu=theta[3];
p=1/(1+exp(-theta[4])); they are useless
matsum(repno, nrepno, typeno, ntypeno, newrepno) ;
cssofmatrix (y, (nrowy*ncoly), newrepno, ntypeno, ybar, css);
logmarg1tmG(sigma2, tau2eta, tau2theta, mu, ybar, css, (ntypeno*ncoly ) , newrepno , ntypeno, result);
}


/*void matsum (double X[],int nX, double N[], int nN, double *result)
void covmaker(double sigma2, double tau2eta, double tau2theta, double r[], int nr, double *Cov)
void dmvnorm(double X[], int nX, double Mu[], double Sigma[],double *logdensity)*/

void logmarg1vectorG(double y[], int ny, double repno[] , 
int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, tau2theta, mu, logdensity[1], submean[ny],subcov[ny*ny],
	subrepno[nrepno],rtsum[ntypeno],sumrtsum,sumtypeno=0,suby[ny];
int t,v,i,counter=0;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2theta=exp(theta[2]); 
mu=theta[3];
matsum(repno,nrepno,typeno,ntypeno,rtsum);
for (t=0;t<ntypeno;t++)
	{
	counter=0;
	for (i=sumtypeno;i<(sumtypeno+typeno[t]);i++){subrepno[(int)(i-sumtypeno)]=repno[i];}
	for (i=0;i<rtsum[t];i++){submean[i]=mu;}
	covmaker(sigma2,tau2eta,tau2theta,subrepno,typeno[t],subcov);

	for (i=sumrtsum;i<(sumrtsum+rtsum[t]);i++){suby[(int)(i-sumrtsum)]=y[i];}
	dmvnorm(suby,rtsum[t],submean,subcov,logdensity);
	result[t]=logdensity[0];
	sumrtsum+=rtsum[t];
	sumtypeno+=typeno[t];
	}
}

void logmarg1vectorAL(double y[], double CholA[], double logdetA, double repno[] , double typeno, double theta[],double  *result)
{
double logk0, logkL, logiL, logkR, logiR,sumys[(int)typeno],sumy2=0,sum,sigma2,tau2eta,tauthetaL,tauthetaR,
	mu,bb[(int)typeno],bd[(int)typeno],dd[(int)typeno],bL[(int)typeno],bR[(int)typeno],dL[(int)typeno],dR[(int)typeno],cL,cR,sumrepno=0;
int i,s,counter;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tauthetaL=sqrt(exp(theta[2])); 
tauthetaR=sqrt(exp(theta[3])); 
mu=theta[4];
counter=0;
	for (s=0;s<typeno;s++)
	{
		sum=0;
			for (i=0;i<repno[s];i++)
			{
				sum+=y[counter];
				sumy2+=y[counter]*y[counter];
				counter++;
			}
		sumys[s]=sum;
		bL[s]=sumys[s]/sigma2+1/(typeno*tauthetaL);
		dL[s]=-1/sqrt(tau2eta*typeno);
		sumrepno+=repno[s];
		bR[s]=sumys[s]/sigma2-1/(typeno*tauthetaR);
		dR[s]=1/sqrt(tau2eta*typeno);
	}	
cL=(mu-tau2eta/(typeno*tauthetaL))/sqrt(tau2eta/typeno);
cR=(-mu-tau2eta/(typeno*tauthetaR))/sqrt(tau2eta/typeno);

	logk0= -0.5*sumrepno*log(2*M_PI*sigma2)-0.5*typeno*log(2*M_PI*tau2eta)+0.5*log(2*M_PI*tau2eta/typeno)-sumy2/(2*sigma2);
	productmaker (bL, typeno ,dL, CholA, bb, bd, dd);
	logkL= -log(2*tauthetaL)-mu/tauthetaL+tau2eta/((2*typeno)*(tauthetaL*tauthetaL));
	logiL= 0.5*bb[0]+0.5*typeno*log(2*M_PI)-0.5*logdetA+pnorm( (cL+bd[0])/sqrt(1+dd[0]),0,1,1,1);				
	productmaker (bR, typeno ,dR, CholA, bb, bd, dd);
	logkR=  -log(2*tauthetaR)+mu/tauthetaR+tau2eta/((2*typeno)*(tauthetaR*tauthetaR));
	logiR=0.5*bb[0]+0.5*typeno*log(2*M_PI)-0.5*logdetA+pnorm( (cR+bd[0])/sqrt(1+dd[0]),0,1,1,1);
if((logkL+logiL)>(logkR+logiR)) 
{*result=logk0+logkL+logiL+log(1+exp(logkR+logiR-logkL-logiL));} else 
{*result=logk0+logkR+logiR+log(exp(logkL+logiL-logkR-logiR)+1);}
}

void Rlogmarg1vectorAL(double *y, double *CholA, double *logdetA, double *repno , int *typeno, double *theta,double  *result)
{
	logmarg1vectorAL(y, CholA, *logdetA, repno , *typeno, theta, result);
}



void logmarg1dataG(double y[], int nrowy,int ncoly, double repno[] , 
int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, tau2theta, mu, logdensity[1], submean[nrowy],subcov[nrowy*nrowy],
	subrepno[nrepno],rtsum[ntypeno],sumrtsum=0,sumtypeno=0,suby[nrowy],tresult[ntypeno*ncoly];
int t,m,i,counter=0;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2theta=exp(theta[2]); 
mu=theta[3];
matsum(repno,nrepno,typeno,ntypeno,rtsum);
counter=0;
for (t=0;t<ntypeno;t++)
	{
	for (i=sumtypeno;i<(sumtypeno+typeno[t]);i++){subrepno[(int)(i-sumtypeno)]=repno[i];}
	for (i=0;i<rtsum[t];i++){submean[i]=mu;}
	covmaker(sigma2,tau2eta,tau2theta,subrepno,typeno[t],subcov);
	for (m=1;m<(ncoly+1);m++)
		{
		for (i=sumrtsum;i<(sumrtsum+rtsum[t]);i++){suby[(int)(i-sumrtsum)]=y[((m-1)*nrowy+i)];}
		dmvnorm(suby,rtsum[t],submean,subcov,logdensity);
		tresult[counter]=logdensity[0];
		counter++;
		}
	sumrtsum+=rtsum[t];
	sumtypeno+=typeno[t];
	}
transpose(tresult,ntypeno,ncoly,result);
}

void logmarg1dataAL(double y[], int nrowy,int ncoly, double repno[] , 
int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double logdetA,sigma2,tau2eta,logdensity[1], subrepno[nrepno],rtsum[ntypeno],sumrtsum=0,
	sumtypeno=0,suby[nrowy],tresult[ntypeno*ncoly],sum,A[nrowy*nrowy];
int t,m,i,counter=0,s;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
matsum(repno,nrepno,typeno,ntypeno,rtsum);

char myL='L';
for (t=0;t<ntypeno;t++)
	{
	logdetA=0;
	int KD=typeno[t],LD=typeno[t]+1,N=typeno[t],info=0;
	for (i=sumtypeno;i<(sumtypeno+typeno[t]);i++){subrepno[(int)(i-sumtypeno)]=repno[i];}
	for (i=0;i<(typeno[t]*typeno[t]);i++)
		{
			A[i]=-1/(typeno[t]*tau2eta);
		}	
	for (i=0;i<typeno[t];i++){A[i*(int)typeno[t]+i]=subrepno[i]/sigma2+1/tau2eta-1/(typeno[t]*tau2eta); }
	F77_CALL(dpbtrf)(&myL, &N, &KD , A, &LD,&info);/* Now A is Cholesky factor of A */
	for (i=0;i<typeno[t];i++){logdetA+=log(A[(int)typeno[t]*i+i]);}/*sums log of main diagonal of cholesky factor which is half of sum of log determinant*/
	logdetA=2*logdetA; /*translates sum log of diagonals of cholesky factor to determinant of A*/
	for (m=1;m<(ncoly+1);m++)
		{
			for (i=sumrtsum;i<(sumrtsum+rtsum[t]);i++){suby[(int)(i-sumrtsum)]=y[((m-1)*nrowy+i)];}
			logmarg1vectorAL(suby, A, logdetA, subrepno , typeno[t], theta, logdensity);
			tresult[counter]=logdensity[0];
			counter++;
		}
		sumrtsum+=rtsum[t];
		sumtypeno+=typeno[t];
	}
transpose(tresult,ntypeno,ncoly,result);
}

void Rlogmarg1dataAL(double *y, int *nrowy, int *ncoly, double *repno , 
int *nrepno, double *typeno, int *ntypeno, double *theta,double *result)
{
logmarg1dataAL(y, *nrowy, *ncoly, repno , *nrepno, typeno, *ntypeno, theta,result);
}


void Rlogmarg1dataG(double *y, int *nrowy, int *ncoly, double *repno , 
int *nrepno, double *typeno, int *ntypeno, double *theta,double *result)
{
logmarg1dataG(y, *nrowy, *ncoly, repno , *nrepno, typeno, *ntypeno, theta,result);

}



void logmarg1datavsG (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, tau2theta, mu, ybar[ntypeno*ncoly],css[ntypeno*ncoly], newrepno[ntypeno],p,l1[ntypeno*ncoly],l0[ntypeno*ncoly] ;
int i;
p=1/(1+exp(-theta[4]));
logmarg1dataG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta,l1);
logmarg0datavsG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta,l0);
for (i=0; i<(ncoly*ntypeno); i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	result[i]=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	result[i]=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
}



void Rlogmarg1datavsG (double *y, int *nrowy, int *ncoly, double *repno ,
int *nrepno, double *typeno, int *ntypeno, double *theta,double  *result)
{
logmarg1datavsG (y, *nrowy, *ncoly, repno , *nrepno, typeno, *ntypeno, theta,result);
}

void Rlogmarg0datavsG (double *y, int *nrowy, int *ncoly, double *repno ,
int *nrepno, double *typeno, int *ntypeno, double *theta,double  *result)
{
logmarg0datavsG (y, *nrowy, *ncoly, repno , *nrepno, typeno, *ntypeno, theta,result);
}




void logmarg1datavsAL (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double sigma2,tau2eta, tau2thetaL,tau2thetaR, mu, ybar[ntypeno*ncoly],css[ntypeno*ncoly], newrepno[ntypeno],p,l1[ntypeno*ncoly],l0[ntypeno*ncoly] ;
int i;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2thetaL=exp(theta[2]); 
tau2thetaR=exp(theta[3]); 
mu=theta[4];
p=1/(1+exp(-theta[5]));
logmarg1dataAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta,l1);
logmarg0datavsAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta,l0);
for (i=0; i<(ncoly*ntypeno); i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	result[i]=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	result[i]=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
}

void Rlogmarg0dataG (double *y, int *nrowy, int *ncoly, double *repno , int *nrepno, double *typeno, int *ntypeno, double *theta,double  *result)
{
logmarg0dataG (y, *nrowy, *ncoly, repno ,*nrepno, typeno, *ntypeno, theta,result);
}





void sumlogmargdataG (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double l0[ntypeno*ncoly], l1[ntypeno*ncoly],p, suml=0 ;
int i;
p=1/(1+exp(-theta[4]));
logmarg0dataG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l0);
logmarg1dataG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l1);
for (i=0; i<(ncoly*ntypeno); i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	suml+=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	suml+=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
*result=suml;
}

void logmargdataG (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double l0[ntypeno*ncoly], l1[ntypeno*ncoly],p, suml=0 ;
int i;
p=1/(1+exp(-theta[4]));
logmarg0dataG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l0);
logmarg1dataG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l1);
for (i=0; i<(ncoly*ntypeno); i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	result[i]=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	result[i]=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
}



void sumlogmargdatavsG (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double q, suml=0,sumlvs0[ncoly],sumlvs1[ncoly], nosumlvs1[(ncoly*ntypeno)],nosumlvs0[(ncoly*ntypeno)],ntypenopointer[1];
ntypenopointer[0]=(double)ntypeno;
q=1/(1+exp(-theta[5]));
int i;
logmarg0datavsG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, nosumlvs0);
logmarg1datavsG (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta,nosumlvs1);

matsum(nosumlvs0,(ntypeno*ncoly), ntypenopointer ,1,sumlvs0);
matsum(nosumlvs1,(ntypeno*ncoly), ntypenopointer ,1,sumlvs1);

for (i=0; i<(ncoly); i++)
	{
	if ( (sumlvs1[i]-sumlvs0[i])>0) 
	suml+=sumlvs1[i] + log(q + (1-q) * exp(sumlvs0[i]-sumlvs1[i]));
	else 
	suml+=sumlvs0[i] + log((1-q) + q * exp(sumlvs1[i]-sumlvs0[i])); 
	}
*result=suml;
}


void sumlogmargdatavsAL (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double q, suml=0,sumlvs0[ncoly],sumlvs1[ncoly], nosumlvs1[(ncoly*ntypeno)],nosumlvs0[(ncoly*ntypeno)],ntypenopointer[1];
ntypenopointer[0]=(double)ntypeno;
q=1/(1+exp(-theta[6]));
int i;
logmarg0datavsAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, nosumlvs0);
logmarg1datavsAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta,nosumlvs1);

matsum(nosumlvs0,(ntypeno*ncoly), ntypenopointer ,1,sumlvs0);
matsum(nosumlvs1,(ntypeno*ncoly), ntypenopointer ,1,sumlvs1);

for (i=0; i<(ncoly); i++)
	{
	if ( (sumlvs1[i]-sumlvs0[i])>0) 
	suml+=sumlvs1[i] + log(q + (1-q) * exp(sumlvs0[i]-sumlvs1[i]));
	else 
	suml+=sumlvs0[i] + log((1-q) + q * exp(sumlvs1[i]-sumlvs0[i])); 
	}
*result=suml;
}



void Rlogmarg0dataAL (double *y, int *nrowy, int *ncoly, double *repno , int *nrepno, double *typeno, int *ntypeno, double *theta,double  *result)
{
logmarg0dataAL (y, *nrowy, *ncoly, repno ,*nrepno, typeno, *ntypeno, theta,result);
}


void Rlogmarg1datavsAL (double *y, int *nrowy, int *ncoly, double *repno , int *nrepno, double *typeno, int *ntypeno, double *theta,double  *result)
{
logmarg1datavsAL (y, *nrowy, *ncoly, repno ,*nrepno, typeno, *ntypeno, theta,result);
}

void RsumlogmargdatavsG (double *y,int *nrowy,int  *ncoly,double  *repno,int *nrepno,double  *typeno,int *ntypeno,double *theta,double *result)
{
sumlogmargdatavsG (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, result);
}

void RsumlogmargdatavsAL (double *y,int *nrowy,int  *ncoly,double  *repno,int *nrepno,double  *typeno,int *ntypeno,double *theta,double *result)
{
sumlogmargdatavsAL (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, result);
}


void RlogmargdataG (double *y,int *nrowy,int  *ncoly,double  *repno,int *nrepno,double  *typeno,int *ntypeno,double *theta,double *result)
{
logmargdataG (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, result);
}

void logmargdataAL (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double l0[ntypeno*ncoly], l1[ntypeno*ncoly],p, suml=0 ;
int i;
p=1/(1+exp(-theta[5]));
logmarg0dataAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l0);
logmarg1dataAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l1);
for (i=0; i<(ncoly*ntypeno); i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	result[i]=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	result[i]=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
}

void RlogmargdataAL (double *y,int *nrowy,int  *ncoly,double  *repno,int *nrepno,double  *typeno,int *ntypeno,double *theta,double *result)
{
logmargdataAL (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, result);
}


void sumlogmargdataAL (double y[], int nrowy, int ncoly, double repno[] , int nrepno, double typeno[], int ntypeno, double theta[],double  *result)
{
double l0[ntypeno*ncoly], l1[ntypeno*ncoly],p, suml=0 ;
int i;
p=1/(1+exp(-theta[5]));
logmarg0dataAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l0);
logmarg1dataAL (y, nrowy, ncoly, repno , nrepno, typeno, ntypeno, theta, l1);
for (i=0; i<(ncoly*ntypeno); i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	suml +=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	suml +=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
*result=suml;
}

void RsumlogmargdataG (double *y,int *nrowy,int  *ncoly,double  *repno,int *nrepno,double  *typeno,int *ntypeno,double *theta,double *result)
{
sumlogmargdataG (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, result);
}


void RsumlogmargdataAL (double *y,int *nrowy,int  *ncoly,double  *repno,int *nrepno,double  *typeno,int *ntypeno,double *theta,double *result)
{
sumlogmargdataAL (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, result);
}



void ldirich(double N[], int nN,double *result)
{
int i;
double sumN=0, sumlagmmaN1=0 ;
	for (i=0;i<nN;i++)
	{
		sumN += N[i];
		sumlagmmaN1 += lgammafn(N[i]+1);
	}

*result= lgammafn(asdouble(nN))+sumlagmmaN1-log(sumN)-lgammafn(sumN+asdouble(nN));
}

void Rldirich(double *N, int *nN,double *result)
{
 ldirich(N, *nN,result);
}



void sumof2rows ( double mat[],int nrow,int ncol,int noi,int noj,double *result) 
	{
	int k;
	for ( k=0;k<ncol;k++)
		{
		result[k]=mat[noi+k*nrow-1]+mat[noj+k*nrow-1];
		}
	}



	
void select2rows(double mat[],int nrow,int ncol,int noi,int noj,double *result) 
{
	int k;
	for (k=0;k<ncol;k++)
		{
		result[k]=mat[noi+k*nrow-1];
		result[(k+ncol)]=mat[noj+k*nrow-1];
		}
}

void Rsumof2rows ( double *mat,int *nrow,int *ncol,int *i,int *j,double *result) 
	{
	sumof2rows ( mat, *nrow,*ncol,*i,*j,result);
	}

void Rselect2rows ( double *mat,int *nrow,int *ncol,int *i,int *j,double *result) 
	{
	select2rows ( mat, *nrow,*ncol,*i,*j,result);
	}
	


void select2individofy ( double y[],int nrowy, int ncoly, double repno[], int nrepno, double typeno[], int ntypeno, int noi,int noj,double *result) 
{
double repnosum[ntypeno]; 
int irowstart=1,startfill=0, i;
matsum (repno, nrepno, typeno, ntypeno, repnosum);
double tresult[(asint(repnosum[noi-1]+repnosum[noj-1])*ncoly) ];

for (i=0;i<(noi-1);i++)
	{
	irowstart += repnosum[i];
	}

pickuprowsofmat(y,  nrowy, ncoly, asint(irowstart), asint(repnosum[(noi-1)] ) , startfill, tresult);

/*it picks up the first part related to noi  and put it in the first line of "result" matrix*/
irowstart=1;
for (i=0;i<(noj-1);i++)
	{
	irowstart += repnosum[i];
	}
startfill=asint(repnosum[(noi-1)] ) * ncoly;
pickuprowsofmat(y,  nrowy, ncoly, asint(irowstart), asint(repnosum[(noj-1)] ) , startfill, tresult);
/*it picks up the second part related to noj  and put it in "result" immidiately after first part taken in the previous step*/

/*we have arranged by rows but we read the matrix by columns so we need to transpose the result matrix*/
transpose(tresult, asint(repnosum[noi-1]+repnosum[noj-1]), ncoly, result);
}


void Rselect2individofy ( double *y,int *nrowy, int *ncoly, double *repno, 
	int *nrepno, double *typeno, int *ntypeno, int *noi,int *noj,double *result) 
	{
	select2individofy (y,*nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, 
		*noi, *noj, result);
	}
	


	
	
	
	
void distmatrixvsG(double y[], int nrowy, int ncoly, double repno[], int nrepno, double typeno[], int ntypeno, double theta[], double *minvalue,double  *minindex)
{
int i,j,l=0;
double newy[nrowy*ncoly], newrepno[nrepno], newtypeno[ntypeno], initlogmarg[1], initdirich[1], distmat;

newtypeno[0]=typeno[0]+typeno[1];
eliminaterowsofmat(typeno,  ntypeno, 1 , 1, 1, 2, 1, 1, newtypeno);
sumlogmargdatavsG (y, nrowy, ncoly, repno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
ldirich(newtypeno,(ntypeno-1),initdirich);
*minvalue=-initlogmarg[0]-initdirich[0];

minindex[0]=1; minindex[1]=2;

for (i=1; i<=(ntypeno-1);i++)
	{
	for (j=(i+1); j<=ntypeno; j++)
		{
		yreorder(y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, newy);
		repnoreorder(repno, nrepno, typeno, ntypeno, i, j, newrepno);
		newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
		eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
		sumlogmargdatavsG (newy, nrowy, ncoly, newrepno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
		ldirich(newtypeno,(ntypeno-1),initdirich);
		distmat=-initlogmarg[0]-initdirich[0];
		if (distmat<*minvalue)
			{
			*minvalue=distmat;
			minindex[0]=i; minindex[1]=j;
			}
		l++;
		}
	}	
}



void distmatrixmatvsG(double y[], int nrowy, int ncoly, double repno[], int nrepno, double typeno[], int ntypeno, double theta[], double *minvalue,double  *minindex,double *distance)
{
int i,j,l=0;
double newy[nrowy*ncoly], newrepno[nrepno], newtypeno[ntypeno], initlogmarg[1], initdirich[1], distmat;

newtypeno[0]=typeno[0]+typeno[1];
eliminaterowsofmat(typeno,  ntypeno, 1 , 1, 1, 2, 1, 1, newtypeno);
sumlogmargdatavsG (y, nrowy, ncoly, repno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
ldirich(newtypeno,(ntypeno-1),initdirich);
*minvalue=-initlogmarg[0]-initdirich[0];

minindex[0]=1; minindex[1]=2;

for (i=1; i<=(ntypeno-1);i++)
	{
	for (j=(i+1); j<=ntypeno; j++)
		{
		yreorder(y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, newy);
		repnoreorder(repno, nrepno, typeno, ntypeno, i, j, newrepno);
		newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
		eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
		sumlogmargdatavsG (newy, nrowy, ncoly, newrepno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
		ldirich(newtypeno,(ntypeno-1),initdirich);
		distmat=-initlogmarg[0]-initdirich[0];
		distance[l]=distmat;
		if (distmat<*minvalue)
			{
			*minvalue=distmat;
			minindex[0]=i; minindex[1]=j;
			}
		l++;
		}
	}	
}





void distmatrixvsAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, double typeno[], int ntypeno, double theta[], double *minvalue,double  *minindex)
{
int i,j,l=0;
double newy[nrowy*ncoly], newrepno[nrepno], newtypeno[ntypeno], initlogmarg[1], initdirich[1], distmat;

newtypeno[0]=typeno[0]+typeno[1];
eliminaterowsofmat(typeno,  ntypeno, 1 , 1, 1, 2, 1, 1, newtypeno);
sumlogmargdatavsAL (y, nrowy, ncoly, repno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
ldirich(newtypeno,(ntypeno-1),initdirich);
*minvalue=-initlogmarg[0]-initdirich[0];

minindex[0]=1; minindex[1]=2;

for (i=1; i<=(ntypeno-1);i++)
	{
	for (j=(i+1); j<=ntypeno; j++)
		{
		yreorder(y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, newy);
		repnoreorder(repno, nrepno, typeno, ntypeno, i, j, newrepno);
		newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
		eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
		sumlogmargdatavsAL (newy, nrowy, ncoly, newrepno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
		ldirich(newtypeno,(ntypeno-1),initdirich);
		distmat=-initlogmarg[0]-initdirich[0];
		if (distmat<*minvalue)
			{
			*minvalue=distmat;
			minindex[0]=i; minindex[1]=j;
			}
		l++;
		}
	}	
}



void distmatrixmatvsAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, double typeno[], int ntypeno, double theta[], double *minvalue,double  *minindex,double *distance)
{
int i,j,l=0;
double newy[nrowy*ncoly], newrepno[nrepno], newtypeno[ntypeno], initlogmarg[1], initdirich[1], distmat;

newtypeno[0]=typeno[0]+typeno[1];
eliminaterowsofmat(typeno,  ntypeno, 1 , 1, 1, 2, 1, 1, newtypeno);
sumlogmargdatavsAL (y, nrowy, ncoly, repno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
ldirich(newtypeno,(ntypeno-1),initdirich);
*minvalue=-initlogmarg[0]-initdirich[0];

minindex[0]=1; minindex[1]=2;

for (i=1; i<=(ntypeno-1);i++)
	{
	for (j=(i+1); j<=ntypeno; j++)
		{
		yreorder(y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, newy);
		repnoreorder(repno, nrepno, typeno, ntypeno, i, j, newrepno);
		newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
		eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
		sumlogmargdatavsAL (newy, nrowy, ncoly, newrepno, nrepno, newtypeno, (ntypeno-1), theta, initlogmarg);
		ldirich(newtypeno,(ntypeno-1),initdirich);
		distmat=-initlogmarg[0]-initdirich[0];
		distance[l]=distmat;
		if (distmat<*minvalue)
			{
			*minvalue=distmat;
			minindex[0]=i; minindex[1]=j;
			}
		l++;
		}
	}	
}

void RdistmatrixmatvsG(double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
double *typeno, int *ntypeno, double *theta, double *minvalue,double  *minindex,double *distance)
{
distmatrixmatvsG(y, *nrowy, *ncoly, repno, *nrepno, 
typeno, *ntypeno, theta, minvalue,minindex,distance);
}

void RdistmatrixmatvsAL(double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
double *typeno, int *ntypeno, double *theta, double *minvalue,double  *minindex,double *distance)
{
distmatrixmatvsAL(y, *nrowy, *ncoly, repno, *nrepno, 
typeno, *ntypeno, theta, minvalue,minindex,distance);
}


void bclustvsG(double y[], int nrowy, int ncoly, double repno[], int nrepno, double theta[],double  *merge,double *height)
{
double typeno[nrepno], typelabel[nrepno], yclust[ncoly*nrowy], repnoclust[nrepno], typenoclust[nrepno], minvalue[1], yresult[nrowy*ncoly],repnoresult[nrepno],typenoresult[nrepno],minindex[2];
int i,j, nrepnoclust=nrepno, ntypenoclust=nrepno, ntypeno=nrepno ;
/*initialization*/
for 	(i=0; i<nrepno;i++)
	{
		typeno[i]=1;
		typelabel[i]=-(i+1);
		repnoclust[i]=repno[i];
		typenoclust [i]=1;
	}
for ( i=0;i<(nrowy*ncoly);i++)
	{
	yclust[i]=y[i];
	}
for (i=1;i<ntypeno;i++)
	{
	distmatrixvsG( yclust, nrowy,  ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta, minvalue, minindex);
	merge[2*i-2]=typelabel[asint(minindex[0]-1)];
	merge[2*i-1]=typelabel[asint(minindex[1]-1)];
	height[i-1]=minvalue[0];
	yreorder(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]), yresult);
		for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
	repnoreorder(repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]),repnoresult);
		for (j=0;j<nrepnoclust;j++){repnoclust[j]=repnoresult[j];}
		
	typenoresult[0]=typenoclust[asint(minindex[0]-1)]+typenoclust[asint(minindex[1]-1)];
	eliminaterowsofmat(typenoclust,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typenoresult);
		for (j=0;j<ntypenoclust;j++){typenoclust[j]=typenoresult[j];}
		ntypenoclust--;
	}
}









void bclustvsAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, double theta[],double  *merge,double *height)
{
double typeno[nrepno], typelabel[nrepno], yclust[ncoly*nrowy], repnoclust[nrepno], typenoclust[nrepno], minvalue[1], yresult[nrowy*ncoly],repnoresult[nrepno],typenoresult[nrepno],minindex[2];
int i,j, nrepnoclust=nrepno, ntypenoclust=nrepno, ntypeno=nrepno ;
/*initialization*/
for 	(i=0; i<nrepno;i++)
	{
		typeno[i]=1;
		typelabel[i]=-(i+1);
		repnoclust[i]=repno[i];
		typenoclust [i]=1;
	}
for ( i=0;i<(nrowy*ncoly);i++)
	{
	yclust[i]=y[i];
	}
for (i=1;i<ntypeno;i++)
	{
	distmatrixvsAL( yclust, nrowy,  ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta, minvalue, minindex);
	merge[2*i-2]=typelabel[asint(minindex[0]-1)];
	merge[2*i-1]=typelabel[asint(minindex[1]-1)];
	height[i-1]=minvalue[0];
	yreorder(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]), yresult);
		for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
	repnoreorder(repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]),repnoresult);
		for (j=0;j<nrepnoclust;j++){repnoclust[j]=repnoresult[j];}
		
	typenoresult[0]=typenoclust[asint(minindex[0]-1)]+typenoclust[asint(minindex[1]-1)];
	eliminaterowsofmat(typenoclust,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typenoresult);
		for (j=0;j<ntypenoclust;j++){typenoclust[j]=typenoresult[j];}
		ntypenoclust--;
	}
}









void RbclustvsG(double *y, int *nrowy, int *ncoly, 
double *repno, int *nrepno, double *theta,double  *merge,double *height)
{
	bclustvsG(y, *nrowy, *ncoly, repno, *nrepno,  theta, merge, height);
}

void RbclustvsAL(double *y, int *nrowy, int *ncoly, 
double *repno, int *nrepno, double *theta,double  *merge,double *height)
{
	bclustvsAL(y, *nrowy, *ncoly, repno, *nrepno,  theta, merge, height);
}

	
void RdistmatrixvsG(double *y, int *nrowy, int *ncoly, double *repno, 
	int *nrepno, double *typeno, int *ntypeno, double *theta, double *minvalue,double  *minindex)
	{
        distmatrixvsG(y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, minvalue, minindex);
	}


void RdistmatrixvsAL(double *y, int *nrowy, int *ncoly, double *repno, 
	int *nrepno, double *typeno, int *ntypeno, double *theta, double *minvalue,double  *minindex)
	{
        distmatrixvsAL(y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, theta, minvalue, minindex);
	}
	
	
	
	
	
	
	
	
	
	
	
	
void fastdistmatrixG(double y[], int nrowy, int ncoly, double repno[], int nrepno, 
double typeno[], int ntypeno, double l0[], double l[], double suml, double theta[],
double *minvalue, double  *minindex)
	{
int i=1,j=2; 
double repnosum[ntypeno],l0ij[ncoly],l1ij[ncoly],p,lijseparate[2*ncoly],newtypeno[ntypeno],newldirich[1];
p=1/(1+exp(-theta[4]));
matsum (repno, nrepno, typeno, ntypeno, repnosum);

/*initialization		*/

			int k, ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			double yij[ijsize ],combinedtypeno[]={1},combinedrepno[]={(repnosum[i-1]+repnosum[j-1])}, sumlij=0;
			double sumlminusij=0,yminusij[ncoly*nrowy-ijsize],sumlijseparate=0, newsuml;
			sumof2rows (l0, ntypeno,ncoly ,i,j,l0ij);
			select2individofy ( y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);
			logmarg1dataG (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			for (k=0; k<ncoly; k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
				sumlij+=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
					sumlij +=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}
			pickuprowsofmat(l,  ntypeno, ncoly, i, 1 ,0, lijseparate);
			pickuprowsofmat(l,  ntypeno, ncoly, j, 1 ,ncoly, lijseparate);
			for (k=0; k<(2*ncoly );k++)
			{
			sumlijseparate += lijseparate[k];
			}
			newsuml=suml+sumlij-sumlijseparate;
			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			*minvalue=-newsuml-newldirich[0];
			minindex[0]=1; minindex[1]=2;
int m=0;
double distmat;
for (i=1; i<=ntypeno;i++)
		{
		for (j=i+1; j<=ntypeno;j++)
			{
			R_CheckUserInterrupt();
			ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			combinedrepno[0]=(repnosum[i-1]+repnosum[j-1]);
			sumlij=0;sumlminusij=0;sumlijseparate=0;
			double yij[ijsize];
/*			double yminusij[ncoly*nrowy-ijsize];*/
			sumof2rows (l0, ntypeno,ncoly ,i,j,l0ij);
			select2individofy ( y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);
			logmarg1dataG (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			for (k=0; k<ncoly; k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
				sumlij+=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
					sumlij +=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}
			pickuprowsofmat(l,  ntypeno, ncoly, i, 1 ,0, lijseparate);
			pickuprowsofmat(l,  ntypeno, ncoly, j, 1 ,ncoly, lijseparate);
/*			printf("\nthe lsep[2*ncol] is %f \n",lijseparate[2*ncoly-1]);*/
			for (k=0; k<(2*ncoly );k++)
			{
			sumlijseparate += lijseparate[k];
			}
			newsuml=suml+sumlij-sumlijseparate;
			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			distmat=-newsuml-newldirich[0];
			if (distmat<*minvalue)
				{
				*minvalue=distmat;
				minindex[0]=i; minindex[1]=j;
				}
			m++;
		}
		

	}
}



void fastdistmatrixvsG(double y[], int nrowy, int ncoly, double repno[], int nrepno, 
double typeno[], int ntypeno, double lvs0[], double sumlvs0[], double lvs1[], double sumlvs1[], double theta[],
double *minvalue, double  *minindex)
	{
int i=1,j=2; 
double repnosum[ntypeno],l1ij[ncoly],l0ij[ncoly],p,lijseparate[2*ncoly],newtypeno[ntypeno],newldirich[1],q;
p=1/(1+exp(-theta[4]));
q=1/(1+exp(-theta[5]));
matsum (repno, nrepno, typeno, ntypeno, repnosum);


			int k, ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			double yij[ijsize ],combinedtypeno[]={1},combinedrepno[]={repnosum[i-1]+repnosum[j-1]}, sumlvs1ijseparate[ncoly];
			double newsuml=0, newsumlvs1[ncoly];
			sumof2rows (lvs0, ntypeno,ncoly ,i,j,l0ij);
			sumof2rows (lvs1, ntypeno,ncoly ,i,j,sumlvs1ijseparate);
			select2individofy ( y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);
			logmarg1dataG (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			
			for (k=0; k<(ncoly); k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) 
				l1ij[k]=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
				l1ij[k]=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}

			for (k=0; k<ncoly; k++)
				{
				newsumlvs1[k]=sumlvs1[k]-sumlvs1ijseparate[k]+l1ij[k];
				}
			

			for (k=0; k<ncoly; k++)
				{
				if ( (newsumlvs1[k]-sumlvs0[k])>0)
				newsuml+=newsumlvs1[k] + log(q + (1-q) * exp(sumlvs0[k]-newsumlvs1[k]));
				else 
				newsuml+=sumlvs0[k] + log((1-q) + q * exp(newsumlvs1[k]-sumlvs0[k])); 
				}
			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			*minvalue=-newsuml-newldirich[0];
			minindex[0]=1; minindex[1]=2;
int m=0;
double distmat;
for (i=1; i<=ntypeno;i++)
		{
		for (j=i+1; j<=ntypeno;j++)
			{
			R_CheckUserInterrupt();
			ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			combinedrepno[0]=(repnosum[i-1]+repnosum[j-1]);
			newsuml=0;
			double yij[ijsize];
			sumof2rows (lvs0, ntypeno,ncoly ,i,j,l0ij);
			sumof2rows (lvs1, ntypeno,ncoly ,i,j,sumlvs1ijseparate);
			select2individofy (y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);

			logmarg1dataG (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			for (k=0; k<(ncoly); k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) 
				l1ij[k]=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
				l1ij[k]=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}

			for (k=0; k<ncoly; k++)
				{
				newsumlvs1[k]=sumlvs1[k]-sumlvs1ijseparate[k]+l1ij[k];
				}
			
			for (k=0; k<ncoly; k++)
				{
				if ( (newsumlvs1[k]-sumlvs0[k])>0)
				newsuml+=newsumlvs1[k] + log(q + (1-q) * exp(sumlvs0[k]-newsumlvs1[k]));
				else 
				newsuml+=sumlvs0[k] + log((1-q) + q * exp(newsumlvs1[k]-sumlvs0[k])); 
				}

			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			distmat=-newsuml-newldirich[0];
			if (distmat<*minvalue)
				{
				*minvalue=distmat;
				minindex[0]=i; minindex[1]=j;
				}
			m++;
		}
		

	}
}

void RfastdistmatrixvsG (double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
double *typeno, int *ntypeno, double *lvs0, double *sumlvs0, double *lvs1, double *sumlvs1, double *theta,
double *minvalue, double  *minindex,double *result)

	{
        fastdistmatrixvsG(y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, lvs0, sumlvs0, lvs1, sumlvs1, theta,minvalue,minindex);	
	}






	
void fastdistmatrixvsAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, 
double typeno[], int ntypeno, double lvs0[], double sumlvs0[], double lvs1[], double sumlvs1[], double theta[],
double *minvalue, double  *minindex)
	{
int i=1,j=2; 
double repnosum[ntypeno],l1ij[ncoly],l0ij[ncoly],p,lijseparate[2*ncoly],newtypeno[ntypeno],newldirich[1],q;
p=1/(1+exp(-theta[5]));
q=1/(1+exp(-theta[6]));
matsum (repno, nrepno, typeno, ntypeno, repnosum);


			int k, ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			double yij[ijsize ],combinedtypeno[]={1},combinedrepno[]={repnosum[i-1]+repnosum[j-1]}, sumlvs1ijseparate[ncoly];
			double newsuml=0, newsumlvs1[ncoly];
			sumof2rows (lvs0, ntypeno,ncoly ,i,j,l0ij);
			sumof2rows (lvs1, ntypeno,ncoly ,i,j,sumlvs1ijseparate);
			select2individofy ( y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);
			logmarg1dataAL (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			
			for (k=0; k<(ncoly); k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) 
				l1ij[k]=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
				l1ij[k]=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}

			for (k=0; k<ncoly; k++)
				{
				newsumlvs1[k]=sumlvs1[k]-sumlvs1ijseparate[k]+l1ij[k];
				}
			

			for (k=0; k<ncoly; k++)
				{
				if ( (newsumlvs1[k]-sumlvs0[k])>0)
				newsuml+=newsumlvs1[k] + log(q + (1-q) * exp(sumlvs0[k]-newsumlvs1[k]));
				else 
				newsuml+=sumlvs0[k] + log((1-q) + q * exp(newsumlvs1[k]-sumlvs0[k])); 
				}
			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			*minvalue=-newsuml-newldirich[0];
			minindex[0]=1; minindex[1]=2;
int m=0;
double distmat;
for (i=1; i<=ntypeno;i++)
		{
		for (j=i+1; j<=ntypeno;j++)
			{
			R_CheckUserInterrupt();
			ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			combinedrepno[0]=(repnosum[i-1]+repnosum[j-1]);
			newsuml=0;
			double yij[ijsize];
			sumof2rows (lvs0, ntypeno,ncoly ,i,j,l0ij);
			sumof2rows (lvs1, ntypeno,ncoly ,i,j,sumlvs1ijseparate);
			select2individofy (y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);

			logmarg1dataAL (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			for (k=0; k<(ncoly); k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) 
				l1ij[k]=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
				l1ij[k]=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}

			for (k=0; k<ncoly; k++)
				{
				newsumlvs1[k]=sumlvs1[k]-sumlvs1ijseparate[k]+l1ij[k];
				}
			
			for (k=0; k<ncoly; k++)
				{
				if ( (newsumlvs1[k]-sumlvs0[k])>0)
				newsuml+=newsumlvs1[k] + log(q + (1-q) * exp(sumlvs0[k]-newsumlvs1[k]));
				else 
				newsuml+=sumlvs0[k] + log((1-q) + q * exp(newsumlvs1[k]-sumlvs0[k])); 
				}

			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			distmat=-newsuml-newldirich[0];
			if (distmat<*minvalue)
				{
				*minvalue=distmat;
				minindex[0]=i; minindex[1]=j;
				}
			m++;
		}
		

	}
}


void RfastdistmatrixvsAL (double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
double *typeno, int *ntypeno, double *lvs0, double *sumlvs0, double *lvs1, double *sumlvs1, double *theta,
double *minvalue, double  *minindex,double *result)

	{
        fastdistmatrixvsAL(y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, lvs0, sumlvs0, lvs1, sumlvs1, theta,minvalue,minindex);	
	}





void fastdistmatrixAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, 
double typeno[], int ntypeno, double l0[], double l[], double suml, double theta[],
double *minvalue, double  *minindex)
	{
int i=1,j=2; 
double repnosum[ntypeno],l0ij[ncoly],l1ij[ncoly],p,lijseparate[2*ncoly],newtypeno[ntypeno],newldirich[1];
p=1/(1+exp(-theta[5]));
matsum (repno, nrepno, typeno, ntypeno, repnosum);

/*initialization		*/

			int k, ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			double yij[ijsize ],combinedtypeno[]={1},combinedrepno[]={(repnosum[i-1]+repnosum[j-1])}, sumlij=0;
			double sumlminusij=0,yminusij[ncoly*nrowy-ijsize],sumlijseparate=0, newsuml;
			sumof2rows (l0, ntypeno,ncoly ,i,j,l0ij);
			select2individofy ( y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);
			logmarg1dataAL (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			for (k=0; k<ncoly; k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
				sumlij+=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
					sumlij +=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}
			pickuprowsofmat(l,  ntypeno, ncoly, i, 1 ,0, lijseparate);
			pickuprowsofmat(l,  ntypeno, ncoly, j, 1 ,ncoly, lijseparate);
			for (k=0; k<(2*ncoly );k++)
			{
			sumlijseparate += lijseparate[k];
			}
			newsuml=suml+sumlij-sumlijseparate;
			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			*minvalue=-newsuml-newldirich[0];
			minindex[0]=1; minindex[1]=2;
int m=0;
double distmat;
for (i=1; i<=ntypeno;i++)
		{
		for (j=i+1; j<=ntypeno;j++)
			{
			R_CheckUserInterrupt();
			ijsize=(asint(repnosum[i-1]+repnosum[j-1])*ncoly);
			combinedrepno[0]=(repnosum[i-1]+repnosum[j-1]);
			sumlij=0;sumlminusij=0;sumlijseparate=0;
			double yij[ijsize];
			double yminusij[ncoly*nrowy-ijsize];
			sumof2rows (l0, ntypeno,ncoly ,i,j,l0ij);
			select2individofy ( y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, i, j, yij);
			logmarg1dataAL (yij,asint(repnosum[i-1]+repnosum[j-1]) , ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);
			for (k=0; k<ncoly; k++)
				{
				if ( (l1ij[k]-l0ij[k])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
				sumlij+=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
				else 
					sumlij +=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
				}
			pickuprowsofmat(l,  ntypeno, ncoly, i, 1 ,0, lijseparate);
			pickuprowsofmat(l,  ntypeno, ncoly, j, 1 ,ncoly, lijseparate);
/*			printf("\nthe lsep[2*ncol] is %f \n",lijseparate[2*ncoly-1]);*/
			for (k=0; k<(2*ncoly );k++)
			{
			sumlijseparate += lijseparate[k];
			}
			newsuml=suml+sumlij-sumlijseparate;
			eliminaterowsofmat(typeno,  ntypeno, 1 , i, 1, j, 1, 1, newtypeno);
			newtypeno[0]=typeno[(i-1)]+typeno[(j-1)];
			ldirich(newtypeno,(ntypeno-1),newldirich);
			distmat=-newsuml-newldirich[0];
			if (distmat<*minvalue)
				{
				*minvalue=distmat;
				minindex[0]=i; minindex[1]=j;
				}
			m++;
		}
	}
}















void RfastdistmatrixG (double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
double *typeno, int *ntypeno, double *l0, double *l, double *suml, double *theta, double *minvalue, double *minindex)
	{
  fastdistmatrixG (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, l0, l, *suml, theta,minvalue,minindex);
	}

void RfastdistmatrixAL (double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
double *typeno, int *ntypeno, double *l0, double *l, double *suml, double *theta, double *minvalue, double *minindex)
	{
  fastdistmatrixAL (y, *nrowy, *ncoly, repno, *nrepno, typeno, *ntypeno, l0, l, *suml, theta,minvalue,minindex);
	}

	

	
void l0arrange(double l0[],int nrow, int ncol,int noi,int noj,double *result)
	{
	double tresult[(nrow-1)*ncol];
	sumof2rows (l0, nrow,ncol ,noi,noj,tresult);
	eliminaterowsofmat(l0,  nrow, ncol, noi,1,noj,1, ncol, tresult);
	transpose(tresult, (nrow-1), ncol,result);
	}

	

void Rl0arrange(double *l0,int *nrow, int *ncol,int *noi,int *noj,double *result)
	{
	l0arrange(l0,*nrow, *ncol,*noi,*noj,result);
	}
	
	
void l1arrangeG(double l1[],int nrowl1, int ncoll1, 
	double y[], int nrowy, int ncoly, double repno[], int nrepno, 
	double typeno[], int ntypeno,double theta[], int noi,int noj,double *result)
	{
	double repnosum[ntypeno],tresult[(nrowl1-1)*ncoll1];
	matsum (repno, nrepno, typeno, ntypeno, repnosum);

	double yij[(asint(repnosum[noi-1]+repnosum[noj-1])*ncoly)],
		     combinedtypeno[]={1},combinedrepno[]={(repnosum[noi-1]+repnosum[noj-1])};
		     


	select2individofy (y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, 
		noi, noj, yij);

     logmarg1dataG (yij,asint(repnosum[noi-1]+repnosum[noj-1]) , 
	     ncoly, combinedrepno, 1, combinedtypeno, 1, theta, tresult);
	eliminaterowsofmat(l1,  nrowl1, ncoll1, noi,1,noj,1, ncoly, tresult);
	transpose(tresult, (nrowl1-1), ncoll1,result);
	}

	
	
void l1arrangevsG(double lvs1[],double lvs0[],int nrowlvs1, int ncollvs1, 
	double y[], int nrowy, int ncoly, double repno[], int nrepno, 
	double typeno[], int ntypeno,double theta[], int noi,int noj,double *result)
	{
	double repnosum[ntypeno],tresult[(nrowlvs1-1)*ncollvs1],l0ij[ncollvs1],l1ij[ncollvs1],p;
	int k;
	p=1/(1+exp(-theta[4]));
	matsum (repno, nrepno, typeno, ntypeno, repnosum);

	double yij[(asint(repnosum[noi-1]+repnosum[noj-1])*ncoly)],
		     combinedtypeno[]={1},combinedrepno[]={(repnosum[noi-1]+repnosum[noj-1])};
		     

	sumof2rows (lvs0, ntypeno,ncoly ,noi,noj,l0ij);

	select2individofy (y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, 
		noi, noj, yij);

	
     logmarg1dataG (yij,asint(repnosum[noi-1]+repnosum[noj-1]) , 
	     ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);

			
		for (k=0; k<(ncoly); k++)
			{
			if ( (l1ij[k]-l0ij[k])>0) 
			tresult[k]=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
			else 
			tresult[k]=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
			}
		     
	
	eliminaterowsofmat(lvs1,  nrowlvs1, ncollvs1, noi,1,noj,1, ncoly, tresult);
	transpose(tresult, (nrowlvs1-1), ncollvs1,result);
	}


void l1arrangevsAL(double lvs1[],double lvs0[],int nrowlvs1, int ncollvs1, 
	double y[], int nrowy, int ncoly, double repno[], int nrepno, 
	double typeno[], int ntypeno,double theta[], int noi,int noj,double *result)
	{
	double repnosum[ntypeno],tresult[(nrowlvs1-1)*ncollvs1],l0ij[ncollvs1],l1ij[ncollvs1],p;
	int k;
	p=1/(1+exp(-theta[5]));
	matsum (repno, nrepno, typeno, ntypeno, repnosum);

	double yij[(asint(repnosum[noi-1]+repnosum[noj-1])*ncoly)],
		     combinedtypeno[]={1},combinedrepno[]={(repnosum[noi-1]+repnosum[noj-1])};
		     

	sumof2rows (lvs0, ntypeno,ncoly ,noi,noj,l0ij);

	select2individofy (y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, 
		noi, noj, yij);

	
     logmarg1dataAL (yij,asint(repnosum[noi-1]+repnosum[noj-1]) , 
	     ncoly, combinedrepno, 1, combinedtypeno, 1, theta, l1ij);

			
		for (k=0; k<(ncoly); k++)
			{
			if ( (l1ij[k]-l0ij[k])>0) 
			tresult[k]=l1ij[k] + log(p + (1-p) * exp(l0ij[k]-l1ij[k]));
			else 
			tresult[k]=l0ij[k] + log((1-p) + p * exp(l1ij[k]-l0ij[k])); 
			}
		     
	
	eliminaterowsofmat(lvs1,  nrowlvs1, ncollvs1, noi,1,noj,1, ncoly, tresult);
	transpose(tresult, (nrowlvs1-1), ncollvs1,result);
	}
	

void Rl1arrangevsG(double *lvs1,double *lvs0,int *nrowlvs1, int *ncollvs1, 
	double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
	double *typeno, int *ntypeno,double *theta, int *noi,int *noj,double *result)
	{
 l1arrangevsG(lvs1,lvs0,*nrowlvs1, *ncollvs1, y, *nrowy, *ncoly, repno, *nrepno, 
	typeno, *ntypeno,theta, *noi,*noj,result);
	}
	
	

void l1arrangeAL(double l1[],int nrowl1, int ncoll1, 
	double y[], int nrowy, int ncoly, double repno[], int nrepno, 
	double typeno[], int ntypeno,double theta[], int noi,int noj,double *result)
	{
	double repnosum[ntypeno],tresult[(nrowl1-1)*ncoll1];
	matsum (repno, nrepno, typeno, ntypeno, repnosum);

	double yij[(asint(repnosum[noi-1]+repnosum[noj-1])*ncoly)],
		     combinedtypeno[]={1},combinedrepno[]={(repnosum[noi-1]+repnosum[noj-1])};
		     


	select2individofy (y, nrowy, ncoly, repno, nrepno, typeno, ntypeno, 
		noi, noj, yij);

     logmarg1dataAL (yij,asint(repnosum[noi-1]+repnosum[noj-1]) , 
	     ncoly, combinedrepno, 1, combinedtypeno, 1, theta, tresult);
	eliminaterowsofmat(l1,  nrowl1, ncoll1, noi,1,noj,1, ncoly, tresult);
	transpose(tresult, (nrowl1-1), ncoll1,result);
	}
	
	
void Rl1arrangeG(double *l1,int *nrowl1, int *ncoll1, 
	double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, 
	double *typeno, int *ntypeno,double *theta, int *noi,int *noj,double *result)
	{
	l1arrangeG(l1,*nrowl1, *ncoll1, y, *nrowy, *ncoly, repno, *nrepno, 
	typeno, *ntypeno,theta, *noi,*noj,result);
	}

	
void fastbclustG(double y[], int nrowy, int ncoly, double repno[], int nrepno, double theta[],double  *merge,double *height)
{

double typeno[nrepno], typelabel[nrepno], yclust[ncoly*nrowy], 
	repnoclust[nrepno], typenoclust[nrepno], minvalue[1], yresult[nrowy*ncoly],
	repnoresult[nrepno],typenoresult[nrepno],minindex[2], sumlclust=0,
	l0clust[nrepno*ncoly],l1clust[nrepno*ncoly],lclust[nrepno*ncoly],l0result[nrepno*ncoly],
	l1result[nrepno*ncoly],	p=1/(1+exp(-theta[4])),typelabelresult[nrepno];
int i,j, nrepnoclust=nrepno, ntypenoclust=nrepno, ntypeno=nrepno;

/*initialization*/
for 	(i=0; i<nrepno;i++)
	{
		typeno[i]=1;
		typelabel[i]=-(i+1);
		repnoclust[i]=repno[i];
		typenoclust [i]=1;
	}
for ( i=0;i<(nrowy*ncoly);i++)
	{
	yclust[i]=y[i];
	}

logmarg0dataG (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l0clust);
logmarg1dataG (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l1clust);
		sumlclust=0;
		for (j=0;j<((ntypenoclust)*ncoly); j++)
		{
		if ( (l1clust[j]-l0clust[j])>0) 
			{
			lclust[j]=l1clust[j] + log(p + (1-p) * exp(l0clust[j]-l1clust[j]));
			sumlclust += lclust[j];
			}
		else 
			{
			lclust[j]= l0clust[j] + log((1-p) + p * exp(l1clust[j]-l0clust[j])); 
			sumlclust += lclust[j];
			}
		}
for (i=1;i<ntypeno;i++)
	{
	Rprintf("%d/%d\n",(i+1),(ntypeno));
		fastdistmatrixG( yclust, nrowy,  ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, 
		l0clust,lclust,sumlclust,theta, minvalue, minindex);
		l0arrange(l0clust,ntypenoclust,ncoly,asint(minindex[0]),asint(minindex[1]),l0result);
		
		l1arrangeG(l1clust, ntypenoclust, ncoly, yclust, nrowy, ncoly, repnoclust, nrepnoclust, 
                typenoclust, ntypenoclust, theta, asint(minindex[0]), asint(minindex[1]),l1result);
		
		for (j=0;j<((ntypenoclust-1)*ncoly);j++)
		{
			l1clust[j]=l1result[j];
			l0clust[j]=l0result[j];
		}

		sumlclust=0;
		for (j=0;j<((ntypenoclust-1)*ncoly); j++)
		{
		if ( (l1clust[j]-l0clust[j])>0) 
			{
			lclust[j]=l1clust[j] + log(p + (1-p) * exp(l0clust[j]-l1clust[j]));
			sumlclust += lclust[j];
			}
		else 
			{
			lclust[j]= l0clust[j] + log((1-p) + p * exp(l1clust[j]-l0clust[j])); 
			sumlclust += lclust[j];
			}
		}






	merge[2*i-2]=typelabel[asint(minindex[0]-1)];
	merge[2*i-1]=typelabel[asint(minindex[1]-1)];
	height[i-1]=minvalue[0];
	typelabelresult[0]=i;
	eliminaterowsofmat(typelabel,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typelabelresult);
	 for (j=0;j<ntypenoclust;j++) {typelabel[j]=typelabelresult[j];}
		
	yreorder(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]), yresult);
		for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
	repnoreorder(repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]),repnoresult);
		for (j=0;j<nrepnoclust;j++){repnoclust[j]=repnoresult[j];}

		
		
	typenoresult[0]=typenoclust[asint(minindex[0]-1)]+typenoclust[asint(minindex[1]-1)];
	eliminaterowsofmat(typenoclust,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typenoresult);
		for (j=0;j<ntypenoclust;j++){typenoclust[j]=typenoresult[j];}
		ntypenoclust--;

	}

}









void fastbclustvsG(double y[], int nrowy, int ncoly, double repno[], int nrepno, double theta[],double  *merge,double *height)
{
double typeno[nrepno], typelabel[nrepno], yclust[ncoly*nrowy], 
	repnoclust[nrepno], typenoclust[nrepno], minvalue[1], yresult[nrowy*ncoly],
	repnoresult[nrepno],typenoresult[nrepno],minindex[2], sumlclust=0,
	l0clust[nrepno*ncoly],l1clust[nrepno*ncoly],lclust[nrepno*ncoly],l0result[nrepno*ncoly],
	l1result[nrepno*ncoly], typelabelresult[nrepno],suml0clust[ncoly],suml1clust[ncoly];
int i,j, nrepnoclust=nrepno, ntypenoclust=nrepno, ntypeno=nrepno;

for 	(i=0; i<nrepno;i++)
	{
		typeno[i]=1;
		typelabel[i]=-(i+1);
		repnoclust[i]=repno[i];
		typenoclust [i]=1;
	}
for ( i=0;i<(nrowy*ncoly);i++)
	{
	yclust[i]=y[i];
	}

logmarg0datavsG (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l0clust);
logmarg1datavsG (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l1clust);

double ntypenoclustpointer[]={ntypenoclust};
matsum(l0clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml0clust);
matsum(l1clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml1clust);


for (i=1;i<ntypeno;i++)
	{
		Rprintf("%d/%d\n",(i+1),(ntypeno));

               fastdistmatrixvsG(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust,
		l0clust, suml0clust, l1clust, suml1clust, theta,minvalue,minindex);	

		l0arrange(l0clust,ntypenoclust,ncoly,asint(minindex[0]),asint(minindex[1]),l0result);
		
		l1arrangevsG(l1clust,l0clust,ntypenoclust, ncoly, yclust, nrowy, ncoly, repnoclust, nrepnoclust, 
		typenoclust, ntypenoclust,theta, asint(minindex[0]), asint(minindex[1]),l1result);
		
		for (j=0;j<((ntypenoclust-1)*ncoly);j++)
		{
			l1clust[j]=l1result[j];
			l0clust[j]=l0result[j];
		}






	merge[2*i-2]=typelabel[asint(minindex[0]-1)];
	merge[2*i-1]=typelabel[asint(minindex[1]-1)];
	height[i-1]=minvalue[0];
	typelabelresult[0]=i;
	eliminaterowsofmat(typelabel,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typelabelresult);
	 for (j=0;j<ntypenoclust;j++) {typelabel[j]=typelabelresult[j];}
		
	yreorder(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]), yresult);
		for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
	repnoreorder(repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]),repnoresult);
		for (j=0;j<nrepnoclust;j++){repnoclust[j]=repnoresult[j];}

		
		
	typenoresult[0]=typenoclust[asint(minindex[0]-1)]+typenoclust[asint(minindex[1]-1)];
	eliminaterowsofmat(typenoclust,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typenoresult);
		for (j=0;j<ntypenoclust;j++){typenoclust[j]=typenoresult[j];}
		ntypenoclust--;
		ntypenoclustpointer[0]=ntypenoclust;
		matsum(l0clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml0clust);
		matsum(l1clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml1clust);

	}

}



void fastbclustvsAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, double theta[],double  *merge,double *height)
{
double typeno[nrepno], typelabel[nrepno], yclust[ncoly*nrowy], 
	repnoclust[nrepno], typenoclust[nrepno], minvalue[1], yresult[nrowy*ncoly],
	repnoresult[nrepno],typenoresult[nrepno],minindex[2], sumlclust=0,
	l0clust[nrepno*ncoly],l1clust[nrepno*ncoly],lclust[nrepno*ncoly],l0result[nrepno*ncoly],
	l1result[nrepno*ncoly], typelabelresult[nrepno],suml0clust[ncoly],suml1clust[ncoly];
int i,j, nrepnoclust=nrepno, ntypenoclust=nrepno, ntypeno=nrepno;

for 	(i=0; i<nrepno;i++)
	{
		typeno[i]=1;
		typelabel[i]=-(i+1);
		repnoclust[i]=repno[i];
		typenoclust [i]=1;
	}
for ( i=0;i<(nrowy*ncoly);i++)
	{
	yclust[i]=y[i];
	}

logmarg0datavsAL (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l0clust);
logmarg1datavsAL (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l1clust);

double ntypenoclustpointer[]={ntypenoclust};
matsum(l0clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml0clust);
matsum(l1clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml1clust);


for (i=1;i<ntypeno;i++)
	{
		Rprintf("%d/%d\n",(i+1),(ntypeno));

               fastdistmatrixvsAL(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust,
		l0clust, suml0clust, l1clust, suml1clust, theta,minvalue,minindex);	

		l0arrange(l0clust,ntypenoclust,ncoly,asint(minindex[0]),asint(minindex[1]),l0result);
		
		l1arrangevsAL(l1clust,l0clust,ntypenoclust, ncoly, yclust, nrowy, ncoly, repnoclust, nrepnoclust, 
		typenoclust, ntypenoclust,theta, asint(minindex[0]), asint(minindex[1]),l1result);
		
		for (j=0;j<((ntypenoclust-1)*ncoly);j++)
		{
			l1clust[j]=l1result[j];
			l0clust[j]=l0result[j];
		}






	merge[2*i-2]=typelabel[asint(minindex[0]-1)];
	merge[2*i-1]=typelabel[asint(minindex[1]-1)];
	height[i-1]=minvalue[0];
	typelabelresult[0]=i;
	eliminaterowsofmat(typelabel,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typelabelresult);
	 for (j=0;j<ntypenoclust;j++) {typelabel[j]=typelabelresult[j];}
		
	yreorder(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]), yresult);
		for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
	repnoreorder(repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]),repnoresult);
		for (j=0;j<nrepnoclust;j++){repnoclust[j]=repnoresult[j];}

		
		
	typenoresult[0]=typenoclust[asint(minindex[0]-1)]+typenoclust[asint(minindex[1]-1)];
	eliminaterowsofmat(typenoclust,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typenoresult);
		for (j=0;j<ntypenoclust;j++){typenoclust[j]=typenoresult[j];}
		ntypenoclust--;
		ntypenoclustpointer[0]=ntypenoclust;
		matsum(l0clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml0clust);
		matsum(l1clust,(ntypenoclust*ncoly), ntypenoclustpointer ,1,suml1clust);

	}

}

void RfastbclustvsG(double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, double *theta,double  *merge,double *height)
{
 fastbclustvsG(y, *nrowy,*ncoly, repno, *nrepno, theta,merge,height);
}

void RfastbclustvsAL(double *y, int *nrowy, int *ncoly, double *repno, int *nrepno, double *theta,double  *merge,double *height)
{
 fastbclustvsAL(y, *nrowy,*ncoly, repno, *nrepno, theta,merge,height);
}














void fastbclustAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, double theta[],double  *merge,double *height)
{

/*	double *typeno = (double *) R_alloc(nrepno * sizeof(double));*/
	
	
double typeno[nrepno], typelabel[nrepno], yclust[ncoly*nrowy], 
	repnoclust[nrepno], typenoclust[nrepno], minvalue[1], yresult[nrowy*ncoly],
	repnoresult[nrepno],typenoresult[nrepno],minindex[2], sumlclust=0,
	l0clust[nrepno*ncoly],l1clust[nrepno*ncoly],lclust[nrepno*ncoly],l0result[nrepno*ncoly],
	l1result[nrepno*ncoly],	p=1/(1+exp(-theta[5])),typelabelresult[nrepno];
int i,j, nrepnoclust=nrepno, ntypenoclust=nrepno, ntypeno=nrepno;

/*initialization*/
for 	(i=0; i<nrepno;i++)
	{
		typeno[i]=1;
		typelabel[i]=-(i+1);
		repnoclust[i]=repno[i];
		typenoclust [i]=1;
	}
for ( i=0;i<(nrowy*ncoly);i++)
	{
	yclust[i]=y[i];
	}

logmarg0dataAL (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l0clust);
logmarg1dataAL (yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, theta,l1clust);
		sumlclust=0;
		for (j=0;j<((ntypenoclust)*ncoly); j++)
		{
		if ( (l1clust[j]-l0clust[j])>0) 
			{
			lclust[j]=l1clust[j] + log(p + (1-p) * exp(l0clust[j]-l1clust[j]));
			sumlclust += lclust[j];
			}
		else 
			{
			lclust[j]= l0clust[j] + log((1-p) + p * exp(l1clust[j]-l0clust[j])); 
			sumlclust += lclust[j];
			}
		}
for (i=1;i<ntypeno;i++)
	{
		Rprintf("%d/%d\n",(i+1),(ntypeno));

		fastdistmatrixAL( yclust, nrowy,  ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, 
		l0clust,lclust,sumlclust,theta, minvalue, minindex);
		l0arrange(l0clust,ntypenoclust,ncoly,asint(minindex[0]),asint(minindex[1]),l0result);
		
		l1arrangeAL(l1clust, ntypenoclust, ncoly, yclust, nrowy, ncoly, repnoclust, nrepnoclust, 
                typenoclust, ntypenoclust, theta, asint(minindex[0]), asint(minindex[1]),l1result);
		
		for (j=0;j<((ntypenoclust-1)*ncoly);j++)
		{
			l1clust[j]=l1result[j];
			l0clust[j]=l0result[j];
		}

		sumlclust=0;
		for (j=0;j<((ntypenoclust-1)*ncoly); j++)
		{
		if ( (l1clust[j]-l0clust[j])>0) 
			{
			lclust[j]=l1clust[j] + log(p + (1-p) * exp(l0clust[j]-l1clust[j]));
			sumlclust += lclust[j];
			}
		else 
			{
			lclust[j]= l0clust[j] + log((1-p) + p * exp(l1clust[j]-l0clust[j])); 
			sumlclust += lclust[j];
			}
		}






	merge[2*i-2]=typelabel[asint(minindex[0]-1)];
	merge[2*i-1]=typelabel[asint(minindex[1]-1)];
	height[i-1]=minvalue[0];
	typelabelresult[0]=i;
	eliminaterowsofmat(typelabel,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typelabelresult);
	 for (j=0;j<ntypenoclust;j++) {typelabel[j]=typelabelresult[j];}
		
	yreorder(yclust, nrowy, ncoly, repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]), yresult);
		for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
	repnoreorder(repnoclust, nrepnoclust, typenoclust, ntypenoclust, asint(minindex[0]), asint(minindex[1]),repnoresult);
		for (j=0;j<nrepnoclust;j++){repnoclust[j]=repnoresult[j];}

		
		
	typenoresult[0]=typenoclust[asint(minindex[0]-1)]+typenoclust[asint(minindex[1]-1)];
	eliminaterowsofmat(typenoclust,  ntypenoclust, 1 , asint(minindex[0]), 1, asint(minindex[1]), 1, 1, typenoresult);
		for (j=0;j<ntypenoclust;j++){typenoclust[j]=typenoresult[j];}
		ntypenoclust--;

	}

}


void RfastbclustG(double *y, int *nrowy, int *ncoly, 
double *repno, int *nrepno, double *theta,double  *merge,double *height)
{
	fastbclustG(y, *nrowy, *ncoly, repno, *nrepno,  theta, merge, height);
}


void RfastbclustAL(double *y, int *nrowy, int *ncoly, 
double *repno, int *nrepno, double *theta,double  *merge,double *height)
{
	fastbclustAL(y, *nrowy, *ncoly, repno, *nrepno,  theta, merge, height);
}




void eliminateindexofvector(double vector[],int nvector, double index[], int nindex, double *result)
	{
		int i,j=0,k=0;
		for (i=0;i<nvector;i++)
		{
		if ((i+1)==asint(index[j])) 
			{
			j++;
			}
			else 
			{
			result[k]=vector[i];
			k++;
			}
		}
		
	}
	
void Reliminateindexofvector(double *vector,int *nvector, double *index, int *nindex, double *result)
	{
		eliminateindexofvector(vector,*nvector, index, *nindex, result);
	}
	
void findindexoflabel (double label[],int nlabel,int desiredlabel,double *result, int *nresult)
	{
		int i,j=0;
		for (i=0;i<nlabel;i++)
		{
			if (asint(label[i])==desiredlabel)
			{
				result[j]=i+1;
				j++;
			}
		}
		*nresult=j;
	}

void Rfindindexoflabel (double *label,int *nlabel,int *desiredlabel,double *result, int *nresult)
	{
	findindexoflabel (label, *nlabel,*desiredlabel,result,nresult);
	}

double maxofvector(double vector[], int nvector)
	{
		double result;
		int i;
		result=vector[0];
		for (i=1;i<nvector;i++)
		{
		if (result<vector[i]) 
			{
				result=vector[i];
			}
		}
	return result;
	}



void yarrangebylabel(double y[], int nrowy, int ncoly, double repno[], 
	int nrepno, double label[], int nlabel,double *ypickedup,double *repnopickedup, double *typenopickedup, int *ntypenopickedup)
	{
	int nindexofdesiredlabel[1],maxlab=asint(maxofvector(label,nlabel));
	double indexofdesiredlabel[nlabel], typickedup[ncoly*nrowy]; 
	int repnostartfill=0,desiredlabel=1,j,i=1,ystartfill=0,irowstarty,k;

		
		for (i=1;i<=maxlab;i++)
		{
			desiredlabel=i;
			findindexoflabel(label, nlabel,desiredlabel,indexofdesiredlabel,nindexofdesiredlabel);
			for  (j=0;j<*nindexofdesiredlabel;j++)
			{
			k=0;irowstarty=1;
			while (k<(asint(indexofdesiredlabel[j])-1))
				{
				irowstarty+=repno[k];
				k++;
				}
			repnopickedup[repnostartfill]=repno[asint(indexofdesiredlabel[j]-1)];
			repnostartfill+=1;
			pickuprowsofmat(y,  nrowy, ncoly, irowstarty, asint(repno[asint(indexofdesiredlabel[j])-1]) , ystartfill, typickedup); 
/*			printf("startfrom %d pick up %d start writing from %d \n",irowstarty,asint(repno[asint(indexofdesiredlabel[j])-1]),ystartfill);*/
			ystartfill+= (asint(repno[asint(indexofdesiredlabel[j])-1])*ncoly);
		}
	        typenopickedup[i-1]=*nindexofdesiredlabel;
		}
transpose(typickedup,nrowy,ncoly,ypickedup);
	*ntypenopickedup=maxlab;
	}
	
void Ryarrangebylabel (double *y, int *nrowy, int *ncoly, double *repno, 
	int *nrepno, double *label, int *nlabel,double *ypickedup, double *repnopickedup, double *typenopickedup,int *ntypenopickedup)
	{
       yarrangebylabel(y, *nrowy, *ncoly, repno, *nrepno, label, *nlabel,ypickedup, repnopickedup,typenopickedup,ntypenopickedup);
	}
	

	

void loglikbylabelG(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
	ldirich(helptypeno, *nhelptypeno, logdirich);
	sumlogmargdataG (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich+*logdirich;
	}

void loglikbylabelvsG(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
	ldirich(helptypeno, *nhelptypeno, logdirich);
	sumlogmargdatavsG (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich+*logdirich;
	}


void loglikbylabelAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
	ldirich(helptypeno, *nhelptypeno, logdirich);
	sumlogmargdataAL (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich+*logdirich;
	}

void loglikbylabelvsAL(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
	ldirich(helptypeno, *nhelptypeno, logdirich);
	sumlogmargdatavsAL (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich+*logdirich;
	}

	
void RloglikbylabelG(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelG(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}

void RloglikbylabelvsG(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelvsG(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}


void RloglikbylabelAL(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelAL(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}

void RloglikbylabelvsAL(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelvsAL(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}
	
void logprobabilitycalculatorG(double y[], int nrowy, int ncoly, double repno[], 
	int nrepno, double label[], int nlabel, int indexpos, double theta[], double *changedlabel, int *nchangedlabel, double *logprob)
	{
	double worklabel[nlabel];
	int i;
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	int maxlab=asint(maxofvector(label,nlabel)),labelround,nhelpvec[1];
	double helpvec[nlabel];
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno],singlelogprob[1], logdirich[1];
	int nhelptypeno[1];

	findindexoflabel (worklabel, nlabel, worklabel[indexpos-1], helpvec, nhelpvec);
	if (nhelpvec[0]>1) {maxlab++;}
	
	for (i=1; i<=maxlab;i++)
		{
		changedlabel[i-1]=i;
		worklabel[indexpos-1]=i;
		yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
		sumlogmargdataG (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, singlelogprob);		
		ldirich(helptypeno, *nhelptypeno, logdirich);
		logprob[i-1]=*singlelogprob+*logdirich;
		}
	
	*nchangedlabel=maxlab;
	}

	
	
	
	
void logprobabilitycalculatorAL(double y[], int nrowy, int ncoly, double repno[], 
	int nrepno, double label[], int nlabel, int indexpos, double theta[], double *changedlabel, int *nchangedlabel, double *logprob)
	{
	double worklabel[nlabel];
	int i;
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	int maxlab=asint(maxofvector(label,nlabel)),labelround,nhelpvec[1];
	double helpvec[nlabel];
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno],singlelogprob[1], logdirich[1];
	int nhelptypeno[1];

	findindexoflabel (worklabel, nlabel, worklabel[indexpos-1], helpvec, nhelpvec);
	if (nhelpvec[0]>1) {maxlab++;}
	
	for (i=1; i<=maxlab;i++)
		{
		changedlabel[i-1]=i;
		worklabel[indexpos-1]=i;
		yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
		sumlogmargdataAL (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, singlelogprob);		
		ldirich(helptypeno, *nhelptypeno, logdirich);
		logprob[i-1]=*singlelogprob+*logdirich;
		}
	
	*nchangedlabel=maxlab;
	}

	

void RlogprobabilitycalculatorG(double *y, int *nrowy, int *ncoly, double *repno, 
	int *nrepno, double *label, int *nlabel, int *indexpos, double *theta, double *changedlabel, int *nchangedlabel, double *logprob)
	{
        logprobabilitycalculatorG(y, *nrowy, *ncoly, repno, 
	*nrepno, label, *nlabel, *indexpos, theta, changedlabel, nchangedlabel,logprob);
	}
	
	void randpermute (int nx, int *result)
	{
	int n,rnd,i, x[nx];
		for (i=1; i<=nx;i++) {x[i-1]=i;}
	for (n=nx; n>0;n--)
		{
		rnd=(int)rand()/(int)(((unsigned)RAND_MAX + 1)/n );
		result[n-1]=x[rnd];
			for (i=rnd+1; i<n;i++)
			{
			x[i-1]=x[i];
			}
		}
	}

	
int rmultinomial (double logprob[], int nlogprob)
	{
		int rnd[nlogprob],i;
		double sumprob=0,prob[nlogprob],logprobadded[nlogprob],maxlogprob;
		maxlogprob=maxofvector(logprob,nlogprob);
		for (i=0;i<nlogprob;i++)
			{
			logprobadded[i]=logprob[i]-maxlogprob;
/*			printf("\n normalized logprob of %d is %d", (i+1), (int)(logprobadded[i]));
*/			}

		for (i=0;i<nlogprob;i++)
			{
			sumprob+=exp(logprobadded[i]);
			}

		
		for (i=0;i<nlogprob;i++)
			{
			prob[i]=exp(logprobadded[i])/sumprob;
/*			printf("\nI sample from %d with percent %d", (i+1), (int)(100*prob[i]));
*/			}
		
		rmultinom(1,prob,nlogprob,rnd);
		for (i=0;i<nlogprob;i++)
		{
		if (rnd[i]>0) {return (i+1);}
		}
	}

		
	
void relabel(double label[], int nlabel,double *result)
	{
		int i, j,k=1, npickupindex[1];
		double pickupindex[nlabel], workinglabel[nlabel];
		for (i=0; i<nlabel;i++){workinglabel[i]=label[i];}
		for (i=0; i<nlabel;i++)
		{
			if(workinglabel[i]>0)
			{
			findindexoflabel (workinglabel, nlabel,(int)workinglabel[i] ,pickupindex,npickupindex);
				for (j=0;j<*npickupindex;j++)
				{
				result[(int)pickupindex[j]-1]=k;
				workinglabel[(int)pickupindex[j]-1]=-1;
				}
			k++;
			}
		}
	}
	

void RmcmcG(double *y, int *nrowy, int *ncoly, double *repno, 
	int *nrepno, double *label, int *nlabel, double *theta, int *mciter,double *result)
	{
	int i=1,j, k,randindex[*nrepno],nchangedlabel[1], mysample,l=0;
	double changedlabel[*nrepno],logprob[*nrepno],helplabel[*nrepno],prob[*nrepno], sumprob;
	for (i=1;i<=*mciter;i++)
		{
		randpermute(*nrepno, randindex);
		for (j=0;j<*nrepno;j++)
			{
/*			printf("\n\ncurrent label is\n");
			for (k=0;k<*nrepno;k++)
				{
				printf("  %d  ",(int)label[k]);
				}
*/			logprobabilitycalculatorG(y, *nrowy, *ncoly, repno, *nrepno, label, *nlabel, randindex[j], 
				theta, changedlabel, nchangedlabel, logprob);
			sumprob=0;
			for (k=1;k<nchangedlabel[0];k++)
				{
					sumprob+=exp(logprob[k-1]);
				}
			for (k=1;k<nchangedlabel[0];k++)
				{
					prob[k-1]=exp(logprob[k-1])/sumprob;
				}
/*			for (k=1;k<=nchangedlabel[0];k++)
				{
				printf("\n position %d can get %d with logprobability %d",randindex[j], k, (int)(logprob[k-1]) );
				}
*/			
			mysample=rmultinomial(logprob,nchangedlabel[0]);
/*			printf("\nmysample for position %d is  %d had logprob %d \n",randindex[j], mysample, (int)logprob[(mysample-1)] );
*/			
			label[randindex[j]-1]=mysample;
			relabel(label,*nlabel ,helplabel);
			for (k=0;k<*nrepno;k++)
				{
				label[k]=helplabel[k];
				result[l]=helplabel[k];
				l++;
				}
			}
		}
	}

	
	
	
	
	
	
	
	
	
	








void RmcmcAL(double *y, int *nrowy, int *ncoly, double *repno, 
	int *nrepno, double *label, int *nlabel, double *theta, int *mciter,double *result)
	{
	int i=1,j, k,randindex[*nrepno],nchangedlabel[1], mysample,l=0;
	double changedlabel[*nrepno],logprob[*nrepno],helplabel[*nrepno],prob[*nrepno], sumprob;
	for (i=1;i<=*mciter;i++)
		{
		randpermute(*nrepno, randindex);
		for (j=0;j<*nrepno;j++)
			{
/*			printf("\n\ncurrent label is\n");
			for (k=0;k<*nrepno;k++)
				{
				printf("  %d  ",(int)label[k]);
				}
*/			logprobabilitycalculatorAL(y, *nrowy, *ncoly, repno, *nrepno, label, *nlabel, randindex[j], 
				theta, changedlabel, nchangedlabel, logprob);
			sumprob=0;
			for (k=1;k<nchangedlabel[0];k++)
				{
					sumprob+=exp(logprob[k-1]);
				}
			for (k=1;k<nchangedlabel[0];k++)
				{
					prob[k-1]=exp(logprob[k-1])/sumprob;
				}
/*			for (k=1;k<=nchangedlabel[0];k++)
				{
				printf("\n position %d can get %d with logprobability %d",randindex[j], k, (int)(logprob[k-1]) );
				}
*/			
			mysample=rmultinomial(logprob,nchangedlabel[0]);
/*			printf("\nmysample for position %d is  %d had logprob %d \n",randindex[j], mysample, (int)logprob[(mysample-1)] );
*/			
			label[randindex[j]-1]=mysample;
			relabel(label,*nlabel ,helplabel);
			for (k=0;k<*nrepno;k++)
				{
				label[k]=helplabel[k];
				result[l]=helplabel[k];
				l++;
				}
			}
		}
	}








	
	
	
	
	/*	printf("\n Label \n");
	for (i=0;i<*nrepno;i++)
		{
		printf(" %d ",(int)label[i]);
		}
		printf("\n ReLabel \n");
	for (i=0;i<*nrepno;i++)
		{
		printf(" %d ",(int)helplabel[i]);
		}*/
/*	printf( "\n");*/




void Rrmultinomial (double *logprob,int *nlogprob,int *result)
	{
	int i;
	GetRNGstate();
		for (i=0;i<10000;i++)
		{
		result[i]=rmultinomial(logprob,*nlogprob);
		}
	}
	
void Rrelabel( double *label, int *nlabel,double *result)
{
	relabel(label,*nlabel,result);
}



void loglikbylabelGunif(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
/*	ldirich(helptypeno, *nhelptypeno, logdirich);*/
	sumlogmargdataG (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich/*+*logdirich*/;
	}

void RloglikbylabelGunif(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelGunif(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}

void loglikbylabelALunif(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
/*	ldirich(helptypeno, *nhelptypeno, logdirich);*/
	sumlogmargdataAL (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich/*+*logdirich*/;
	}
	
	
void RloglikbylabelALunif(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelALunif(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}


void loglikbylabelvsGunif(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
/*	ldirich(helptypeno, *nhelptypeno, logdirich);*/
	sumlogmargdatavsG (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich/*+*logdirich*/;
	}

	
void RloglikbylabelvsGunif(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelvsGunif(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}

void loglikbylabelvsALunif(double y[], int nrowy, int ncoly, double repno[], int nrepno, double label[], int nlabel,
	double theta[], double *result)
	{
	double worklabel[nlabel],logdirich[1],resultnodirich[1];
	int i,nhelpno,nhelptypeno[1];
	for (i=0;i<nlabel;i++){worklabel[i]=label[i];}
	double helpy[ncoly*nrowy], helprepno[nrepno], helptypeno[nrepno];
	yarrangebylabel (y, nrowy, ncoly, repno, nrepno, worklabel, nlabel,helpy, helprepno, helptypeno, nhelptypeno);
/*	ldirich(helptypeno, *nhelptypeno, logdirich);*/
	sumlogmargdatavsAL (helpy, nrowy, ncoly, helprepno, nrepno, helptypeno, *nhelptypeno, theta, resultnodirich);
	*result=*resultnodirich/*+*logdirich*/;
	}

void RloglikbylabelvsALunif(double *y,int *nrowy,int *ncoly,double *repno,int *nrepno,double  *label,int *nlabel,double *theta,double *result)
	{
	loglikbylabelvsALunif(y,*nrowy,*ncoly,repno,*nrepno, label,*nlabel,
	theta,result);
	}
