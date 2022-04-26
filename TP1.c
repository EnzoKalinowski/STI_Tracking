#include <stdio.h>
#include <stdlib.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include <math.h>
#include <dirent.h>
#include<string.h>

#define TRESHOLD 42
#define EROSION_SIZE 8

// #define MEDIAN 
// #define AVERAGE
#define TEST_INTEREST_POINT
// #define MAIN_WHILE

void binarize(byte **I, byte **B, int treshold, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			if (I[i][j] > treshold)
			{
				B[i][j] = 255;
			}
			else
			{
				B[i][j] = 0;
			}
		}
	}
}

void binarize_rgb8(rgb8 **I, byte **B, int treshold, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	int moy;
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			moy= (I[i][j].r +I[i][j].g + I[i][j].b)/3;

			if (moy > treshold)
			{
				B[i][j] = 255;
			}
			else
			{
				B[i][j] = 0;
			}
		}
	}
}

void convert_rgb8_to_byte(rgb8 **I, byte **B, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	int moy;
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			moy= (I[i][j].r +I[i][j].g + I[i][j].b)/3;

				B[i][j] = moy;
			
		}
	}
}

void substract(byte **I, byte**I2, byte **D, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	int moy;
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			D[i][j]=abs(I[i][j]-I2[i][j]);
			
		}
	}
}

byte** erosion(byte **I, int erode_size, long nrl, long nrh, long ncl, long nch)
{
	int i, j, k;
	int x, y;
	byte **R;
	R = bmatrix(nrl, nrh, ncl, nch);

	
	for (i = nrl + 1; i < nrh - 1; i++)
	{
		for (j = ncl + 1; j < nch - 1; j++)
		{
			R[i][j] = I[i][j];
			
		}
	}

	for (i = nrl + 1; i < nrh - 1; i++)
	{
		for (j = ncl + 1; j < nch - 1; j++)
		{
			for (x = -1; x <= 1; x++)
				{
					for (y = -1; y <= 1; y++)
					{
						if (I[i + x][j + y] == 0)
						{
							R[i][j] = 0;
						}
					}
				}
		}
	}

	if(erode_size>1)
	{
		return erosion(R,erode_size-1,nrl,nrh,ncl,nch);
	}else
	{
		return R;
	}
	
}

byte** dilatation(byte **I, int erode_size, long nrl, long nrh, long ncl, long nch)
{
	int i, j, k;
	int x, y;
	byte **R;
	R = bmatrix(nrl, nrh, ncl, nch);

	
	for (i = nrl + 1; i < nrh - 1; i++)
	{
		for (j = ncl + 1; j < nch - 1; j++)
		{
			R[i][j] = I[i][j];
			
		}
	}

	for (i = nrl + 1; i < nrh - 1; i++)
	{
		for (j = ncl + 1; j < nch - 1; j++)
		{
			for (x = -1; x <= 1; x++)
				{
					for (y = -1; y <= 1; y++)
					{
						if (I[i + x][j + y] == 255)
						{
							R[i][j] = 255;
						}
					}
				}
		}
	}

	if(erode_size>1)
	{
		return dilatation(R,erode_size-1,nrl,nrh,ncl,nch);
	}else
	{
		return R;
	}
	
}

byte ** ouverture(byte **I, int erode_size, long nrl, long nrh, long ncl, long nch)
{
	byte **R_erosion;
	byte **R_dilatation;

	R_erosion = erosion(I,erode_size, nrl, nrh, ncl, nch);
	R_dilatation = dilatation(R_erosion,erode_size, nrl, nrh, ncl, nch);

	return R_dilatation;
}

byte ** fermeture(byte **I, int erode_size, long nrl, long nrh, long ncl, long nch)
{
	byte **R_erosion;
	byte **R_dilatation;

	R_dilatation = dilatation(I,erode_size, nrl, nrh, ncl, nch);
	R_erosion = erosion(R_dilatation,erode_size, nrl, nrh, ncl, nch);

	return R_erosion;
}

void reference_by_time_average(char *folderPath, rgb8 **R, int nb_img ,long nrl, long nrh, long ncl, long nch)
{
	char filePath[100];
	int sum[nrh][nch][3];
	rgb8 **I;
		
	for(int k=1; k<nb_img;k++)
	{	
		I=rgb8matrix0(nrl, nrh, ncl, nch);

		//load
		sprintf(filePath,"%s%s%03d.%s",folderPath,"fomd",k,"ppm");
		I=LoadPPM_rgb8matrix(filePath, &nrl, &nrh, &ncl, &nch);

		for (int i = nrl; i <= nrh; i++)
		{
			for (int j = ncl; j <= nch; j++)
			{
				sum[i][j][0]+=I[i][j].r;
				sum[i][j][1]+=I[i][j].g;
				sum[i][j][2]+=I[i][j].b;
			}
		}
	}
		for (int i = nrl; i <= nrh; i++)
		{
			for (int j = ncl; j <= nch; j++)
			{
				R[i][j].r=sum[i][j][0]/nb_img;
				R[i][j].g=sum[i][j][1]/nb_img;
				R[i][j].b=sum[i][j][2]/nb_img;
			}
		}
}

static int intCompare(const void *p1, const void *p2)
{
    int int_a = * ( (int*) p1 );
    int int_b = * ( (int*) p2 );

    if ( int_a == int_b ) return 0;
    else if ( int_a < int_b ) return -1;
    else return 1;
}

void reference_by_median(char *folderPath, rgb8 **R, int nb_img, int start_index,long nrl, long nrh, long ncl, long nch)
{

	char filePath[100];
	rgb8 **I;
	int *pxR=malloc(nb_img*sizeof(int));
	int *pxG=malloc(nb_img*sizeof(int));
	int *pxB=malloc(nb_img*sizeof(int));

		
	for (int i = nrl; i <= nrh; i++)
	{
		printf("%d\n",i);
		for (int j = ncl; j <= nch; j++)
		{
			for(int k=start_index; k<nb_img+start_index;k++)
			{
				//load img
				sprintf(filePath,"%s%s%03d.%s",folderPath,"fomd",k,"ppm");
				I=LoadPPM_rgb8matrix(filePath, &nrl, &nrh, &ncl, &nch);

				pxR[k-start_index]=I[i][j].r;
				pxG[k-start_index]=I[i][j].g;
				pxB[k-start_index]=I[i][j].b;

				free_rgb8matrix(I,nrl, nrh, ncl, nch);	
			}
			//calculate median
			R[i][j].r=median(pxR,nb_img);
			R[i][j].g=median(pxG,nb_img);
			R[i][j].b=median(pxB,nb_img);	
		}
	}
}

int median(int * tab, int tab_size)
{
	int i, j, index_min, tmp;

	// sort
	qsort(tab,tab_size,sizeof(int),intCompare);
	/*for(i=0; i<tab_size-1; i++)
	{
		index_min=i;
		for(j=i+1; j<tab_size; j++)
		{
			if(tab[j]<tab[i])
			{
				index_min=j;
			}
		}
		//swap values
		tmp=tab[i];
		tab[i]=tab[index_min];
		tab[index_min]=tmp;
	}*/
	return tab[(tab_size-1)/2];
}

int ** tags(byte **I, long nrl, long nrh, long ncl, long nch)
{
	int i,j, x,y;
	int A,B,C;
	int newtag=1;
	int **tag=imatrix(nrl,nrh,ncl,nch);
	
	for(i=nrl;i<=nrh;i++)
	{	
		for(j=ncl;j<=nch;j++)
		{
			if(i==nrl)
			{
				B=-2;
			}
			else
			{
				B=I[i-1][j];
			}

			if(j==ncl)
			{
				A=-1;
				
			}
			else
			{
				A=I[i][j-1];
			}
			
			C=I[i][j];

			if(C==0)
			{
				tag[i][j]=0;
			}else
			{
				if(C==A && C!=B)//1
			{
				tag[i][j]=tag[i][j-1];
			}
			if(C==B && C!=A)//2
			{
				tag[i][j]=tag[i-1][j];
			}
			if(C!=B && C!=A)//3
			{
				tag[i][j]=newtag;
				newtag++;
			}
			if(C==B && C==A && tag[i][j-1]==tag[i-1][j])//4
			{
				tag[i][j]=tag[i-1][j];
			}
			if(C==B && C==A && tag[i-1][j]!=tag[i][j-1])//5
			{
				tag[i][j]=tag[i-1][j];
				
				for(x=nrl+1;x<=i;x++)
				{	
					for(y=ncl+1;y<nch;y++)
					{
						if(x==i && y==j)
						{
							break;
						}
						if(tag[x][y]==tag[i][j-1])
						{
							tag[x][y]=tag[i-1][j];
						}	
					}
				}
			}
			}		
		}
	}
	return tag;
}

void write_infoTag(int **tag, rgb8 **C, byte **N, char **filename,long nrl, long nrh, long ncl, long nch)
{
	int i,j,k;
	int nbLabel;
	int maxLabel=0;
	int *sizeLabel;
	int label;
	int imin=0,jmin=0,imax=0,jmax=0;
	int sumi=0,sumj=0;
	int cpt=0;
	int moyi,moyj;
	//label maximum
	for(i=nrl;i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			if(tag[i][j]>maxLabel){
				maxLabel=tag[i][j];
			}
		}
	}
	sizeLabel=malloc(maxLabel*sizeof(int));

	
	for(i=nrl;i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			k=tag[i][j]-1;
			sizeLabel[k]++;	
		}
	}

	for (k=0;k<maxLabel;k++)
	{
		if(sizeLabel[k]>0)
		{
			label=k+1;
			for(i=nrl;i<=nrh;i++)
			{
				for(j=ncl;j<=nch;j++)
				{
					if(imin>i)
				{
					imin=i;
				}else if (imax<i)
				{
					imax=i;
				}

				if(jmin>j)
				{
					jmin=j;
				}else if (imax<j)
				{
					jmax=j;
				}
				
				sumi+=i;
				sumj+=j;
				cpt++;
				}
			}
		}
	}	
		
		

	FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("Error opening the file %s", filename);
    }
    // write to the text file

    fprintf(fp, "%d", nbLabel);

    // close the file
    fclose(fp);

}

//----------------------------TP2------------------------------//
void sobel(byte **I, double **Ix, double **Iy, long nrl, long nrh, long ncl, long nch)
{
	int i,j;
	long total;
	int SobelH[3][3]={{-1,0,1},{-2,0,2},{-1,0,1} };
	int SobelV[3][3]={{-1,-2,-1},{0,0,0},{1,2,1} };

	for(i=nrl+1;i<nrh;i++)
	{	
		for(j=ncl+1;j<nch;j++)
		{
			//convolution Sobel horizontal
			total=I[i-1][j-1]*SobelH[0][0] + I[i][j-1]*SobelH[1][0] + I[i+1][j-1]*SobelH[2][0] + I[i-1][j]*SobelH[0][1] + I[i][j]*SobelH[1][1] + I[i+1][j]*SobelH[2][1] + I[i-1][j+1]*SobelH[0][2] + I[i][j+1]*SobelH[1][2] + I[i+1][j+1]*SobelH[2][2];
			Ix[i][j]=total/4;

			//convolution Sobel horizontal
			total=I[i-1][j-1]*SobelV[0][0] + I[i][j-1]*SobelV[1][0] + I[i+1][j-1]*SobelV[2][0] + I[i-1][j]*SobelV[0][1] + I[i][j]*SobelV[1][1] + I[i+1][j]*SobelV[2][1] + I[i-1][j+1]*SobelV[0][2] + I[i][j+1]*SobelV[1][2] + I[i+1][j+1]*SobelV[2][2];
			Iy[i][j]=total/4;
		}
	}

}

double ** create_gaussian_filter(float sigma, int size)
{
	double **G=dmatrix(0,size-1,0,size-1);
	float x,y;
	double sum=0;
	
	for (int i=0;i<size;i++){
		for(int j=0;j<size;j++){
			x=i-size/2;
			y=j-size/2;
			G[i][j]=exp(- (pow(x,2)+pow(y,2))/(2*pow(sigma,2)) )/(2*PI*pow(sigma,2));
			sum+=G[i][j];
		}
	}
	printf("Gaussian filter :\n");
	for (int i=0;i<size;i++){
		for(int j=0;j<size;j++){
			G[i][j]/=sum;
			printf("%lf ",G[i][j]);

		}
		printf("\n");
	}

	return G;
}

double ** harris(byte **I, double** filter, int filter_size, float lambda, long nrl, long nrh, long ncl, long nch)
{
	double ** C;
	double ** Ix;
	double ** Iy;
	double Ix_square, Iy_square, IxIy;

	int max=filter_size/2;

	C= dmatrix(nrl, nrh, ncl, nch);
	Ix= dmatrix(nrl,nrh,ncl,nch);
	Iy= dmatrix(nrl,nrh,ncl,nch);

	sobel(I,Ix,Iy,nrl,nrh,ncl,nch);

	for (int x=nrl+max; x<=nrh-max; x++)
		{
			for(int y=ncl+max; y<=nch-max; y++)
			{
				Ix_square=0;
				Iy_square=0;
				IxIy=0;
				
				for (int i=-max;i<=max;i++)
				{
					for(int j=-max;j<=max;j++)
					{
						Ix_square+= filter[i+max][j+max]*pow(Ix[x+i][y+j],2);
						Iy_square+= filter[i+max][j+max]*pow(Iy[x+i][y+j],2);
						IxIy+= filter[i+max][j+max]*(Ix[x+i][y+j]*Iy[x+i][y+j]);
					}
				}
				C[x][y]=(Ix_square*Iy_square)-(IxIy)-(lambda*pow(Ix_square+Iy_square,2));
			}
		}

	free_dmatrix(Ix,nrl, nrh, ncl, nch);
	free_dmatrix(Iy,nrl, nrh, ncl, nch);
	return C;
}

double ** gradient_direction_interest_points(byte **I, double** filter, int filter_size, long nrl, long nrh, long ncl, long nch)
{
	double ** C;
	double ** Ix;
	double ** Iy;
	double Ix_square,Ix_square_filtered, Iy_square,Iy_square_filtered, IxIy, IxIy_filtered;

	int max=filter_size/2;
	double GradFilter[3][3]={ {1,1,1},{1,0,1},{1,1,1} };

	C= dmatrix(nrl, nrh, ncl, nch);
	Ix=dmatrix(nrl,nrh,ncl,nch);
	Iy=dmatrix(nrl,nrh,ncl,nch);

	sobel(I,Ix,Iy,nrl,nrh,ncl,nch);

	for (int x=nrl+max; x<nrh-max; x++)
		{
			for(int y=ncl+max; y<nch-max; y++)
			{
				Ix_square=0;
				Ix_square_filtered=0;
				Iy_square=0;
				Iy_square_filtered=0;
				IxIy=0;
				IxIy_filtered=0;
				
				for (int i=-max;i<=max;i++)
				{
					for(int j=-max;j<=max;j++)
					{
						Ix_square+= pow(Ix[x+i][y+j],2);
						Ix_square_filtered+= GradFilter[i+max][j+max]*pow(Ix[x+i][y+j],2);
						Iy_square+= pow(Iy[x+i][y+j],2);
						Iy_square_filtered+= GradFilter[i+max][j+max]*pow(Iy[x+i][y+j],2);
						IxIy+=Ix[x+i][y+j]*Iy[x+i][y+j];
						IxIy_filtered+= GradFilter[i+max][j+max]*Ix[x+i][y+j]*Iy[x+i][y+j];
					}
				}

				C[x][y]=(Ix_square*Iy_square_filtered)+(Iy_square*Ix_square_filtered)-(2*IxIy*IxIy_filtered);
			}
		}


	free_dmatrix(Ix,nrl, nrh, ncl, nch);
	free_dmatrix(Iy,nrl, nrh, ncl, nch);
	return C;
}

byte ** convolve(byte **I, double** filter, int filter_size, long nrl, long nrh, long ncl, long nch)
{	
	byte ** out;
	out=bmatrix(nrl,nrh,ncl,nch);

	double ** C;
	C= dmatrix(nrl, nrh, ncl, nch);
	int max=filter_size/2;

	for (int x=nrl+max; x<nrh-max; x++)
			{
				for(int y=ncl+max; y<nch-max; y++)
				{
					
					for (int i=-max;i<=max;i++)
					{
						for(int j=-max;j<=max;j++)
						{
							C[x][y]+= filter[i+max][j+max]*I[x+i][y+i];					
						}
					}
				}
			}

	double max_C=max_dmatrix(C,nrl, nrh, ncl, nch);
	printf("%f\n",max_C);
	for (int i=nrl; i<=nrh; i++)
	{
		for(int j=ncl; j<=nch; j++)
		{
			out[i][j]=(C[i][j]*255)/max_C;
		}
	}
	free_dmatrix(C, nrl, nrh, ncl, nch);
	return out;
}


void main()
{
	DIR *d;
    struct dirent *dir;
	rgb8 **I;
	rgb8 **I2;
	byte **G;//nuance de gris
	byte **G2;
	byte **D;//difference entre 2 images
	byte **B;//binaire de la diff
	byte **M;//morphoMath
	byte **SobelX;
	byte ** SobelY;
	rgb8 **R;//reference
	int **tag;


	long nrl, nrh, ncl, nch;

	char *folderPath="./img/Sequences/Fomd/ppm/";
	char filePath[100];
	char filePath2[100];
	char savePath[100];
	char infoTagFolder[100]="./InfoTag/Fomd/";
	char infoTagFile[100];

#ifdef TEST_INTEREST_POINT
	byte **J;
	rgb8 **Jcolor;
	byte **out;

	Jcolor=LoadPPM_rgb8matrix("./img/Sequences/Fomd/ppm/fomd001.ppm",&nrl, &nrh, &ncl, &nch);
	J=bmatrix(nrl, nrh, ncl, nch);

	// J=LoadPGM_bmatrix("../Images/Test/carreTrou.pgm",&nrl, &nrh, &ncl, &nch);
	convert_rgb8_to_byte(Jcolor,J ,nrl, nrh, ncl, nch);
	SavePGM_bmatrix(J,nrl, nrh, ncl, nch,"./img/GrayScales.pgm");
	out=bmatrix0(nrl, nrh, ncl, nch);

	double** Gauss= create_gaussian_filter(0.3,3);
	double GradFilter[3][3]={ {1,1,1},{1,0,1},{1,1,1} };
	double **C;
	// double **Interest2;
	C=harris(J,Gauss,3,0.01,nrl,nrh,ncl,nch);
	// Interest2= gradient_direction_interest_points(J,GradFilter,3,nrl,nrh,ncl,nch);

	double max_C=max_dmatrix(C,nrl, nrh, ncl, nch);
	for (int i=nrl; i<=nrh; i++)
	{
		for(int j=ncl; j<=nch; j++)
		{			
			if(C[i][j]>0)
			{
				// out[i][j]=((C[i][j]-min_C)*255)/(max_C-min_C);
				out[i][j]=((C[i][j])*255)/(max_C);
			}
		}
	}

	SavePGM_bmatrix(out,nrl, nrh, ncl, nch,"./img/Interest_points_harris.pgm");
	// SavePGM_bmatrix(Interest2,nrl, nrh, ncl, nch,"./img/Interest_points_gradient.pgm");

	free_bmatrix(out,nrl, nrh, ncl, nch);
	free_dmatrix(C,nrl, nrh, ncl, nch);
	// free_bmatrix(Interest2,nrl, nrh, ncl, nch);
	free_bmatrix(J,nrl, nrh, ncl, nch);

#endif

	I=LoadPPM_rgb8matrix("./img/Sequences/Fomd/ppm/fomd001.ppm",&nrl, &nrh, &ncl, &nch);
	R=rgb8matrix(nrl, nrh, ncl, nch);

#ifdef AVERAGE
	reference_by_time_average(folderPath,R,870,nrl, nrh, ncl, nch);
	SavePPM_rgb8matrix(R,nrl, nrh, ncl, nch,"./img/reference_average.ppm");
#endif
#ifdef MEDIAN
	reference_by_median(folderPath,R,870,1,nrl, nrh, ncl, nch);
	SavePPM_rgb8matrix(R,nrl, nrh, ncl, nch,"./img/reference_median.ppm");
#endif

#ifdef MAIN_WHILE
	for(int i=1; i<870;i++)
	{	
		//load
		sprintf(filePath,"%s%s%03d.%s",folderPath,"fomd",i,"ppm");
		I=LoadPPM_rgb8matrix(filePath, &nrl, &nrh, &ncl, &nch);
		// sprintf(filePath2,"%s%s%03d.%s",fold/erPath,"lbox",i+1,"ppm");
		I2=LoadPPM_rgb8matrix("./img/reference_median.ppm", &nrl, &nrh, &ncl, &nch);
		
		//convert gray shades
		G=bmatrix(nrl, nrh, ncl, nch);
		convert_rgb8_to_byte(I, G, nrl, nrh, ncl, nch);
		G2=bmatrix(nrl, nrh, ncl, nch);
		convert_rgb8_to_byte(I2, G2, nrl, nrh, ncl, nch);

		//substract
		D=bmatrix(nrl, nrh, ncl, nch);
		substract(G,G2,D,nrl, nrh, ncl, nch);

		//binarize
		B=bmatrix(nrl, nrh, ncl, nch);
		binarize(D,B,TRESHOLD,nrl, nrh, ncl, nch);

		//save
		sprintf(savePath,"%s%s%3d.%s","./img/Resultats/Fomd/","fomd",i,"pgm");	
		SavePGM_bmatrix(B,nrl, nrh, ncl, nch,savePath);

		//morphomath
		M=bmatrix(nrl, nrh, ncl, nch);
		M=fermeture(B,EROSION_SIZE,nrl, nrh, ncl, nch);

		//tags
		tag=tags(M,nrl, nrh, ncl, nch);
		sprintf(infoTagFile,"%s%s%3d.%s",infoTagFolder,"fomd",i,"txt");	

	

		sprintf(savePath,"%s%s%3d.%s","./img/ResultatsMorpho/Fomd/","fomd",i,"pgm");	
		SavePGM_bmatrix(M,nrl, nrh, ncl, nch,savePath);

		
	}
	
    
	free_rgb8matrix(I,nrl, nrh, ncl, nch);
	free_rgb8matrix(R,nrl, nrh, ncl, nch);
	free_rgb8matrix(I2,nrl, nrh, ncl, nch);
	free_bmatrix(G,nrl, nrh, ncl, nch);
	free_bmatrix(G2,nrl, nrh, ncl, nch);
	free_bmatrix(D,nrl, nrh, ncl, nch);
	free_bmatrix(B,nrl, nrh, ncl, nch);
#endif	

}
