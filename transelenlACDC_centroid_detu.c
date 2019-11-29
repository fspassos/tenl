/*
propagacao de eletrons em cadeias unidmensionais nao lineares com efeitos de campo externo
FREDERICO SALGUEIRO PASSOS
fspassos@ifal.edu.br
Última modificação: 11/11/2019
*/
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*=========================BIBLIOTECAS======================================*/
#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<stdlib.h>
#include<time.h>
/*=========================================================================*/


/*=========PARAMETROS GLOBAIS=============*/
#define rk  1e-4 //Parametro Runge-Kutta
#define tmax 10000 //Tempo de simulaca
#define sitiosiniciais 100
#define dnl 0.5 //variacao do nlparameter
#define newline 100 //Acrescimo de sitios /*apenas usar numeros pares*/
#define tmedio 100
#define omega 0.6
#define F0 0.6

double Fw=0.8*F0;
double phase =0.0;
double Sigma=1;
double domega=0.06;

int data = 500;//alocacao de memoria
int n = sitiosiniciais;



/*======================FUNCOES====================*/
double modulo(double complex fn) //funcao do modulo quadrado da funcao de onda
{
	double probn;
	probn = creal(fn)*creal(fn)+cimag(fn)*cimag(fn);
	return probn;

}

double complex rungefunc (double complex fn,double complex fa,double complex fp, double qui,int j, double t) //aplicacao da DNLS
{
	double complex func;
	double complex scho = -1i;
	double E,j0=data/2;
	E=F0+Fw*sin((omega+domega)*t+phase);
	func = scho*(fa+fp-fn*qui*(modulo(fn))+E*(j-j0)*fn);
	return func;
}








int main()
{

	double timecount1=time(NULL),timecount2; //Contador de tempo de calculo

	/*=============ENTRADA DE VALORES MIN E MAX DA NAO LINEARIDADE ================*/
	float nlparametermin;
	float nlparametermax;
	printf("Entre com o nlparameter min\n");
	scanf("%f",&nlparametermin);
	nlparametermax=nlparametermin;
	float nlparameter=nlparametermin;// inicializacao do parametro nao linear


/*====================CRIACAO DO ARQUIVO DE SAIDA====================*/
	char filename1[250];
	FILE *fil1;
	sprintf(filename1,"X S%.2f t%d nli %.2f Fw=%.2f F0=%.2f dw=%.2f.dat",Sigma,tmax,nlparametermin,Fw,F0,domega);
	fil1=fopen(filename1,"w");

/*	char filename2[250]; //USAR CASO PRECISE DE OUTRO ARQUIVO
	FILE *fil2;
	sprintf(filename2,"tccPxt S%.2f t%d nli %.2f nlf %.2f Fw=%.2f F0=%.2f.dat",Sigma,tmax,nlparametermin,nlparametermax,Fw,F0);
	fil2=fopen(filename2,"w");*/

/*=======INTEIROS DE CONTAGEM=========*/
	int i,j,write=0;
	int countpar=0;
	double countfn=0;


/*===================VARIAVEIS DO SISTEMA====================*/
	double complex f[data];
	double t,qui;
	double par=0.0,partotal=0.0;
	double x=0;

	/*==============PARAMETROS DE RUNGE-KUTTA 8TH ORDER================*/
	double complex k1[data],k2[data],k3[data],k4[data],k5[data],k6[data],k7[data],k8[data],k9[data],k10[data],k11[data];
	double s21=sqrt(21.0);
	double k2a=rk*1.0/2.0;
	double sk2 = k2a;
	double k3a=rk*1.0/4.0;
	double sk3 = k3a*2;
	double k4a=rk*1.0/7.0;
	double k4b=rk*(-7.0-3.0*s21)/98.0;
	double k4c=rk*(21.0+5.0*s21)/49.0;
	double sk4 = k4a+k4b+k4c;
	double k5a=rk*(11.0+s21)/84.0;
	double k5b=rk*(18.0+4.0*s21)/63.0;
	double k5c=rk*(21.0-s21)/252.0;
	double sk5 = k5a+k5b+k5c;
	double k6a=rk*(5.0+s21)/48.0;
	double k6b=rk*(9.0+s21)/36.0;
	double k6c=rk*(-231.0+14.0*s21)/360.0;
	double k6d=rk*(63.0-7.0*s21)/80.0;
	double sk6 = k6a+k6b+k6c+k6d;
	double k7a=rk*(10.0-s21)/42.0;
	double k7b=rk*(-432.0+92.0*s21)/315.0;
	double k7c=rk*(633.0-145.0*s21)/90.0;
	double k7d=rk*(-504.0+115.0*s21)/70.0;
	double k7e=rk*(63.0-13.0*s21)/35.0;
	double sk7 = k7a+k7b+k7c+k7d+k7e;
	double k8a=rk*1.0/14.0;
	double k8b=rk*(14.0-3.0*s21)/126.0;
	double k8c=rk*(13.0-3.0*s21)/63.0;
	double k8d=rk*1.0/9.0;
	double sk8 = k8a+k8b+k8c+k8d;
	double k9a=rk*1.0/32.0;
	double k9b=rk*(91.0-21.0*s21)/576.0;
	double k9c=rk*11.0/72.0;
	double k9d=rk*(-385.0-75.0*s21)/1152.0;
	double k9e=rk*(63.0+13.0*s21)/128.0;
	double sk9 = k9a+k9b+k9c+k9d+k9e;
	double k10a=rk*1.0/14.0;
	double k10b=rk*1.0/9.0;
	double k10c=rk*(-733.0-147.0*s21)/2205.0;
	double k10d=rk*(515.0+111.0*s21)/504.0;
	double k10e=rk*(-51.0-11.0*s21)/56.0;
	double k10f=rk*(132.0+28.0*s21)/245.0;
	double sk10 = k10a+k10b+k10c+k10d+k10e+k10f;
	double k11a=rk*(-42.0+7.0*s21)/18.0;
	double k11b=rk*(-18.0+28.0*s21)/45.0;
	double k11c=rk*(-273.0-53.0*s21)/72.0;
	double k11d=rk*(301.0+53.0*s21)/72.0;
	double k11e=rk*(28.0-28.0*s21)/45.0;
	double k11f=rk*(49.0-7.0*s21)/18.0;
	double sk11 = k11a+k11b+k11c+k11d+k11e+k11f;
	int n1,n2;


/*=============LOOP EXTERNO DE VARIACAO DO domega*/
	do {
n1=(data/2)-n;
n2=(data/2)+n;



/*=======reiniciando a cadeia para cada nlparameter=========*/
		t=0;
		countfn=0;
		n=sitiosiniciais;
		i=0;
		qui=nlparameter;
		for (i=0;i<=data;i=i+1) // loop de limpeza da memoria
		{
			f[i]=0;
		}
/*===================CONDICAO INICIAL DO PACOTE============*/
			if (Sigma==0)//PACOTE TIPO DELTA
				{
					f[data/2]=1;
				}
			else //PACOTE TIPO GAUSSIANA
			{
				double mu, ASigma=0;
				for (j = 0; j < data; j++) {
					mu=j-(data/2);
					f[j]=exp(-mu*mu/Sigma);
					ASigma = ASigma+modulo(f[j]);
				}
				for (j = n1; j < n2; j++) {
					f[j]=f[j]/sqrt(ASigma);
				}
			}



/*======================LOOP TEMPORAL=====================*/

			do
			{
				for ( j = n1; j < n2; j++) {
					k1[j]=rungefunc(f[j],f[j-1],f[j+1],qui,j,t);
				}

				for ( j = n1; j < n2; j++) {
					k2[j]=rungefunc(f[j]+k2a*k1[j],f[j-1]+k2a*k1[j-1],f[j+1]+k2a*k1[j+1],qui,j,t+sk2);
				}

				for ( j = n1; j < n2; j++) {
					k3[j]=rungefunc(f[j]+ k3a*( k1[j] + k2[j]),f[j-1]+ k3a*( k1[j-1] + k2[j-1]),f[j+1]+ k3a*( k1[j+1] + k2[j+1]),qui,j,t+sk3);
				}

				for ( j = n1; j < n2; j++) {
					k4[j]=rungefunc(f[j]+ k1[j]*k4a + k2[j]*k4b + k3[j]*k4c,f[j-1]+ k1[j-1]*k4a + k2[j-1]*k4b + k3[j-1]*k4c,f[j+1]+ k1[j+1]*k4a + k2[j+1]*k4b + k3[j+1]*k4c,qui,j,t+sk4);
				}

				for ( j = n1; j < n2; j++) {
					k5[j]=rungefunc(f[j]+ k1[j]*k5a + k3[j]*k5b + k4[j]*k5c,f[j-1]+ k1[j-1]*k5a + k3[j-1]*k5b + k4[j-1]*k5c,f[j+1]+ k1[j+1]*k5a + k3[j+1]*k5b + k4[j+1]*k5c,qui,j,t+sk5);
				}

				for ( j = n1; j < n2; j++) {
					k6[j]=rungefunc(f[j]+ k1[j]*k6a + k3[j]*k6b + k4[j]*k6c + k5[j]*k6d,f[j-1]+ k1[j-1]*k6a + k3[j-1]*k6b + k4[j-1]*k6c + k5[j-1]*k6d,f[j+1]+ k1[j+1]*k6a + k3[j+1]*k6b + k4[j+1]*k6c + k5[j+1]*k6d,qui,j,t+sk6);
				}

				for ( j = n1; j < n2; j++) {
					k7[j]=rungefunc(f[j]+ k1[j]*k7a + k3[j]*k7b + k4[j]*k7c + k5[j]*k7d+ k6[j]*k7e,f[j-1]+ k1[j-1]*k7a + k3[j-1]*k7b + k4[j-1]*k7c + k5[j-1]*k7d+ k6[j-1]*k7e,f[j+1]+ k1[j+1]*k7a + k3[j+1]*k7b + k4[j+1]*k7c + k5[j+1]*k7d+ k6[j+1]*k7e,qui,j,t+sk7);
				}

				for ( j = n1; j < n2; j++) {
					k8[j]=rungefunc(f[j]+ k1[j]*k8a + k5[j]*k8b + k6[j]*k8c + k7[j]*k8d,f[j-1]+ k1[j-1]*k8a + k5[j-1]*k8b + k6[j-1]*k8c + k7[j-1]*k8d,f[j+1]+ k1[j+1]*k8a + k5[j+1]*k8b + k6[j+1]*k8c + k7[j+1]*k8d,qui,j,t+sk8);
				}

				for ( j = n1; j < n2; j++) {
					k9[j]=rungefunc(f[j]+ k1[j]*k9a + k5[j]*k9b + k6[j]*k9c + k7[j]*k9d + k8[j]*k9e,f[j-1]+ k1[j-1]*k9a + k5[j-1]*k9b + k6[j-1]*k9c + k7[j-1]*k9d + k8[j-1]*k9e,f[j+1]+ k1[j+1]*k9a + k5[j+1]*k9b + k6[j+1]*k9c + k7[j+1]*k9d + k8[j+1]*k9e,qui,j,t+sk9);
				}

				for ( j = n1; j < n2; j++) {
					k10[j]=rungefunc(f[j]+ k1[j]*k10a + k5[j]*k10b + k6[j]*k10c + k7[j]*k10d +k8[j]*k10e + k9[j]*k10f,f[j-1]+ k1[j-1]*k10a + k5[j-1]*k10b + k6[j-1]*k10c + k7[j-1]*k10d +k8[j-1]*k10e + k9[j-1]*k10f,f[j+1]+ k1[j+1]*k10a + k5[j+1]*k10b + k6[j+1]*k10c + k7[j+1]*k10d +k8[j+1]*k10e + k9[j+1]*k10f,qui,j,t+sk10);
				}

				for ( j = n1; j < n2; j++) {
					k11[j]=rungefunc(f[j]+ k5[j]*k11a + k6[j]*k11b + k7[j]*k11c + k8[j]*k11d +k9[j]*k11e + k10[j]*k11f,f[j-1]+ k5[j-1]*k11a + k6[j-1]*k11b + k7[j-1]*k11c + k8[j-1]*k11d +k9[j-1]*k11e + k10[j-1]*k11f,f[j+1]+ k5[j+1]*k11a + k6[j+1]*k11b + k7[j+1]*k11c + k8[j+1]*k11d +k9[j+1]*k11e + k10[j+1]*k11f,qui,j,t+sk11);
				}

				for ( j = n1; j < n2; j++) {
					f[j] = f[j] + rk*(9.0*(k1[j] + k11[j]) + 49.0*(k8[j] + k10[j]) + 64.0*k9[j])/180.0;
				}



				countfn=countfn+modulo(f[data/2]); // contador para a media temporal da probabilidade de retorno a origem
				i=i+1;

				t=t+rk; // evolucao do tempo



/*=======================================================================*/
				write = write+1;//contador para escrita em arquivo
				/*printf("%f %.2f\r",t,nlparameter); //acompanhamento de tela*/

/*========PARTICIPACAO EM FUNCAO DO TEMPO==============*/
if(write%1000==0){
					for (j = n1; j < n2; j++) {
						x=x+j*modulo(f[j]);
					}
								//////////////////////// Salvando cada ponto no arquivo
							//	fprintf(fil1,"%f %f %f\n",t,x,modulo(f[data/2]));
							fprintf(fil1,"%f %f %d\n",t,x,n);
								fclose(fil1);//abrir e fechar o arquivo para cada loop afim de salvar os dados sempre
								fil1=fopen(filename1,"a");
								x=0;
							}


/*================ACRESCIMO DE SITIOS ==================*/
		if(modulo(f[n1+10])>rk||modulo(f[n2-10])>rk)
		{
			n=n+newline;
			n1=(data/2)-n;
			n2=(data/2)+n;
		}

	}while(t<tmax);



/*====================Salvando a probabilidade de retorno media.=====================*/
			/*fprintf(fil1,"%f %f\n",nlparameter,countfn/i);
			fclose(fil1);//abrir e fechar o arquivo para cada loop afim de salvar os dados sempre
			fil1=fopen(filename1,"a");*/

	nlparameter = nlparameter+dnl;	//acrescimo no parametro nao linear
	}while(nlparameter<=nlparametermax);

	fclose(fil1);

	timecount2=time(NULL);
	printf("\n O tempo de execucao foi de %.2f s\n",timecount2-timecount1);




}
