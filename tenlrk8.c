/*Transmissao de Eletrons em cadeias unidimensionais com não linearidade de
Hoinstein e campos externos
Autor: Frederico Passos
Orientador: Wandearley Dias
Versao = 6
n 15/11/19*/

/*===========BIBLIOTECAS=======*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>


#include "dados.h" //DADOS INICIAIS




/*============/FUNCOES\================*/

/*=============MODULO QUADRADO============*/
double modulo(double complex psin)
{
  double MODULO;
  MODULO = cabs(psin)*cabs(psin);
  return MODULO;
}

/*=============EQ DE SCHRODINGER NL COM CAMPO============*/

double complex dnlse (double complex fn,double complex fa,double complex fp, double chi,int j, int n, double t) //aplicacao da DNLS
{
	double complex func;
	double complex scho = -1i;
	double E,j0=n/2,feq=omega+domega;
//  double Fw=0, omega=0, domega=0, phase=0;
	E=F0+Fw*sin(feq*t+phase);//CAMPO EXTERNO
	func = scho*(fa+fp-fn*chi*(modulo(fn))+E*(j-j0)*fn);//DNLSE
	return func;
}

/*=============METODO RK4============*/

double complex rungekutta4(_Complex double *c,double *chi, int n, double t)
{
  double complex k1[n],k2[n],k3[n],k4[n];
  int j;
  for (j = 1; j < n; j++) {k1[j]=rk*dnlse(c[j],c[j-1],c[j+1],chi[j],j,n,t);}//printf("%d %f %f\n",j,n,creal(k1[j]),cimag(k1[j]));}
  for (j = 1; j < n; j++) {k2[j]=rk*dnlse(c[j]+k1[j]/2,c[j-1]+k1[j-1]/2,c[j+1]+k1[j+1]/2,chi[j],j,n,t+rk/2);}
  for (j = 1; j < n; j++) {k3[j]=rk*dnlse(c[j]+k2[j]/2,c[j-1]+k2[j-1]/2,c[j+1]+k2[j+1]/2,chi[j],j,n,t+rk/2);}
  for (j = 1; j < n; j++) {k4[j]=rk*dnlse(c[j]+k3[j],c[j-1]+k3[j-1],c[j+1]+k3[j+1],chi[j],j,n,t+rk);}
  for (j = 0; j < n; j++) {c[j]=c[j]+(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;}
  }

/*=============METODO RK8============*/

double complex rungekutta8(_Complex double *c,double *chi, int n, double t)
  {
    double complex k1[n],k2[n],k3[n],k4[n],k5[n],k6[n],k7[n],k8[n],k9[n],k10[n],k11[n];
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
    int j,n1=1,n2=n;
    for ( j = n1; j < n2; j++) {
      k1[j]=dnlse(c[j],c[j-1],c[j+1],chi[j],j,n,t);
    }

    for ( j = n1; j < n2; j++) {
      k2[j]=dnlse(c[j]+k2a*k1[j],c[j-1]+k2a*k1[j-1],c[j+1]+k2a*k1[j+1],chi[j],j,n,t+sk2);
    }

    for ( j = n1; j < n2; j++) {
      k3[j]=dnlse(c[j]+ k3a*( k1[j] + k2[j]),c[j-1]+ k3a*( k1[j-1] + k2[j-1]),c[j+1]+ k3a*( k1[j+1] + k2[j+1]),chi[j],j,n,t+sk3);
    }

    for ( j = n1; j < n2; j++) {
      k4[j]=dnlse(c[j]+ k1[j]*k4a + k2[j]*k4b + k3[j]*k4c,c[j-1]+ k1[j-1]*k4a + k2[j-1]*k4b + k3[j-1]*k4c,c[j+1]+ k1[j+1]*k4a + k2[j+1]*k4b + k3[j+1]*k4c,chi[j],j,n,t+sk4);
    }

    for ( j = n1; j < n2; j++) {
      k5[j]=dnlse(c[j]+ k1[j]*k5a + k3[j]*k5b + k4[j]*k5c,c[j-1]+ k1[j-1]*k5a + k3[j-1]*k5b + k4[j-1]*k5c,c[j+1]+ k1[j+1]*k5a + k3[j+1]*k5b + k4[j+1]*k5c,chi[j],j,n,t+sk5);
    }

    for ( j = n1; j < n2; j++) {
      k6[j]=dnlse(c[j]+ k1[j]*k6a + k3[j]*k6b + k4[j]*k6c + k5[j]*k6d,c[j-1]+ k1[j-1]*k6a + k3[j-1]*k6b + k4[j-1]*k6c + k5[j-1]*k6d,c[j+1]+ k1[j+1]*k6a + k3[j+1]*k6b + k4[j+1]*k6c + k5[j+1]*k6d,chi[j],j,n,t+sk6);
    }

    for ( j = n1; j < n2; j++) {
      k7[j]=dnlse(c[j]+ k1[j]*k7a + k3[j]*k7b + k4[j]*k7c + k5[j]*k7d+ k6[j]*k7e,c[j-1]+ k1[j-1]*k7a + k3[j-1]*k7b + k4[j-1]*k7c + k5[j-1]*k7d+ k6[j-1]*k7e,c[j+1]+ k1[j+1]*k7a + k3[j+1]*k7b + k4[j+1]*k7c + k5[j+1]*k7d+ k6[j+1]*k7e,chi[j],j,n,t+sk7);
    }

    for ( j = n1; j < n2; j++) {
      k8[j]=dnlse(c[j]+ k1[j]*k8a + k5[j]*k8b + k6[j]*k8c + k7[j]*k8d,c[j-1]+ k1[j-1]*k8a + k5[j-1]*k8b + k6[j-1]*k8c + k7[j-1]*k8d,c[j+1]+ k1[j+1]*k8a + k5[j+1]*k8b + k6[j+1]*k8c + k7[j+1]*k8d,chi[j],j,n,t+sk8);
    }

    for ( j = n1; j < n2; j++) {
      k9[j]=dnlse(c[j]+ k1[j]*k9a + k5[j]*k9b + k6[j]*k9c + k7[j]*k9d + k8[j]*k9e,c[j-1]+ k1[j-1]*k9a + k5[j-1]*k9b + k6[j-1]*k9c + k7[j-1]*k9d + k8[j-1]*k9e,c[j+1]+ k1[j+1]*k9a + k5[j+1]*k9b + k6[j+1]*k9c + k7[j+1]*k9d + k8[j+1]*k9e,chi[j],j,n,t+sk9);
    }

    for ( j = n1; j < n2; j++) {
      k10[j]=dnlse(c[j]+ k1[j]*k10a + k5[j]*k10b + k6[j]*k10c + k7[j]*k10d +k8[j]*k10e + k9[j]*k10f,c[j-1]+ k1[j-1]*k10a + k5[j-1]*k10b + k6[j-1]*k10c + k7[j-1]*k10d +k8[j-1]*k10e + k9[j-1]*k10f,c[j+1]+ k1[j+1]*k10a + k5[j+1]*k10b + k6[j+1]*k10c + k7[j+1]*k10d +k8[j+1]*k10e + k9[j+1]*k10f,chi[j],j,n,t+sk10);
    }

    for ( j = n1; j < n2; j++) {
      k11[j]=dnlse(c[j]+ k5[j]*k11a + k6[j]*k11b + k7[j]*k11c + k8[j]*k11d +k9[j]*k11e + k10[j]*k11f,c[j-1]+ k5[j-1]*k11a + k6[j-1]*k11b + k7[j-1]*k11c + k8[j-1]*k11d +k9[j-1]*k11e + k10[j-1]*k11f,c[j+1]+ k5[j+1]*k11a + k6[j+1]*k11b + k7[j+1]*k11c + k8[j+1]*k11d +k9[j+1]*k11e + k10[j+1]*k11f,chi[j],j,n,t+sk11);
    }

    for ( j = n1; j < n2; j++) {
      c[j] = c[j] + rk*(9.0*(k1[j] + k11[j]) + 49.0*(k8[j] + k10[j]) + 64.0*k9[j])/180.0;
    }
    }

/*=============MEDIDA DO CENTROIDE============*/
double complex centroide(_Complex double *c, int n)
{
  double x=0;
  int j;
  for (j = 0; j < n; j++) {
    x+=j*modulo(c[j]);
  }
  return x;
}

/*=============MEDIDA DA PARTICIPACAO TOTAL============*/

double complex participacao(_Complex double *c, int n)
{
    double soma2, soma4;
    int j;
    for (j = 0; j < n; j++) {
      soma2+=modulo(c[j]);
      soma4+=modulo(c[j])*modulo(c[j]);
    }
    return soma2/soma4;
}

/*=============MEDIDA DA PARTICIPACAO ESQUERDA DA POSICAO INICIAL============*/

double complex participacaoL(_Complex double *c, int n)
{
    double soma2, soma4;
    int j;
    for (j = 0; j <= n/2; j++) {
      soma2+=modulo(c[j]);
      soma4+=modulo(c[j])*modulo(c[j]);
    }
    return soma2/soma4;
}

/*=============MEDIDA DA PARTICIPACAO DIREITA DA POSICAO INICIAL============*/

double complex participacaoR(_Complex double *c, int n)
{
    double soma2, soma4;
    int j;
    for (j = n/2; j <= n; j++) {
      soma2+=modulo(c[j]);
      soma4+=modulo(c[j])*modulo(c[j]);
    }
    return soma2/soma4;
}


/*=============TIPO DO PACOTE - INICIALIZACAO NO ARQUIVO EXTERNO============*/

double complex condicaoinicial(_Complex double *c, int n)
{
  int j;
	if (Sigma==0)//PACOTE TIPO DELTA
		{
			c[n/2]=1;
		}
	else //PACOTE TIPO GAUSSIANA
	{
		double mu, ASigma=0;
		for (j = 0; j < n; j++) {
			mu=j-(n/2);
			c[j]=exp(-mu*mu/Sigma);
			ASigma = ASigma+modulo(c[j]);
		}
		for (j = 0; j < n; j++) {
			c[j]=c[j]/sqrt(ASigma);
		}
	}
}

/*=========================/FUNCAO PRINCIPAL\=========================*/
int main() {
double timecount1=time(NULL),timecount2;//CONTATOR DE TEMPO DE EXECUCAO
  _Complex double *psi; //PONTEIRO DAS AMPLITUDES DE PROBABLIDADE
  double *chi; //PONTEIRO PARA O PARAMETRO NAO LINEAR
  int n=500,count; //TAMANHO DA REDE E O CONTATOR
  double complex c[n];
  double retorno; //RETORNO

  double nlparameter=nlparametermin; //INICIALIZACAO DO PARAMETRO DE LOOP DA NAO LINEARIDADE


  /*=====ALOCACAO DE MEMORIA========*/
  psi = (_Complex double*) malloc(sizeof(_Complex double)*n);
  chi = (double*) malloc(sizeof(double)*n);

/*========CRIACAO DE ARQUIVO========*/
  char filename1[250];
  FILE *fil1;
  sprintf(filename1,"RK4_centro nli=%.2f F0=%.2f Dw=%.2f S=%.d.dat",nlparametermin,F0, domega, Sigma);
  fil1=fopen(filename1,"w");

  double t=0; //PARAMETRO DE LOOP TEMPORAL
  int j; //PARAMETRO DE LOOP DE SITIOS
  double centr; //CENTROIDE

  do//LOOP DO PARAMETRO NAO LINEAR
  {
    for (j = 0; j < n; j++) {psi[j]=0;chi[j]=nlparameter;}//LIMPEZA DE MEMORIA E ALOCACAO DO PARAMETRO NAO LINEAR

    condicaoinicial(psi,n);//INICIALIZACAO DA FUNCAO DE ONDA

    do//LOOP TEMPORAL
    {
        rungekutta8(psi,chi, n, t);//METODO DE INTEGRACAO NUMERICA

/*===========ESCRITA NO ARQUIVO=============*/
        count++;
        if (count%100==0) {
          centr=centroide(psi,n);
          fprintf(fil1, "%f %f\n",t, centr);
          fclose(fil1);//abrir e fechar o arquivo para cada loop afim de salvar os dados sempre
          fil1=fopen(filename1,"a");
        }
/*===========================================*/


        retorno+=modulo(psi[n/2]); //MEDIDA DO RETORNO

        t=t+rk; //INCREMENTO TEMPORAL
      }while (t<tmax);


/*========REINICIALIZACAO DAS VARIAVEIS=========*/
      t=0; //GARANTE O LOOP TEMPORAL
      count=0;
      retorno=0;

      nlparameter=nlparameter+dl;//INCREMENTO DO PARAMETRO NAO LINEAR
  } while(nlparameter<=nlparametermax);

    fclose(fil1);//FECHAMENTO DO ARQUIVO DE SAIDA
    	timecount2=time(NULL);
    	printf("\n O tempo de execucao foi de %.2f s\n",timecount2-timecount1);//EXIBIR O TEMPO DE EXECUCAO
}
