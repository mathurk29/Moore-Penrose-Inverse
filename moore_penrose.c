#include<stdio.h>
#include<stdlib.h>
#include "svd.c"	

void main()
{
    int m,n;
    printf("\n\t---Moore Penrose Inverse to Solve System of Linear Equations by:---\n");
    printf("\n\t----------------------Anish Agarwal(2012026)-----------------------\n");
    printf("\n\t----------------------Apoorva Gupta(2012040)-----------------------\n");
    printf("\n\t----------------------Hemang Dewan(2012102)------------------------\n");
    printf("\n\t------------------Prashant Shrivastava(2012174)--------------------\n");
    printf("\n\t-------------------------------------------------------------------\n");
    printf("\nEnter size of Matrix: \n");
    scanf("%d %d",&m,&n);
    float **a,**v,*w,**d,*b,**pinv,**temp_pinv,**x,temp,**ut,**vtt,**dt;
    a=malloc(m*sizeof(float*));
    b=malloc(m*sizeof(float));
    v=malloc(n*sizeof(float*));
    w=malloc(n*sizeof(float));
    d=malloc(m*sizeof(float*));
    x=malloc(n*sizeof(float*));
    ut=malloc(m*sizeof(float*));
    vtt=malloc(n*sizeof(float*));
    pinv=malloc(n*sizeof(float*));
    temp_pinv=malloc(n*sizeof(float*));
    dt=malloc(n*sizeof(float*));
    
    int i,j,k;
    for(i=0;i<m;i++)
    {
        a[i]=malloc(n*sizeof(float));

    }
    for(i=0;i<m;i++)
    {
        ut[i]=malloc(m*sizeof(float));

    }
    for(i=0;i<n;i++)
    {
        v[i]=malloc(n*sizeof(float));

    }
    for(i=0;i<n;i++)
    {
        vtt[i]=malloc(n*sizeof(float));

    }
    for(i=0;i<m;i++)
    {
        d[i]=malloc(n*sizeof(float));
    }
    for(i=0;i<m;i++)
    {
        dt[i]=malloc(m*sizeof(float));
    }
  for(i=0;i<n;i++)
    {
        x[i]=malloc(1*sizeof(float));

    }
  for(i=0;i<n;i++)
    {
        pinv[i]=malloc(m*sizeof(float));
    }
  for(i=0;i<n;i++)
    {
        temp_pinv[i]=malloc(m*sizeof(float));
    }

    printf("Equation form:AX = B\n");
    printf("Enter Matrix A %d X %d\n",m,n);
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            scanf("%f",&a[i][j]);
        }
    }
	
	
    printf("Enter Matrix B: \n");
    for(i=0;i<m;i++)
        {
            scanf("%f",&b[i]);
        }

    printf("\n A is\n");
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%f ",a[i][j]);
        }
        printf("\n");
    }
    dsvd(a,m,n,w,v);
    printf("\n Singular Vector is\n");
    for(i=0;i<n;i++)
    {
        printf("%f ",w[i]);
    }
    printf("\n");
    
    printf("\n U is\n");

    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            printf("%f ",a[i][j]);
        }
        printf("\n");
    }

    printf("\n V is\n");

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%f ",v[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for(i=0;i<m;i++){
      for (j=0;j<n;j++){
        if (i==j){
    	    d[i][j]=w[i];
        }
        else d[i][j] = 0;
        }
    }

    //Reciprocal of Sigma

    for(i=0;i<m;i++){
      for(j=0;j<n;j++){
        if(d[i][j]!=0){
            d[i][j]=(1/d[i][j]);
        }
      }
    }
    
    //Transpose Of Sigma
    
    for(i=0;i<n;i++){
       for(j=0;j<m;j++){
          dt[i][j]=d[j][i];
       }
    }

    //U transpose
    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            ut[i][j]=a[j][i];
        }
    }

    // V 
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            vtt[i][j]=v[j][i];
        }
    }
    
    printf("\n Diagonal Matrix is \n");
    //Diagonal Matrix Formed.
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            printf("%f ",dt[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    // Now A-pseudoinverse = V.(diagonal matrix).(U Transpose)
     for (i=0;i<n;i++){
        for(j=0;j<n;j++){
            for(k=0;k<n;k++){
                temp_pinv[i][j]+=vtt[i][k]*dt[k][j];
            }
        }
    }

     for (i=0;i<n;i++){
        for(j=0;j<m;j++){
            for(k=0;k<m;k++){
                pinv[i][j]+=temp_pinv[i][k]*ut[k][j];
            }
        }
    }
    
    //pinv formed of size n X m
    // Our Equation was Ax=b,now x=pinv.b so we take input b matrix also.(constant's matrix) of size m X 1
    //final answer in x matrix of size n X 1

    for (i=0;i<n;i++){
        for(j=0;j<1;j++){
            for(k=0;k<m;k++){
                x[i][j]+=pinv[i][k]*b[k];
            }
        }
    }
    
    printf("\n Solution of System of Linear Equations is: \n");
        for(i=0;i<n;i++)
        {
            for(j=0;j<1;j++)
            {
                printf("%f ",x[i][j]);
            }
            printf("\n");
        }
}

