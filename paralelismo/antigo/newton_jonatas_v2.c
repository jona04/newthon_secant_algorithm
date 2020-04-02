#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "cont_quad_knapsack.h"

#define INFINITO_POSITIVO 999999
#define INFINITO_NEGATIVO -999999
#define MAX_IT 10

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

void allocate_cqk_problem(unsigned n, cqk_problem *restrict p) {
    p->n = n;
    p->d = (double *) malloc(p->n*sizeof(double));
    p->a = (double *) malloc(p->n*sizeof(double));
    p->b = (double *) malloc(p->n*sizeof(double));
    p->low = (double *) malloc(p->n*sizeof(double));
    p->up = (double *) malloc(p->n*sizeof(double));
    if (!p->d || !p->a || !p->b || !p->low || !p->up) {
        fprintf(stderr, "Memory allocation error, line %d, file %s\n",
                __LINE__, __FILE__);
        exit(1);
    }
}

void free_cqk_problem(cqk_problem *restrict p) {
    p->n = 0;
    free(p->d);
    p->d = NULL;
    free(p->a);
    p->a = NULL;
    free(p->b);
    p->b = NULL;
    free(p->low);
    p->low = NULL;
    free(p->up);
    p->up = NULL;
}

void imprime_resultado_x(cqk_problem *restrict p,double *x){
    double total = 0;
    for (int i = 0; i < p->n;i++){
     	total = total + x[i]*p->b[i];
        //  printf("x%d = %f  l = %f   u = %f  \n",i,x[i],p->low[i],p->up[i]);
    }
    printf("r %f\n",p->r);
    printf("total %f \n",total);
}

void initial_lambda(cqk_problem *restrict p,double *lambda,double *slopes)
{
    double s0 = 0.0;
    double q0 = 0.0;
    #pragma omp parallel
    {   
        #pragma omp for reduction(+:s0,q0)
        for (int i = 0; i < p->n;i++)
        {
            slopes[i] = (p->b[i]/p->d[i])*p->b[i];

            s0 = s0 + (p->a[i] * p->b[i])/p->d[i];
            q0 = q0 + (p->b[i] * p->b[i]) / p->d[i];
        }
    }
    
    *lambda = (p->r-s0)/q0;
}

void phi_lambda(cqk_problem *restrict p,double *restrict lambda,double *restrict phi,double *restrict x,int *Ik,double *restrict deriv,double *slopes,double *restrict r)
{
    *deriv = 0.0;
    double soma_deriv = 0.0;
    *phi = *r;
    *phi = *phi * -1;
    double soma_phi = *phi;
    // *phi = 0.0;
    // double soma_phi = *phi;

    #pragma omp parallel
    {   
        #pragma omp for reduction(+:soma_deriv,soma_phi)
        for (int i =0;i<p->n;i++)
        {
            if (Ik[i] != -1){
                x[i] = (p->b[i] * (*lambda) + p->a[i])/p->d[i];

                if (x[i] < p->low[i])
                    x[i] = p->low[i];
                else if (x[i] > p->up[i])
                    x[i] = p->up[i];
                else
                    soma_deriv = soma_deriv + slopes[i];
                soma_phi = soma_phi + p->b[i] * x[i];
            }
        }
    }

    *deriv = soma_deriv;
    *phi = soma_phi;
    // printf("valor lambda no phi = %f \n",*lambda);
    // printf("valor de phi = %f \n",*phi);
    // printf("r = %f \n",*r);

}

void fix_variable_low(cqk_problem *restrict p,double *x,int *Ik,int pos_neg,double *restrict r)
{
    double soma_r = *r;
    #pragma omp parallel
    {   
        #pragma omp for reduction(+:soma_r)
        for (int i =0;i<p->n;i++)
        {
            if (Ik[i] != -1){
                
                    if (x[i] <= p->low[i]){
                        Ik[i] = -1;
                        soma_r = soma_r - x[i]*p->b[i];
                        // fix++;
                    }
                
            }
        }
    }
    *r = soma_r;
    // printf("fix = %d \n",fix);
}

void fix_variable_up(cqk_problem *restrict p,double *x,int *Ik,int pos_neg,double *restrict r)
{
    double soma_r = *r;
    #pragma omp parallel
    {   
        #pragma omp for reduction(+:soma_r)
        for (int i =0;i<p->n;i++)
        {
            if (Ik[i] != -1){
                
                if (x[i] >= p->up[i]){
                    Ik[i] = -1;
                    soma_r = soma_r - x[i]*p->b[i];
                    // fix++;
                }
                
            }
        }
    }
    *r = soma_r;
    // printf("fix = %d \n",fix);
}

double breakpoint_to_the_right(cqk_problem *restrict p,double *lambda,int *Ik)
{
    double next_break = INFINITO_POSITIVO;
    
    for (int i = 0; i < p->n;i++)
    {
        if (Ik[i] != -1){
            double pos_break = (p->d[i]*p->low[i] - p->a[i])/p->b[i];
            if (pos_break > *lambda && pos_break < INFINITO_POSITIVO)
                next_break = pos_break;
        }
    }
    return next_break;
}

double breakpoint_to_the_left(cqk_problem *restrict p,double *lambda,int *Ik)
{
    double next_break = INFINITO_NEGATIVO;
            
    for (int i = 0; i < p->n;i++)
    {
        if (Ik[i] != -1){
            double neg_break = (p->d[i]*p->up[i] - p->a[i])/p->b[i];
            if (neg_break < *lambda && neg_break > INFINITO_NEGATIVO)
                next_break = neg_break;
        }
    }
    return next_break;
}

double secant(cqk_problem *restrict p, double *x,int *Ik,double *restrict alfa,double *restrict beta,double *restrict r,double *restrict phi_alfa, double *restrict phi_beta)
{
    double lambda_secant = 0.0;

    lambda_secant = *alfa - *phi_alfa * ((*beta - *alfa)/(*phi_beta - *phi_alfa));

    if (lambda_secant == *alfa || lambda_secant == *beta)
        lambda_secant = 0.5*(*alfa + *beta);
    
    double soma_phi = 0;

    #pragma omp parallel
    {   
        #pragma omp for reduction(+:soma_phi)
        for (int i =0;i<p->n;i++)
        {
            if (Ik[i] != -1){
                x[i] = (p->b[i] * (lambda_secant) + p->a[i])/p->d[i];
                if (x[i] < p->low[i]){
                    x[i] = p->low[i];
                }else if (x[i] > p->up[i])
                    {
                        x[i] = p->up[i];
                    }
                soma_phi = soma_phi + p->b[i] * x[i];
            }
        }
    }
    if (soma_phi < *r) 
    {
        return MAX(lambda_secant,*alfa);
    }else{
        return MIN(lambda_secant,*beta);
    }

}

int newton_jonatas(cqk_problem *restrict p, double *x)
{
    double phi = 0.0;
    double lambda = 0.0;
    double alfa = INFINITO_NEGATIVO;
    double beta = INFINITO_POSITIVO;
    double phi_alfa = 0.0;
    double phi_beta = 0.0;
    double deriv;         /* Derivative of phi */
    double r = p->r;

    double *slopes = (double *) malloc(p->n*sizeof(double));
    int *Ik = (int *) malloc(p->n*sizeof(int));

    initial_lambda(p,&lambda,slopes);
    
    phi_lambda(p,&lambda,&phi,x,Ik,&deriv,slopes,&r);
    
    int it;
    for (it = 1; it < MAX_IT;it++){

        if (phi > 0)
        {
            alfa = lambda;
            phi_alfa = phi;
            double lambda_n = 0.0;

            if (deriv > 0.0)
            {
                lambda_n = lambda - ( (phi)/deriv);
                
                if (fabs(lambda_n - lambda) <= 0.0000000001) {
                    phi = 0.0;
                    break;
                }

                if( lambda_n < beta)
                {
                    lambda = lambda_n;
                    fix_variable_low(p,x,Ik,1,&r);
                    
                }else
                {
                    lambda = secant(p,x,Ik,&alfa,&beta,&phi_alfa,&phi_beta,&r);
                    printf("seacnt \n");
                } 

            }
            if (deriv == 0.0)
            {
                lambda = breakpoint_to_the_right(p,&lambda,Ik);
                fix_variable_low(p,x,Ik,1,&r); //1 de positivo

                /* Test for infeasibility */
                if (lambda <= INFINITO_NEGATIVO || 
                    lambda >= INFINITO_POSITIVO) {
                    // interval.pos_lambda = interval.neg_lambda = lambda;
                    break;
                }
            }
        }else{
            beta = lambda;
            phi_beta = phi;
            // printf("negativo \n");
            double lambda_n = 0.0;

            if (deriv > 0.0)
            {
                lambda_n = lambda - ( (phi)/deriv);

                if (fabs(lambda_n - lambda) <= 0.0000000001) {
                    phi = 0.0;
                    break;
                }
                if (lambda_n > alfa)
                {
                    lambda = lambda_n;
                    fix_variable_up(p,x,Ik,0,&r);
                    // printf("deriv 22 \n");
                }else
                {
                    lambda = secant(p,x,Ik,&alfa,&beta,&phi_alfa,&phi_beta,&r);
                    printf("seacnt \n");
                }

                

            }
            if (deriv == 0.0){

                lambda = breakpoint_to_the_left(p,&lambda,Ik);
                fix_variable_up(p,x,Ik,0,&r); //1 de positivo

                /* Test for infeasibility */
                if (lambda <= INFINITO_NEGATIVO || 
                    lambda >= INFINITO_POSITIVO) {
                    // interval.pos_lambda = interval.neg_lambda = lambda;
                    break;
                }

            }
        }

        phi_lambda(p,&lambda,&phi,x,Ik,&deriv,slopes,&r);  
        
        // if (fabs(phi - r) <= 0.00001)
        // {   
        //     phi = 0.0;
        //     break;
        // }

        if (phi < 0.0000001 && phi > -0.0000001)
        {
            phi = 0.0;
            break;
        }

    }


    // imprime_resultado_x(p,x);

    // free(slopes);
    // free(Ik);
    if (phi == 0.0)
    {
        return it;
    }else if (alfa == beta){
        printf("problema sem solucao");
        return -1;
    }else
    {
        printf("problema no metodo");
        return -2;
    }
}