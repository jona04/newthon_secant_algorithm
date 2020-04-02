#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "cont_quad_knapsack.h"

#define INFINITO_POSITIVO 999999
#define INFINITO_NEGATIVO -999999
#define MAX_IT 100

#define NUM_THREADS 1

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
    

    double s0;
    double q0;
    // double teste = 0,testelow = 0,testelow2 = 0,testeup = 0,testeup2 = 0;

    // #pragma omp parallel
    // {   
        // #pragma omp for reduction(+:s0,q0)
        for (int i = 0; i < p->n;i++)
        {
            slopes[i] = (p->b[i]/p->d[i])*p->b[i];

            s0 = s0 + (p->a[i] * p->b[i])/p->d[i];
            q0 = q0 + (p->b[i] * p->b[i]) / p->d[i];

            // teste += ((((p->low[i] * p->d[i])-p->a[i] )/p->b[i]) + 
            //         (((p->up[i] * p->d[i])-p->a[i] )/p->b[i]))
            //         /2;
            // testelow2 += ((p->low[i] * p->d[i])-p->a[i] )/p->b[i];
            // testeup2 += ((p->up[i] * p->d[i])-p->a[i] )/p->b[i];
        }
    
    // }
    *lambda = (p->r-s0)/q0;
}

void phi_lambda(cqk_problem *restrict p,double lambda,double *restrict phi,double *restrict x,double *restrict deriv,double *slopes,double r)
{
    *deriv = 0.0;
    // double soma_deriv = 0.0;
    *phi = r * -1;
    // double soma_phi = *phi;
    
    // int qtd_thread = omp_get_max_threads();
    // int porc = p->n / qtd_thread;
    
    

    #pragma omp parallel
    { 
        int i,id,nthrds;
        id = omp_get_thread_num();
        // printf("id = %d \n", id);
        double soma_d=0.0;
        double soma_p=0.0;
        
        omp_set_num_threads(NUM_THREADS);

        nthrds = omp_get_num_threads();
        for (i=id; i< p->n; i=i+nthrds)
        // for (i=id*5; i< p->n/2;i++)  
        {
            // printf("id = %d / i = %d \n",id,i);
            x[i] = (p->b[i] * lambda + p->a[i])/p->d[i];

            if (x[i] < p->low[i])
                x[i] = p->low[i];
            else if (x[i] > p->up[i])
                x[i] = p->up[i];
            else
                soma_d = soma_d + slopes[i];
            soma_p = soma_p + p->b[i] * x[i];

            
        }
        #pragma omp atomic
        *deriv += soma_d;
        *phi += soma_p;


        // #pragma omp for reduction(+:soma_deriv,soma_phi)
        // for (int i =0;i<p->n;i++)
        // {
        //     printf("id = %d / i = %d \n",id,i);
        //     x[i] = (p->b[i] * (*lambda) + p->a[i])/p->d[i];

        //     if (x[i] < p->low[i])
        //         x[i] = p->low[i];
        //     else if (x[i] > p->up[i])
        //         x[i] = p->up[i];
        //     else
        //         soma_deriv = soma_deriv + slopes[i];
        //     soma_phi = soma_phi + p->b[i] * x[i];

        // }

    }

    // *deriv = soma_deriv;
    // *phi = soma_phi;

}

double breakpoint_to_the_right(cqk_problem *restrict p,double *lambda)
{
    double next_break = INFINITO_POSITIVO;
    
    for (int i = 0; i < p->n;i++)
    {
        double pos_break = (p->d[i]*p->low[i] - p->a[i])/p->b[i];
        if (pos_break > *lambda && pos_break < INFINITO_POSITIVO)
            next_break = pos_break;
        
    }
    return next_break;
}

double breakpoint_to_the_left(cqk_problem *restrict p,double *lambda)
{
    double next_break = INFINITO_NEGATIVO;
            
    for (int i = 0; i < p->n;i++)
    {

        double neg_break = (p->d[i]*p->up[i] - p->a[i])/p->b[i];
        if (neg_break < *lambda && neg_break > INFINITO_NEGATIVO)
            next_break = neg_break;

    }
    return next_break;
}

double secant(cqk_problem *restrict p, double *x,double *restrict alfa,double *restrict beta,double *restrict r,double *restrict phi_alfa, double *restrict phi_beta)
{
    double lambda_secant = 0.0;

    lambda_secant = *alfa - *phi_alfa * ((*beta - *alfa)/(*phi_beta - *phi_alfa));

    if (lambda_secant == *alfa || lambda_secant == *beta)
        lambda_secant = 0.5*(*alfa + *beta);
    
    double soma_phi = 0;

    // #pragma omp parallel
    // {   
        // #pragma omp for reduction(+:soma_phi)
        for (int i =0;i<p->n;i++)
        {
            x[i] = (p->b[i] * (lambda_secant) + p->a[i])/p->d[i];
            if (x[i] < p->low[i]){
                x[i] = p->low[i];
            }else if (x[i] > p->up[i])
                {
                    x[i] = p->up[i];
                }
            soma_phi = soma_phi + p->b[i] * x[i];

        }
    // }
    if (soma_phi < *r) 
    {
        return MAX(lambda_secant,*alfa);
    }else{
        return MIN(lambda_secant,*beta);
    }

}

int newton_jonatas(cqk_problem *restrict p, double *x)
{
    double phi;
    double lambda;
    double alfa = INFINITO_NEGATIVO;
    double beta = INFINITO_POSITIVO;
    double phi_alfa = 0.0;
    double phi_beta = 0.0;
    double deriv;         /* Derivative of phi */
    double r = p->r;

    double *slopes = (double *) malloc(p->n*sizeof(double));

    initial_lambda(p,&lambda,slopes);
    
    phi_lambda(p,lambda,&phi,x,&deriv,slopes,r);

    int it = 1;
    while (phi != 0.0 && it <= MAX_IT) {
        // printf("steps %d \n",it);
        if (phi > 0)
        {
            
            beta = lambda;
            
            // printf("positivo \n");
            double lambda_n = 0.0;

            if (deriv > 0.0)
            {
                lambda_n = lambda - ( phi/deriv);

                if (fabs(lambda_n - lambda) <= 0.00000000001) {
                    phi = 0.0;
                    break;
                }
                if (lambda_n > alfa)
                {
                    lambda = lambda_n;
                    // printf("lambda = %f\n",lambda);
                }else
                {
                    phi_beta = phi;
                    lambda = secant(p,x,&alfa,&beta,&phi_alfa,&phi_beta,&r);
                    printf("secant 2 \n");
                }

                

            }
            if (deriv == 0.0){
                
                lambda = breakpoint_to_the_left(p,&lambda);
                printf("break 2 \n");
                /* Test for infeasibility */
                if (lambda <= INFINITO_NEGATIVO || 
                    lambda >= INFINITO_POSITIVO) {
                    // interval.pos_lambda = interval.neg_lambda = lambda;
                    break;
                }

            }

            
        }else{


            // printf("negativo \n");
            alfa = lambda;
            
            double lambda_n = 0.0;

            if (deriv > 0.0)
            {
                lambda_n = lambda - ( phi/deriv);
                
                if (fabs(lambda_n - lambda) <= 0.00000000001) {
                    phi = 0.0;
                    break;
                }

                if( lambda_n < beta)
                {
                    lambda = lambda_n;
                    // printf("lambda = %f\n",lambda);
                }else
                {
                    phi_alfa = phi;
                    lambda = secant(p,x,&alfa,&beta,&phi_alfa,&phi_beta,&r);
                    printf("seacnt 1 \n");
                } 

            }
            if (deriv == 0.0)
            {
                lambda = breakpoint_to_the_right(p,&lambda);

                printf("break 1 \n");
                /* Test for infeasibility */
                if (lambda <= INFINITO_NEGATIVO || 
                    lambda >= INFINITO_POSITIVO) {
                    // interval.pos_lambda = interval.neg_lambda = lambda;
                    break;
                }
            }
            
        }

        phi_lambda(p,lambda,&phi,x,&deriv,slopes,r);

        it++;
    }

    // printf("tempo 2 = %f\n",t_phi);
    // imprime_resultado_x(p,x);

    free(slopes);
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