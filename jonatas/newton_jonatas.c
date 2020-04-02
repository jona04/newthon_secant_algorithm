#include <math.h>
#include <stdio.h>
#include <stdlib.h>


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

void initial_lambda(cqk_problem *restrict p,double *lambda,double *slopes,unsigned *ind)
{
    double s0 = 0.0;
    double q0 = 0.0;
    for (int i = 0; i < p->n;i++)
    {
        slopes[i] = (p->b[i]/p->d[i])*p->b[i];
        ind[i] = i;
        s0 = s0+ (p->a[i] * p->b[i])/p->d[i];
        q0 = q0 + (p->b[i] * p->b[i]) / p->d[i];
    }
    
    *lambda = (p->r-s0)/q0;
}

void phi_lambda(cqk_problem *restrict p,double *restrict lambda,double *restrict phi,double *restrict x,double *restrict deriv,double *slopes,double *restrict r,unsigned *ind,unsigned n)
{
    *deriv = 0.0;
    // *phi = 0.0;
    *phi = *r;
    *phi = *phi * -1;
    // double soma_phi = *phi;
    // double soma_deriv = 0.0;

    
    for (int i =0;i<n;i++)
    {

        unsigned ii = ind[i];


        x[ii] = (p->b[ii] * (*lambda) + p->a[ii])/p->d[ii];

        if (x[ii] < p->low[ii])
            x[ii] = p->low[ii];
        else if (x[ii] > p->up[ii])
            x[ii] = p->up[ii];
        else
            *deriv += slopes[ii];
        
        *phi += p->b[ii] * x[ii];
    
    }
    
    // *deriv = soma_deriv;
    // *phi = soma_phi;
    // printf("*deriv %f \n",*deriv);
    // printf("valor lambda no phi = %f \n",*lambda);
    // printf("valor de phi = %f \n",*phi);
    // printf("r = %f \n",*r);

}

void fix_variable_low(cqk_problem *restrict p,double *x,double *restrict r,unsigned *ind,unsigned *restrict n)
{
    double soma_r = *r;
    unsigned len = 0;

    for (int i =0;i<*n;i++)
    {
        unsigned ii = ind[i];
        if (x[ii] <= p->low[ii]){
            *r -= x[ii]*p->b[ii];
            // fix++;
        }else{
            ind[len] = ii;
            // printf("ind[len] %d %d\n",len,ii);
            ++len;
        }
    
    }
    // printf("soma_r %f \n",soma_r);
    // *r = soma_r;
    *n = len;
}

void fix_variable_up(cqk_problem *restrict p,double *x,double *restrict r,unsigned *ind,unsigned *restrict n)
{
    double soma_r = *r;
    unsigned len = 0;
    for (int i =0;i<*n;i++)
    {
        unsigned ii = ind[i];
        if (x[ii] >= p->up[ii]){
            *r -= x[ii]*p->b[ii];
        }else {
            ind[len] = ii;
            ++len;
        }
        
    }
    // printf("len %d \n",len);
    // *r = soma_r;
    *n = len;
}

double breakpoint_to_the_right(cqk_problem *restrict p,double *lambda,unsigned *restrict ind, unsigned n)
{
    double next_break = INFINITO_POSITIVO;
            
    for (int i = 0; i < n;i++)
    {
        unsigned ii = ind[i];
        double pos_break = (p->d[ii]*p->low[ii] - p->a[ii])/p->b[ii];
        if (pos_break > *lambda && pos_break < INFINITO_POSITIVO)
            next_break = pos_break;
        
    }
    return next_break;
}

double breakpoint_to_the_left(cqk_problem *restrict p,double *lambda,unsigned *restrict ind, unsigned n)
{
    double next_break = INFINITO_NEGATIVO;
            
    for (int i = 0; i < n;i++)
    {
        unsigned ii = ind[i];
        double neg_break = (p->d[ii]*p->up[ii] - p->a[ii])/p->b[ii];
        if (neg_break < *lambda && neg_break > INFINITO_NEGATIVO)
            next_break = neg_break;
        
    }
    return next_break;
}

double secant(cqk_problem *restrict p, double *x,double *restrict alfa,double *restrict beta,double *restrict phi_alfa, double *restrict phi_beta,double *restrict r,unsigned *restrict ind,unsigned n)
{

    double lambda_secant = 0.0;

    // lambda_secant = *alfa - *phi_alfa * ((*beta - *alfa)/(*phi_beta - *phi_alfa));

    // if (lambda_secant == *alfa || lambda_secant == *beta)
        // lambda_secant = 0.5*(*alfa + *beta);
    lambda_secant = 0.5*(*alfa + *beta);
    double soma_phi = 0;
    for (int i =0;i<n;i++)
    {
        unsigned ii = ind[i];
        x[ii] = (p->b[ii] * (lambda_secant) + p->a[ii])/p->d[ii];
        if (x[ii] < p->low[ii]){
            x[ii] = p->low[ii];
        }else if (x[ii] > p->up[ii])
            {
                x[ii] = p->up[ii];
            }
        soma_phi = soma_phi + p->b[ii] * x[ii];
        
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
    unsigned n = p->n;
    double phi = 0.0;
    double lambda = 0.0;
    double alfa = INFINITO_NEGATIVO;
    double beta = INFINITO_POSITIVO;
    double phi_alfa = 0.0;
    double phi_beta = 0.0;
    double deriv;         /* Derivative of phi */
    double r = p->r;
    double *slopes = (double *) malloc(n*sizeof(double));
    unsigned *restrict ind = (unsigned *) malloc(n*sizeof(unsigned));

    initial_lambda(p,&lambda,slopes,ind);
    // printf("lambda init %f \n",lambda);

    phi_lambda(p,&lambda,&phi,x,&deriv,slopes,&r,ind,n);
    
     printf("Initial Lambda=%e, phi - r = %e r= %f\n", lambda, phi, p->r);

    int it = 1;
    // for (it = 1; it < MAX_IT;it++){
    while (phi != 0.0 && it <= MAX_IT) {
        // printf("steps %d \n",it);

        if (phi > 0)
        {

            // printf("positivo \n");
            beta = lambda;
            phi_beta = phi;
            double lambda_n = 0.0;

            if (deriv > 0.0)
            {
                lambda_n = lambda - ( (phi)/deriv);
                // printf("lambda_n = %f\n",lambda_n);
                if (fabs(lambda_n - lambda) <= 0.0000000001) {
                    phi = 0.0;
                    break;
                }

                if (lambda_n > alfa)
                {

                    lambda = lambda_n;
                    fix_variable_up(p,x,&r,ind,&n);
                    // printf("lambda = %f \n",lambda);
                    
                }else
                {
                    lambda = secant(p,x,&alfa,&beta,&phi_alfa,&phi_beta,&r,ind,n);
                    printf("secant 2\n");
                }

            }
            if (deriv == 0.0){

                lambda = breakpoint_to_the_left(p,&lambda,ind,n);
                fix_variable_up(p,x,&r,ind,&n);
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
            phi_alfa = phi;

            double lambda_n = 0.0;

            if (deriv > 0.0)
            {
                lambda_n = lambda - ( (phi)/deriv);
                if (fabs(lambda_n - lambda) <= 0.00000000001) {
                    phi = 0.0;
                    break;
                }
                // printf("lambda_n positivo = %f \n",lambda_n);
                // printf("beta positivo = %f \n",beta);
                if( lambda_n < beta)
                {
                    lambda = lambda_n;
                    fix_variable_low(p,x,&r,ind,&n);
                    // printf("lambda = %f \n",lambda);
                }else
                {
                    lambda = secant(p,x,&alfa,&beta,&phi_alfa,&phi_beta,&r,ind,n);
                    printf("secant 1\n");
                } 
            }
            if (deriv == 0.0)
            {
                lambda = breakpoint_to_the_right(p,&lambda,ind,n);
                fix_variable_low(p,x,&r,ind,&n);
                printf("break 2 \n");
                /* Test for infeasibility */
                if (lambda <= INFINITO_NEGATIVO || 
                    lambda >= INFINITO_POSITIVO) {
                    // interval.pos_lambda = interval.neg_lambda = lambda;
                    break;
                }

            }

        }

        phi_lambda(p,&lambda,&phi,x,&deriv,slopes,&r,ind,n);  
        printf("phi = %f \n",phi);
        if (phi < 0.0000001 && phi > -0.0000001)
        {
            phi = 0.0;
            break;
        }
        it++;
    }

    
    // free(slopes);
    // free(Ik);
    imprime_resultado_x(p,x);

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