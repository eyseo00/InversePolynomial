
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

//#define NTRU_N 509
//#define NTRU_N 677
//#define NTRU_N 701
//#define NTRU_N 821
#define NTRU_N 509
#define iteration 10000


typedef struct{
  uint16_t coeffs[NTRU_N];
} poly;


void poly_R2_inv_variable(poly *r, poly *a);
void poly_R2_inv_variable_combine(poly *r, poly *a);

int main()
{
    poly inv_variable, inv_variable_combine;
    poly gf;

    int i, k, tmp1, tmp2, avg_variable, avg_variable_combine;
    double duration;
    clock_t start, finish;

    tmp1 = 0; tmp2 = 0;
    avg_variable = 0; avg_variable_combine = 0;

    for(k=0; k< iteration; k++){
        //For binary random coefficient generation
        for(i=0; i<NTRU_N-1; i++) gf.coeffs[i] = rand()&1;
        //For ternary random coefficient generation
        //for(i=0; i<NTRU_N-1; i++)gf.coeffs[i] = rand()&3;

        gf.coeffs[NTRU_N-1] = 0;
    
        /////////////////////////////////////////////////////////////////////////////
        //Run NTRU R2_variable inverse
        start =clock();
        poly_R2_inv_variable(&inv_variable,&gf);
        finish = clock();
        duration = (double)(finish-start);
        avg_variable = avg_variable + duration;
        ////////////////////////////////////////////////////////////////////////////
        ///
        ////////////////////////////////////////////////////////////////////////////
        //Run NTRU R2_variable_combine inverse
        start =clock();
        poly_R2_inv_variable_combine(&inv_variable_combine,&gf);
        finish = clock();
        duration = (double)(finish-start);
        avg_variable_combine = avg_variable_combine + duration;
        ////////////////////////////////////////////////////////////////////////////
        ///
        ////////////////////////////////////////////////////////////////////////////
        //Comparing results
        //Checking that the result is correct
        for(i=0; i<NTRU_N-1; i++)
        {
            if(inv_variable.coeffs[i] != inv_variable_combine.coeffs[i]) tmp2 = tmp2 + 1;
        }
        if(tmp2 != 0)
        {
            printf("Alarm k: %d\n", k);
            for(i=0; i<NTRU_N; i++)printf("%d",inv_variable_combine.coeffs[i]);
            printf("\n");
        }
        tmp2 = 0;
        ////////////////////////////////////////////////////////////////////////////
    }
    avg_variable = avg_variable/iteration;
    avg_variable_combine = avg_variable_combine/iteration;
    printf("avg_variable: %d, avg_variable_combine: %d\n", avg_variable, avg_variable_combine);
    return 0;
}
    
void poly_R2_inv_variable(poly *r, poly *a)
{
    poly f, g, f_row, s_row, tmp_1,tmp_2;
    int i;
    int16_t delta, f_deg, g_deg, B_k, initial_B;
    
    delta = 0;
    B_k = 0;
    initial_B = 0;

    //Fill the inital coefficients of g(x), f_row(x), s_row(x)
    for (i = 0;i < NTRU_N;++i)
    {
      g.coeffs[i] = 1;
      f_row.coeffs[i] = 0;
      s_row.coeffs[i] = 0;
    }
    s_row.coeffs[0] = 1;
    
    //Reduce the degree of f(x) by using mod g(x)
    for (i = 0;i < NTRU_N-1;++i) f.coeffs[i] = (a->coeffs[i] ^ a->coeffs[NTRU_N-1]) & 1;
    f.coeffs[NTRU_N-1] = 0;
    
    //For the case of f(x) has zero constant term
    if(f.coeffs[0] == 0)
    {
        initial_B = 0;
        while( f.coeffs[initial_B] == 0) initial_B = initial_B + 1;
        for(i = 0;i < NTRU_N-initial_B;++i) f.coeffs[i] = f.coeffs[i+initial_B];
        for(i = 0;i < initial_B;++i) f.coeffs[NTRU_N-1-i] = 0;
        
        s_row.coeffs[0] = 0;
        s_row.coeffs[NTRU_N-initial_B] = 1;
    }
    
    //Calculate the degree of f(x) and set the initial degree of g(x)
    i = NTRU_N-1;
    while(f.coeffs[i] == 0) i = i - 1;
    f_deg = i;
    g_deg = NTRU_N-1;
    
    //Reducing the degrees while g_deg is greater than zero
    while(g_deg>0)
    {
        //delta is zero(0000) when g_deg is greater than or equal to f_deg.
        //delta is minus one(1111) when f_deg is greater than g_deg.
        delta = (g_deg - f_deg)>>15;
        
        //Calculate B_k
        i = 0;
        while(f.coeffs[i] == g.coeffs[i]) i = i + 1;
        B_k = i;
              
        //Calculate x^{-B_k}*(f(x)+g(x))
        for (i = 0;i < NTRU_N-B_k;++i)
        {
            tmp_1.coeffs[i] = f.coeffs[i+B_k] ^ g.coeffs[i+B_k];
            tmp_2.coeffs[i] = f_row.coeffs[i+B_k]^s_row.coeffs[i+B_k];
        }
        for (i = 0;i < B_k;++i)
        {
            tmp_1.coeffs[NTRU_N-1-i] = 0;
            tmp_2.coeffs[NTRU_N-B_k+i] = f_row.coeffs[i]^s_row.coeffs[i];
        }
        
        for (i = 0;i < NTRU_N;++i)
        {
            f.coeffs[i] = ((~delta)&f.coeffs[i])^(delta&g.coeffs[i]);
            s_row.coeffs[i] = ((~delta)&s_row.coeffs[i])^(delta&f_row.coeffs[i]);
        }
        g = tmp_1;
        f_row = tmp_2;
        
        //Update g_deg and f_deg.
        //The order of i and f_deg should not be changed.
        i = ((~delta)&g_deg) + (delta&f_deg) - B_k;
        f_deg = ((~delta)&f_deg) + (delta&g_deg);
        
        while(g.coeffs[i] == 0) i = i - 1;
        g_deg = i;
        
        //Initialize B_k
        B_k = 0;
    }
    
    //Save the final inverse result into poly r
    for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] = (f_row.coeffs[i] ^ f_row.coeffs[NTRU_N-1]) & 1;
    r->coeffs[NTRU_N-1] = 0;
}

void poly_R2_inv_variable_combine(poly *r, poly *a)
{
    poly f, g, f_row, s_row, tmp_1, tmp_2, tmp_3, tmp_4;
    int16_t i, j, k, tmp1, delta1, delta2, f_deg, g_deg, diff_deg, B_k1, B_k2, initial_B, sum;
    
    delta1 = 0; delta2 = 0;
    B_k1 = 0; B_k2 = 0;
    initial_B = 0;

    //Fill the inital coefficients of g(x), f_row(x), s_row(x)
    for (i = 0;i < NTRU_N;++i)
    {
      g.coeffs[i] = 1;
      f_row.coeffs[i] = 0;
      s_row.coeffs[i] = 0;
    }
    s_row.coeffs[0] = 1;
    
    //Reduce the degree of f(x) by using mod g(x)
    for (i = 0;i < NTRU_N-1;++i) f.coeffs[i] = (a->coeffs[i] ^ a->coeffs[NTRU_N-1]);
    f.coeffs[NTRU_N-1] = 0;
    
    //For the case of f(x) has zero constant term
    if(f.coeffs[0] == 0)
    {
        initial_B = 0;
        while( f.coeffs[initial_B] == 0) initial_B = initial_B + 1;
        for(i = 0;i < NTRU_N-initial_B;++i) f.coeffs[i] = f.coeffs[i+initial_B];
        for(i = 0;i < initial_B;++i) f.coeffs[NTRU_N-1-i] = 0;
        
        s_row.coeffs[0] = 0;
        s_row.coeffs[NTRU_N-initial_B] = 1;
    }
    
    //Calculate the degree of f(x) and set the initial degree of g(x)
    i = NTRU_N-1;
    while(f.coeffs[i] == 0) i = i - 1;
    f_deg = i;
    g_deg = NTRU_N-1;
    
    while(g_deg>0 && f_deg>0)
    {
        //delta1 is zero(0000) when g_deg is greater than or equal to f_deg.
        //delta1 is minus one(1111) when f_deg is greater than g_deg.
        diff_deg = g_deg - f_deg;
        delta1 = diff_deg>>15 ;
        
        //Calculate B_k1
        i = 0;
        while(f.coeffs[i] == g.coeffs[i]) i = i + 1;
        B_k1 = i;
        
        //Calculate B_k2
        if(delta1 == 0)
        {
            if(diff_deg >= B_k1) delta2 = 0;
            else delta2 = -1;

            i = 0;
            while((f.coeffs[i+B_k1]^g.coeffs[i+B_k1]^f.coeffs[i]) == 0) i = i + 1;
            B_k2 = i;
        }
        else
        {
            if(-diff_deg >= B_k1) delta2 = 0;
            else delta2 = -1;

            i = 0;
            while((f.coeffs[i+B_k1]^g.coeffs[i+B_k1]^g.coeffs[i]) == 0) i = i + 1;
            B_k2 = i;
        }
        
       sum = B_k1 + B_k2;
        
       if( delta1 == -1 && delta2 == -1)
       {
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = g.coeffs[i+B_k2]^g.coeffs[i+sum]^f.coeffs[i+sum];
               tmp_2.coeffs[i] = f_row.coeffs[i+B_k2]^f_row.coeffs[i+sum]^s_row.coeffs[i+sum];
               tmp_3.coeffs[i] = g.coeffs[i+B_k1]^f.coeffs[i+B_k1];
               tmp_4.coeffs[i] = f_row.coeffs[i+B_k1]^s_row.coeffs[i+B_k1];
            }
            for (j = 0;j < B_k1;++j)
            {
                tmp_1.coeffs[NTRU_N-sum+j] = g.coeffs[NTRU_N-B_k1+j];
                tmp_2.coeffs[NTRU_N-sum+j] = f_row.coeffs[NTRU_N-B_k1+j]^f_row.coeffs[j]^s_row.coeffs[j];
                tmp_3.coeffs[NTRU_N-1-j] = 0;
                tmp_4.coeffs[NTRU_N-B_k1+j] = f_row.coeffs[j]^s_row.coeffs[j];
            }
            for (k = 0;k < B_k2;++k)
            {
                tmp_1.coeffs[NTRU_N-1-k] = 0;
                tmp_2.coeffs[NTRU_N-B_k2+k] = f_row.coeffs[k]^f_row.coeffs[B_k1+k]^s_row.coeffs[B_k1+k];
                tmp_3.coeffs[NTRU_N-sum+k] = g.coeffs[NTRU_N-B_k2+k]^f.coeffs[NTRU_N-B_k2+k];
                tmp_4.coeffs[NTRU_N-sum+k] = f_row.coeffs[NTRU_N-B_k2+k]^s_row.coeffs[NTRU_N-B_k2+k];
            }
            g = tmp_1;
            f_row = tmp_2;
            f = tmp_3;
            s_row = tmp_4;
           
            i = g_deg - B_k2;
            while(g.coeffs[i] == 0) i = i - 1;
            g_deg = i;
           
            j = f_deg - B_k1;
            while(f.coeffs[j] == 0) j = j - 1;
            f_deg = j;
       }
       else if(delta1 == -1 && delta2 == 0)
       {
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = g.coeffs[i+B_k2]^g.coeffs[i+sum]^f.coeffs[i+sum];
               tmp_2.coeffs[i] = f_row.coeffs[i+B_k2]^f_row.coeffs[i+sum]^s_row.coeffs[i+sum];
           }
           for (i = 0;i < B_k1;++i)
           {
               tmp_1.coeffs[NTRU_N-sum+i] = g.coeffs[NTRU_N-B_k1+i];
               tmp_2.coeffs[NTRU_N-sum+i] = f_row.coeffs[NTRU_N-B_k1+i]^f_row.coeffs[i]^s_row.coeffs[i];
           }
           for (i = 0;i < B_k2;++i)
           {
               tmp_1.coeffs[NTRU_N-1-i] = 0;
               tmp_2.coeffs[NTRU_N-B_k2+i] = f_row.coeffs[i]^f_row.coeffs[B_k1+i]^s_row.coeffs[B_k1+i];
           }
           f = g;
           g = tmp_1;
           s_row = f_row;
           f_row = tmp_2;

           tmp1 = f_deg;
           f_deg = g_deg;
           
           i = tmp1 - sum;
           while(g.coeffs[i] == 0) i = i - 1;
           g_deg = i;
       }
       else if(delta1 == 0 && delta2 == -1)
       {
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = f.coeffs[i+B_k2]^g.coeffs[i+sum]^f.coeffs[i+sum];
               tmp_2.coeffs[i] = s_row.coeffs[i+B_k2]^f_row.coeffs[i+sum]^s_row.coeffs[i+sum];
               tmp_3.coeffs[i] = g.coeffs[i+B_k1]^f.coeffs[i+B_k1];
               tmp_4.coeffs[i] = f_row.coeffs[i+B_k1]^s_row.coeffs[i+B_k1];
           }
           for (i = 0;i < B_k1;++i)
           {
               tmp_1.coeffs[NTRU_N-sum+i] = f.coeffs[NTRU_N-B_k1+i];
               tmp_2.coeffs[NTRU_N-sum+i] = s_row.coeffs[NTRU_N-B_k1+i]^f_row.coeffs[i]^s_row.coeffs[i];
               tmp_3.coeffs[NTRU_N-1-i] = 0;
               tmp_4.coeffs[NTRU_N-B_k1+i] = f_row.coeffs[i]^s_row.coeffs[i];
           }
           for (i = 0;i < B_k2;++i)
           {
               tmp_1.coeffs[NTRU_N-1-i] = 0;
               tmp_2.coeffs[NTRU_N-B_k2+i] = s_row.coeffs[i]^f_row.coeffs[B_k1+i]^s_row.coeffs[B_k1+i];
               tmp_3.coeffs[NTRU_N-sum+i] = g.coeffs[NTRU_N-B_k2+i]^f.coeffs[NTRU_N-B_k2+i];
               tmp_4.coeffs[NTRU_N-sum+i] = f_row.coeffs[NTRU_N-B_k2+i]^s_row.coeffs[NTRU_N-B_k2+i];
           }
           g = tmp_1;
           f_row = tmp_2;
           f = tmp_3;
           s_row = tmp_4;

           tmp1 = f_deg;
           i = g_deg - B_k1;
           while(f.coeffs[i] == 0) i = i - 1;
           f_deg = i;
           
           j = tmp1 - B_k2;
           while(g.coeffs[j] == 0) j = j - 1;
           g_deg = j;
       }
       else if(delta1 == 0 && delta2 == 0)
       {
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = f.coeffs[i+B_k2]^g.coeffs[i+sum]^f.coeffs[i+sum];
               tmp_2.coeffs[i] = s_row.coeffs[i+B_k2]^f_row.coeffs[i+sum]^s_row.coeffs[i+sum];
           }
           for (i = 0;i < B_k1;++i)
           {
               tmp_1.coeffs[NTRU_N-sum+i] = f.coeffs[NTRU_N-B_k1+i];
               tmp_2.coeffs[NTRU_N-sum+i] = s_row.coeffs[NTRU_N-B_k1+i]^f_row.coeffs[i]^s_row.coeffs[i];
           }
           for (i = 0;i < B_k2;++i)
           {
               tmp_1.coeffs[NTRU_N-1-i] = 0;
               tmp_2.coeffs[NTRU_N-B_k2+i] = s_row.coeffs[i]^f_row.coeffs[B_k1+i]^s_row.coeffs[B_k1+i];
           }
           g = tmp_1;
           f_row = tmp_2;

           i = g_deg - sum;
           while(g.coeffs[i] == 0) i = i - 1;
           g_deg = i;
       }
       B_k1 = 0; B_k2 = 0;
    }
    if(g_deg == 0)
    {
        for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] = (f_row.coeffs[i] ^ f_row.coeffs[NTRU_N-1]);
        r->coeffs[NTRU_N-1] = 0;
    }
    else if(f_deg == 0)
    {
        for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] = (s_row.coeffs[i] ^ s_row.coeffs[NTRU_N-1]);
        r->coeffs[NTRU_N-1] = 0;
    }
}



