


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


void poly_S3_inv_variable(poly *r, poly *a);
void poly_S3_inv_variable_combine(poly *r, poly *a);


int main()
{
    poly inv_variable, inv_variable_combine;
    poly gf;

    int i, k, tmp1, tmp2, avg_variable, avg_variable_combine;
    double duration;
    clock_t start, finish;

    tmp1 = 0; tmp2 = 0;
    avg_variable = 0; avg_variable_combine = 0;

    for(k=0; k< iteration; k++)
    {
        //For ternary random coefficient generation
        for(i=0; i<NTRU_N-1; i++)gf.coeffs[i] = rand()&3;

        gf.coeffs[NTRU_N-1] = 0;
    
        /////////////////////////////////////////////////////////////////////////////
        //Run NTRU S3_variable inverse
        start =clock();
        poly_S3_inv_variable(&inv_variable,&gf);
        finish = clock();
        duration = (double)(finish-start);
        avg_variable = avg_variable + duration;
        ////////////////////////////////////////////////////////////////////////////
        ///
        ////////////////////////////////////////////////////////////////////////////
        //Run NTRU S3_variable_combine inverse
        start =clock();
        poly_S3_inv_variable_combine(&inv_variable_combine,&gf);
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
            for(i=0; i<NTRU_N; i++)printf("%d",inv_variable.coeffs[i]);
            printf("\n");
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

static inline uint8_t mod3(uint8_t a) /* a between 0 and 9 (I think max is 15)*/
{
    int16_t t, c;
    a = (a >> 2) + (a & 3); /* between 0 and 4 */
    t = a - 3;
    c = t >> 5;
    return (uint8_t) (t^(c&(a^t)));
}

static inline uint8_t mod3h(uint8_t a) /* a between 0 and 20 */
{
    int16_t t, c, k, l;
    a = (a >> 4) + (a & 3) + ((a>>2)&3);
    a = (a >> 2) + (a & 3); /* between 0 and 4 */
    t = a - 3;
    c = t >> 5;
    return (uint8_t) (t^(c&(a^t)));
}

//mod3 multiplication funcion
//However, it is slower than ordinary multiplication.
/*static inline uint16_t mul3(uint16_t a, uint16_t b)
{
    return (uint16_t) (((a&b)|((a&b)>>1))&1)|(((a^b)&((a^b)<<1))&2);
}*/

void poly_S3_inv_variable(poly *r, poly *a)
{
    poly f, g, f_row, s_row, tmp_1,tmp_2;
    int i;
    int16_t delta, f_deg, g_deg, B_k, initial_B, tmp;
  
    delta = 0;
    B_k = 0;
    initial_B = 0;

    for (i = 0;i < NTRU_N;++i)
    {
        g.coeffs[i] = 1;
        f_row.coeffs[i] = 0;
        s_row.coeffs[i] = 0;
    }
    s_row.coeffs[0] = 1;

    for (i = 0;i < NTRU_N-1;++i) f.coeffs[i] = mod3((uint8_t)(a->coeffs[i] + 2*a->coeffs[NTRU_N-1]));
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
  
    while(g_deg>0)
    {
        //delta is zero(0000) when g_deg is greater than or equal to f_deg.
        //delta is minus one(1111) when f_deg is greater than g_deg.
        delta = (g_deg - f_deg)>>15;
      
        //Calculate B_k
        i = 0;
        tmp = mod3((uint8_t)((g.coeffs[0]*f.coeffs[0]) + ((f.coeffs[0]^3)*g.coeffs[0])));
        while(tmp == 0)
        {
            i = i + 1;
            tmp = mod3((uint8_t)((g.coeffs[0]*f.coeffs[i]) + ((f.coeffs[0]^3)*g.coeffs[i])));
        }
        B_k = i;
        
        for (i = 0;i < NTRU_N-B_k;++i)
        {
            tmp_1.coeffs[i] = mod3((uint8_t)((g.coeffs[0]*f.coeffs[i+B_k]) + ((f.coeffs[0]^3)*g.coeffs[i+B_k])));
            tmp_2.coeffs[i] = mod3((uint8_t)((g.coeffs[0]*s_row.coeffs[i+B_k]) + ((f.coeffs[0]^3)*f_row.coeffs[i+B_k])));
        }
      
        for (i = 0;i < B_k;++i)
        {
            tmp_1.coeffs[NTRU_N-1-i] = 0;
            tmp_2.coeffs[NTRU_N-B_k+i] = mod3((uint8_t)((g.coeffs[0]*s_row.coeffs[i]) + ((f.coeffs[0]^3)*f_row.coeffs[i])));
        }

        B_k = 0;
        
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
        
    }
    for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] =  mod3((uint8_t)(g.coeffs[0]*(f_row.coeffs[i] + 2*f_row.coeffs[NTRU_N-1])));
    r->coeffs[NTRU_N-1] = 0;
}

void poly_S3_inv_variable_combine(poly *r, poly *a)
{
    poly f, g, f_row, s_row, tmp_1, tmp_2, tmp_3, tmp_4;
    int16_t i, j, tmp1, delta1, delta2, f_deg, g_deg, diff_deg, B_k1, B_k2, initial_B, sum, t_con1, t_con2;
    int16_t t1, t2, t3;
    
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
    for (i = 0;i < NTRU_N-1;++i) f.coeffs[i] = mod3((uint8_t)(a->coeffs[i] + 2*a->coeffs[NTRU_N-1]));
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
        t_con1= mod3((uint8_t)((g.coeffs[0]*f.coeffs[0]) + ((f.coeffs[0]^3)*g.coeffs[0])));
        while(t_con1 == 0)
        {
            i = i + 1;
            t_con1 = mod3((uint8_t)((g.coeffs[0]*f.coeffs[i]) + ((f.coeffs[0]^3)*g.coeffs[i])));
        }
        B_k1 = i;
        
        //Calculate B_k2
        if(delta1 == 0)
        {
            if(diff_deg >= B_k1) delta2 = 0;
            else delta2 = -1;

            i = 0;
            t1 = mod3((uint8_t)(f.coeffs[0]*g.coeffs[0]));
            t2 = (f.coeffs[0]|(f.coeffs[0]<<1))&2;//mod3(f.coeffs[0]*(-f.coeffs.[0]))
            t_con2 = mod3((uint8_t)(t1*f.coeffs[B_k1]+t2*g.coeffs[B_k1]+(t_con1^3)*f.coeffs[0]));
            while(t_con2 == 0)
            {
                i = i + 1;
                t_con2 = mod3((uint8_t)(t1*f.coeffs[i+B_k1]+t2*g.coeffs[i+B_k1]+(t_con1^3)*f.coeffs[i]));
            }
            B_k2 = i;
        }
        else
        {
            if(-diff_deg >= B_k1) delta2 = 0;
            else delta2 = -1;
            
            i = 0;
            t1 = (g.coeffs[0]|(g.coeffs[0]>>1))&1;//g.coeffs[0]*g.coeffs[0]
            t2 = mod3((uint8_t)(g.coeffs[0] * (f.coeffs[0]^3)));
            t_con2 = mod3((uint8_t)(t1*f.coeffs[B_k1]+t2*g.coeffs[B_k1]+(t_con1^3)*g.coeffs[0]));
            while(t_con2 == 0)
            {
                i = i + 1;
                t_con2 = mod3((uint8_t)(t1*f.coeffs[i+B_k1]+t2*g.coeffs[i+B_k1]+(t_con1^3)*g.coeffs[i]));
            }
            B_k2 = i;
        }
        
       sum = B_k1 + B_k2;
        
       if( delta1 == -1 && delta2 == -1)
       {
           t1 = (g.coeffs[0]|(g.coeffs[0]>>1))&1;//g.coeffs[0]*g.coeffs[0]
           t2 = mod3((uint8_t)(g.coeffs[0] * (f.coeffs[0]^3)));
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = mod3((uint8_t)(t1*f.coeffs[i+sum]+t2*g.coeffs[i+sum]+(t_con1^3)*g.coeffs[i+B_k2]));
               tmp_2.coeffs[i] = mod3((uint8_t)(t1*s_row.coeffs[i+sum]+t2*f_row.coeffs[i+sum]+(t_con1^3)*f_row.coeffs[i+B_k2]));
               tmp_3.coeffs[i] = mod3((uint8_t)(g.coeffs[0]*f.coeffs[i+B_k1]+(f.coeffs[0]^3)*g.coeffs[i+B_k1]));
               tmp_4.coeffs[i] = mod3((uint8_t)(g.coeffs[0]*s_row.coeffs[i+B_k1]+(f.coeffs[0]^3)*f_row.coeffs[i+B_k1]));
           }
           for (i = 0;i < B_k1;++i)
           {
               tmp_1.coeffs[NTRU_N-sum+i] = mod3((uint8_t)((t_con1^3)*g.coeffs[NTRU_N-B_k1+i]));
               tmp_2.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(t1*s_row.coeffs[i]+t2*f_row.coeffs[i]+(t_con1^3)*f_row.coeffs[NTRU_N-B_k1+i]));
               tmp_3.coeffs[NTRU_N-1-i] = 0;
               tmp_4.coeffs[NTRU_N-B_k1+i] = mod3((uint8_t)(g.coeffs[0]*s_row.coeffs[i]+(f.coeffs[0]^3)*f_row.coeffs[i]));
           }
           for (i = 0;i < B_k2;++i)
           {
               tmp_1.coeffs[NTRU_N-1-i] = 0;
               tmp_2.coeffs[NTRU_N-B_k2+i] = mod3((uint8_t)(t1*s_row.coeffs[i+B_k1]+t2*f_row.coeffs[i+B_k1]+(t_con1^3)*f_row.coeffs[i]));
               tmp_3.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(g.coeffs[0]*f.coeffs[NTRU_N-B_k2+i]+(f.coeffs[0]^3)*g.coeffs[NTRU_N-B_k2+i]));
               tmp_4.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(g.coeffs[0]*s_row.coeffs[NTRU_N-B_k2+i]+(f.coeffs[0]^3)*f_row.coeffs[NTRU_N-B_k2+i]));
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
           t1 = (g.coeffs[0]|(g.coeffs[0]>>1))&1;//g.coeffs[0]*g.coeffs[0]
           t2 = mod3((uint8_t)(g.coeffs[0] * (f.coeffs[0]^3)));
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = mod3((uint8_t)(t1*f.coeffs[i+sum]+t2*g.coeffs[i+sum]+(t_con1^3)*g.coeffs[i+B_k2]));
               tmp_2.coeffs[i] = mod3((uint8_t)(t1*s_row.coeffs[i+sum]+t2*f_row.coeffs[i+sum]+(t_con1^3)*f_row.coeffs[i+B_k2]));
           }
           for (i = 0;i < B_k1;++i)
           {
               tmp_1.coeffs[NTRU_N-sum+i] = mod3((uint8_t)((t_con1^3)*g.coeffs[NTRU_N-B_k1+i]));
               tmp_2.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(t1*s_row.coeffs[i]+t2*f_row.coeffs[i]+(t_con1^3)*f_row.coeffs[NTRU_N-B_k1+i]));
           }
           for (i = 0;i < B_k2;++i)
           {
               tmp_1.coeffs[NTRU_N-1-i] = 0;
               tmp_2.coeffs[NTRU_N-B_k2+i] = mod3((uint8_t)(t1*s_row.coeffs[i+B_k1]+t2*f_row.coeffs[i+B_k1]+(t_con1^3)*f_row.coeffs[i]));
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
           t1 = mod3((uint8_t)(f.coeffs[0]*g.coeffs[0]));
           t2 = (f.coeffs[0]|(f.coeffs[0]<<1))&2;//mod3(f.coeffs[0]*(-f.coeffs.[0]))
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = mod3((uint8_t)(t1*f.coeffs[i+sum]+t2*g.coeffs[i+sum]+(t_con1^3)*f.coeffs[i+B_k2]));
               tmp_2.coeffs[i] = mod3((uint8_t)(t1*s_row.coeffs[i+sum]+t2*f_row.coeffs[i+sum]+(t_con1^3)*s_row.coeffs[i+B_k2]));
               tmp_3.coeffs[i] = mod3((uint8_t)(g.coeffs[0]*f.coeffs[i+B_k1]+(f.coeffs[0]^3)*g.coeffs[i+B_k1]));
               tmp_4.coeffs[i] = mod3((uint8_t)(g.coeffs[0]*s_row.coeffs[i+B_k1]+(f.coeffs[0]^3)*f_row.coeffs[i+B_k1]));
           }
           for (i = 0;i < B_k1;++i)
           {
               tmp_1.coeffs[NTRU_N-sum+i] = mod3((uint8_t)((t_con1^3)*f.coeffs[NTRU_N-B_k1+i]));
               tmp_2.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(t1*s_row.coeffs[i]+t2*f_row.coeffs[i]+(t_con1^3)*s_row.coeffs[NTRU_N-B_k1+i]));
               tmp_3.coeffs[NTRU_N-1-i] = 0;
               tmp_4.coeffs[NTRU_N-B_k1+i] = mod3((uint8_t)(g.coeffs[0]*s_row.coeffs[i]+(f.coeffs[0]^3)*f_row.coeffs[i]));
           }
           for (i = 0;i < B_k2;++i)
           {
               tmp_1.coeffs[NTRU_N-1-i] = 0;
               tmp_2.coeffs[NTRU_N-B_k2+i] = mod3((uint8_t)(t1*s_row.coeffs[i+B_k1]+t2*f_row.coeffs[i+B_k1]+(t_con1^3)*s_row.coeffs[i]));
               tmp_3.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(g.coeffs[0]*f.coeffs[NTRU_N-B_k2+i]+(f.coeffs[0]^3)*g.coeffs[NTRU_N-B_k2+i]));
               tmp_4.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(g.coeffs[0]*s_row.coeffs[NTRU_N-B_k2+i]+(f.coeffs[0]^3)*f_row.coeffs[NTRU_N-B_k2+i]));
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
           t1 = mod3((uint8_t)(f.coeffs[0]*g.coeffs[0]));
           t2 = (f.coeffs[0]|(f.coeffs[0]<<1))&2;//mod3(f.coeffs[0]*(-f.coeffs.[0]))
           for (i = 0;i < NTRU_N-sum;++i)
           {
               tmp_1.coeffs[i] = mod3((uint8_t)(t1*f.coeffs[i+sum]+t2*g.coeffs[i+sum]+(t_con1^3)*f.coeffs[i+B_k2]));
               tmp_2.coeffs[i] = mod3((uint8_t)(t1*s_row.coeffs[i+sum]+t2*f_row.coeffs[i+sum]+(t_con1^3)*s_row.coeffs[i+B_k2]));
           }
           for (i = 0;i < B_k1;++i)
           {
               tmp_1.coeffs[NTRU_N-sum+i] = mod3((uint8_t)((t_con1^3)*f.coeffs[NTRU_N-B_k1+i]));
               tmp_2.coeffs[NTRU_N-sum+i] = mod3((uint8_t)(t1*s_row.coeffs[i]+t2*f_row.coeffs[i]+(t_con1^3)*s_row.coeffs[NTRU_N-B_k1+i]));
           }
           for (i = 0;i < B_k2;++i)
           {
               tmp_1.coeffs[NTRU_N-1-i] = 0;
               tmp_2.coeffs[NTRU_N-B_k2+i] = mod3((uint8_t)(t1*s_row.coeffs[i+B_k1]+t2*f_row.coeffs[i+B_k1]+(t_con1^3)*s_row.coeffs[i]));
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
        for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] =  mod3((uint8_t)(g.coeffs[0]*(f_row.coeffs[i] + 2*f_row.coeffs[NTRU_N-1])));
        r->coeffs[NTRU_N-1] = 0;
    }
    else if(f_deg == 0)
    {
        for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] =  mod3((uint8_t)(f.coeffs[0]*(s_row.coeffs[i] + 2*s_row.coeffs[NTRU_N-1])));
        r->coeffs[NTRU_N-1] = 0;
    }
}



