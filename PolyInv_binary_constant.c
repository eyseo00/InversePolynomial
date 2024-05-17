
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


void poly_R2_inv(poly *r, poly *a);
void poly_R2_inv_constant(poly *r, poly *a);

int main()
{
    poly inv_org, inv_constant;
    poly gf;

    int i, k, tmp1, tmp2, avg_org, avg_constant;
    double duration;
    clock_t start, finish;

    tmp1 = 0; tmp2 = 0;
    avg_org = 0; avg_constant = 0;

    for(k=0; k< iteration; k++){
        //For binary random coefficient generation
        for(i=0; i<NTRU_N-1; i++) gf.coeffs[i] = rand()&1;
        //For ternary random coefficient generation
        //for(i=0; i<NTRU_N-1; i++)gf.coeffs[i] = rand()&3;

        gf.coeffs[NTRU_N-1] = 0;
    
        /////////////////////////////////////////////////////////////////////////////
        //Run NTRU R2 original inverse
        start =clock();
        poly_R2_inv(&inv_org,&gf);
        finish = clock();
        duration = (double)(finish-start);
        avg_org = avg_org + duration;
        ////////////////////////////////////////////////////////////////////////////
        ///
        ////////////////////////////////////////////////////////////////////////////
        //Run NTRU R2_constant inverse
        start =clock();
        poly_R2_inv_constant(&inv_constant,&gf);
        finish = clock();
        duration = (double)(finish-start);
        avg_constant = avg_constant + duration;
        ////////////////////////////////////////////////////////////////////////////
        ///
        ////////////////////////////////////////////////////////////////////////////
        //Comparing results
        //Checking that the result is correct
        for(i=0; i<NTRU_N-1; i++)
        {
            if(inv_org.coeffs[i] != inv_constant.coeffs[i]) tmp2 = tmp2 + 1;
        }
        if(tmp2 != 0)
        {
            printf("Alarm k: %d\n", k);
            for(i=0; i<NTRU_N; i++)printf("%d",inv_org.coeffs[i]);
            printf("\n");
            for(i=0; i<NTRU_N; i++)printf("%d",inv_constant.coeffs[i]);
            printf("\n");
        }
        tmp2 = 0;
        ////////////////////////////////////////////////////////////////////////////
    }
    avg_org = avg_org/iteration;
    avg_constant = avg_constant/iteration;
    printf("avg_org: %d, avg_constant: %d\n", avg_org, avg_constant);
    return 0;
}
    
static inline int16_t both_negative_mask(int16_t x,int16_t y)
{
  return (x & y) >> 15;
}

void poly_R2_inv(poly *r, poly *a)
{
  poly f, g, v, w;
  size_t i, loop;
  int16_t delta,sign,swap,t;

  for (i = 0;i < NTRU_N;++i) v.coeffs[i] = 0;
  for (i = 0;i < NTRU_N;++i) w.coeffs[i] = 0;
  w.coeffs[0] = 1;

  for (i = 0;i < NTRU_N;++i) f.coeffs[i] = 1;
  for (i = 0;i < NTRU_N-1;++i) g.coeffs[NTRU_N-2-i] = (a->coeffs[i] ^ a->coeffs[NTRU_N-1]) & 1;
  g.coeffs[NTRU_N-1] = 0;

  delta = 1;

  for (loop = 0;loop < 2*(NTRU_N-1)-1;++loop) {
    for (i = NTRU_N-1;i > 0;--i) v.coeffs[i] = v.coeffs[i-1];
    v.coeffs[0] = 0;

    sign = g.coeffs[0] & f.coeffs[0];
    swap = both_negative_mask(-delta,-(int16_t) g.coeffs[0]);
    delta ^= swap & (delta ^ -delta);
    delta += 1;

    for (i = 0;i < NTRU_N;++i) {
      t = swap&(f.coeffs[i]^g.coeffs[i]); f.coeffs[i] ^= t; g.coeffs[i] ^= t;
      t = swap&(v.coeffs[i]^w.coeffs[i]); v.coeffs[i] ^= t; w.coeffs[i] ^= t;
    }

    for (i = 0;i < NTRU_N;++i) g.coeffs[i] = g.coeffs[i]^(sign&f.coeffs[i]);
    for (i = 0;i < NTRU_N;++i) w.coeffs[i] = w.coeffs[i]^(sign&v.coeffs[i]);
    for (i = 0;i < NTRU_N-1;++i) g.coeffs[i] = g.coeffs[i+1];
    g.coeffs[NTRU_N-1] = 0;
  }

  for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] = v.coeffs[NTRU_N-2-i];
  r->coeffs[NTRU_N-1] = 0;
}

void poly_R2_inv_constant(poly *r, poly *a)
{
    poly f, g, f_row, s_row, tmp_g, tmp_f, tmp_f_row, tmp_s_row;
    int i, j;
    int16_t a0, a1, b0, b1, g_deg, f_deg, delta, tmp_g_deg, tmp_f_deg, mu, n_mu, n_delta, con1, con2;
    int16_t t1, t2, t3, t4, t5, t6, t8, t9, t10, t11, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11;
    int16_t gG, gF, g1, g2, fG, fF, f1, f2;

    for (i = 0;i < NTRU_N;++i)
    {
        f_row.coeffs[i] = 0;
        s_row.coeffs[i] = 0;
    }
    s_row.coeffs[0] = 1;
  
    for (i = 0;i < NTRU_N;++i) g.coeffs[i] = 1;
    for (i = 0;i < NTRU_N-1;++i)
    {
        f.coeffs[i] = (a->coeffs[i] ^ a->coeffs[NTRU_N-1]) & 1;
    }
    f.coeffs[NTRU_N-1] = 0;
  
    g_deg = NTRU_N;
    f_deg = NTRU_N;
    delta = 0;//g_deg>f_deg-->delta=0-->reduce first row polynomial in the summation

    for(j = 0;j < NTRU_N;++j)
    {
        a0 = g.coeffs[0];a1 = g.coeffs[1];
        b0 = f.coeffs[0];b1 = f.coeffs[1];

        mu = (a0&b1)^(a1&b0);
        n_mu = mu^1;

        n_delta = delta^1;
        t1 = (a0|a1);
        t2 = (b0|b1);

        t3 = t1&n_delta&b0;//(a0+a1)*n_delta*b0
        t4 = t1&n_delta&(b0^1)&b1;//(a0+a1)*n_delta*n_b0*b1
        t5 = t2&delta&a0&mu;//(b0+b1)*delta*a0*mu
        t6 = t2&delta&a0&n_mu;//(b0+b1)*delta*a0*n_mu
        t8 = t2&delta&(a0^1)&a1;//(b0+b1)&delta*n_a0*a1
        t9 = (t1^1)&n_delta;//n_a0*n_a1*n_delta
        t10 = (t1^1)&delta;//n_a0*n_a1*delta
        t11 = t1&(t2^1);//(a0+a1)*n_b0*n_b1
    
        m1 = t6 + t11;
        m2 = t4 + t5 + t8;
        m3 = t3 + (t5&b0) + t9;

        m4 = t10;
        m5 = t3&mu;
        m6 = (t3&a0) + (t4&a0) + t5;

        m7 = t5&b0;
        m8 = ((t6 + t8)&b0) + t10;

        m9 = t3 + t9;
        m10 = t4 + t5 + t8;
        m11 = t6 + t11;

        gG = (n_delta|t1); gF = t10;
        g1 = t4 + t5 + t8;
        g2 = t3 + t9;

        fG = gF; fF = gG; f1 = g1;
        f2 = t6 + t10 + t11;

        tmp_g_deg = gG*g_deg + gF*f_deg - g1 - 2*g2;
        tmp_f_deg = fG*g_deg + fF*f_deg - f1 - 2*f2;
        g_deg = tmp_g_deg;
        f_deg = tmp_f_deg;
        
        con1 = NTRU_N - 1; con2 = NTRU_N - 2;

        for (i = 0;i < NTRU_N-2;++i)
        {
            tmp_g.coeffs[i] = ((m1&g.coeffs[i]) ^ (m4&f.coeffs[i])) ^ ((m2&g.coeffs[i+1]) ^ (m5&f.coeffs[i+1])) ^ ((m3&g.coeffs[i+2]) ^ (m6&f.coeffs[i+2]));
            tmp_f.coeffs[i] = ((m9&f.coeffs[i]) ^ (m7&g.coeffs[i+1])) ^ ((m10&f.coeffs[i+1]) ^ (m8&g.coeffs[i+2])) ^ (m11&f.coeffs[i+2]);
            tmp_f_row.coeffs[i] = ((m1&f_row.coeffs[i]) ^ (m4&s_row.coeffs[i])) ^ ((m2&f_row.coeffs[i+1]) ^ (m5&s_row.coeffs[i+1])) ^ ((m3&f_row.coeffs[i+2]) ^ (m6&s_row.coeffs[i+2]));
            tmp_s_row.coeffs[i] = ((m9&s_row.coeffs[i]) ^ (m7&f_row.coeffs[i+1])) ^ ((m10&s_row.coeffs[i+1]) ^ (m8&f_row.coeffs[i+2])) ^ (m11&s_row.coeffs[i+2]);
        }
        
        tmp_g.coeffs[con2] = ((m1&g.coeffs[con2]) ^ (m4&f.coeffs[con2])) ^ ((m2&g.coeffs[con1]) ^ (m5&f.coeffs[con1]));
        tmp_g.coeffs[con1] = ((m1&g.coeffs[con1]) ^ (m4&f.coeffs[con1]));
        tmp_f.coeffs[con2] = (m9&f.coeffs[con2]) ^ ((m7&g.coeffs[con1]) ^ (m10&f.coeffs[con1]));
        tmp_f.coeffs[con1] = (m9&f.coeffs[con1]);
        tmp_f_row.coeffs[con2] = ((m1&f_row.coeffs[con2]) ^ (m4&s_row.coeffs[con2])) ^ ((m2&f_row.coeffs[con1]) ^(m5&s_row.coeffs[con1])) ^ ((m3&f_row.coeffs[0]) ^ (m6&s_row.coeffs[0]));
        tmp_f_row.coeffs[con1] = ((m1&f_row.coeffs[con1]) ^ (m4&s_row.coeffs[con1])) ^ ((m2&f_row.coeffs[0]) ^ (m5&s_row.coeffs[0])) ^((m3&f_row.coeffs[1]) ^ (m6&s_row.coeffs[1]));
        tmp_s_row.coeffs[con2] = (m9&s_row.coeffs[con2]) ^ ((m7&f_row.coeffs[con1]) ^ (m10&s_row.coeffs[con1])) ^((m8&f_row.coeffs[0]) ^ (m11&s_row.coeffs[0]));
        tmp_s_row.coeffs[con1] = (m9&s_row.coeffs[con1]) ^ ((m7&f_row.coeffs[0]) ^ (m10&s_row.coeffs[0])) ^ ((m8&f_row.coeffs[1]) ^(m11&s_row.coeffs[1]));

        g = tmp_g;
        f = tmp_f;
        f_row = tmp_f_row;
        s_row = tmp_s_row;
      
        delta = ((g_deg - f_deg - 1)>>15)&1;
    }
    for (i = 0;i < NTRU_N-1;++i) r->coeffs[i] = (f_row.coeffs[i] ^ f_row.coeffs[NTRU_N-1]);
    r->coeffs[NTRU_N-1] = 0;
}



