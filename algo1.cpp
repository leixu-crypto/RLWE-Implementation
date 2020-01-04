#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/RR.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace NTL;

ZZ r = ZZ(10);
ZZ q = ZZ(1009);
int l, d, k;

void merge(mat_ZZ &a, mat_ZZ &b, int dir=0)
{//0 for row extending, otherwise column extending.
    mat_ZZ tmpa = a, tmpb=b;

    if (dir == 0)
    {
        long n = a.NumRows();
        long m1 = a.NumCols();
        long m2 = b.NumCols();
        a.SetDims(n, m1+m2);
        for (int i = 0 ; i < n; i++){
            for (int j = 0 ; j < m1; j++)
                a[i][j] = tmpa[i][j];
            for (int j = 0; j < m2; j++)
                a[i][j+m1] = tmpb[i][j];
        }
    }
    else{
        long n1 = a.NumRows();
        long n2 = b.NumRows();
        long m = b.NumCols();
        a.SetDims(n1+n2, m);

        for (int i = 0; i < n1; i++)
            for (int j = 0 ; j < m; j++)
                a[i][j] = tmpa[i][j];

        for (int i = 0 ; i < n2; i++)
            for (int j = 0; j < m; j++)
                a[i+n1][j] = tmpb[i][j];
    }
}

void merge(mat_ZZ_p &a, mat_ZZ_p &b, int dir=0)
{//0 for row extending, otherwise column extending.
    mat_ZZ_p tmpa = a, tmpb=b;

    if (dir == 0)
    {
        long n = a.NumRows();
        long m1 = a.NumCols();
        long m2 = b.NumCols();
        a.SetDims(n, m1+m2);
        for (int i = 0 ; i < n; i++){
            for (int j = 0 ; j < m1; j++)
                a[i][j] = tmpa[i][j];
            for (int j = 0; j < m2; j++)
                a[i][j+m1] = tmpb[i][j];
        }
    }
    else{
        long n1 = a.NumRows();
        long n2 = b.NumRows();
        long m = b.NumCols();
        a.SetDims(n1+n2, m);

        for (int i = 0; i < n1; i++)
            for (int j = 0 ; j < m; j++)
                a[i][j] = tmpa[i][j];

        for (int i = 0 ; i < n2; i++)
            for (int j = 0; j < m; j++)
                a[i+n1][j] = tmpb[i][j];
    }
}



int algo1(mat_ZZ_p& A1, long m2, mat_ZZ_p& A2, mat_ZZ& S)
{
    long n = A1.NumRows();
    long m1 = A1.NumCols();
    long m = m1+m2;

    mat_ZZ H, Hp, G, P, U, T, R, C;
    RR tmp, t2;
    ZZ tmpz;

    P.SetDims(m2, m1);
    H.SetDims(m1, m1);
    Hp.SetDims(m1, m1);
    G.SetDims(m1, m2);
    R.SetDims(m1, m2);

    ident(U, m2);
    ident(C, m1);
    clear(H);
    clear(Hp);

    for (int i = 0 ; i < m1; i++)
        H[i][i] = q, Hp[i][i] = q-1;


    ident(T, l);

    for (int i = 0; i < l-1; i++)
    {
        T[i][i+1] = -r;
    }

    for (int i = 0; i < m1; i++)
    {
        for (int j = 0; j < l-1; j++)
        {
            U[i*l+j][i*l+j+1] = T[j][j+1];
        }

        for (int k = 0; k < m1; k++){
            G[k][i*l+(l-1)] = Hp[k][i];
        }

        P[i*l+(l-1)][i] = 1; // e_jl
        for (int j = l-2; j>=0; j--)
            for (int k = 0; k < m1; k++)
                G[k][i*l+j] = G[k][i*l + j+1] / r;
    }



    srand(time(NULL));
    int rr;
    clear(R);


    for (int i = 0; i < d; i++)
        for (int j = 0; j < m2; j++)
        {
            rr = rand()%4;
            if (rr == 0) R[i][j] = 1;
            else if (rr == 1) R[i][j] = -1;
            else R[i][j] = 0;
        }



    mat_ZZ tmp_A1, tmp_A2;
    conv(tmp_A1, A1);

    tmp_A2 = -tmp_A1*(R+G);
    merge(tmp_A1, tmp_A2);

    G = (G+R)*U;
    R = R*P - C;

    merge(G, R);

    merge(U, P);

    mat_ZZ_p res;
    merge(G, U, 1);
    conv(A2, tmp_A2);
    S = G;
    conv(res, tmp_A1*G);

    return 0;
}


int main(){
    //n=512
    clock_t endTime, startTime;
    startTime = clock();
    long n=32, m, m1, m2;

    RR t1, t2;
    ZZ tmp_z;
    mat_ZZ_p A1, A2;
    mat_ZZ S;

    ZZ_p::init(q);
    conv(t1, n);
    conv(t2, 2);
    t1 = log10(t1);
    t2 = log10(t2);
    t1 = 3*t1/t2; // k=3logn
    CeilToZZ(tmp_z, t1);
    conv(k, tmp_z);

    conv(t1, q);
    conv(t2, 2);
    log10(t1, t1); // l = ceil(logr_q), r=10
    log10(t2, t2);
    t2 = t1/t2;
    t2 *= n*(1+1/3);
    CeilToZZ(tmp_z, t1);
    conv(l, tmp_z);
    CeilToZZ(tmp_z, t2);
    conv(d, tmp_z);

    ZZ_pX p;
    ZZ_pX ch;
    BuildIrred(p, n);
    ZZ_pE::init(p);
    ZZ_pEX f, g, h;
    random(f, n+10);
    SetCoeff(f, n+10);

    int m1bd = (d/n) + 1;
    k = k>m1bd?k:m1bd;

    m1 = n*k;
    m2 = m1*l;
    cout<<"n m1 m2: "<<n<<" "<<m1<<" "<<m2<<endl;
    cout<<"d k l: "<<d<<" "<<k<<" "<<l<<endl;

    A1.SetDims(n, m1);
    A2.SetDims(n, m2);
    for (int i = 0; i < k; i++)
    {
        random(h, n+10);
        g = MinPolyMod(h, f);
        conv(ch, g[0]);

        if (CompMod(g, h, f) != 0) // check that g(h) = 0 mod f
            Error("g(h) != 0 mod f");
        for (int j = 0; j < n; j++)
            GetCoeff(A1[0][i*n+j],ch, j);

        for (int k = 1; k < n; k++)
            for (int j = 0; j < n; j++){
                A1[k][i*n+j] = A1[k-1][i*n+((j-1+n)%n)];
            }
    }

    endTime = clock();
	cout << "Setup Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;


    startTime = clock();
    algo1(A1, m2, A2, S);
    endTime = clock();
	cout << "Key generation Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    startTime = clock();
    merge(A1, A2);
    mat_ZZ tmp_mat;
    mat_ZZ_p res;
    conv(tmp_mat, A1);
    tmp_mat = tmp_mat*S;
    conv(res, tmp_mat);
    cout<<IsZero(res)<<endl;
    endTime = clock();
	cout << "Verification Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    return 0;
}