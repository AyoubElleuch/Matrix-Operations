/*
    Mohamed Ayoub Eleuch
    CS-04
    m.eleuch@innopolis.university
*/

#include<iostream>
#include<cstdio>
#include<vector>
#include<math.h>

using namespace std;


class ColumnVector{
    public:
        int n;
        vector<double> v;
        ColumnVector(int m){
            v = vector<double>(m, 0);
            n = m;
        }

        void operator = (ColumnVector a){
            for(int i=0; i<n; i++){
                v[i] = a.v[i];
            }
        }

        ColumnVector operator + (const ColumnVector b){
            ColumnVector res(n);
            for(int i=0; i<n; i++){
                res.v[i] = v[i] + b.v[i];
            }
            return res;
        }

        ColumnVector operator - (const ColumnVector b){
            ColumnVector res(n);
            for(int i=0; i<n; i++){
                res.v[i] = v[i] - b.v[i];
            }
            return res;
        }
};

class Matrix{
    public:
        int m, n;
        vector<vector<double>> v;
        Matrix(int a, int b){
            m = a;
            n = b;
            v = vector<vector<double>>(m, vector<double>(n, 0));
        }

        void operator = (Matrix a){
            for(int i=0; i<m; i++){
                for(int j=0; j<n; j++){
                    v[i][j] = a.v[i][j];
                }
            }
        }

        void operator = (ColumnVector a){
            for(int i=0; i<m; i++){
                v[i][0] = a.v[i];
            }
        }

        Matrix operator + (const Matrix b){
            Matrix res(m, b.n);
            for(int i=0; i<m; i++){
                for(int j=0; j<n; j++){
                    res.v[i][j] = v[i][j] + b.v[i][j];
                }
            }
            return res;
        }

        Matrix operator - (const Matrix b){
            Matrix res(m, b.n);
            for(int i=0; i<m; i++){
                for(int j=0; j<n; j++){
                    res.v[i][j] = v[i][j]-b.v[i][j];
                }
            }
            return res;
        }

        Matrix operator * (const Matrix b){
            Matrix res(m, b.n);
            for(int i=0; i<res.m; i++){
                for(int j=0; j<b.n; j++){
                    double sum = 0;
                    for(int k=0; k<n; k++){
                        sum = sum + (v[i][k]*b.v[k][j]);
                    }
                    if(fabs(sum)<=1e-4) sum = 0;
                    res.v[i][j] = sum;
                }
            }
            return res;
        }

        ColumnVector operator * (const ColumnVector b){
            ColumnVector res(m);
            for(int i=0; i<m; i++){
                double sum = 0;
                for(int j=0; j<n; j++){
                    sum = sum + (v[i][j]*b.v[j]);
                }
                if(fabs(sum)<=1e-4) sum = 0;
                res.v[i] = sum;
            }
            return res;
        }
        
        Matrix T(){
            Matrix res(n, m);
            for(int i=0; i<m; i++){
                for(int j=0; j<n; j++){
                    res.v[j][i] = v[i][j];
                }
            }
            return res;
        }
};

class SquareMatrix : public Matrix{
    public:
        SquareMatrix(int a):Matrix(a, a){}
        double determinant(){
            double result = 1;
            for(int i=0; i<n; i++){
                result*=v[i][i];
            }
            return result;
        }
};

class IdentityMatrix : public SquareMatrix{
    public:
        IdentityMatrix(int a):SquareMatrix(a){
            for(int i=0; i<a; i++){
                v[i][i] = 1;
            }
        }
};

class EliminationMatrix : public IdentityMatrix{
    public:
        EliminationMatrix(int a, int x, int y, Matrix matrix):IdentityMatrix(a){
            if(matrix.v[y-1][y-1]!=0){
                v[x-1][y-1] = -(double(matrix.v[x-1][y-1]/matrix.v[y-1][y-1]));
            }
        }
};

class PermutationMatrix : public IdentityMatrix{
    public:
        PermutationMatrix(int a, int r1, int r2):IdentityMatrix(a){
            v[r1-1][r1-1] = 0;
            v[r1-1][r2-1] = 1;
            v[r2-1][r2-1] = 0;
            v[r2-1][r1-1] = 1;
        }
};

void inverse(SquareMatrix* A){
    int n = A->n;
    Matrix* I = new IdentityMatrix(n);
    for(int i=1; i<n; i++){
        int maxPivot = i-1;
        for(int j=i-1; j<n; j++){
            if(fabs(A->v[j][i-1])>fabs(A->v[maxPivot][i-1])){
                maxPivot = j;
            }
        }
        if(maxPivot+1!=i){
            Matrix* P = new PermutationMatrix(n, maxPivot+1, i);
            *(Matrix*)A = *P * *A;
            (*I) = *P* *I;
        }
        for(int j=i+1; j<=n; j++){
            if(A->v[j-1][i-1]!=0){
                Matrix* E = new EliminationMatrix(n, j, i, *((Matrix*)A));
                *(Matrix*)A = *E * *A;
                (*I) = *E* *I;
            }
        }
    }
    for(int i=n; i>1; i--){
        for(int j=i-1; j>=1; j--){
            if(A->v[j-1][i-1]!=0){
                double div = -(A->v[j-1][i-1]/A->v[i-1][i-1]);
                Matrix* temp = new IdentityMatrix(n);
                temp->v[j-1][i-1] = div;
                *(Matrix*)A = *temp * *A;
                (*I) = *temp* *I;
            }
        }
    }
    for(int i=0; i<n; i++){
        Matrix* temp = new IdentityMatrix(n);
        if(A->v[i][i]!=0){
            temp->v[i][i]=(1/A->v[i][i]);
            *(Matrix*)A = *temp * *A;
            (*I) = *temp* *I;
        }
    }
    *(Matrix*)A = *I;
}


ostream& operator<<(ostream& output, Matrix matrix){
    for(int i=0; i<matrix.m; i++){
        for(int j=0; j<matrix.n; j++){
            output << matrix.v[i][j] << " ";
        }
        output << endl;
    }
    return output;
}

istream& operator>>(istream& input, Matrix& matrix){
    for(int i=0; i<matrix.m; i++){
        for(int j=0; j<matrix.n; j++){
            input >> matrix.v[i][j];
        }
    }
    return input;
}

ostream& operator<<(ostream& output, ColumnVector columnVector){
    for(int i=0; i<columnVector.n; i++){
        output << columnVector.v[i] << endl;
    }
    return output;
}

istream& operator>>(istream& input, ColumnVector& columnVector){
    for(int i=0; i<columnVector.n; i++){
        input >> columnVector.v[i];
    }
    return input;
}

void error(){
    cout << "Error: the dimensional problem occurred" << endl;
}

ColumnVector LeastSquareApproximation(int n, int deg, vector<double> ts, vector<double> bs){
    Matrix* A = new Matrix(n, deg+1);
    ColumnVector B(n);
    for(int i=0; i<deg+1; i++){
        for(int j=0; j<n; j++){
            if(i==0) A->v[j][i]=1;
            else if(i==1) A->v[j][i] = ts[j];
            else{
                A->v[j][i] = pow(A->v[j][1], i);
            }
        }
    }
    for(int i=0; i<n; i++){
        B.v[i] = bs[i];
    }
    Matrix* ATA = new Matrix(deg+1, deg+1);
    *ATA = A->T()* *A;
    Matrix* Ainv = new SquareMatrix(deg+1);
    (*Ainv) = *ATA;
    inverse((SquareMatrix*)Ainv);
    ColumnVector ATb(deg+1);
    ATb = A->T()*B;
    ColumnVector x(deg+1);
    x = (*Ainv)*ATb;
    return x;
}

// The template for GNUplot program was taken for Professor Ivan Konyukhov's tutorial
// The link for the tutorial: https://web.microsoftstream.com/video/9969a851-2920-4dca-9fb5-d775750cf1bc

#ifdef WIN32
#define GNUPLOT_NAME "c:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

int main(){
    #ifdef WIN32
        FILE* pipe = _popen(GNUPLOT_NAME, "w");
    #else
        FILE* pipe = _popen(GNUPLOT_NAME, "w");
    #endif

    if(pipe!=NULL){
        // Set of points
        vector<double> ts = {-4, -2, 0, 2};
        vector<double> bs = {0, 3.7, -1, 8};

        // Degree of the polynomial (3rd degree in this case)
        int deg = 3;

        // Number of points to be plotted and the x interval
        const double nbpoints = ts.size()*25;
        const double step = ((ts[ts.size()-1])-ts[0]) / nbpoints;

        // Solve it using LeastSquareApproximation
        ColumnVector solution = LeastSquareApproximation(4, deg, ts, bs);

        // Get the coefficients of the curve
        vector<double> coefficients = solution.v;


        fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'Least Square Approximation' with lines");
    
        // For each point calculate the x and y coordinates and plot it
        for(int i=0; i<nbpoints+1; i++){
            double x = ts[0] + i * step;
            double y = 0;
            // Go over the coefficients of and calculate the value of y
            for(int j = int(coefficients.size())-1; j>=0; j--){
                y = y + coefficients[j]*pow(x, j);
            }
            fprintf(pipe, "%f\t%f\n", x, y);
        }

        fflush(pipe);

        #ifdef WIN32
            _pclose(pipe);
        #else
            pclose(pipe);
        #endif
    }
    else{
        cout << "Could not open pipe" << endl;
    }

    return 0;
}
