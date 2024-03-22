void matrix_product(int n, int m, int k, double* a, double* b, double* c) {
    
    for (int i = 0; i < n*k; ++i) {
        c[i] = 0;
    }
    
    double sum00, sum01, sum02, sum10, sum11, sum12, sum20, sum21, sum22;
    
    int v3 = n%3;
    int h3 = k%3;
    
    for (int i = 0; i < v3; ++i) {
        for (int j = 0; j < h3; ++j) {
            sum00 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
            }
            
            c[k*i + j] = sum00;
        }
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                double factor = a[m*i + p];
                sum00 += factor * b[k*p + j];
                sum01 += factor * b[k*p + j + 1];
                sum02 += factor * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01;
            c[k*i + j + 2] = sum02;
        }
    }
    
    for (int i = v3; i < n; i+=3) {   
        for (int j = 0; j < h3; ++j) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                double factor = b[k*p + j];
                sum00 += a[m*i + p] * factor;
                sum01 += a[m*(i + 1) + p] * factor;
                sum02 += a[m*(i + 2) + p] * factor;
            }
            c[k*i + j] = sum00;
            c[k*(i + 1) + j] = sum01;
            c[k*(i + 2) + j] = sum02;
        }
        
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0; sum10 = 0; sum11 = 0; sum12 = 0;
            sum20 = 0; sum21 = 0; sum22 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*i + p] * b[k*p + j + 1];
                sum02 += a[m*i + p] * b[k*p + j + 2];
                sum10 += a[m*(i + 1) + p] * b[k*p + j];
                sum11 += a[m*(i + 1) + p] * b[k*p + j + 1];
                sum12 += a[m*(i + 1) + p] * b[k*p + j + 2];
                sum20 += a[m*(i + 2) + p] * b[k*p + j];
                sum21 += a[m*(i + 2) + p] * b[k*p + j + 1];
                sum22 += a[m*(i + 2) + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01; 
            c[k*i + j + 2] = sum02;
            c[k*(i + 1) + j] = sum10;
            c[k*(i + 1)  + j + 1] = sum11;
            c[k*(i + 1) + j + 2] = sum12;
            c[k*(i + 2) + j] = sum20;
            c[k*(i + 2) + j + 1] = sum21;
            c[k*(i + 2) + j + 2] = sum22;
        }
    }    
}
