float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,4);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
}

void calculateLambda(Matrix &C, mesh m, int i){

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    zeroes(C,4,3);
    C.at(0).at(0) = -2*x1+x3+x4;    C.at(0).at(1) = -2*x1+x2+x4;    C.at(0).at(2) = -2*x1+x2+x3; 
    C.at(1).at(0) = 2*x2+x3+x4;     C.at(1).at(1) = -x1+x3;         C.at(1).at(2) = -x1+x4; 
    C.at(2).at(0) = -x1+x2;         C.at(2).at(1) = x2+2*x3+x4;     C.at(2).at(2) = -x1+x4;
    C.at(3).at(0) = -x1+x2;         C.at(3).at(1) = -x1+x3;         C.at(3).at(2) = x2+x3+2*x4;

}

void calculateMu(Matrix &C, mesh m, int i){

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    zeroes(C,4,3);
    C.at(0).at(0) = -(10*x1+5*x3+5*x4-2*(3*pow(y1,2)+2*y1*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2)));      C.at(0).at(1) = -(10*x1+5*x2+5*x4-2*(3*pow(y1,2)+2*y1*(y2+y4)+pow(y2,2)+y2*y4+pow(y4,2)));      C.at(0).at(2) = -(10*x1+5*x2+5*x3-2*(3*pow(y1,2)+2*y1*(y2+y3)+pow(y2,2)+y2*y3+pow(y3,2))); 
    C.at(1).at(0) = 10*x2+5*x3+5*x4-2*(3*pow(y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2));         C.at(1).at(1) = -(5*x1-5*x3-2*(y1-y3)*(y1+2*y2+y3+y4));                                         C.at(1).at(2) = -(5*x1-5*x4-2*(y1-y4)*(y1+2*y2+y3+y4)); 
    C.at(2).at(0) = -(5*x1-5*x2-2*(y1-y2)*(y1+y2+2*y3+y4));                                         C.at(2).at(1) = 5*x2+10*x3+5*x4-2*(pow(y2,2)+y2*(2*y3+y4)+3*pow(y3,2)+2*y3*y4+pow(y4,2));       C.at(2).at(2) = -(5*x1-5*x4-2*(y1-y4)*(y1+y2+2*y3+y4));
    C.at(3).at(0) = -(5*x1-5*x2-2*(y1-y2)*(y1+y2+y3+2*y4));                                         C.at(3).at(1) = -(5*x1-5*x3-2*(y1-y3)*(y1+y2+y3+2*y4));                                         C.at(3).at(2) = 5*x2+5*x3+2*(5*x4-pow(y2,2)-y2*(y3+2*y4)-pow(y3,2)-2*y3*y4-3*pow(y4,2));

}

void calculateOmega(Matrix &C, mesh m, int i){

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    zeroes(C,4,3);
    C.at(0).at(0) = 2*x1+x2+x3+x4+2*y1+y2+y3+y4;    C.at(0).at(1) = 2*x1+x2+x3+x4+2*y1+y2+y3+y4;    C.at(0).at(2) = 2*x1+x2+x3+x4+2*y1+y2+y3+y4; 
    C.at(1).at(0) = x1+2*x2+x3+x4+y1+2*y2+y3+y4;    C.at(1).at(1) = x1+2*x2+x3+x4+y1+2*y2+y3+y4;    C.at(1).at(2) = x1+2*x2+x3+x4+y1+2*y2+y3+y4; 
    C.at(2).at(0) = x1+x2+2*x3+x4+y1+y2+2*y3+y4;    C.at(2).at(1) = x1+x2+2*x3+x4+y1+y2+2*y3+y4;    C.at(2).at(2) = x1+x2+2*x3+x4+y1+y2+2*y3+y4;
    C.at(3).at(0) = x1+x2+x3+2*x4+y1+y2+y3+2*y4;    C.at(3).at(1) = x1+x2+x3+2*x4+y1+y2+y3+2*y4;    C.at(3).at(2) = x1+x2+x3+2*x4+y1+y2+y3+2*y4;

}

void calculatePhi(Matrix &C, mesh m, int i){

    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    zeroes(C,4,3);
    C.at(0).at(0) = 3*pow(y1,2)+2*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2);     C.at(0).at(1) = 3*pow(y1,2)+2*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2);     C.at(0).at(2) = 3*pow(y1,2)+2*y1*(y2+y3+y4)+pow(y2,2)+y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2); 
    C.at(1).at(0) = pow(y1,2)+y1*(2*y2+y3+y4)+3*pow(y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2);   C.at(1).at(1) = pow(y1,2)+y1*(2*y2+y3+y4)+3*pow(y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2);   C.at(1).at(2) = pow(y1,2)+y1*(2*y2+y3+y4)+3*pow(y2,2)+2*y2*(y3+y4)+pow(y3,2)+y3*y4+pow(y4,2);
    C.at(2).at(0) = pow(y1,2)+y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+3*pow(y3,2)+2*y3*y4+pow(y4,2); C.at(2).at(1) = pow(y1,2)+y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+3*pow(y3,2)+2*y3*y4+pow(y4,2); C.at(2).at(2) = pow(y1,2)+y1*(y2+2*y3+y4)+pow(y2,2)+y2*(2*y3+y4)+3*pow(y3,2)+2*y3*y4+pow(y4,2);
    C.at(3).at(0) = pow(y1,2)+y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*pow(y4,2); C.at(3).at(1) = pow(y1,2)+y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*pow(y4,2); C.at(3).at(2) = pow(y1,2)+y1*(y2+y3+2*y4)+pow(y2,2)+y2*(y3+2*y4)+pow(y3,2)+2*y3*y4+3*pow(y4,2);

}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixB,matrixC,matrixK,matrixL;
    float swr,J,Determinant;
    Matrix lambda_matrix, mu_matrix, omega_matrix, phi_matrix, Alpha, Beta, Alpha_t,Beta_t;

    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
        
    
    /* [ B+C  K ]
       [  0   L ]
    */


    //Matrix K
    float real_k = (float) J/(24*pow(Determinant,2));
    calculateLambda(lambda_matrix, m, e);
    productRealMatrix(real_k,productMatrixMatrix(lambda_matrix,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,4),3,3,4),4,3,4),matrixK);

    //Matrix L
    float real_l = (float) (J/(120*pow(Determinant,2)));
    calculateMu(mu_matrix, m, e);
    productRealMatrix(real_l,productMatrixMatrix(mu_matrix,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,4),3,3,4),4,3,4),matrixL);
    
    //Matrix B
    float real_b = (float) (J/(120*Determinant));
    calculateOmega(omega_matrix, m, e);
    productRealMatrix(real_b,productMatrixMatrix(omega_matrix,productMatrixMatrix(Alpha,Beta,3,3,4),4,3,4),matrixB);

    //Matrix C
    float real_c = (float)(J/(360*Determinant));
    calculatePhi(phi_matrix, m, e);
    productRealMatrix(real_c,productMatrixMatrix(phi_matrix,productMatrixMatrix(Alpha,Beta,3,3,4),4,3,4),matrixC);

    //Matrix M
    Matrix M;
    zeroes(M,8);

    ubicarSubMatriz(M,0,3,0,3, sumMatrix(matrixB,matrixC,4,4));
    ubicarSubMatriz(M,0,3,4,7,matrixK);
    ubicarSubMatriz(M,4,7,4,7,matrixL);
    ubicarSubMatriz(M,4,7,0,3,matrixL);

    return M;
}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux;

    float swr = m.getParameter(SHOWERS_RATIO);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,8);
    b_aux.at(0)=1;
    b_aux.at(1)=1;
    b_aux.at(2)=1;
    b_aux.at(3)=1;
    productRealVector((swr*J)/24,b_aux,b);
    
    return b;
}
