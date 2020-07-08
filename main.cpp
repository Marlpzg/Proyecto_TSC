#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "display_tools.h"
#include "sel_tools.h"
#include "components.h"
#include "sel.h"
#include "assembly.h"

int main(int argc, char *argv[])
{
    char filename[150];
    strcpy(filename,argv[1]);

    vector<Matrix> localKs;
    vector<Vector> localbs;
    Matrix K;
    Vector b;
    Vector T;

    cout << "IMPLEMENTACION DEL METODO DE LOS ELEMENTOS FINITOS\n"
         << "ECUACIONES DE MARIO LOPEZ - 00046317\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallaYCondiciones(m,filename);
    
    crearSistemasLocales(m,localKs,localbs);
    
    zeroes(K,2*m.getSize(NODES));
    zeroes(b,2*m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);
    
    applyDirichlet(m,K,b);
    
    cout << "M Global: " << endl;
    showMatrix(K);
    cout << endl;

    cout << "b Global: " << endl;
    showVector(b);
    cout << endl;

    zeroes(T,b.size());
    calculate(K,b,T);

    cout << "La respuesta es: \n";
    showVector(T);

    writeResults(m,T,filename);

    return 0;
}
