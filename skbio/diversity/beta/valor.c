/*
 Valor - funciones de phylogenetic_diversity, unweighted_unifrac y unweighted_unifrac_all
 
 Avila Cortes Karina 
 Ramirez Garcia Juan Carlos
 
 - En este modulo se obtiene el valor de Pdiv                          
 - Entradas-  tabla de especies, arol
 - Salidas:  resultados de Pdiv
 
 - El programa corre de la siguiente forma: ./valor
 - Version: 15.2                                                                          
*/

//bibliotecas

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//estructura de arbol
typedef struct arbol
{
    float valor;
    char nombre;
}ARBOL;

int main (void)
{
    ARBOL ini[4][9],
    conjunto[4][9],
    ini_1[4][9],
    ini_2[4][9],
    ini_3[4][9],
    muestra_1[4][9],
    muestra_2[4][9],
    muestra_3[4][9];
    
    float conjunto_1  =  0;
    float conjunto_2  =  0;
    float conjunto_3  =  0;
    float conjunto_7  =  0;
    float conjunto_8  =  0;
    float conjunto_9  =  0;
    
    float conjunto_4[9];
    float conjunto_5[9];
    float conjunto_6[9];

    float contador;
    
    float valor_1,valor_2,valor_3;
    float total[3][3];
    float conjunto_total[3][3];
    
    int f_arreglo[9][3];
    int a,b,c,d,e,f  =  0;
    
    //valores designados a los arreglos
    
    f_arreglo[0][0] = 1;
    f_arreglo[0][1] = 1;
    f_arreglo[0][2] = 5;
    
    f_arreglo[1][0] = 1;
    f_arreglo[1][1] = 2;
    f_arreglo[1][2] = 0;
    
    f_arreglo[2][0] = 3;
    f_arreglo[2][1] = 1;
    f_arreglo[2][2] = 0;
    
    f_arreglo[3][0] = 0;
    f_arreglo[3][1] = 2;
    f_arreglo[3][2] = 0;
    
    f_arreglo[4][0] = 0;
    f_arreglo[4][1] = 0;
    f_arreglo[4][2] = 0;
    
    f_arreglo[5][0] = 0;
    f_arreglo[5][1] = 0;
    f_arreglo[5][2] = 3;
    
    f_arreglo[6][0] = 0;
    f_arreglo[6][1] = 0;
    f_arreglo[6][2] = 1;
    
    f_arreglo[7][0] = 1;
    f_arreglo[7][1] = 0;
    f_arreglo[7][2] = 1;
    
    f_arreglo[8][0] = 0;
    f_arreglo[8][1] = 0;
    f_arreglo[8][2] = 0;

    //valores designados al arbol
    conjunto[0][0].valor  =  0.35;
    conjunto[0][1].valor  =  0.35;
    conjunto[0][2].valor  =  0.35;
    conjunto[0][3].valor  =  0.35;
    conjunto[0][4].valor  =  0.35;
    
    conjunto[0][5].valor  =  0.2;
    conjunto[0][6].valor  =  0.2;
    conjunto[0][7].valor  =  0.2;
    conjunto[0][8].valor  =  0.2;
    
    conjunto[1][0].valor  =  0.3;
    conjunto[1][1].valor  =  0.3;
    conjunto[1][2].valor  =  0.3;
    conjunto[1][3].valor  =  0.3;
    conjunto[1][4].valor  =  0.3;
    
    conjunto[1][5].valor  =  0.3;
    conjunto[1][6].valor  =  0.3;
    conjunto[1][7].valor  =  0.7;
    conjunto[1][8].valor  =  0.7;
    
    conjunto[2][0].valor  =  0.2;
    conjunto[2][1].valor  =  0.3;
    conjunto[2][2].valor  =  0.2;
    conjunto[2][3].valor  =  0.2;
    conjunto[2][4].valor  =  0.9;
    
    conjunto[2][5].valor  =  0.2;
    conjunto[2][6].valor  =  0.3;
    conjunto[2][7].valor  =  0.3;
    conjunto[2][8].valor  =  0.4;
    
    conjunto[3][0].valor  =  0;
    conjunto[3][1].valor  =  0;
    conjunto[3][2].valor  =  0.5;
    conjunto[3][3].valor  =  0.3;
    conjunto[3][4].valor  =  0;
    
    conjunto[3][5].valor  =  0;
    conjunto[3][6].valor  =  0;
    conjunto[3][7].valor  =  0;
    conjunto[3][8].valor  =  0;
    

    for(a = 0; a < 4; a++)
    {
        for(b = 0; b < 9; b++)
        {
            ini[a][b].valor  =  0;
            ini[a][b].nombre  =  '\0';
            
            muestra_1[a][b].valor  =  0;
            muestra_1[a][b].nombre  =  '\0';
            
            muestra_2[a][b].valor  =  0;
            muestra_2[a][b].nombre  =  '\0';
        
            muestra_3[a][b].valor  =  0;
            muestra_3[a][b].nombre  =  '\0';
        }
    }
        
    //arbitrariamente se asignan nombres al arbol
    
    conjunto[1][0].nombre  =  'a';
    conjunto[1][1].nombre  =  'a';
    conjunto[1][2].nombre  =  'c';
    conjunto[1][3].nombre  =  'c';
    conjunto[1][4].nombre  =  'c';
    conjunto[1][5].nombre  =  'f';
    conjunto[1][6].nombre  =  'f';
    conjunto[1][7].nombre  =  'h';
    conjunto[1][8].nombre  =  'h';
    
    conjunto[2][0].nombre  =  'a';
    conjunto[2][1].nombre  =  'b';
    conjunto[2][2].nombre  =  'c';
    conjunto[2][3].nombre  =  'c';
    conjunto[2][4].nombre  =  'e';
    conjunto[2][5].nombre  =  'f';
    conjunto[2][6].nombre  =  'g';
    conjunto[2][7].nombre  =  'h';
    conjunto[2][8].nombre  =  'i';
    
    conjunto[3][0].nombre  =  '\0';
    conjunto[3][1].nombre  =  '\0';
    conjunto[3][2].nombre  =  'c';
    conjunto[3][3].nombre  =  'd';
    conjunto[3][4].nombre  =  '\0';
    conjunto[3][5].nombre  =  '\0';
    conjunto[3][6].nombre  =  '\0';
    conjunto[3][7].nombre  =  '\0';
    conjunto[3][8].nombre  =  '\0';
    
    conjunto[0][0].nombre  =  'a';
    conjunto[0][1].nombre  =  'a';
    conjunto[0][2].nombre  =  'a';
    conjunto[0][3].nombre  =  'a';
    conjunto[0][4].nombre  =  'a';
    
    conjunto[0][5].nombre  =  'f';
    conjunto[0][6].nombre  =  'f';
    conjunto[0][7].nombre  =  'f';
    conjunto[0][8].nombre  =  'f';

    for(b = 0; b < 10; b++)
    {
        conjunto_4[b] = 0;
        conjunto_5[b] = 0;
        conjunto_6[b] = 0;
    }
    
    //MUESTRA_1
    
    a = 0;
    b = 0;
    c = 0;
    e = 0;
    
    //obtiene valores de las muestras - arbol
    
    while(c < 9)
    {
        if (f_arreglo[c][a] !=  0)
        {
            if(c == 0)
            {
                for(d = 0; d < 4; d++)
                {
                    ini[d][b].valor = conjunto[d][b].valor;
                    ini[d][b].nombre = conjunto[d][b].nombre;
                    
                    if (ini[d][b].valor > 0)
                    {
                        muestra_1[d][b].valor = conjunto[d][b].valor;
                        muestra_1[d][b].nombre = conjunto[d][b].nombre;
                        
                        conjunto_1  =  muestra_1[d][b].valor + conjunto_1;
                    }
                }
            }
            else
            {
                f = 0;
                //valores de arbol
                for(d = 0; d < 4; d++)
                {
                    ini[d][b].valor = conjunto[d][b].valor;
                    ini[d][b].nombre = conjunto[d][b].nombre;
                    
                    if(ini[f][e].nombre !=  conjunto[d][b].nombre)
                    {
                        if (ini[d][b].valor > 0)
                        {
                            muestra_1[d][b].valor = conjunto[d][b].valor;
                            muestra_1[d][b].nombre = conjunto[d][b].nombre;
                            conjunto_1  =  muestra_1[d][b].valor + conjunto_1;
                        }
                    }
                    f++;
                }
                e = b;
            }
            conjunto_4[c] = 1; //verfica partes del arbol
        }
        
        b++;
        c++;
    }
    
    for(a = 0; a < 4; a++)
    {
        for(b = 0; b < 9; b++)
        {
            ini[a][b].valor  =  0;
            ini[a][b].nombre  =  '\0';
        }
    }
    
    //MUESTRA_2
    a = 1;
    b = 0;
    c = 0;
    e = 0;
    
    while(c < 9)
    { 
        if (f_arreglo[c][a] != 0)
        {
            conjunto_5[c] = 1;
            if(c == 0)
            {
                for(d = 0; d<4; d++)
                {
                    ini[d][b].valor = conjunto[d][b].valor;
                    ini[d][b].nombre = conjunto[d][b].nombre;
                                                
                    if (ini[d][b].valor > 0)
                    {
                        muestra_2[d][b].valor = conjunto[d][b].valor;
                        muestra_2[d][b].nombre = conjunto[d][b].nombre;
                        
                        conjunto_2  =  muestra_2[d][b].valor + conjunto_2;
                    }
                }
            }
            
            else
            {
                f = 0;
                for(d = 0; d < 4; d++)
                {
                    ini[d][b].valor = conjunto[d][b].valor;
                    ini[d][b].nombre = conjunto[d][b].nombre;
                    
                    if(ini[f][e].nombre !=  conjunto[d][b].nombre)
                    {
                        if (ini[d][b].valor > 0)
                        {
                            muestra_2[d][b].valor = conjunto[d][b].valor;
                            muestra_2[d][b].nombre = conjunto[d][b].nombre;
                            
                            conjunto_2  =  muestra_2[d][b].valor + conjunto_2;
                        }
                    }
                    
                    f++;
                }
                e = b;
            }
        }
        b++;
        c++;
    }
    
    for(a = 0; a < 4; a++)
    {
        for(b = 0; b < 9; b++)
        {
            ini[a][b].nombre  =  '\0';
            ini[a][b].valor  =  0;
        }
    }
    
    //MUESTRA_3
    a = 2;
    b = 0;
    c = 0;
    e = 0;
    
    while(c < 9)
    {
        if (f_arreglo[c][a] != 0)
        {
            conjunto_6[c] = 1;
            if(c == 0)
            {
                for(d = 0; d < 4;d++)
                {
                    ini[d][b].valor = conjunto[d][b].valor;
                    ini[d][b].nombre = conjunto[d][b].nombre;
                    if (ini[d][b].valor > 0)
                    {
                        muestra_3[d][b].valor = conjunto[d][b].valor;
                        muestra_3[d][b].nombre = conjunto[d][b].nombre;
                        
                        conjunto_3  =  muestra_3[d][b].valor + conjunto_3;
                    }
                }
            }
            else
            {
                f = 0;
                for(d = 0; d < 4; d++)
                {
                    ini[d][b].valor = conjunto[d][b].valor;
                    ini[d][b].nombre = conjunto[d][b].nombre;
                    
                    if(ini[f][e].nombre !=  conjunto[d][b].nombre)
                    {
                        if (ini[d][b].valor > 0)
                        {
                            muestra_3[d][b].valor = conjunto[d][b].valor;
                            muestra_3[d][b].nombre = conjunto[d][b].nombre;
                            
                            conjunto_3  =  muestra_3[d][b].valor + conjunto_3;
                        }
                    }
                    f++;
                }
                e = b;
            }
        }
        b++;
        c++;
    }
       
    //separa muestras
    
    for(a = 0; a < 4; a++)
    {
        for(b = 0; b < 9; b++)
        {
            ini_1[a][b].valor  =  0;
            ini_1[a][b].nombre  =  '\0';
            
            ini_2[a][b].valor  =  0;
            ini_2[a][b].nombre  =  '\0';
            
            ini_3[a][b].valor  =  0;
            ini_3[a][b].nombre  =  '\0';
        }
    }
    
    //inicia con "nodos" de muestra_1 y muestra_2, luego separa las muestras
    for(a = 3; a >= 0; a--)
    {
        for(b = 8; b >= 0; b--)
        {
            if(muestra_1[a][b].nombre !=  '\0')
            {
                ini_1[a][b].valor = muestra_1[a][b].valor;
                ini_1[a][b].nombre = muestra_1[a][b].nombre;
            }
            
            else if(muestra_2[a][b].nombre !=  '\0')
            {
                ini_1[a][b].valor = muestra_2[a][b].valor;
                ini_1[a][b].nombre = muestra_2[a][b].nombre;
            }
        }
    }
    
    for(a = 0; a < 9; a++)
    {
        if(conjunto_4[a] != conjunto_5[a])
        {
            for(b = 0; b < 4; b++)
            {
                if(ini_1[b][a].nombre !=  '\0')
                {
                    for(e = 0; e < 9; e++)
                    {
                        if(ini_1[b][e].nombre  ==  ini_1[b][a].nombre && e != a)
                        {
                           ini_1[b][a].nombre = '\0';
                           ini_1[b][a].valor = 0;
                       }
                    }
                }
            }
        }
    }
    //identifica muestras_1 y muestras_3
    for(a = 3; a >= 0; a--)
    {
        for(b = 8; b >= 0; b--)
        {
            if(muestra_1[a][b].nombre !=  '\0')
            {
                ini_2[a][b].valor = muestra_1[a][b].valor;
                ini_2[a][b].nombre = muestra_1[a][b].nombre;
            }
            else if(muestra_3[a][b].nombre !=  '\0')
            {
                ini_2[a][b].valor = muestra_3[a][b].valor;
                ini_2[a][b].nombre = muestra_3[a][b].nombre;
            }
        }
    }
    
    for(a = 0; a < 9; a++)
    {
        if(conjunto_4[a] !=  conjunto_6[a])
        {
            for(b = 0; b < 4; b++)
            {
                if(ini_2[b][a].nombre !=  '\0')
                {
                    for(e = 0; e < 9; e++)
                    {
                        if(ini_2[b][e].nombre  ==  ini_2[b][a].nombre && e != a)
                        {
                            ini_2[b][a].nombre = '\0';
                            ini_2[b][a].valor = 0;
                        }
                    }
                }
            }
        }
    }
    
    //inicia con "nodos" de muestra_3 y muestra_1, luego separa las muestras
    for(a = 3; a >= 0; a--)
    {
        for(b = 8; b >= 0; b--)
        {
            if(muestra_2[a][b].nombre !=  '\0')
            {
                ini_3[a][b].valor = muestra_2[a][b].valor;
                ini_3[a][b].nombre = muestra_2[a][b].nombre;
            }
            else if(muestra_3[a][b].nombre !=  '\0')
            {
                ini_3[a][b].nombre = muestra_3[a][b].nombre;
                ini_3[a][b].valor = muestra_3[a][b].valor;
            }
        }
    }

    for(a = 0; a < 9; a++)
    {
        if(conjunto_5[a] != conjunto_6[a])
        {
            for(b = 0; b < 4; b++)
            {
                if(ini_3[b][a].nombre != '\0')
                {
                    for(e = 0; e < 9; e++)
                    {
                        if(ini_3[b][e].nombre  ==  ini_3[b][a].nombre && e != a)
                        {
                            ini_3[b][a].nombre = '\0';
                            ini_3[b][a].valor = 0;
                        }
                    }
                }
            }
        }
    }
    
    //obtendra valores de unweighted_unifrac
    conjunto_1 = 0; 
    for (a = 0; a < 9; a++)
    {
        if(conjunto_4[a] != conjunto_6[a])
        {
            for(b = 0; b < 4; b++)
            {
                conjunto_7 = conjunto_7 + ini_2[b][a].valor;
            }
        }
        else
        {
            for(b = 0; b < 4; b++)
            {
                conjunto_1 = conjunto_1 + ini_2[b][a].valor;
            }
        }
    }
    
    conjunto_1 = conjunto_1 + conjunto_7; //valor total de las muestras
    
    conjunto_2 = 0;
    
    for (a = 0; a < 9; a++)
    {
        if(conjunto_4[a] != conjunto_5[a])
        {
            for(b = 0; b < 4; b++)
            {
                conjunto_8 = conjunto_8 + ini_1[b][a].valor;
            }
        }
        
        else
        
        {
            for(b = 0; b < 4; b++)
            {
                conjunto_2 = conjunto_2 + ini_1[b][a].valor;
            }
        }
    }
    
    conjunto_2 = conjunto_2 + conjunto_8;
    conjunto_3 = 0;

    for (a = 0; a < 9; a++)
    {
        if(conjunto_5[a] != conjunto_6[a])
        {
            for(b = 0; b < 4; b++)
            {
                conjunto_9 = conjunto_9 + ini_3[b][a].valor;
            }
        }
        else
        {
            for(b = 0; b < 4; b++)
            {
                conjunto_3 = conjunto_3 + ini_3[b][a].valor;
            }
        }
    }
    
    conjunto_3 = conjunto_3 + conjunto_9;
    
    //Por ultimo toma los valores de unweighted_unifrac

    conjunto_total[2][0] = conjunto_7 / conjunto_1;
    conjunto_total[2][1] = conjunto_9 / conjunto_3;
    conjunto_total[1][0] = conjunto_8 / conjunto_2;
    
    printf("%f\n", conjunto_total[2][0]);
    printf("%f\n", conjunto_total[2][1]);
    printf("%f\n", conjunto_total[1][0]);
    
    // unweighted_unifrac_all
    for(a = 0; a < 3; a++)
    {
        for(b = 0; b < 3; b++)
        {
            total[a][b] = 0;
        }
    }
    
    //disenio de f_arreglo
    for(a = 0; b < 3; a++)
    {
        for(b = 0; b < a; b++)
        {
            total[a][b] = total[b][a] = conjunto_total[a][b];
        }
    }

//resultados
    for (a = 0; a < 3; a++)
    {
        for(b = 0; b < 3; b++)
        {
            printf("%f ",total[a][b]);
        }
        printf("\n");
    }

}
