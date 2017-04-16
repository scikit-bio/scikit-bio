
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//Identificadores de los Samples
char sample_ids[20] = {'A','B','C','D','E','F','G','H','I','J','K','L'};
char tree[100];
char separador[5] = "0";
char separador2[5] = "00";
char separador3[5] = "000";

//Dimensiones de la Matriz
int samples, observaciones;

//Indicadores
int i,j,opc,bandera = 0;

//Funcion Principal
void main()
{
    //Mensajes de Bienvenida
    printf("Phylogenetic Diversity (PD)\n");
    printf("\n");
    
    //Le Pedimos al Usuario las Dimensiones de la Tabla
    printf("Samples Deseados: ");
    scanf("%d",&samples);
    printf("Observaciones Deseadas (OTU): ");
    scanf("%d",&observaciones);
    printf("\n");
    
    //Creamos la Tabla
    double matriz [observaciones][samples];
    
    //Creamos el Arreglo de IDs de los Samples
    double OTU[samples];
    
    //Ciclo para Ingresar Valores a la Tabla (Por OTU)
    for (i = 0; i < observaciones; i++) {
        for (j = 0; j < samples; j++) {
            printf("Dato en Posicion [%d][%d]: ",i,j);
            scanf("%lf",&matriz[i][j]);
        }
    }
    
    printf("\n");
    
    //Imprimimos el Arreglo de Ids de los Samples
    for (j = 0; j < samples; j++) {
        printf(" %c    ",sample_ids[j]);
    }
    
    printf("\n");
    
    //Imprimimos la Tabla
    for (i = 0; i < observaciones; i++) {
        for (j = 0; j < samples; j++) {
             printf("[%1.1f]",matriz[i][j]);
            }
        printf("\n");
    }
    
     while(bandera != 1)
     {
         
         for (i = 0; i < samples; i++)
         {
             OTU[i] = 0;
         }
    
         //Opciones
         printf("\nOpciones: \n");
         printf("1. Observed_OTU\n");
         printf("2. Generar TreeNode\n");
         printf("3. Calcular Phylogenetic_Diversity\n");
         printf("4. Salir\n");
         printf("Opcion: ");
         scanf("%d",&opc);
    
         //Calcular observed_OTU
         if(opc == 1)
         {
        
             printf("\n");
        
             for (i = 0; i < samples; i++)
             {
                 for (j = 0; j < observaciones; j++)
                 {
                     if(matriz[j][i] != 0)
                     {
                         OTU[i]++;
                     }

                 }
                 printf("Observed OTU (Sample [%c]): [%1.1f]\n",sample_ids[i],OTU[i]);
             }
             
             bandera = 0;
         }
         
         //Generar TreeNode
         if(opc == 2)
         {
             
             printf("\n");
             printf("Ingresa el TreeNode en el Formato Newick (Formato Especial): ");
             scanf("%s",tree);
             printf("\n");
             printf("tree: %s\n",tree);
             bandera = 0;
         }
         
         //Calcular PD
         if(opc == 3)
         {
             
             for (i = 0; i < samples; i++)
             {
                 for (j = 0; j < observaciones; j++)
                 {
                     if(matriz[j][i] != 0)
                     {
                         OTU[i]++;
                     }
                     
                 }
             }
             
             printf("\nEn el sample [%c] hay [%1.1f] OTUs \n",sample_ids[0],OTU[0]);
             
             float num1;
             float num2;
             float num3;
             float PD1;
             char nums[20];
             char nums2[20];
             char nums3[20];
             
             strcpy(nums,strstr(tree,separador));
             strcpy(nums2,strstr(tree,separador2));
             strcpy(nums3,strstr(tree,separador3));
             
             num1 = atof(nums);
             num2 = atof(nums2);
             num3 = atof(nums3);
             
             PD1 = num1 + num2 + num3;
             
             printf( "\ntree = %s\n\n", tree);
             printf( "Ramas del nodo a la Raiz : %1.2f | %1.2f | %1.2f \n",num1,num2, num3);
             printf( "Phylogenetic Diversity: %1.3f", PD1);
             printf("\n");
             bandera = 0;
             
         }
         
         //Salir
         if(opc == 4)
         {
             
             printf("Saliendo...\n");
             bandera = 1;
         }
     }
}