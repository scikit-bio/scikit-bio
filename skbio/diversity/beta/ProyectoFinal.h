/*****************************************************************************************
* Function aplicaciones_final                                                            * 
* -This program obtain the value of the philogenetic diversity                           * 
* Syntax:                                                                                *
*  ./funcion1                                                                            *
* Inputs:                                                                                *
*  int arr[9][3] - table of the samples vs species observed                              *
*  BASE todo[4][9] - phylogenetic tree                                                   *
* Outputs:                                                                               *
*  float final[3][3] - array of the resultsets that philogenetic diversity               *
* Functions:                                                                             *
*  int aplicaciones_final() - this function realize the operation of the functions       *
*   phylogenetic_diversity, unweighted_unifrac and unweighted_unifrac_all. The functions *
*   are implied into aplicaciones_final                                                  *
* Structures:                                                                            *
*  struct base:                                                                          *
*   -char ident - identifies a different branches in the tree                            *
*   -float valor - value of the branches in the tree                                     *
* Authors:                                                                               *
*  Gustavo Adolfo Avendaño Martinez                                                      *
*  Amalinalli Lopez Mejia                                                                *
* Version: 10.5                                                                          *
* ***************************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//create structure for tree
typedef struct base{
    char ident;
    float valor;
}BASE;

int aplicaciones_final (){
    BASE todo[4][9],temp1[4][9],temp2[4][9],temp3[4][9],temp[4][9],sample1[4][9],sample2[4][9],sample3[4][9]; //representation of samples
    float suma1=0,suma2=0,suma3=0,cont,dif1[9],dif2[9],dif3[9],sum1=0,sum2=0,sum3=0,val1,val2,val3,final[3][3],cfinal[3][3]; //final matriz, sum of the values of shared nodes and no shared
    int arr[9][3]; //table of samples vs species
    int i,j,k,l,n,m=0;  //counters
    //array of table samples and species
    arr[0][0]=1;
    arr[0][1]=1;
    arr[0][2]=5;
    arr[1][0]=1;
    arr[1][1]=2;
    arr[1][2]=0;
    arr[2][0]=3;
    arr[2][1]=1;
    arr[2][2]=0;
    arr[3][0]=0;
    arr[3][1]=2;
    arr[3][2]=0;
    arr[4][0]=0;
    arr[4][1]=0;
    arr[4][2]=0;
    arr[5][0]=0;
    arr[5][1]=0;
    arr[5][2]=3;
    arr[6][0]=0;
    arr[6][1]=0;
    arr[6][2]=1;
    arr[7][0]=1;
    arr[7][1]=0;
    arr[7][2]=1;
    arr[8][0]=0;
    arr[8][1]=0;
    arr[8][2]=0;

    //array with values of tree's branches
    todo[0][0].valor = 0.35;
    todo[0][1].valor = 0.35;
    todo[0][2].valor = 0.35;
    todo[0][3].valor = 0.35;
    todo[0][4].valor = 0.35;
    todo[0][5].valor = 0.2;
    todo[0][6].valor = 0.2;
    todo[0][7].valor = 0.2;
    todo[0][8].valor = 0.2;
    
    todo[1][0].valor = 0.3;
    todo[1][1].valor = 0.3;
    todo[1][2].valor = 0.3;
    todo[1][3].valor = 0.3;
    todo[1][4].valor = 0.3;
    todo[1][5].valor = 0.3;
    todo[1][6].valor = 0.3;
    todo[1][7].valor = 0.7;
    todo[1][8].valor = 0.7;
    
    todo[2][0].valor = 0.2;
    todo[2][1].valor = 0.3;
    todo[2][2].valor = 0.2;
    todo[2][3].valor = 0.2;
    todo[2][4].valor = 0.9;
    todo[2][5].valor = 0.2;
    todo[2][6].valor = 0.3;
    todo[2][7].valor = 0.3;
    todo[2][8].valor = 0.4;
    
    todo[3][0].valor = 0;
    todo[3][1].valor = 0;
    todo[3][2].valor = 0.5;
    todo[3][3].valor = 0.3;
    todo[3][4].valor = 0;
    todo[3][5].valor = 0;
    todo[3][6].valor = 0;
    todo[3][7].valor = 0;
    todo[3][8].valor = 0;
    
    //initialize array of samples
    for(i=0;i<4;i++){
        for(j=0;j<9;j++){
            sample1[i][j].valor = 0;
            sample1[i][j].ident = '\0';
            sample2[i][j].valor = 0;
            sample2[i][j].ident = '\0';
            sample3[i][j].valor = 0;
            sample3[i][j].ident = '\0';
            temp[i][j].valor = 0;
            temp[i][j].ident = '\0';
        }
    }
        
    //array with identifiers for the tree's branches
    //each letter represents a different branches of the tree
    todo[0][0].ident = 'a';
    todo[0][1].ident = 'a';
    todo[0][2].ident = 'a';
    todo[0][3].ident = 'a';
    todo[0][4].ident = 'a';
    todo[0][5].ident = 'f';
    todo[0][6].ident = 'f';
    todo[0][7].ident = 'f';
    todo[0][8].ident = 'f';
    todo[1][0].ident = 'a';
    todo[1][1].ident = 'a';
    todo[1][2].ident = 'c';
    todo[1][3].ident = 'c';
    todo[1][4].ident = 'c';
    todo[1][5].ident = 'f';
    todo[1][6].ident = 'f';
    todo[1][7].ident = 'h';
    todo[1][8].ident = 'h';
    todo[2][0].ident = 'a';
    todo[2][1].ident = 'b';
    todo[2][2].ident = 'c';
    todo[2][3].ident = 'c';
    todo[2][4].ident = 'e';
    todo[2][5].ident = 'f';
    todo[2][6].ident = 'g';
    todo[2][7].ident = 'h';
    todo[2][8].ident = 'i';
    todo[3][0].ident = '\0';
    todo[3][1].ident = '\0';
    todo[3][2].ident = 'c';
    todo[3][3].ident = 'd';
    todo[3][4].ident = '\0';
    todo[3][5].ident = '\0';
    todo[3][6].ident = '\0';
    todo[3][7].ident = '\0';
    todo[3][8].ident = '\0';

    
    //inicialize array of dif
    //dif is a array that representing the different of shared branches or no shared
    for(j=0;j<10;j++){
           dif1[j]=0;
           dif2[j]=0;
           dif3[j]=0;
        }
    
    //Sample 1
    k=0;
    i=0;
    j=0;
    //this loop obtain the values of the sample1 according with the tree and the table
    while(k<9){
        //check that the species has seen in this sample
        if (arr[k][i] != 0){
            //if is the first iteration
            if(k==0){
                //full the first column in the array of sample1
                for(l=0;l<4;l++){
                    //saved the repeated branches
                    temp[l][j].valor=todo[l][j].valor;
                    temp[l][j].ident=todo[l][j].ident;
                    if (temp[l][j].valor>0){
                        sample1[l][j].valor=todo[l][j].valor;
                        sample1[l][j].ident=todo[l][j].ident;
                        suma1 = sample1[l][j].valor + suma1;
                    }
                }
            //if is the second hereinafter iteration
            }else{
                n=0;
                //full the second hereinafter column in the array of sample1
                for(l=0;l<4;l++){
                    //saved the repeated branches
                    temp[l][j].valor=todo[l][j].valor;
                    temp[l][j].ident=todo[l][j].ident;
                    if(temp[n][m].ident != todo[l][j].ident){
                        if (temp[l][j].valor>0){
                            sample1[l][j].valor=todo[l][j].valor;
                            sample1[l][j].ident=todo[l][j].ident;
                            suma1 = sample1[l][j].valor + suma1;
                        }
                    }
                    n++;
                }
                m=j;
            }
            //identifies the branches that has had seen
            dif1[k]=1;
        }
        j++;
        k++;
    }
    //reebot temp
    for(i=0;i<4;i++){
        for(j=0;j<9;j++){
            temp[i][j].valor = 0;
            temp[i][j].ident = '\0';
        }
    }
     //Sample 2
    k=0;
    i=1;
    j=0;
    m=0;
    //this loop obtain the values of the sample2 according with the tree and the table
    while(k<9){
        //check that the species has seen in this sample
        if (arr[k][i]!= 0){
            //identifies the branches that has had seen
            dif2[k]=1;
            //if is the first iteration
            if(k==0){
                //full the first column in the array of sample2
                for(l=0;l<4;l++){
                    //saved the repeated branches
                    temp[l][j].valor=todo[l][j].valor;
                    temp[l][j].ident=todo[l][j].ident;
                    if (temp[l][j].valor>0){
                        sample2[l][j].valor=todo[l][j].valor;
                        sample2[l][j].ident=todo[l][j].ident;
                        suma2 = sample2[l][j].valor + suma2;
                    }
                }
            //if is the second hereinafter iteration
            }else{
                n=0;
                //full the second hereinafter column in the array of sample2
                for(l=0;l<4;l++){
                    //saved the repeated branches
                    temp[l][j].valor=todo[l][j].valor;
                    temp[l][j].ident=todo[l][j].ident;
                    if(temp[n][m].ident != todo[l][j].ident){
                        if (temp[l][j].valor>0){
                            sample2[l][j].valor=todo[l][j].valor;
                            sample2[l][j].ident=todo[l][j].ident;
                            suma2 = sample2[l][j].valor + suma2;
                        }
                    }
                    n++;
                }
                m=j;
            }
        }
        j++;
        k++;
    }
    
    //reebot temp
    for(i=0;i<4;i++){
        for(j=0;j<9;j++){
            temp[i][j].valor = 0;
            temp[i][j].ident = '\0';
        }
    }
    
    //Sample 3
    k=0;
    i=2;
    j=0;
    m=0;
    //this loop obtain the values of the sample3 according with the tree and the table
    while(k<9){
        //check that the species has seen in this sample
        if (arr[k][i]!= 0){
            //identifies the branches that has had seen
            dif3[k]=1;
            //if is the first iteration
            if(k==0){
                //full the first column in the array of sample3
                for(l=0;l<4;l++){
                    //saved the repeated branches
                    temp[l][j].valor=todo[l][j].valor;
                    temp[l][j].ident=todo[l][j].ident;
                    if (temp[l][j].valor>0){
                        sample3[l][j].valor=todo[l][j].valor;
                        sample3[l][j].ident=todo[l][j].ident;
                        suma3 = sample3[l][j].valor + suma3;
                    }
                }
            //if is the second hereinafter iteration
            }else{
                n=0;
                //full the second hereinafter column in the array of sample3
                for(l=0;l<4;l++){
                    //saved the repeated branches
                    temp[l][j].valor=todo[l][j].valor;
                    temp[l][j].ident=todo[l][j].ident;
                    if(temp[n][m].ident != todo[l][j].ident){
                        if (temp[l][j].valor>0){
                            sample3[l][j].valor=todo[l][j].valor;
                            sample3[l][j].ident=todo[l][j].ident;
                            suma3 = sample3[l][j].valor + suma3;
                        }
                    }
                    n++;
                }
                m=j;
            }
        }
        j++;
        k++;
    }
       
    //combines arrangements of sample1, sample2 and sample3 for show the shared branches and not shared   
  
    //reebot temp
    for(i=0;i<4;i++){
        for(j=0;j<9;j++){
            temp1[i][j].valor = 0;
            temp1[i][j].ident = '\0';
            temp2[i][j].valor = 0;
            temp2[i][j].ident = '\0';
            temp3[i][j].valor = 0;
            temp3[i][j].ident = '\0';
        }
    }
    
    //combines arrangements of sample1 and sample2
    for(i=3;i>=0;i--){
        for(j=8;j>=0;j--){
            if(sample1[i][j].ident != '\0'){
                temp1[i][j].valor=sample1[i][j].valor;
                temp1[i][j].ident=sample1[i][j].ident;
            }else if(sample2[i][j].ident != '\0'){
                temp1[i][j].valor=sample2[i][j].valor;
                temp1[i][j].ident=sample2[i][j].ident;
            }
        }
    }
    
    //clean the new array of sample1sample2 for see the shared branches and not shared
    for(i=0;i<9;i++){
        //if exists one branches unique in a sample, clean this branches with shared branches
        if(dif1[i]!=dif2[i]){
            //check all nivel of the tree
            for(j=0;j<4;j++){
                //if a nodes exists
                if(temp1[j][i].ident != '\0'){
                    //compared all the nodes of the same nivel and clean the not shared branches
                    for(m=0;m<9;m++){
                        //jump the same node and check the shared branches
                        if(temp1[j][m].ident == temp1[j][i].ident && m!=i){
                            temp1[j][i].valor=0;
                            temp1[j][i].ident='\0';
                        }
                    }
                }
            }
        }
    }
    
    //combines arrangements of sample1 and sample3
    for(i=3;i>=0;i--){
        for(j=8;j>=0;j--){
            if(sample1[i][j].ident != '\0'){
                temp2[i][j].valor=sample1[i][j].valor;
                temp2[i][j].ident=sample1[i][j].ident;
            }else if(sample3[i][j].ident != '\0'){
                temp2[i][j].valor=sample3[i][j].valor;
                temp2[i][j].ident=sample3[i][j].ident;
            }
        }
    }
    
    
    //clean the new array of sample1sample3 for see the shared branches and not shared
    for(i=0;i<9;i++){
        //if exists one branches unique in a sample, clean this branches with shared branches
        if(dif1[i]!= dif3[i]){
            //check all nivel of the tree
            for(j=0;j<4;j++){
                //if a nodes exists
                if(temp2[j][i].ident != '\0'){
                    //compared all the nodes of the same nivel and clean the not shared branches
                    for(m=0;m<9;m++){
                        //jump the same node and check the shared branches
                        if(temp2[j][m].ident == temp2[j][i].ident && m!=i){
                            temp2[j][i].valor=0;
                            temp2[j][i].ident='\0';
                        }
                    }
                }
            }
        }
    }
    
  
    //combines arrangements of sample2 and sample3
    for(i=3;i>=0;i--){
        for(j=8;j>=0;j--){
            if(sample2[i][j].ident != '\0'){
                temp3[i][j].valor=sample2[i][j].valor;
                temp3[i][j].ident=sample2[i][j].ident;
            }else if(sample3[i][j].ident != '\0'){
                temp3[i][j].valor=sample3[i][j].valor;
                temp3[i][j].ident=sample3[i][j].ident;
            }
        }
    }
    
    //clean the new array of sample2sample3 for see the shared branches and not shared
    for(i=0;i<9;i++){
        //if exists one branches unique in a sample, clean this branches with shared branches
        if(dif2[i]!=dif3[i]){
            //check all nivel of the tree
            for(j=0;j<4;j++){
                //if a nodes exists
                if(temp3[j][i].ident != '\0'){
                    //compared all the nodes of the same nivel and clean the not shared branches
                    for(m=0;m<9;m++){
                        //jump the same node and check the shared branches
                        if(temp3[j][m].ident == temp3[j][i].ident && m!=i){
                            temp3[j][i].valor=0;
                            temp3[j][i].ident='\0';
                        }
                    }
                }
            }
        }
    }
    
    //obtain the value of the unifrac of sample1 and sample2
    suma1=0;
    //check the view of the sample1 and sample 2
    for (i=0;i<9;i++){
        //check the not shared branches
        if(dif1[i]!=dif3[i]){
            //sum all the values of the species
            for(j=0;j<4;j++){
                sum1=sum1+temp2[j][i].valor;
            }
        //check the shared branches
        }else {
            //sum all the values of the species
            for(j=0;j<4;j++){
                suma1=suma1+temp2[j][i].valor;
            }
        }
    }
    //sum all the values of the observed nodes in sample1 and sample2
    suma1=suma1+sum1;
    
    //obtain the value of the unifrac of sample1 and sample3
    suma2=0;
    //check the view of the sample1 and sample 3
    for (i=0;i<9;i++){
        //check the not shared branches
        if(dif1[i]!=dif2[i]){
            //sum all the values of the species
            for(j=0;j<4;j++){
                sum2=sum2+temp1[j][i].valor;
            }
        //check the shared branches
        }else {
             //sum all the values of the species
            for(j=0;j<4;j++){
                suma2=suma2+temp1[j][i].valor;
            }
        }
    }
    //sum all the values of the observed nodes in sample1 and sample3
    suma2=suma2+sum2;
    
    //obtain the value of the unifrac of sample2 and sample3
    suma3=0;
    //check the view of the sample2 and sample 3
    for (i=0;i<9;i++){
        //check the not shared branches
        if(dif2[i]!=dif3[i]){
            //sum all the values of the species
            for(j=0;j<4;j++){
                sum3=sum3+temp3[j][i].valor;
            }
        //check the shared branches
        }else {
            //sum all the values of the species
            for(j=0;j<4;j++){
                suma3=suma3+temp3[j][i].valor;
            }
        }
    }
    //sum all the values of the observed nodes in sample2 and sample3
    suma3=suma3+sum3;
    
    //obtain unifrac value
    cfinal[1][0]=sum2/suma2;
    cfinal[2][0]=sum1/suma1;
    cfinal[2][1]=sum3/suma3;
    
    //unifrac_all
    //initialize final
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            final[i][j]=0;
        }
    }
    //accommodate the final array according with samples combined
    for(i=0;i<3;i++){
        for(j=0;j<i;j++){
            final[i][j]=final[j][i]=cfinal[i][j];
        }
    }
    //print the final array
    for (i=0;i<3;i++){
        for(j=0;j<3;j++){
            printf("%f ",final[i][j]);
        }
        printf("\n");
    }

}