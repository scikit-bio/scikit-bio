/**
 * @file scikit-bio.c
 *
 * @brief This project is able to calculate the Unweighted Unifrac between two samples, the Phylogenetic Diversity of a specific sample and Unweighted Unifrac of all possible combination of samples. The data is read from a file given by the user, also the user has to give the file name of the fail containing the tree information with Newick method.
 
 * @author Rafael Alvarez, Alonso De Cosio and Bruno Ramos Berho.
 * @date 5/12/2014
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

float unweighted_unifrac(char tree[], float **table, char **otus, char **samples, int sample1, int sample2);
float uniqueLenght(char otu[], char tree[]);
float sumAllLengths(char tree[]);
float phylogenetic_diversity(char tree[], float **table, char **otus, int sample1);
void rowsColumns(char fileName[],int *rows, int *columns);
void readFile(char fileName[],float **table,char **samples, char **otus, int columns, int rows);
int getSample(char sample[], char **samples, int columns);

/**
 *  main: This function is the root of the proyect. It creates the structure of the program. First it reads the information and create the necesary memory. The function will show a menu on screen where the user can select 4 different options, 1.Unweighted Unifrac, 2.Phylogenetic Diversity, 3.Unweighted Unifrac All, 4. Exit. At the end all the memory used is released.
 *
 *  @param argc The function doesn't recieve any arguments.
 *  @param argv The function doesn't recieve any arguments.
 *
 *  @return The function doesn't return any important data.
 */
int main(int argc, const char * argv[])
{
    int i,k,sample1, sample2, rows=0, columns=0, option = 0;
    char tree[1000]; //"(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0),(OTU4:0.75,OTU5:0.75):1.25))root;";
    char **samples, **otus, fileName[30], sample[30], treeFile[30];
    float **table, result;
    FILE *fp;
    
    system("clear");
    system("more title.txt");
    printf("\nName of biom file containing the table with type (csv): ");
    scanf("%s", fileName);
    printf("\nName of the file containing the Newick Tree: ");
    scanf("%s", treeFile);
    fp=fopen(treeFile,"r");
    if (fp == NULL) {
        printf("\n\n\t ERROR: The file doesn't exist.\n\n");
        exit(1);
    }
    fgets(tree,1000,fp);
    fclose(fp);
    
    rowsColumns(fileName, &rows, &columns);
    
    table = malloc(rows*sizeof(int *));//Matrix with the data from the testing file
    for (i=0; i < rows; i++)
        table[i] = malloc(columns*sizeof(int));
    
    samples = malloc(columns*sizeof(char *));
    for (i=0; i < columns; ++i)
        samples[i]=malloc(30*sizeof(char));
    
    otus = malloc(rows*sizeof(char *));
    for (i=0; i < rows; ++i)
        otus[i]=malloc(30*sizeof(char));
    
    while(option != 4){
        readFile(fileName,table,samples,otus, columns, rows);
        
        printf("\n\nWhat do you want to know?\n");
        printf("\t- (1) Unweighted Unifrac.\n");
        printf("\t- (2) Phylogenetic Diversity.\n");
        printf("\t- (3) Unweighted Unifrac All.\n");
        printf("\t- (4) Exit.\n\n");
        printf("Please select an option (1,2,3,4) and then press <enter> :  ");
        scanf("%d", &option);
        
        switch (option) {
            case 1:
                printf("Name of Sample 1: ");
                scanf("%s", sample);
                sample1 = getSample(sample,samples, columns);
                if (sample1 > columns) {
                    printf("\n*** The sample doesn't not exist.***\n");
                    break;
                }
                printf("Name of Sample 2: ");
                scanf("%s", sample);
                sample2 = getSample(sample,samples, columns);
                if (sample2 > columns) {
                    printf("\n*** The sample doesn't not exist.***\n");
                    break;
                }
                result = unweighted_unifrac(tree, table, otus, samples, sample1, sample2);
                printf("\nUnweighted Unifrac %s/%s = %f\n\n", samples[sample1],samples[sample2],result);
                break;
            case 2:
                printf("Name of Sample: ");
                scanf("%s", sample);
                sample1 = getSample(sample,samples, columns);
                if (sample1 > columns) {
                    printf("\n*** The sample doesn't not exist.***\n");
                    break;
                }
                result = phylogenetic_diversity(tree, table, otus, sample1);
                printf("\nPhylogenetic Diversity for %s = %f\n\n", samples[sample1],result);
                break;
            case 3:
                //unewightet unifrac all
                printf("\n\nUnweighted Unifrac All: \n\n\t");
                for(i=0 ; i<columns ; i++)
                    printf("%s\t", otus[i]);
                printf("\n");
                for(i=0 ; i<columns ; i++){
                    printf("%s\t", otus[i]);
                    for(k=0 ; k<columns ; k++)
                        printf("%.3f\t",unweighted_unifrac(tree, table, otus, samples, i, k));
                    printf("\n");
                }
                printf("\n");
                break;
            case 4:
                break;
            default:
                printf("\n**** I'm sorry, but I don't recognize that option. ****\n\n");
                break;
        }
        
        if(option != 4){
            printf("Press enter to continue... ");
            char enter = 0;
            while (enter != '\r' && enter != '\n'){
                enter = getchar();
                enter = getchar();
            }
        }
        system("clear");
        
    }
    printf("\n**** Thank you! Goodbye.****\n\n");
    free(samples);
    free(table);
    free(otus);
    
}//main

/**
 *  getSample: This function gets a string with the sample name and returns the index where the sample is on the samples array.
 *
 *  @param sample  The sample string
 *  @param samples The samples array
 *  @param columns The number of samples
 *
 *  @return An int which is the index on samples that contains the sample.
 */
int getSample(char sample[], char **samples, int columns){
    int i;
    for (i=0; i<columns; i++) {
        if(strcmp(sample, samples[i]) == 0)
            return i;
    }
    return i;
}

/**
 *  rowsColumns: This function calculates the number of samples and otus that are on the file given by the user.
 *
 *  @param fileName The name of the file containing the data
 *  @param rows     The number of rows on the file
 *  @param columns  The number of columns on the file
 */
void rowsColumns(char fileName[],int *rows, int *columns){
    FILE *fp;
    int i, x=0, y=0;
    fp=fopen(fileName,"r");
    if (fp == NULL) {
        printf("\n\n\t ERROR: The file doesn't exist.\n\n");
        exit(1);
    }
    char line[1000];
    fgets(line,300,fp);
    for(i=0; i < strlen(line);i++)
        if(line[i]== ',')
            x++;
    while(fgets(line,300,fp)!=NULL)// haremos una lectura linea por linea del archivo
        y++;
    
    *rows = y;
    *columns = x;
    fclose(fp);
}


/**
 *  unweighted_unifrac: This function calculates the unweighted unifrac between two samples given by the user.
 *
 *  @param tree    String containing the information of the tree.
 *  @param table   A matrix of floats containing the information of the data observed
 *  @param otus    An array that contains all the OTUs
 *  @param samples An array that contains the name of all samples
 *  @param sample1 The sample to be compared (as an index of the samples array)
 *  @param sample2 The other sample to be compared (as an index of the samples array)
 *
 *  @return The function returns a float which represents the unweighted_unifrac.
 */
float unweighted_unifrac(char tree[], float **table, char **otus, char **samples, int sample1, int sample2){
    int i;
    char otu[10];
    float unique = 0, except = 0, observed, uniFrac;
    
    for(i=0; i<5; i++)
    {
        if (table[i][sample1] > 0 && table[i][sample2] == 0 ) {
            strcpy(otu,otus[i]);
            unique += uniqueLenght(otu,tree);
        }
        else if (table[i][sample1] == 0 && table[i][sample2] > 0 ) {
            strcpy(otu,otus[i]);
            unique += uniqueLenght(otu,tree);
        }
        else if (table[i][sample1] == 0 && table[i][sample2] == 0 ) {
            strcpy(otu,otus[i]);
            except += uniqueLenght(otu,tree);
        }
    }
    
    observed = sumAllLengths(tree) - except;//The sum of the lenght of al the branches minus the unique lenght of all the samples not being compared.
    uniFrac = unique/observed;
    return uniFrac;
}

/**
 *  uniqueLenght: This functions goes through the tree string to obtain the lenght of the branch of a certain OTU.
 *
 *  @param otu  The OTU wich branch length wants to be known.
 *  @param tree The string that contains the tree information.
 *
 *  @return The function returns a floar representing the length of the branch o a certain OTU.
 */
float uniqueLenght(char otu[], char tree[]){
    char * pch, lengthString[15];
    int i;
    pch = strstr (tree,otu) + strlen(otu) + 1;

    for(i=0; *pch != ')'&& *pch != ','; i++){
        lengthString[i] = *pch;
        pch ++;
    }
    lengthString[i] = '\0';

    return atof (lengthString);
}

/**
 *  sumAllLengths: This function sums all the lengths of all the branches of a tree.
 *
 *  @param tree String containing the tree
 *
 *  @return The function returns a float that represents the sum of the lengths of all the branches.
 */
float sumAllLengths(char tree[]){
    char * pch, lengthString[15];
    int i;
    float sum = 0;
    
    while (* pch != ';') {
        if(*pch == ':'){
            pch ++;
            for(i=0; *pch != ')'&& *pch != ','; i++){
                lengthString[i] = *pch;
                pch ++;
            }
            lengthString[i] = '\0';
            sum += atof (lengthString);
        }
        pch ++;
    }
    return sum;
}

/**
 *  phylogenetic_diversity: This function calculates the Phylogenetic Diversity of a certain sample, it goes through the string containing the tree information. The function returns the sum of the lengths of all the branches that lead to a OTU where the sample was observed.
 *
 *  @param tree    The string containing the tree.
 *  @param table   The matrix of floats containing the observed samples data.
 *  @param otus    The array containing the name of all OTUs
 *  @param sample1 The sample which the PD is going to be calculated with.
 *
 *  @return The function returns a float representing the Phylogenetic Diversity of sample1 given by the user.
 */
float phylogenetic_diversity(char tree[], float **table, char **otus, int sample1){
    int i,j, numOtu = 0, nodo = 0, parentesis = -2;
    char otu[30],lengthString[15];
    float length, diversity = 0;
    
    for(i=0 ; i< strlen(tree); i++){
        if(tree[i] == '(')
            parentesis++;
        if(tree[i] == ')')
            parentesis--;
        if(tree[i-1]==':'){
            if(tree[i-2]==')'){
                //parentesis--;
                if(nodo == 1){
                    for(j=0; tree[i] != ')'&& tree[i] != ','; j++){
                        lengthString[j] = tree[i];
                        i++;
                    }
                    i--;
                    lengthString[j] = '\0';
                    diversity += atof (lengthString);
                }
                if(parentesis < 1)
                    nodo = 0;
            }
            else{
                if(table[numOtu][sample1] > 0){//The OTU contains the sample
                    for(j=0; tree[i] != ')'&& tree[i] != ','; j++){
                        lengthString[j] = tree[i];
                        i++;
                    }
                    i--;
                    lengthString[j] = '\0';
                    diversity += atof (lengthString);
                    nodo = 1;
                }
                numOtu++;
            }
        }
    }

    return diversity;
}

/**
 *  readFile: This function reads the file given by the user that contains the observed samples information. The file needs to be a cvs format with the first row having all the samples, and the following rows with the first element being the name of the OTU following with the data. All the information is stored on the arrays and matrix on the program for manipulation.
 *
 *  @param fileName The file that contains the data.
 *  @param table    The matrix of floats where the data will be stored
 *  @param samples  The array where the sample names will be stored
 *  @param otus     The array where the OTUs names will be stored
 *  @param columns  The number of samples
 *  @param rows     The number of OTUs.
 */
void readFile(char fileName[],float **table,char **samples, char **otus, int columns, int rows){
    int i=0,j=0,flag=0;
    char datos[1000];
    char *token;
    FILE *fp;
    fp=fopen(fileName,"r");

    flag=0;
    while(fgets(datos,1000,fp)!=NULL)// haremos una lectura linea por linea del archivo
    {
        i=0;
        if(flag!=0)
        {
            token = strtok(datos,"," );
            strcpy(otus[j],token);
            while( (token=strtok(NULL,",")) != NULL )
            {
                table[j][i]=atoi(token);
                i++;
                
            }
            j++;
        }
        if(flag==0)
        {
            token = strtok(datos,"," );
            strcpy(samples[j],token);
            
            while( (token=strtok(NULL,",")) != NULL )
            {
                j++;
                strcpy(samples[j],token);
            }
            flag=1;
            j=0;
        }
    }
    fclose(fp);
    i = strlen (samples[columns-1]);
    samples[columns-1][i-2] = '\0';
    
    printf("\nConstructed from biom file:\n\n#OTU ID\t\t");
    for (j = 0; j < columns; ++j)
    {
        printf("%s\t\t",samples[j]);
    }
 printf("\n");
    for (i = 0; i < rows; ++i)
    {
        printf("%s\t\t",otus[i]);
        for ( j = 0; j < columns; ++j)
        {
            printf("%.2f\t\t",table[i][j]);
        }
        printf("\n");
    }
}




