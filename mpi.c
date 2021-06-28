
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "mpi.h"

#define MAX    256
#define MASTER   0  // ranque mestre:  não modificar
//#define COMMENT1     // se não quer os printf para os desc_query colocar esta linha em comentário
//#define COMMENT2     // se não quer os printf para os seq_dna com resultado colocar esta linha em comentário
//#define COMMENT3     // se não quer os printf para os seq_dna sem resultado  colocar esta linha em comentário

// Cor de texto no stdout
#define KNRM  "\x1B[0m"   // normal
#define KRED  "\x1B[31m"  // red
#define KGRN  "\x1B[32m"  // green
#define KYEL  "\x1B[33m"  // yellow
#define KBLU  "\x1B[34m"  // blue
#define KMAG  "\x1B[35m"  // magenta
#define KCYN  "\x1B[36m"  // cyan 
#define KWHT  "\x1B[37m"  // white

#define COR KYEL // cor para os printfs dos seq_dna com resultado

typedef struct Sequencia
 {
  int  lseq;    // comprimento da sequência de DNA
  char dna[70000];    // sequencia de DNA    
  int result;   // resultado da busca
 } sequencia;

MPI_File fdatabase,fp,fquery,fp2,fout;
sequencia tab[22000];
char *seq,*seq_query;
char desc[22000][100];
int rank,nprocs;  // ranque e quantia de processadores

void programa_sequencial();

// Boyers-Moore-Hospool-Sunday algorithm for string matching
int bmhs(char *string, int n, char *substr, int m)
 {
  int d[MAX];
  int i, j, k;

  // pre-processing
  for ( j = 0; j < MAX; j++) d[j] = m + 1;
  
  for (j = 0; j < m; j++) d[(int) substr[j]] = m - j;

  // searching
  i = m - 1;

  while ( i < n )
   {
    k = i;
    j = m - 1;

    while ( (j >= 0) && (string[k] == substr[j]) ) { j--; k--; }

    if ( j < 0 ) return k + 1;

    i = i + d[(int) string[i + 1]];
   }

  return -1;
 }

int ler_sequencias() // leitura e indexação das sequencias de DNA
 {
  char line[100];
  int i,num,count,test,dest,etiq = 0;
  long int pos,lenght,lg;
  MPI_Status status;

  pos = 0; // offset no arquivo dna.in
  num = 0; // contador do número de sequencia

  dest = 1;
  test = 0;

  do
   {
    MPI_File_read_at(fdatabase,pos,line,100,MPI_CHAR,&status);
    MPI_Get_count(&status,MPI_CHAR,&count);

    if ( count != 0 ) // count == 0 => End of File
     {
      for ( i = 0 ; i < 100 ; i++ ) if ( line[i] == '\r' ) { line[i] = '\0' ; i = 100; } 
      lenght = strlen(line);
      line[lenght] = '\0';

      if ( line[0] == '>' ) // leitura da descrição
       {
        if ( test )
         {
          lg = strlen(seq);
          MPI_Send(&lg,1,MPI_INT,dest,etiq,MPI_COMM_WORLD);    // envio do comprimento da sequencia
          MPI_Send(seq,lg,MPI_CHAR,dest,etiq,MPI_COMM_WORLD);  // envio da sequencia 
          dest++;
          if ( dest == nprocs ) dest = 1;
         }
        else test = 1; 

        seq[0] = '\0'; // string nula. Nova sequência
        strcpy(desc[num],line); // salvar a descrição da sequencia
        num++;
       }
      else  strcat(seq,line); // acrescentar o pedaço da sequência ao buffer

      pos += lenght+2; // comprimento + \r \n : novo offset para a leitura de dna.in
     }
   }
  while ( count != 0 ); // count == 0 => End of File

  lg = strlen(seq);
  MPI_Send(&lg,1,MPI_INT,dest,etiq,MPI_COMM_WORLD);    // envio do comprimento da sequencia
  MPI_Send(seq,lg,MPI_CHAR,dest,etiq,MPI_COMM_WORLD);  // envio da sequencia 

  lenght = 0;  // finalizar a leitura das sequencias
  for ( dest = 1 ; dest < nprocs ; dest++ ) MPI_Send(&lenght,1,MPI_INT,dest,etiq,MPI_COMM_WORLD);

  return num; // número de sequência de DNA no arquivo
 }

void ler_querys(int num)
 {
  long int pos;
  int proc,ind,i,dest,test,lenght,count,result,etiq = 0;
  char desc_query[100],str[50000],line[101];
  MPI_Status status;
 
  pos = 0;
 
  do
   {
    // leitura da descrição do query
    MPI_File_read_at(fquery,pos,desc_query,100,MPI_CHAR,&status);
    MPI_Get_count(&status,MPI_CHAR,&count);

    if ( count != 0 ) // count == 0 => End of File
     {
      // eliminar o \n
      for ( i = 0 ; i < 100 ; i++ ) if ( desc_query[i] == '\n' ) { desc_query[i] = '\0' ; break; } 

      lenght = strlen(desc_query);
      pos += lenght + 1; // comprimento + \n : novo offset para a leitura de query.in

      seq_query[0] = '\0';

      do
       {
        MPI_File_read_at(fquery,pos,line,100,MPI_CHAR,&status);
        MPI_Get_count(&status,MPI_CHAR,&count);

        if ( (count != 0) && (line[0] != '>') )
         {
          for ( i = 0 ; i < 100 ; i++ ) if ( line[i] == '\n' ) { line[i] = '\0' ; break; } 
          lenght = strlen(line);
          if ( i == 100 ) pos += lenght; // + comprimento : novo offset para a leitura de query.in
          else pos += lenght+1;          // + comprimento + \n : novo offset para a leitura de query.in
          strcat(seq_query,line);
         }
       }
      while ( (count != 0) && (line[0] != '>') );

      sprintf(str,"\n%s\n",desc_query); // descrição do query
      lenght = strlen(str);
      MPI_File_write(fout,str,lenght,MPI_CHAR,&status); // escrever a descrição do query no arquivo dna.out
      #ifdef COMMENT1
      printf("\n%s\n",desc_query);
      fflush(stdout);
      #endif

      for ( dest = 1 ; dest < nprocs ; dest++ ) // envio de seq_query para cada processo
       {
        lenght = strlen(seq_query);
        MPI_Send(&lenght,1,MPI_INT,dest,etiq,MPI_COMM_WORLD);
        MPI_Send(seq_query,lenght,MPI_CHAR,dest,etiq,MPI_COMM_WORLD);
       }

      test = 0;

      for ( ind = 0 , proc = 1 ; ind < num ; ind++ )
       {
        MPI_Recv(&result,1,MPI_INT,proc,etiq,MPI_COMM_WORLD,&status);    // recepção do resultado

        if ( result > 0 ) 
         {
          sprintf(str,"%s\n%d\n",desc[ind],result);  // impressão do resultado da busca
          lenght = strlen(str);
          MPI_File_write(fout,str,lenght,MPI_CHAR,&status); // escrever no arquivo dna.out
          #ifdef COMMENT2
          printf("%s%s %d%s\n",COR,desc[ind],result,KNRM);
          fflush(stdout);
          #endif
          test = 1;
         } 
        #ifdef COMMENT3
        else
         {
          printf("%s\n",desc[ind]);
          fflush(stdout);
         }
        #endif

        proc++;
        if ( proc == nprocs ) proc = 1;
       }

      if ( !test )
       {
        sprintf(str,"NOT FOUND\n");
        lenght = strlen(str);
        MPI_File_write(fout,str,lenght,MPI_CHAR,&status); // escrever no arquivo dna.out
       }
     } // end   if ( count != 0 )
   }
  while ( count != 0 ); // count == 0 => End of File 

  lenght = 0;  // finalizar todos os processos
  for ( dest = 1 ; dest < nprocs ; dest++ ) MPI_Send(&lenght,1,MPI_INT,dest,etiq,MPI_COMM_WORLD);
 }

void rank_master()
 {
  int num;

  num = ler_sequencias(); // leitura das sequencias, retorna o número de sequencia
  ler_querys(num);        // leitura dos querys e busca nas sequencias de DNA
 }

void rank_not_master()
 { 
  int i,ind,lseq,lquery,dest,etiq = 0;
  MPI_Status status;

  ind = 0;
  MPI_Recv(&lseq,1,MPI_INT,0,etiq,MPI_COMM_WORLD,&status);    // recepção do comprimento da sequencia

  while ( lseq != 0 ) // lseq == 0  => finalizar o processo
   {
    MPI_Recv(seq,lseq,MPI_CHAR,0,etiq,MPI_COMM_WORLD,&status);  // recepção da sequencia
    strcpy(tab[ind].dna,seq);  // salvar a sequencia na memória
    tab[ind].lseq = lseq;      // salvar o comprimento da sequencia na memória

    ind++;
    MPI_Recv(&lseq,1,MPI_INT,0,etiq,MPI_COMM_WORLD,&status);    // recepção do comprimento da sequencia
   }

  do
   {
    MPI_Recv(&lquery,1,MPI_INT,0,etiq,MPI_COMM_WORLD,&status);          // recepção do comprimento do query

    if ( lquery != 0 )
     {
      MPI_Recv(seq_query,lquery,MPI_CHAR,0,etiq,MPI_COMM_WORLD,&status);  // recepção do query

      // aplicação da função bmhs para todas as sequencias que foram atribuidas ao processo
      for ( i = 0 ; i < ind ; i++ ) tab[i].result = bmhs(tab[i].dna,tab[i].lseq,seq_query,lquery);

      for ( i = 0 , dest = 0 ; i < ind ; i++ ) // para todas as sequencias atribuidas
        MPI_Send(&(tab[i].result),1,MPI_INT,dest,etiq,MPI_COMM_WORLD); // envio do resultado para a sequencia
     }
   }
  while ( lquery != 0 );
 }

int main(int argc, char *argv[])
 {
  int i,etiq = 0;
  double dt,t,t1,t2;
  MPI_Status status;

  /* Chamada inicial para o MPI */
  MPI_Init(&argc,&argv);

  t1 = MPI_Wtime(); // tempo inicial

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs); // Determina o número do processos

  if ( nprocs > 1 )
   {
    seq = (char*) malloc(sizeof(char) * 1000001);

    if ( seq == NULL )
     {
      printf("erro de alocação de memória : seq\n");
      MPI_Finalize();
      exit (-1);
     }

    seq_query = (char*) malloc(sizeof(char) * 1000001);

    if ( seq_query == NULL )
     {
      printf("erro de alocação de memória : seq_query\n");
      MPI_Finalize();
      exit (-1);
     }

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);   // Determina o rank

    MPI_File_open(MPI_COMM_WORLD,"query.in",MPI_MODE_RDONLY,MPI_INFO_NULL,&fquery);
    MPI_File_open(MPI_COMM_WORLD,"dna.in",MPI_MODE_RDONLY,MPI_INFO_NULL,&fdatabase);
    MPI_File_delete( "dna.out", MPI_INFO_NULL );
    MPI_File_open(MPI_COMM_WORLD,"dna.out",MPI_MODE_CREATE|MPI_MODE_RDWR,MPI_INFO_NULL,&fout);

    if  ( rank == MASTER ) rank_master(); //  processo mestre
    else rank_not_master(); // processo não mestre

    free(seq);
    free(seq_query);

    MPI_File_close(&fout);
    MPI_File_close(&fquery);
    MPI_File_close(&fdatabase);
   }
  else programa_sequencial();

  t2 = MPI_Wtime();
  dt = (t2-t1);

  if ( rank == MASTER )
   {
    for ( i = 1 ; i < nprocs ; i++ )
     {
      MPI_Recv(&t,1,MPI_INT,i,etiq,MPI_COMM_WORLD,&status); // recepção do tempo de execução dos processos não mestres
      if ( t > dt ) dt = t;
     }
 
    printf("np = %d tempo de execução: %.1f s\n",nprocs,dt);   // tempo em segundos
   }
  else MPI_Send(&dt,1,MPI_INT,0,etiq,MPI_COMM_WORLD); // envio do tempo de execução do processo

  MPI_Finalize();

  return 0;
 }

/***********************************************************************************************************************/
/*             PROGRAMA SEQUENCIAL                                                                                     */
/***********************************************************************************************************************/

static inline void remove_eol(char *line)
 {
  int i = strlen(line) - 1;

  while (line[i] == '\n' || line[i] == '\r')
   {
    line[i] = 0;
    i--;
   }
 }

FILE *fdatabase2, *fquery2, *fout2;
char *bases,*str;

void openfiles()
 {
  fdatabase2 = fopen("dna.in", "r+");

  if (fdatabase2 == NULL)
   {
    perror("dna.in");
    exit(EXIT_FAILURE);
   }

  fquery2 = fopen("query.in", "r");

  if (fquery2 == NULL)
   {
    perror("query.in");
    exit(EXIT_FAILURE);
   }

  fout2 = fopen("dna.out", "w");

  if (fout2 == NULL)
   {
    perror("fout");
    exit(EXIT_FAILURE);
   }
 }

void closefiles()
 {
  fflush(fdatabase2);
  fclose(fdatabase2);

  fflush(fquery2);
  fclose(fquery2);

  fflush(fout2);
  fclose(fout2);
 }

void programa_sequencial()
 {
  char desc_dna[100], desc_query[100];
  char line[100];
  int i,test,found,result,max;

  bases = (char*) malloc(sizeof(char) * 4000001);

  if (bases == NULL)
   {
    perror("malloc");
    exit(EXIT_FAILURE);
   }

  str = (char*) malloc(sizeof(char) * 4000001);

  if (str == NULL)
   {
    perror("malloc str");
    exit(EXIT_FAILURE);
   }

  openfiles();

  fgets(desc_query, 100, fquery2); // leitura da descrição do query
  remove_eol(desc_query);

  test = 1;

  while (!feof(fquery2)) // até o final do arquivo query
   {
    #ifdef COMMENT1
    printf("\n%s\n",desc_query);
    fflush(stdout);
    #endif
    fprintf(fout2, "\n%s\n", desc_query); // escitura da descrição do query no arquivo dna.out
    fflush(fout2);

    fgets(line, 100, fquery2);  // leitura de um pedaço do query
    remove_eol(line);
    str[0] = 0;
    i = 0;

     do
     {
      strcat(str+i,line); // juntar os pedaços do query
      if (fgets(line,100,fquery2) == NULL) break;  // leitura de um pedaço do query
      remove_eol(line);
      i += 80;
     }
    while (line[0] != '>'); // até chegar a uma nova descrição 

    strcpy(desc_query, line); // copiar a nova descrição do query

    // read database and search
    found = 0;
    fseek(fdatabase2, 0, SEEK_SET);  // colocar o ponteiro do arquivo dna.in no inicio do arquivo
    fgets(line, 100, fdatabase2);    // ler a descrição da sequencia de DNA
    remove_eol(line);

    while (!feof(fdatabase2)) 
     {
      strcpy(desc_dna, line);  // copiar a descrição da sequencia de DNA
      bases[0] = 0;
      i = 0;
      fgets(line, 100, fdatabase2);  // ler um pedaço da sequencia de DNA
      remove_eol(line);

      do
       {
        strcat(bases + i, line); // juntar os pedaços da sequencia
        if (fgets(line, 100, fdatabase2) == NULL) break; // ler um pedaço da sequencia de DNA

        if ( test )
         {
          test = 0;
          if ( line[51] == '\n' ) max = 50; // linha de 50 caracteres
          else max = 80; // linha de 80 caracteres
         }

        remove_eol(line);
        i += max;
       }
      while (line[0] != '>'); // até chegar na próxima descrição de sequência

      result = bmhs(bases, strlen(bases), str, strlen(str));

      if (result > 0)
       {
        #ifdef COMMENT2
        printf("%s%s %d%s\n",COR,desc_dna,result,KNRM);
        fflush(stdout);
        #endif
        fprintf(fout2, "%s\n%d\n", desc_dna, result);  // imprimir o resultado
        fflush(fout2); // esvaziar o buffer. Forçar a impressão                    
        found++;
       }
      #ifdef COMMENT3
      else
       {
        printf("%s\n",desc_dna);
        fflush(stdout);
       }
      #endif
     }

    if (!found) fprintf(fout2, "NOT FOUND\n");
   }

  closefiles();
  free(str);
  free(bases);
 }

