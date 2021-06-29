# include <stdio.h>
# include <stdlib.h>
# include "el1.h"

#define LEN 50
#define TOTALITER 10
#define OFFSET 0
#define LEN_GAUSS 100000
#define LEN_IDX 1500
 
int main()
{
   int sparsityVect[LEN];
   int counterVect[TOTALITER];
   int counterVect2[TOTALITER];
   double *sol = NULL;
  
   FILE* file = fopen("rng_gauss.csv", "r");
   int i = 0;
   int j = 0;
   double num = 0;
   double rng_gauss[LEN_GAUSS];
   int d2 = 142;

   while(fscanf(file, "%lf,", &num) > 0)
   {
       rng_gauss[i++] = num;
      if (i>=LEN_GAUSS) break;
   }

   fclose(file);

   int num2 = 0;
   int rng_idx[LEN_IDX];
   i = 0;
   file = fopen("rng_idx.csv", "r");
   while(fscanf(file, "%d,", &num2) > 0){
      rng_idx[i] = num2-1;
      i++;
      if (i>=LEN_IDX) break;
   }

   fclose(file);


   for(i=0;i<LEN;i++){
      sparsityVect[i] = OFFSET + 2*(i+1);
   }

   for(i=0;i<TOTALITER;i++){
      counterVect[i] = 10000*i;
      counterVect2[i] = 150*i;
   }
   
   file = fopen("res_el1.csv", "w");
   if (file == NULL)
   {
      printf("Error opening file!\n");
      exit(1);
   }


   for(i=0;i<LEN;i++){
      for(j=0;j<TOTALITER;j++){
         printf("NEW PROBLEM: Sparsity %d, Iteration %d\n", sparsityVect[i], j);
         sol = solveEl1sparse(sparsityVect[i], counterVect[j], counterVect2[j], rng_gauss, rng_idx);
         fprintf(file, "%4.2f,%4.2f,%4.2f,%4.2f\n", sol[0],sol[1],sol[2],sol[3]);
      }
   }

   fclose(file);
   return 0;
}
