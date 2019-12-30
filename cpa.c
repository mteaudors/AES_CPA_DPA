#include "cpa.h"

u_int8_t intermediate(u_int8_t plain , u_int8_t key) {
    return sbox[plain^key];
}

int read_data(char *path, void *buffer, int size_element, int nb_element) {
    FILE *f = fopen(path, "r");
    fread(buffer , size_element, nb_element , f);
    fclose(f);
}


int main(int argc, char *argv[]) {

    // define the number of threads to use 
    int nb_thread = (argc==2) ? atoi(argv[1]) : NB_THREAD;
    omp_set_num_threads(nb_thread);


    // To measure the time
    clock_t start, end;
    double cpu_time_used;

    // byte of the key to attack
    u_int8_t bnum;
    
    // Plaintext
    u_int8_t PT[16*NB_TRACES]; 

    // Correct key
    u_int8_t key[16];  

    // Guess key
    u_int8_t guess[16]; 

    // Raw values of power consuption
    double *traces = (double*)malloc(sizeof(double)*N);    

    // vector for the mean of the power consuption at each instant
    double T_mean[NB_SAMPLES];                        

    // read data from files
    read_data("traces.raw",traces,sizeof(double),N);
    read_data("plain.raw",PT,16*sizeof(u_int8_t),NB_TRACES);
    read_data("key.raw",key,16*sizeof(u_int8_t),1);

    // Compute iteratively the mean of the power consuption at each instant
    for(int i = 0 ; i < NB_SAMPLES ; ++i) {
        T_mean[i] = 0;
        for(int l = 0 ; l < NB_TRACES ; ++l) {
            T_mean[i] = T_mean[i] + (traces[l*NB_SAMPLES+i]-T_mean[i])/(l+1); 
        } 
    }   

    // start timer
    start = clock();

    // loop over the 16 bytes of the key
    #pragma omp parallel for
    for(bnum = 1 ; bnum<=16; ++bnum) {

        // matrix for model values
        u_int8_t *P = (u_int8_t*)malloc(sizeof(u_int8_t*)*NB_TRACES*256);         

        // vector for the mean of the HW according to each key value (mean on each colum of P)
        double H_mean[256];  
        
        // Generate data using the model ( construct P )
        for(int i = 0 ; i < NB_TRACES ; ++i) {
            for(int k = 0 ; k < 256 ; ++k) {
                P[i*256 + k] = HW[intermediate(PT[i*16 + bnum - 1],k)];
            }
        }

        // Compute iteratively the mean of the HW for each key value
        for(int k = 0 ; k < 256 ; ++k) {
            H_mean[k] = 0;
            for(int i = 0 ; i < NB_TRACES ; ++i) {
                H_mean[k] = H_mean[k] + (P[i*256+k]-H_mean[k])/(i+1); 
            } 
        }

        // Compute Pearson correlation coefficients ( matrix R )

        double max = 0;
        double A,B,C;
    
        for(int i = 0 ; i < NB_SAMPLES ; ++i) {
            for(int k = 0 ; k < 256 ; ++k) {

                A = 0;
                B = 0; 
                C = 0;
                
                for(int l = 0 ; l < NB_TRACES ; ++l) {
                    A += (P[l*256+k] - H_mean[k])*(traces[l*NB_SAMPLES+i] - T_mean[i]);
                    B += (P[l*256+k] - H_mean[k])*(P[l*256+k] - H_mean[k]);
                    C += (traces[l*NB_SAMPLES+i] - T_mean[i])*(traces[l*NB_SAMPLES+i] - T_mean[i]);
                }
                double r = fabs(A/(sqrt(B*C)));
                if(r > max) {
                    max = r;
                    guess[bnum-1] = k;
                }
            }
        }
        free(P);
    }


    
    // end timer
    end = clock();

    printf("Correct key : \t");
    for(int i = 0 ; i < 16 ; ++i) {
        printf("%02x ",key[i]);
    }
    printf("\n");
    
    printf("Guess : \t");
    for(int i = 0 ; i < 16 ; ++i) {
        printf("%02x ",guess[i]);
    }
    printf("\n");

    double success;
    for(int i = 0 ; i < 16 ; ++i) {
        success += (key[i] == guess[i]);
    }

    success = (success/16)*100;
    printf("Success attack : %0.1f %%\n",success);

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Total Computation time : %0.2f s\n",cpu_time_used);

    printf("\n");
    printf("Number of thread : %d\n",nb_thread);
    printf("Average time per thread : %0.2f s\n",cpu_time_used/nb_thread);
    free(traces);
    
    return 0;
}