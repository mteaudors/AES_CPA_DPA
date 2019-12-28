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

    // byte of the key to attack
    u_int8_t bnum = 1;
    // *** TODO *** Add a loop for bnum = 1 to 16
    
    // Plaintext
    u_int8_t PT[16*NB_TRACES]; 

    // Key
    u_int8_t key[16];   

    double *traces = (double*)malloc(sizeof(double)*N);    

    // matrix for model values
    u_int8_t *P = (u_int8_t*)malloc(sizeof(u_int8_t*)*NB_TRACES*256);         

    // matrix for Pearson coefficients         
    double *R = (double*)malloc(sizeof(double)*NB_SAMPLES*256);                  

    // vector for the mean of the HW according to each key value (mean on each colum of P)
    double H_mean[256];                          

    // read data from files
    read_data("traces.raw",traces,sizeof(double),N);
    read_data("plain.raw",PT,16*sizeof(u_int8_t),NB_TRACES);
    read_data("key.raw",key,16*sizeof(u_int8_t),1);

    // Generate data using the model ( construct P )
    for(int i = 0 ; i < NB_TRACES ; ++i) {
        for(int k = 0 ; k < 256 ; ++k) {
            P[i*256 + k] = HW[intermediate(PT[i*16 + bnum - 1],k)];
        }
    }

    // for(int i = 0 ; i < 2 ; ++i) {
    //     for(int k = 0 ; k < 256 ; ++k) {
    //         printf("%d ",P[i*256+k]);
    //     }
    //     printf("\n");
    // }

    // Compute iteratively the mean of the HW for each key value
    for(int k = 0 ; k < 256 ; ++k) {
        H_mean[k] = 0;
        for(int i = 0 ; i < NB_TRACES ; ++i) {
            H_mean[k] = H_mean[k] + (P[i*256+k]-H_mean[k])/(i+1); 
        } 
    }

    for(int i = 0 ; i < 256 ; ++i) {
        printf("%f ",H_mean[i]);
    }

    printf("\n");
    printf("Correct key : ");
    for(int i = 0 ; i < 16 ; ++i) {
        printf("%02x ",key[i]);
    }
    printf("\n");
    
    free(traces);
    free(P);
    free(R);

    return 0;
}