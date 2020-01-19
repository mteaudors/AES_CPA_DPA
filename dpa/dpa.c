#include "dpa.h"

u_int8_t intermediate(u_int8_t plain , u_int8_t key) {
    return sbox[plain^key];
}

int read_data(char *path, void *buffer, int size_element, int nb_element) {
    FILE *f = fopen(path, "r");
    
    if(f == NULL) {
        fprintf(stderr , "Error opening file %s : %s\n",path,strerror(errno));   
        exit(EXIT_FAILURE);     
    }

    if(fread(buffer , size_element, nb_element , f) == 0) {
        fprintf(stderr , "Error reading file %s : %s\n",path,strerror(errno));   
        exit(EXIT_FAILURE);     
    }
    fclose(f);
}

void init(double *group, int length){
    for(int i = 0; i < length; ++i)
        group[i] = 0;
}

void add(double *res, double *group, double *traces, int length){
    for(int i = 0; i < length; ++i)
        res[i] = group[i] + traces[i];
}

void subtract(double *res, double *group, double *traces, int start, int length){
    for(int i = start; i < start+length; ++i)
        res[i - start] = traces[i] - group[i - start];
}

void divide(double *res, double *group, double denumerator, int length){
    for(int i = 0; i < length; ++i)
        res[i] = group[i] / denumerator;
}

// Calculate the absolute value of the difference of two double buffers
void abs_sub(double *res, double *grp1, double *grp2, int length){
    for(int i = 0; i < length; ++i)
        res[i] = fabs(grp1[i] - grp2[i]);
}

double max(double *group, int length){
    double max = 0;
    for(int i = 0; i < length; ++i)
        if(group[i] > max)
            max = group[i];

    return max;
}

void print_arr(double *arr, int length){
    for(int i = 0; i < length; ++i)
        printf("%f\n", arr[i]);
    printf("\n");
}

int arg_max(double *group, int length){
    double max = 0;
    int argmax = 0;
    for(int i = 0; i < length; ++i){
        if(group[i] > max){
            max = group[i];
            argmax = i;
        }
    }
    return argmax;
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

    // read data from files
    read_data("traces.raw",traces,sizeof(double),N);
    read_data("plain.raw",PT,16*sizeof(u_int8_t),NB_TRACES);
    read_data("key.raw",key,16*sizeof(u_int8_t),1);

    // start timer
    start = clock();

    // loop over the 16 bytes of the key
    #pragma omp parallel for
    for(bnum = 0; bnum < 16; ++bnum) {
        //init a buffer of length 256, max value of a byte, and initialize it
        double group[256];
        init(group, 256);

        // Temporary buffer that will contain a trace
        double tmp[NB_SAMPLES];

        for(int k = 0; k < 256; ++k){
            // Two separate groups to discriminate the traces
            double grp1[NB_SAMPLES];
            double grp2[NB_SAMPLES];
            
            init(grp1, NB_SAMPLES);
            init(grp2, NB_SAMPLES);

            int nb_traces_g1 = 0;
            int nb_traces_g2 = 0;

            /*
             * The next part calculate the means iterativaly. However,
             * it requires a first loop to calculate the length of the
             * two sub-groups, as we don't know know it firsthand.
             * The second loop then calculates the means iteratively of
             * the two sub-groups.
             *
             * Using this method reduces the cancellation phenomena, but
             * increases the time taken overall by a slight amount.
             */
            for(int i = 0; i < NB_TRACES; ++i){
                int hw_of_byte = intermediate(PT[i*16 + bnum], k);
                // Discriminate over the second to last bit
                int first_byte = hw_of_byte & 0b0000001;

                if(first_byte == 1)
                    nb_traces_g1++;
                else
                    nb_traces_g2++;
            }

            for(int i = 0; i < NB_TRACES; ++i){
                int hw_of_byte = intermediate(PT[i*16 + bnum], k);
                int first_byte = hw_of_byte & 0b0000001;

                // Calculate each mean iteratively
                if(first_byte == 1){
                    subtract(tmp, grp1, traces, i*NB_SAMPLES, NB_SAMPLES);
                    divide(tmp, tmp, (double)nb_traces_g1, NB_SAMPLES);
                    add(grp1, grp1, tmp, NB_SAMPLES);
                }else{
                    subtract(tmp, grp2, traces, i*NB_SAMPLES, NB_SAMPLES);
                    divide(tmp, tmp, (double)nb_traces_g2, NB_SAMPLES);
                    add(grp2, grp2, tmp, NB_SAMPLES);
                }
            }
            // Calculate the absolute value of the difference of means
            abs_sub(tmp, grp1, grp2, NB_SAMPLES);
            // Select the maximum for each potential sub-key
            group[k] = max(tmp, NB_SAMPLES);
        }
        // Select the maximum argument, our guess for a specific key byte
        int result = arg_max(group, 256);
        guess[bnum] = result;
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
