//> Some macro definitions

//> General Settings
#define FIX_RANDOMNESS                  (true)
#define RUN_CERES_SOLVER_ON             (false)

//> Thresholds Settings
#define CERTIFY_THRESH                  (-1e-09)   //> -1e-09 accroding to the paper. Negative because of numerical error.
#define DIFF_CONSECUTIVE_SOLS_THRESH    (3e-10)    //> 3e-10 according to the paper

#define LOG_INFOR_MESG(info_msg)        printf("\033[1;32m[INFO] %s\033[0m\n", std::string(info_msg).c_str() );
#define LOG_FILE_ERROR(err_msg)         printf("\033[1;31m[ERROR] File %s not found!\033[0m\n", std::string(err_msg).c_str() );
#define LOG_ERROR(err_msg)              printf("\033[1;31m[ERROR] %s\033[0m\n", std::string(err_msg).c_str() );
#define LOG_DATA_LOAD_ERROR(err_msg)    printf("\033[1;31m[DATA LOAD ERROR] %s not loaded successfully!\033[0m\n", std::string(err_msg).c_str() );

#define PRINT_VECTOR3D(name, vec)       printf("%s = [%f, %f, %f]\n", std::string(name).c_str(), vec(0), vec(1), vec(2));
#define PRINT_VECTORXD(name, vec)       printf("%s = [", std::string(name).c_str()); for(int i = 0; i < vec.size(); i++) {printf("%f ", vec(i));} printf("]\n");

#define PRINT_ESSENTIAL(id1, id2, mat)  printf("E%d%d = [\n", id1, id2); \
                                        for(int i = 0; i < 3; i++) { \
                                            for(int j = 0; j < 3; j++) printf("%10.7f ", mat(i,j)); \
                                            printf(";\n"); \
                                        } printf("]\n");

