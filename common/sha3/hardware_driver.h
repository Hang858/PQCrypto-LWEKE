#ifndef HARDWARE_DRIVER_H
#define HARDWARE_DRIVER_H

#include <stdint.h>
#include <stddef.h>

#define OP_ALG_SHA3_256   0
#define OP_ALG_SHA3_512   1
#define OP_ALG_SM3      2
#define OP_ALG_SHAKE128 3    
#define OP_ALG_SHAKE256 4
#define OP_MODE_NORMAL 0
#define OP_MODE_LMS 1
#define OP_MODE_SPHINCS 2
#define INITIAL_HW_BUFFER_SIZE 4096
#define INITIAL_HW_OUTPUT_SIZE 1024

int OP_hash(uint8_t alg, uint8_t mode, int n, void *input, int input_len, void *output);

#endif