#include <stdint.h>
#include <stddef.h>
#include "vector_op.h"

int OP_vector_mul(uint16_t *z_out, const uint16_t *x_in, const uint16_t *y_in, uint16_t length, uint16_t q) {
    if (z_out == NULL || x_in == NULL || y_in == NULL) {
        return -1; 
    }
    if (q == 0) {
        return -1;
    }
    uint64_t accumulator = 0;
    for (int i = 0; i < length; ++i) {
        accumulator += (uint64_t)x_in[i] * y_in[i];
    }

    *z_out = (uint16_t)(accumulator % q);
    return 0;
}

void matrix_transpose(const uint16_t *in, uint16_t *out, int rows, int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            out[j * rows + i] = in[i * cols + j];
        }
    }
}