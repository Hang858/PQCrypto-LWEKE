#include "matrix_op.h"

int OP_matrix_mul_8x8(uint16_t Z[8][8], const uint16_t X[8][8], const uint16_t Y[8][8], uint16_t q) {
    // 这里是参考实现，需要替换为硬件指令
    for (int i = 0; i < OP_DIM; i++) {
        for (int j = 0; j < OP_DIM; j++) {
            uint32_t acc = 0;
            for (int k = 0; k < OP_DIM; k++) {
                acc += (uint32_t)X[i][k] * (uint32_t)Y[k][j];
            }
            if (q == 0) {
                Z[i][j] = (uint16_t)acc;
            } else {
                Z[i][j] = (uint16_t)(acc % q);
            }
        }
    }
    return 0;
}

// 标准加载
void OP_load_8x8(matrix8x8_t *dst, const uint16_t *src, int stride) {
    for (int r = 0; r < OP_DIM; r++) {
        for (int c = 0; c < OP_DIM; c++) {
            dst->val[r][c] = src[r * stride + c];
        }
    }
}
 
// 安全加载，当矩阵不是 8 的倍数时
void OP_load_8x8_safe(matrix8x8_t *dst, const uint16_t *src_base, 
                                 int r_offset, int c_offset, int rows, int cols) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if ((r_offset + i) < rows && (c_offset + j) < cols) {
                dst->val[i][j] = src_base[(r_offset + i) * cols + (c_offset + j)];
            } else {
                dst->val[i][j] = 0;
            }
        }
    }
}

// 安全转置加载：处理边界，不足补0
void OP_load_transposed_8x8_safe(matrix8x8_t *dst, const uint16_t *src_base, 
                                            int r_offset, int c_offset, int rows, int cols) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            int src_r = r_offset + j;
            int src_c = c_offset + i;
            
            if (src_r < rows && src_c < cols) {
                dst->val[i][j] = src_base[src_r * cols + src_c];
            } else {
                dst->val[i][j] = 0;
            }
        }
    }
}

// 安全累加回写：处理边界
void OP_store_accumulate_8x8_safe(uint16_t *dst_base, const matrix8x8_t *src, 
                                             int r_offset, int c_offset, int rows, int cols) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if ((r_offset + i) < rows && (c_offset + j) < cols) {
                uint32_t val = (uint32_t)dst_base[(r_offset + i) * cols + (c_offset + j)] + src->val[i][j];
                dst_base[(r_offset + i) * cols + (c_offset + j)] = (uint16_t)(val & 0xFFF);
            }
        }
    }
}

// 转置加载 (将 src 的列变为 dst 的行)
void OP_load_transposed_8x8(matrix8x8_t *dst, const uint16_t *src, int stride) {
    for (int r = 0; r < OP_DIM; r++) {
        for (int c = 0; c < OP_DIM; c++) {
            dst->val[r][c] = src[c * stride + r];
        }
    }
}

// 累加回写
void OP_store_accumulate_8x8(uint16_t *dst, const matrix8x8_t *src, int stride) {
    for (int r = 0; r < OP_DIM; r++) {
        for (int c = 0; c < OP_DIM; c++) {
            dst[r * stride + c] += src->val[r][c];
        }
    }
}

// 掩码覆盖回写
void OP_store_write_mask_8x8(uint16_t *dst, const matrix8x8_t *src, int stride, uint16_t mask) {
    for (int r = 0; r < OP_DIM; r++) {
        for (int c = 0; c < OP_DIM; c++) {
            dst[r * stride + c] = src->val[r][c] & mask;
        }
    }
}

// 8x8 矩阵加法：Acc += Res
void OP_add_8x8(matrix8x8_t *acc, const matrix8x8_t *res) {
    for (int r = 0; r < 8; r++) {
        for (int c = 0; c < 8; c++) {
            acc->val[r][c] += res->val[r][c];
        }
    }
}

// 清零
void OP_clear_8x8(matrix8x8_t *dst) {
    memset(dst, 0, sizeof(matrix8x8_t));
}