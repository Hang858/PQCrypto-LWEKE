#ifndef MATRIX_OP_H
#define MATRIX_OP_H

#include <stdint.h>
#include <string.h>

#define OP_DIM 8

typedef struct {
    uint16_t val[OP_DIM][OP_DIM];
} __attribute__((aligned(32))) matrix8x8_t;


/**
 * 8x8 矩阵乘法：Z = X * Y (Mod q)
 */
int OP_matrix_mul_8x8(uint16_t Z[8][8], const uint16_t X[8][8], const uint16_t Y[8][8], uint16_t q);


// --- 数据搬运指令---

/**
 * [OP_load] 标准加载
 * 从大矩阵 src 中加载一个 8x8 块到寄存器 dst
 * dst[r][c] = src[r * stride + c]
 */
void OP_load_8x8(matrix8x8_t *dst, const uint16_t *src, int stride);

/**
 * [OP_load_transposed] 转置加载
 * 从大矩阵 src 中加载并转置一个 8x8 块到寄存器 dst
 * dst[r][c] = src[c * stride + r]
 */
void OP_load_transposed_8x8(matrix8x8_t *dst, const uint16_t *src, int stride);

/**
 * [OP_store_accumulate] 累加回写
 * 将寄存器 src 的值累加回大矩阵 dst
 * dst[r * stride + c] += src[r][c]
 */
void OP_store_accumulate_8x8(uint16_t *dst, const matrix8x8_t *src, int stride);

/**
 * [OP_store_write_mask] 掩码覆盖回写
 * 将寄存器 src 的值截断后写入大矩阵 dst
 * dst[r * stride + c] = src[r][c] & mask
 */
void OP_store_write_mask_8x8(uint16_t *dst, const matrix8x8_t *src, int stride, uint16_t mask);

/**
 * [OP_add] 矩阵加法
 * acc[r][c] += res[r][c]
 */
void OP_add_8x8(matrix8x8_t *acc, const matrix8x8_t *res);

/**
 * [OP_clear] 清零寄存器
 */
void OP_clear_8x8(matrix8x8_t *dst);

void OP_clear_8x8(matrix8x8_t *dst);

/**
 * [load_8x8_safe] 安全加载指令
 * 从大矩阵(线性内存)中加载一个 8x8 块到寄存器。
 * 如果 (r_offset + i, c_offset + j) 在 (rows, cols) 范围内，加载内存值。
 * 否则，寄存器对应位置填充 0。
 * * @param dst        目标 8x8 寄存器
 * @param src_base   源矩阵的首地址
 * @param r_offset   当前块的起始行坐标 (相对于 src_base)
 * @param c_offset   当前块的起始列坐标 (相对于 src_base)
 * @param rows       源矩阵的总行数 (边界上限)
 * @param cols       源矩阵的总列数 (边界上限，同时也作为内存跨度 stride)
 */
void OP_load_8x8_safe(matrix8x8_t *dst, const uint16_t *src_base, 
                   int r_offset, int c_offset, int rows, int cols);

/**
 * [load_transposed_8x8_safe] 安全转置加载指令
 * 从大矩阵中加载一个 8x8 块，并在加载过程中进行转置。
 * dst[i][j] 读取自 src_base 的 Row (r_offset + j), Col (c_offset + i)。
 * 如果源坐标越界，填充 0。
 * @param dst        目标 8x8 寄存器
 * @param src_base   源矩阵的首地址
 * @param r_offset   当前块在源矩阵中的起始行坐标
 * @param c_offset   当前块在源矩阵中的起始列坐标
 * @param rows       源矩阵的总行数
 * @param cols       源矩阵的总列数
 */
void OP_load_transposed_8x8_safe(matrix8x8_t *dst, const uint16_t *src_base, 
                              int r_offset, int c_offset, int rows, int cols);

/**
 * [store_accumulate_8x8_safe] 安全累加回写指令
 * 将 8x8 寄存器中的值累加回大矩阵内存中。
 * dst_base[pos] += src[i][j]。
 * 只有当目标坐标 (r_offset + i, c_offset + j) 在 (rows, cols) 范围内时才写入。
 * @param dst_base   目标矩阵的首地址
 * @param src        源 8x8 寄存器
 * @param r_offset   写入位置的起始行坐标
 * @param c_offset   写入位置的起始列坐标
 * @param rows       目标矩阵的总行数
 * @param cols       目标矩阵的总列数
 */
void OP_store_accumulate_8x8_safe(uint16_t *dst_base, const matrix8x8_t *src, 
                               int r_offset, int c_offset, int rows, int cols);

#endif