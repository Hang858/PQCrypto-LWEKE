#include <stdint.h>

/**
 * @brief 算子2. 向量乘法  (元素数量可变)
 * * Z = (X * Y) mod q
 * 
 * @param z_out     输出：一个 16-bit 元素。
 * @param x_in      输入：向量  (length个 16-bit 元素)。
 * @param y_in      输入：向量  (length个 16-bit 元素)。
 * @param length    输入：向量长度。（最大长度：Frodo: 1344, scloudplus: 1136）
 * @param q         输入：模数 (16-bit) 2的次幂。
 * @return          0: 成功; -1: 失败。
 * 待定：数据无符号
 */
int OP_vector_mul(uint16_t *z_out, const uint16_t *x_in, const uint16_t *y_in, uint16_t length, uint16_t q);


/**
 * @brief 矩阵转置辅助函数
 * 将输入矩阵 (rows x cols) 转置为 (cols x rows)
 * @param in  输入矩阵指针
 * @param out 输出矩阵指针
 * @param rows 输入矩阵的行数
 * @param cols 输入矩阵的列数
 */
void matrix_transpose(const uint16_t *in, uint16_t *out, int rows, int cols);