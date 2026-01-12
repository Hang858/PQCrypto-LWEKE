/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: matrix arithmetic functions used by the KEM
*********************************************************************************************/

#if defined(USE_AES128_FOR_A)
#if !defined(USE_OPENSSL)
    #include "../../common/aes/aes.h"
#else
    #include "../../common/aes/aes_openssl.h"
#endif
#elif defined (USE_SHAKE128_FOR_A)
#if !defined(USE_AVX2)
    #include "../../common/sha3/fips202.h"
#else
    #include "../../common/sha3/fips202x4.h"
#endif
#endif    
#if defined(USE_AVX2)
    #include <immintrin.h>
#endif
#include "vector_op.h"

int frodo_mul_add_as_plus_e(uint16_t *out, const uint16_t *s, const uint16_t *e, const uint8_t *seed_A) 
{
    int i, j, k, p;
    // 矩阵 A 分配内存
    size_t size_A = (size_t)PARAMS_N * PARAMS_N * sizeof(uint16_t);
    uint16_t *MatrixA = malloc(size_A);

    for (i = 0; i < (PARAMS_N*PARAMS_NBAR); i += 2) {    
        *((uint32_t*)&out[i]) = *((uint32_t*)&e[i]);
    }    
    
#if defined(USE_AES128_FOR_A)
    int16_t a_row_temp[8*PARAMS_N] = {0};                              
#if !defined(USE_OPENSSL)
    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);   
#else
    EVP_CIPHER_CTX *aes_key_schedule;    
    int len;
    if (!(aes_key_schedule = EVP_CIPHER_CTX_new())) handleErrors();    
    if (1 != EVP_EncryptInit_ex(aes_key_schedule, EVP_aes_128_ecb(), NULL, seed_A, NULL)) handleErrors();    
#endif
#endif

    for (i = 0; i < PARAMS_N; i += 8) {

#if defined(USE_AES128_FOR_A)
        for (j = 0; j < PARAMS_N; j += PARAMS_STRIPE_STEP) {
            for(p=0; p<8; p++) a_row_temp[j + 1 + p*PARAMS_N] = UINT16_TO_LE(j);
        }
        for (j = 0; j < PARAMS_N; j += PARAMS_STRIPE_STEP) {
            for(p=0; p<8; p++) a_row_temp[j + p*PARAMS_N] = UINT16_TO_LE(i+p);
        }
#if !defined(USE_OPENSSL)
        AES128_ECB_enc_sch((uint8_t*)a_row_temp, 8*PARAMS_N*sizeof(int16_t), aes_key_schedule, (uint8_t*)(MatrixA + i*PARAMS_N));
#else   
        if (1 != EVP_EncryptUpdate(aes_key_schedule, (uint8_t*)(MatrixA + i*PARAMS_N), &len, (uint8_t*)a_row_temp, 8*PARAMS_N*sizeof(int16_t))) handleErrors();
#endif

#elif defined (USE_SHAKE128_FOR_A)       
#if !defined(USE_AVX2)
        uint8_t seed_A_separated[2 + BYTES_SEED_A];
        uint16_t* seed_A_origin = (uint16_t*)&seed_A_separated;
        memcpy(&seed_A_separated[2], seed_A, BYTES_SEED_A);
        for (p = 0; p < 8; p++) {
            seed_A_origin[0] = UINT16_TO_LE(i + p);
            shake128((unsigned char*)(MatrixA + (i+p)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        }
#else
        uint8_t seed_A_sep[4][2 + BYTES_SEED_A];
        uint16_t* seed_A_origin[4];
        for(int x=0; x<4; x++) { 
             seed_A_origin[x] = (uint16_t*)&seed_A_sep[x];
             memcpy(&seed_A_sep[x][2], seed_A, BYTES_SEED_A);
        }

        for(int x=0; x<4; x++) seed_A_origin[x][0] = UINT16_TO_LE(i + x);
        shake128_4x((unsigned char*)(MatrixA + i*PARAMS_N), (unsigned char*)(MatrixA + (i+1)*PARAMS_N), (unsigned char*)(MatrixA + (i+2)*PARAMS_N), (unsigned char*)(MatrixA + (i+3)*PARAMS_N), 
                    (unsigned long long)(2*PARAMS_N), seed_A_sep[0], seed_A_sep[1], seed_A_sep[2], seed_A_sep[3], 2 + BYTES_SEED_A);
        
        for(int x=0; x<4; x++) seed_A_origin[x][0] = UINT16_TO_LE(i + 4 + x);
        shake128_4x((unsigned char*)(MatrixA + (i+4)*PARAMS_N), (unsigned char*)(MatrixA + (i+5)*PARAMS_N), (unsigned char*)(MatrixA + (i+6)*PARAMS_N), (unsigned char*)(MatrixA + (i+7)*PARAMS_N), 
                    (unsigned long long)(2*PARAMS_N), seed_A_sep[0], seed_A_sep[1], seed_A_sep[2], seed_A_sep[3], 2 + BYTES_SEED_A);
#endif
#endif
        for (k = 0; k < 8 * PARAMS_N; k++) {
            MatrixA[k] = LE_TO_UINT16(MatrixA[k]);
        }
    }
    uint16_t q = (1 << PARAMS_LOGQ);
    uint16_t res;
    for (k = 0; k < PARAMS_NBAR; k++) {
        const uint16_t *s_col_vec = &s[k * PARAMS_N];
        
        for (i = 0; i < PARAMS_N; i++) {
            const uint16_t *a_row_vec = &MatrixA[i * PARAMS_N];
            OP_vector_mul(&res, a_row_vec, s_col_vec, PARAMS_N, q);
            out[i * PARAMS_NBAR + k] += res;
        }
    }
    free(MatrixA);
    
#if defined(USE_AES128_FOR_A)
    AES128_free_schedule(aes_key_schedule);
#endif
    return 1;
}


int frodo_mul_add_sa_plus_e(uint16_t *out, const uint16_t *s, uint16_t *e, const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
  // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in: J.W. Bos, M. Ofner, J. Renes, 
  // T. Schneider, C. van Vredendaal, "The Matrix Reloaded: Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j, q, p; 
    size_t size_A = (size_t)PARAMS_N * PARAMS_N * sizeof(uint16_t);
    uint16_t *MatrixA = (uint16_t *)malloc(size_A);

#if defined(USE_AES128_FOR_A)
#if !defined(USE_OPENSSL)
    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);
#else
    EVP_CIPHER_CTX *aes_key_schedule;
    int len;
    if (!(aes_key_schedule = EVP_CIPHER_CTX_new())) handleErrors();
    if (1 != EVP_EncryptInit_ex(aes_key_schedule, EVP_aes_128_ecb(), NULL, seed_A, NULL)) handleErrors();
#endif
    // Initialize matrix used for encryption
    ALIGN_HEADER(32) uint16_t Ainit[PARAMS_N*8] ALIGN_FOOTER(32) = {0};
       
    for(j = 0; j < PARAMS_N; j+=8) {
        Ainit[0*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[1*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[2*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[3*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[4*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[5*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[6*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[7*PARAMS_N + j + 1] = UINT16_TO_LE(j);
    }

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i+=8) {
        // Generate 8 rows of A on-the-fly using AES
        for (q = 0; q < 8; q++) {
            for (p = 0; p < PARAMS_N; p+=8) {
                Ainit[q*PARAMS_N + p] = UINT16_TO_LE(i+q);
            }
        }

        size_t A_len = 8 * PARAMS_N * sizeof(uint16_t);
#if !defined(USE_OPENSSL)
        AES128_ECB_enc_sch((uint8_t*)Ainit, A_len, aes_key_schedule, (uint8_t*)(MatrixA + i*PARAMS_N));
#else   
        if (1 != EVP_EncryptUpdate(aes_key_schedule, (uint8_t*)(MatrixA + i*PARAMS_N), &len, (uint8_t*)Ainit, A_len)) handleErrors();
#endif

#elif defined (USE_SHAKE128_FOR_A)  // SHAKE128
#if !defined(USE_AVX2)
    uint8_t seed_A_separated[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin = (uint16_t*)&seed_A_separated;
    memcpy(&seed_A_separated[2], seed_A, BYTES_SEED_A);

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i+=8) {
        seed_A_origin[0] = UINT16_TO_LE(i + 0);
        shake128((unsigned char*)(MatrixA + (i+0)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 1);
        shake128((unsigned char*)(MatrixA + (i+1)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 2);
        shake128((unsigned char*)(MatrixA + (i+2)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 3);
        shake128((unsigned char*)(MatrixA + (i+3)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 4);
        shake128((unsigned char*)(MatrixA + (i+4)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 5);
        shake128((unsigned char*)(MatrixA + (i+5)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 6);
        shake128((unsigned char*)(MatrixA + (i+6)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 7);
        shake128((unsigned char*)(MatrixA + (i+7)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A); 
#else  // Using vector intrinsics
    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_2[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_3[2 + BYTES_SEED_A];
    uint16_t *seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t *seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    uint16_t *seed_A_origin_2 = (uint16_t*)&seed_A_separated_2;
    uint16_t *seed_A_origin_3 = (uint16_t*)&seed_A_separated_3;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_2[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_3[2], seed_A, BYTES_SEED_A);

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i+=8) {
        // Generate hash output
        // First 4 rows
        seed_A_origin_0[0] = UINT16_TO_LE(i + 0);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 1);
        seed_A_origin_2[0] = UINT16_TO_LE(i + 2);
        seed_A_origin_3[0] = UINT16_TO_LE(i + 3);
        shake128_4x((unsigned char*)(MatrixA + (i+0)*PARAMS_N), (unsigned char*)(MatrixA + (i+1)*PARAMS_N), (unsigned char*)(MatrixA + (i+2)*PARAMS_N), (unsigned char*)(MatrixA + (i+3)*PARAMS_N),
                    (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, seed_A_separated_2, seed_A_separated_3, 2 + BYTES_SEED_A);
        // Second 4 rows
        seed_A_origin_0[0] = UINT16_TO_LE(i + 4);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 5);
        seed_A_origin_2[0] = UINT16_TO_LE(i + 6);
        seed_A_origin_3[0] = UINT16_TO_LE(i + 7);
        shake128_4x((unsigned char*)(MatrixA + (i+4)*PARAMS_N), (unsigned char*)(MatrixA + (i+5)*PARAMS_N), (unsigned char*)(MatrixA + (i+6)*PARAMS_N), (unsigned char*)(MatrixA + (i+7)*PARAMS_N),
                   (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, seed_A_separated_2, seed_A_separated_3, 2 + BYTES_SEED_A);
    }
#endif
    }
#endif
    for(size_t k=0; k < (size_t)PARAMS_N * PARAMS_N; k++) {
        MatrixA[k] = LE_TO_UINT16(MatrixA[k]);
    }
    uint16_t *MatrixA_Trans = (uint16_t *)malloc(size_A);
    matrix_transpose(MatrixA, MatrixA_Trans, PARAMS_N, PARAMS_N);
    free(MatrixA);
    uint16_t q_mod = (1 << PARAMS_LOGQ);
    uint16_t res;

    for (j = 0; j < PARAMS_NBAR; j++) {
        const uint16_t *s_row = &s[j * PARAMS_N];
        
        for (i = 0; i < PARAMS_N; i++) {
            const uint16_t *at_row = &MatrixA_Trans[i * PARAMS_N];
            OP_vector_mul(&res, s_row, at_row, PARAMS_N, q_mod);
            e[j * PARAMS_N + i] = (e[j * PARAMS_N + i] + res) & (q_mod - 1);
        }
    }

    memcpy((unsigned char*)out, (unsigned char*)e, 2*PARAMS_N*PARAMS_NBAR);
    free(MatrixA_Trans);
#if defined(USE_AES128_FOR_A)
    AES128_free_schedule(aes_key_schedule);
#endif
    return 1;
}


void frodo_mul_bs(uint16_t *out, const uint16_t *b, const uint16_t *s) 
{ // Multiply by s on the right
  // Inputs: b (N_BAR x N), s (N x N_BAR)
  // Output: out = b*s (N_BAR x N_BAR)
    int i, j, k;
    uint16_t q_mod = (1 << PARAMS_LOGQ);
    uint16_t res;

    for (i = 0; i < PARAMS_NBAR; ++i) {
        for (j = 0; j < PARAMS_NBAR; ++j) {
            
            OP_vector_mul(&res, 
                          &b[i * PARAMS_N], 
                          &s[j * PARAMS_N], 
                          PARAMS_N, 
                          q_mod);
            
            out[i * PARAMS_NBAR + j] = res & (q_mod - 1);
        }
    }
}


void frodo_mul_add_sb_plus_e(uint16_t *out, const uint16_t *b, const uint16_t *s, const uint16_t *e) 
{ // Multiply by s on the left
  // Inputs: b (N x N_BAR), s (N_BAR x N), e (N_BAR x N_BAR)
  // Output: out = s*b + e (N_BAR x N_BAR)
    int i, j, k;
    size_t size_BT = (size_t)PARAMS_NBAR * PARAMS_N * sizeof(uint16_t);
    uint16_t *B_Trans = (uint16_t *)malloc(size_BT);
    matrix_transpose(b, B_Trans, PARAMS_N, PARAMS_NBAR);
    uint16_t q_mod = (1 << PARAMS_LOGQ);
    uint16_t res;
    for (k = 0; k < PARAMS_NBAR; k++) {
        for (i = 0; i < PARAMS_NBAR; i++) {
            OP_vector_mul(&res, 
                          &s[k * PARAMS_N], 
                          &B_Trans[i * PARAMS_N], 
                          PARAMS_N, 
                          q_mod);
            uint16_t val = e[k * PARAMS_NBAR + i] + res;
            out[k * PARAMS_NBAR + i] = val & (q_mod - 1);
        }
    }
    free(B_Trans);
}


void frodo_add(uint16_t *out, const uint16_t *a, const uint16_t *b) 
{ // Add a and b
  // Inputs: a, b (N_BAR x N_BAR)
  // Output: c = a + b

    for (int i = 0; i < (PARAMS_NBAR*PARAMS_NBAR); i++) {
        out[i] = (a[i] + b[i]) & ((1<<PARAMS_LOGQ)-1);
    }
}


void frodo_sub(uint16_t *out, const uint16_t *a, const uint16_t *b) 
{ // Subtract a and b
  // Inputs: a, b (N_BAR x N_BAR)
  // Output: c = a - b

    for (int i = 0; i < (PARAMS_NBAR*PARAMS_NBAR); i++) {
        out[i] = (a[i] - b[i]) & ((1<<PARAMS_LOGQ)-1);
    }
}


void frodo_key_encode(uint16_t *out, const uint16_t *in) 
{ // Encoding
    unsigned int i, j, npieces_word = 8;
    unsigned int nwords = (PARAMS_NBAR*PARAMS_NBAR)/8;
    uint64_t temp, mask = ((uint64_t)1 << PARAMS_EXTRACTED_BITS) - 1;
    uint16_t* pos = out;

    for (i = 0; i < nwords; i++) {
        temp = 0;
        for(j = 0; j < PARAMS_EXTRACTED_BITS; j++) 
            temp |= ((uint64_t)((uint8_t*)in)[i*PARAMS_EXTRACTED_BITS + j]) << (8*j);
        for (j = 0; j < npieces_word; j++) { 
            *pos = (uint16_t)((temp & mask) << (PARAMS_LOGQ - PARAMS_EXTRACTED_BITS));  
            temp >>= PARAMS_EXTRACTED_BITS;
            pos++;
        }
    }
}


void frodo_key_decode(uint16_t *out, const uint16_t *in)
{ // Decoding
    unsigned int i, j, index = 0, npieces_word = 8;
    unsigned int nwords = (PARAMS_NBAR * PARAMS_NBAR) / 8;
    uint16_t temp, maskex=((uint16_t)1 << PARAMS_EXTRACTED_BITS) -1, maskq =((uint16_t)1 << PARAMS_LOGQ) -1;
    uint8_t  *pos = (uint8_t*)out;
    uint64_t templong;

    for (i = 0; i < nwords; i++) {
        templong = 0;
        for (j = 0; j < npieces_word; j++) {  // temp = floor(in*2^{-11}+0.5)
            temp = ((in[index] & maskq) + (1 << (PARAMS_LOGQ - PARAMS_EXTRACTED_BITS - 1))) >> (PARAMS_LOGQ - PARAMS_EXTRACTED_BITS);
            templong |= ((uint64_t)(temp & maskex)) << (PARAMS_EXTRACTED_BITS * j);
            index++;
        }
	for(j = 0; j < PARAMS_EXTRACTED_BITS; j++) 
	    pos[i*PARAMS_EXTRACTED_BITS + j] = (templong >> (8*j)) & 0xFF;
    }
}
