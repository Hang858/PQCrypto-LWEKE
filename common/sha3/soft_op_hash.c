#include <stdint.h>
#include <stddef.h>
#include <string.h> 
#include "hardware_driver.h" 

#define OPH_NROUNDS 24
#define OPH_ROL(a, offset) ((a << offset) ^ (a >> (64 - offset)))

#define OPH_SHAKE128_RATE 168
#define OPH_SHAKE256_RATE 136
#define OPH_SHA3_256_RATE 136
#define OPH_SHA3_512_RATE 72

// --- 内部辅助函数 ---

static uint64_t oph_load64(const uint8_t x[8]) {
    unsigned int i;
    uint64_t r = 0;
    for (i = 0; i < 8; i++)
        r |= (uint64_t)x[i] << 8 * i;
    return r;
}

static void oph_store64(uint8_t x[8], uint64_t u) {
    unsigned int i;
    for (i = 0; i < 8; i++)
        x[i] = u >> 8 * i;
}

static const uint64_t OPH_KeccakF_RoundConstants[OPH_NROUNDS] = {
    (uint64_t)0x0000000000000001ULL, (uint64_t)0x0000000000008082ULL,
    (uint64_t)0x800000000000808aULL, (uint64_t)0x8000000080008000ULL,
    (uint64_t)0x000000000000808bULL, (uint64_t)0x0000000080000001ULL,
    (uint64_t)0x8000000080008081ULL, (uint64_t)0x8000000000008009ULL,
    (uint64_t)0x000000000000008aULL, (uint64_t)0x0000000000000088ULL,
    (uint64_t)0x0000000080008009ULL, (uint64_t)0x000000008000000aULL,
    (uint64_t)0x000000008000808bULL, (uint64_t)0x800000000000008bULL,
    (uint64_t)0x8000000000008089ULL, (uint64_t)0x8000000000008003ULL,
    (uint64_t)0x8000000000008002ULL, (uint64_t)0x8000000000000080ULL,
    (uint64_t)0x000000000000800aULL, (uint64_t)0x800000008000000aULL,
    (uint64_t)0x8000000080008081ULL, (uint64_t)0x8000000000008080ULL,
    (uint64_t)0x0000000080000001ULL, (uint64_t)0x8000000080008008ULL};

static void oph_permute(uint64_t state[25]) {
    int round;
    uint64_t Aba, Abe, Abi, Abo, Abu;
    uint64_t Aga, Age, Agi, Ago, Agu;
    uint64_t Aka, Ake, Aki, Ako, Aku;
    uint64_t Ama, Ame, Ami, Amo, Amu;
    uint64_t Asa, Ase, Asi, Aso, Asu;
    uint64_t BCa, BCe, BCi, BCo, BCu;
    uint64_t Da, De, Di, Do, Du;
    uint64_t Eba, Ebe, Ebi, Ebo, Ebu;
    uint64_t Ega, Ege, Egi, Ego, Egu;
    uint64_t Eka, Eke, Eki, Eko, Eku;
    uint64_t Ema, Eme, Emi, Emo, Emu;
    uint64_t Esa, Ese, Esi, Eso, Esu;

    Aba = state[0];  Abe = state[1];  Abi = state[2];  Abo = state[3];  Abu = state[4];
    Aga = state[5];  Age = state[6];  Agi = state[7];  Ago = state[8];  Agu = state[9];
    Aka = state[10]; Ake = state[11]; Aki = state[12]; Ako = state[13]; Aku = state[14];
    Ama = state[15]; Ame = state[16]; Ami = state[17]; Amo = state[18]; Amu = state[19];
    Asa = state[20]; Ase = state[21]; Asi = state[22]; Aso = state[23]; Asu = state[24];

    for (round = 0; round < OPH_NROUNDS; round += 2) {
        BCa = Aba ^ Aga ^ Aka ^ Ama ^ Asa;
        BCe = Abe ^ Age ^ Ake ^ Ame ^ Ase;
        BCi = Abi ^ Agi ^ Aki ^ Ami ^ Asi;
        BCo = Abo ^ Ago ^ Ako ^ Amo ^ Aso;
        BCu = Abu ^ Agu ^ Aku ^ Amu ^ Asu;
        Da = BCu ^ OPH_ROL(BCe, 1);
        De = BCa ^ OPH_ROL(BCi, 1);
        Di = BCe ^ OPH_ROL(BCo, 1);
        Do = BCi ^ OPH_ROL(BCu, 1);
        Du = BCo ^ OPH_ROL(BCa, 1);

        Aba ^= Da; BCa = Aba;
        Age ^= De; BCe = OPH_ROL(Age, 44);
        Aki ^= Di; BCi = OPH_ROL(Aki, 43);
        Amo ^= Do; BCo = OPH_ROL(Amo, 21);
        Asu ^= Du; BCu = OPH_ROL(Asu, 14);
        Eba = BCa ^ ((~BCe) & BCi);
        Eba ^= (uint64_t)OPH_KeccakF_RoundConstants[round];
        Ebe = BCe ^ ((~BCi) & BCo);
        Ebi = BCi ^ ((~BCo) & BCu);
        Ebo = BCo ^ ((~BCu) & BCa);
        Ebu = BCu ^ ((~BCa) & BCe);

        Abo ^= Do; BCa = OPH_ROL(Abo, 28);
        Agu ^= Du; BCe = OPH_ROL(Agu, 20);
        Aka ^= Da; BCi = OPH_ROL(Aka, 3);
        Ame ^= De; BCo = OPH_ROL(Ame, 45);
        Asi ^= Di; BCu = OPH_ROL(Asi, 61);
        Ega = BCa ^ ((~BCe) & BCi);
        Ege = BCe ^ ((~BCi) & BCo);
        Egi = BCi ^ ((~BCo) & BCu);
        Ego = BCo ^ ((~BCu) & BCa);
        Egu = BCu ^ ((~BCa) & BCe);

        Abe ^= De; BCa = OPH_ROL(Abe, 1);
        Agi ^= Di; BCe = OPH_ROL(Agi, 6);
        Ako ^= Do; BCi = OPH_ROL(Ako, 25);
        Amu ^= Du; BCo = OPH_ROL(Amu, 8);
        Asa ^= Da; BCu = OPH_ROL(Asa, 18);
        Eka = BCa ^ ((~BCe) & BCi);
        Eke = BCe ^ ((~BCi) & BCo);
        Eki = BCi ^ ((~BCo) & BCu);
        Eko = BCo ^ ((~BCu) & BCa);
        Eku = BCu ^ ((~BCa) & BCe);

        Abu ^= Du; BCa = OPH_ROL(Abu, 27);
        Aga ^= Da; BCe = OPH_ROL(Aga, 36);
        Ake ^= De; BCi = OPH_ROL(Ake, 10);
        Ami ^= Di; BCo = OPH_ROL(Ami, 15);
        Aso ^= Do; BCu = OPH_ROL(Aso, 56);
        Ema = BCa ^ ((~BCe) & BCi);
        Eme = BCe ^ ((~BCi) & BCo);
        Emi = BCi ^ ((~BCo) & BCu);
        Emo = BCo ^ ((~BCu) & BCa);
        Emu = BCu ^ ((~BCa) & BCe);

        Abi ^= Di; BCa = OPH_ROL(Abi, 62);
        Ago ^= Do; BCe = OPH_ROL(Ago, 55);
        Aku ^= Du; BCi = OPH_ROL(Aku, 39);
        Ama ^= Da; BCo = OPH_ROL(Ama, 41);
        Ase ^= De; BCu = OPH_ROL(Ase, 2);
        Esa = BCa ^ ((~BCe) & BCi);
        Ese = BCe ^ ((~BCi) & BCo);
        Esi = BCi ^ ((~BCo) & BCu);
        Eso = BCo ^ ((~BCu) & BCa);
        Esu = BCu ^ ((~BCa) & BCe);

        BCa = Eba ^ Ega ^ Eka ^ Ema ^ Esa;
        BCe = Ebe ^ Ege ^ Eke ^ Eme ^ Ese;
        BCi = Ebi ^ Egi ^ Eki ^ Emi ^ Esi;
        BCo = Ebo ^ Ego ^ Eko ^ Emo ^ Eso;
        BCu = Ebu ^ Egu ^ Eku ^ Emu ^ Esu;

        Da = BCu ^ OPH_ROL(BCe, 1);
        De = BCa ^ OPH_ROL(BCi, 1);
        Di = BCe ^ OPH_ROL(BCo, 1);
        Do = BCi ^ OPH_ROL(BCu, 1);
        Du = BCo ^ OPH_ROL(BCa, 1);

        Eba ^= Da; BCa = Eba;
        Ege ^= De; BCe = OPH_ROL(Ege, 44);
        Eki ^= Di; BCi = OPH_ROL(Eki, 43);
        Emo ^= Do; BCo = OPH_ROL(Emo, 21);
        Esu ^= Du; BCu = OPH_ROL(Esu, 14);
        Aba = BCa ^ ((~BCe) & BCi);
        Aba ^= (uint64_t)OPH_KeccakF_RoundConstants[round + 1];
        Abe = BCe ^ ((~BCi) & BCo);
        Abi = BCi ^ ((~BCo) & BCu);
        Abo = BCo ^ ((~BCu) & BCa);
        Abu = BCu ^ ((~BCa) & BCe);

        Ebo ^= Do; BCa = OPH_ROL(Ebo, 28);
        Egu ^= Du; BCe = OPH_ROL(Egu, 20);
        Eka ^= Da; BCi = OPH_ROL(Eka, 3);
        Eme ^= De; BCo = OPH_ROL(Eme, 45);
        Esi ^= Di; BCu = OPH_ROL(Esi, 61);
        Aga = BCa ^ ((~BCe) & BCi);
        Age = BCe ^ ((~BCi) & BCo);
        Agi = BCi ^ ((~BCo) & BCu);
        Ago = BCo ^ ((~BCu) & BCa);
        Agu = BCu ^ ((~BCa) & BCe);

        Ebe ^= De; BCa = OPH_ROL(Ebe, 1);
        Egi ^= Di; BCe = OPH_ROL(Egi, 6);
        Eko ^= Do; BCi = OPH_ROL(Eko, 25);
        Emu ^= Du; BCo = OPH_ROL(Emu, 8);
        Esa ^= Da; BCu = OPH_ROL(Esa, 18);
        Aka = BCa ^ ((~BCe) & BCi);
        Ake = BCe ^ ((~BCi) & BCo);
        Aki = BCi ^ ((~BCo) & BCu);
        Ako = BCo ^ ((~BCu) & BCa);
        Aku = BCu ^ ((~BCa) & BCe);

        Ebu ^= Du; BCa = OPH_ROL(Ebu, 27);
        Ega ^= Da; BCe = OPH_ROL(Ega, 36);
        Eke ^= De; BCi = OPH_ROL(Eke, 10);
        Emi ^= Di; BCo = OPH_ROL(Emi, 15);
        Eso ^= Do; BCu = OPH_ROL(Eso, 56);
        Ama = BCa ^ ((~BCe) & BCi);
        Ame = BCe ^ ((~BCi) & BCo);
        Ami = BCi ^ ((~BCo) & BCu);
        Amo = BCo ^ ((~BCu) & BCa);
        Amu = BCu ^ ((~BCa) & BCe);

        Ebi ^= Di; BCa = OPH_ROL(Ebi, 62);
        Ego ^= Do; BCe = OPH_ROL(Ego, 55);
        Eku ^= Du; BCi = OPH_ROL(Eku, 39);
        Ema ^= Da; BCo = OPH_ROL(Ema, 41);
        Ese ^= De; BCu = OPH_ROL(Ese, 2);
        Asa = BCa ^ ((~BCe) & BCi);
        Ase = BCe ^ ((~BCi) & BCo);
        Asi = BCi ^ ((~BCo) & BCu);
        Aso = BCo ^ ((~BCu) & BCa);
        Asu = BCu ^ ((~BCa) & BCe);
    }

    state[0] = Aba;  state[1] = Abe;  state[2] = Abi;  state[3] = Abo;  state[4] = Abu;
    state[5] = Aga;  state[6] = Age;  state[7] = Agi;  state[8] = Ago;  state[9] = Agu;
    state[10] = Aka; state[11] = Ake; state[12] = Aki; state[13] = Ako; state[14] = Aku;
    state[15] = Ama; state[16] = Ame; state[17] = Ami; state[18] = Amo; state[19] = Amu;
    state[20] = Asa; state[21] = Ase; state[22] = Asi; state[23] = Aso; state[24] = Asu;
}

static void oph_init(uint64_t s[25]) {
    unsigned int i;
    for (i = 0; i < 25; i++)
        s[i] = 0;
}

static unsigned int oph_absorb(uint64_t s[25], unsigned int pos,
                                  unsigned int r, const uint8_t *in,
                                  size_t inlen) {
    unsigned int i;
    while (pos + inlen >= r) {
        for (i = pos; i < r; i++)
            s[i / 8] ^= (uint64_t)*in++ << 8 * (i % 8);
        inlen -= r - pos;
        oph_permute(s);
        pos = 0;
    }
    for (i = pos; i < pos + inlen; i++)
        s[i / 8] ^= (uint64_t)*in++ << 8 * (i % 8);
    return i;
}

static void oph_finalize(uint64_t s[25], unsigned int pos, unsigned int r, uint8_t p) {
    s[pos / 8] ^= (uint64_t)p << 8 * (pos % 8);
    s[r / 8 - 1] ^= 1ULL << 63;
}

static unsigned int oph_squeeze(uint8_t *out, size_t outlen, uint64_t s[25],
                                   unsigned int pos, unsigned int r) {
    unsigned int i;
    while (outlen) {
        if (pos == r) {
            oph_permute(s);
            pos = 0;
        }
        for (i = pos; i < r && i < pos + outlen; i++)
            *out++ = s[i / 8] >> 8 * (i % 8);
        outlen -= i - pos;
        pos = i;
    }
    return pos;
}

static void oph_absorb_once(uint64_t s[25], unsigned int r,
                               const uint8_t *in, size_t inlen, uint8_t p) {
    unsigned int i;
    oph_init(s);

    while (inlen >= r) {
        for (i = 0; i < r / 8; i++)
            s[i] ^= oph_load64(in + 8 * i);
        in += r;
        inlen -= r;
        oph_permute(s);
    }
    for (i = 0; i < inlen; i++)
        s[i / 8] ^= (uint64_t)in[i] << 8 * (i % 8);
    
    s[i / 8] ^= (uint64_t)p << 8 * (i % 8);
    s[(r - 1) / 8] ^= 1ULL << 63;
}


/**
 * @brief OP_hash 软件参考
 * 
 */
int OP_hash(uint8_t alg, uint8_t mode, int n, void *input, int input_len, void *output) {
    // 基础参数检查
    if (input == NULL || output == NULL || input_len < 0 || n <= 0) {
        return -1;
    }

    uint64_t s[25];
    uint8_t *out_ptr = (uint8_t *)output;
    const uint8_t *in_ptr = (const uint8_t *)input;
    unsigned int rate;
    uint8_t padding;

    // 1. 根据算法选择参数
    switch (alg) {
        case OP_ALG_SHAKE128:
            rate = OPH_SHAKE128_RATE;
            padding = 0x1F;
            break;
        case OP_ALG_SHAKE256:
            rate = OPH_SHAKE256_RATE;
            padding = 0x1F;
            break;
        case OP_ALG_SHA3_256:
            if (n != 32) return -1;
            rate = OPH_SHA3_256_RATE;
            padding = 0x06;
            break;
        case OP_ALG_SHA3_512:
            if (n != 64) return -1;
            rate = OPH_SHA3_512_RATE;
            padding = 0x06;
            break;
        default:
            return -1; // 不支持的算法
    }

    // 2. 核心逻辑执行
    
    if (alg == OP_ALG_SHA3_256 || alg == OP_ALG_SHA3_512) {
        oph_absorb_once(s, rate, in_ptr, (size_t)input_len, padding);
        oph_permute(s);
        
        unsigned int i;
        for (i = 0; i < (unsigned int)n / 8; i++) {
            oph_store64(out_ptr + 8 * i, s[i]);
        }
        return 0;
    }

    oph_init(s);
    unsigned int pos = 0;
    pos = oph_absorb(s, pos, rate, in_ptr, (size_t)input_len);
    oph_finalize(s, pos, rate, padding);
    
    pos = rate; 
    oph_squeeze(out_ptr, (size_t)n, s, pos, rate);

    return 0;
}