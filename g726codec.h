/*
 * Cut from ffmpeg, only suitable for 16 bits wide
 *
 * CARE: It's different from CCITT's algorithm
 * */


#ifndef __G726_FFMPEG_H__
#define __G726_FFMPEG_H__
#include <stdint.h>
#if 0
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long int uint64_t;
#endif
#define INT_MAX 2147483647
#define FFMAX(a,b) ((a) > (b) ? (a) : (b))
#define av_log2_16bit(x) (31 - __builtin_clz((x)|1))
#define BUF_BITS (8 * sizeof(uint32_t))

#define AV_WL32(p, val) do {                 \
        uint32_t d = (val);                     \
        ((uint8_t*)(p))[0] = (d);               \
        ((uint8_t*)(p))[1] = (d)>>8;            \
        ((uint8_t*)(p))[2] = (d)>>16;           \
        ((uint8_t*)(p))[3] = (d)>>24;           \
    } while(0)
#define AV_RL32(x)                                \
    (((uint32_t)((const uint8_t*)(x))[3] << 24) |    \
               (((const uint8_t*)(x))[2] << 16) |    \
               (((const uint8_t*)(x))[1] <<  8) |    \
                ((const uint8_t*)(x))[0])

/**
 * G.726 11-bit float.
 * G.726 Standard uses rather odd 11-bit floating point arithmetic for
 * numerous occasions. It's a mystery to me why they did it this way
 * instead of simply using 32-bit integer arithmetic.
 */
typedef struct Float11 {
    uint8_t sign;   /**< 1 bit sign */
    uint8_t exp;    /**< 4 bits exponent */
    uint8_t mant;   /**< 6 bits mantissa */
} Float11;

typedef struct G726Tables {
    const int* quant;         /**< quantization table */
    const int16_t* iquant;    /**< inverse quantization table */
    const int16_t* W;         /**< special table #1 ;-) */
    const uint8_t* F;         /**< special table #2 */
} G726Tables;

static const int quant_tbl16[] =                  /**< 16kbit/s 2 bits per sample */
           { 260, INT_MAX };
static const int16_t iquant_tbl16[] =
           { 116, 365, 365, 116 };
static const int16_t W_tbl16[] =
           { -22, 439, 439, -22 };
static const uint8_t F_tbl16[] =
           { 0, 7, 7, 0 };
static const G726Tables G726Tables_pool[] =
           {{ quant_tbl16, iquant_tbl16, W_tbl16, F_tbl16 }};

typedef struct G726Context {
    G726Tables tbls;    /**< static tables needed for computation */

    Float11 sr[2];      /**< prev. reconstructed samples */
    Float11 dq[6];      /**< prev. difference */
    int a[2];           /**< second order predictor coeffs */
    int b[6];           /**< sixth order predictor coeffs */
    int pk[2];          /**< signs of prev. 2 sez + dq */

    int ap;             /**< scale factor control */
    int yu;             /**< fast scale factor */
    int yl;             /**< slow scale factor */
    int dms;            /**< short average magnitude of F[i] */
    int dml;            /**< long average magnitude of F[i] */
    int td;             /**< tone detect */

    int se;             /**< estimated signal for the next iteration */
    int sez;            /**< estimated second order prediction */
    int y;              /**< quantizer scaling factor for the next iteration */
    int code_size;
    int little_endian;  /**< little-endian bitstream as used in aiff and Sun AU */
} G726Context;

typedef struct PutBitContext {
    uint32_t bit_buf;
    int bit_left;
    uint8_t *buf, *buf_ptr, *buf_end;
} PutBitContext;

typedef struct GetBitContext {
    const uint8_t *buffer, *buffer_end;
    int index;
    int size_in_bits;
    int size_in_bits_plus8;
} GetBitContext;


static inline unsigned int zero_extend(unsigned int val, unsigned int bits)
{
    return (val << ((8 * sizeof(int)) - bits)) >> ((8 * sizeof(int)) - bits);
}

static inline Float11* i2f(int i, Float11* f)
{
	f->sign = (i < 0);
	if (f->sign)
		i = -i;
	f->exp = av_log2_16bit(i) + !!i;
	f->mant = i? (i<<6) >> f->exp : 1<<5;
	return f;
}

static inline int16_t mult(Float11* f1, Float11* f2)
{
	int res, exp;

	exp = f1->exp + f2->exp;
	res = (((f1->mant * f2->mant) + 0x30) >> 4);
	res = exp > 19 ? res << (exp - 19) : res >> (19 - exp);
	return (f1->sign ^ f2->sign) ? -res : res;
}

static inline int sgn(int value)
{
	return (value < 0) ? -1 : 1;
}
/**
 * Clip a signed integer value into the amin-amax range.
 * @param a value to clip
 * @param amin minimum value of the clip range
 * @param amax maximum value of the clip range
 * @return clipped value
 */
static inline int av_clip(int a, int amin, int amax)
{
    if      (a < amin) return amin;
    else if (a > amax) return amax;
    else               return a;
}
/**
 * Clear high bits from an unsigned integer starting with specific bit position
 * @param  a value to clip
 * @param  p bit position to clip at
 * @return clipped value
 */
static inline uint8_t av_mod_uintp2(unsigned int a, unsigned int p)
{
    return a & ((1U << p) - 1);
}




int g726_encode_init(G726Context **c);

/**
 * Encode some frames
 * @param  out 		Data after encoding
 * @param  inS16Len	Data len before encoding
 * @param  in		Data  before encoding
 * @return clipped value
 */
int g726_encode_frame(G726Context **c, uint8_t *out, int inS16Len, const int16_t *in);

int g726_encode_destroy(G726Context **c);


int g726_decode_init(G726Context **c);
/**
 * Decode some frames
 * @param  in		Data  before decoding
 * @param  inU8Len	Data len before decoding
 * @param  out 		Data after decoding
 * @return clipped value
 */
int g726_decode_frame(G726Context **c, const uint8_t *in, int inU8Len, int16_t *out);

int g726_decode_destroy(G726Context **c);

void g726_decode_flush(G726Context **c);
#endif
