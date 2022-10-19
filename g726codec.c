#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "g726codec.h"


/**
 * Paragraph 4.2.2 page 18: Adaptive quantizer.
 */
static inline uint8_t quant(G726Context* c, int d)
{
	int sign, exp, i, dln;

	sign = i = 0;
	if (d < 0) {
		sign = 1;
		d = -d;
	}
    exp = av_log2_16bit(d);
    dln = ((exp<<7) + (((d<<7)>>exp)&0x7f)) - (c->y>>2);

    while (c->tbls.quant[i] < INT_MAX && c->tbls.quant[i] < dln)
        ++i;

    if (sign)
        i = ~i;
    if (c->code_size != 2 && i == 0) /* I'm not sure this is a good idea */
        i = 0xff;

    return i;
}

/**
 * Paragraph 4.2.3 page 22: Inverse adaptive quantizer.
 */
static inline int16_t inverse_quant(G726Context* c, int i)
{
    int dql, dex, dqt;

    dql = c->tbls.iquant[i] + (c->y >> 2);
    dex = (dql>>7) & 0xf;        /* 4-bit exponent */
    dqt = (1<<7) + (dql & 0x7f); /* log2 -> linear */
    return (dql < 0) ? 0 : ((dqt<<dex) >> 7);
}

/**
 * Clip a signed integer into the -(2^p),(2^p-1) range.
 * @param  a value to clip
 * @param  p bit position to clip at
 * @return clipped value
 */
static inline int av_clip_intp2(int a, int p)
{
    if (((unsigned int)a + (1 << p)) & ~((2 << p) - 1))
        return (a >> 31) ^ ((1 << p) - 1);
    else
        return a;
}

static int16_t g726_decode(G726Context* c, int I)
{
    int dq = 0, re_signal = 0, pk0 = 0, fa1 = 0, i = 0, tr = 0, ylint = 0, ylfrac = 0, thr2 = 0, al = 0, dq0 = 0;
    Float11 f;
    int I_sig = I >> (c->code_size - 1);

	dq = inverse_quant(c, I);

    /* Transition detect */
    ylint = (c->yl >> 15);
    ylfrac = (c->yl >> 10) & 0x1f;
    thr2 = (ylint > 9) ? 0x1f << 10 : (0x20 + ylfrac) << ylint;
    tr= (c->td == 1 && dq > ((3*thr2)>>2));

    if (I_sig)  /* get the sign */
        dq = -dq;
    re_signal = (int16_t)(c->se + dq);

    /* Update second order predictor coefficient A2 and A1 */
    pk0 = (c->sez + dq) ? sgn(c->sez + dq) : 0;
    dq0 = dq ? sgn(dq) : 0;
    if (tr) {
        c->a[0] = 0;
        c->a[1] = 0;
        for (i=0; i<6; i++)
            c->b[i] = 0;
    } else {
        /* This is a bit crazy, but it really is +255 not +256 */
        fa1 = av_clip_intp2((-c->a[0]*c->pk[0]*pk0)>>5, 8);
        c->a[1] += 128*pk0*c->pk[1] + fa1 - (c->a[1]>>7);
        c->a[1] = av_clip(c->a[1], -12288, 12288);
        c->a[0] += 64*3*pk0*c->pk[0] - (c->a[0] >> 8);
        c->a[0] = av_clip(c->a[0], -(15360 - c->a[1]), 15360 - c->a[1]);

        for (i=0; i<6; i++)
            c->b[i] += 128*dq0*sgn(-c->dq[i].sign) - (c->b[i]>>8);
    }

    /* Update Dq and Sr and Pk */
    c->pk[1] = c->pk[0];
    c->pk[0] = pk0 ? pk0 : 1;
    c->sr[1] = c->sr[0];
    i2f(re_signal, &c->sr[0]);
    for (i=5; i>0; i--)
        c->dq[i] = c->dq[i-1];
    i2f(dq, &c->dq[0]);
    c->dq[0].sign = I_sig; /* Isn't it crazy ?!?! */

    c->td = c->a[1] < -11776;

    /* Update Ap */
    c->dms += (c->tbls.F[I]<<4) + ((- c->dms) >> 5);
    c->dml += (c->tbls.F[I]<<4) + ((- c->dml) >> 7);
    if (tr)
        c->ap = 256;
    else {
        c->ap += (-c->ap) >> 4;
        if (c->y <= 1535 || c->td || abs((c->dms << 2) - c->dml) >= (c->dml >> 3))
            c->ap += 0x20;
    }

    /* Update Yu and Yl */
    c->yu = av_clip(c->y + c->tbls.W[I] + ((-c->y)>>5), 544, 5120);
    c->yl += c->yu + ((-c->yl)>>6);

    /* Next iteration for Y */
    al = (c->ap >= 256) ? 1<<6 : c->ap >> 2;
    c->y = (c->yl + (c->yu - (c->yl>>6))*al) >> 6;

    /* Next iteration for SE and SEZ */
    c->se = 0;
    for (i=0; i<6; i++)
        c->se += mult(i2f(c->b[i] >> 2, &f), &c->dq[i]);
    c->sez = c->se >> 1;
    for (i=0; i<2; i++)
        c->se += mult(i2f(c->a[i] >> 2, &f), &c->sr[i]);
    c->se >>= 1;

	return av_clip(re_signal * 4, -0xffff, 0xffff);
}

static int g726_reset(G726Context *c)
{
    int i;

    c->tbls = G726Tables_pool[c->code_size - 2];
    for (i=0; i<2; i++) {
        c->sr[i].mant = 1<<5;
        c->pk[i] = 1;

		c->a[i] = 0;
    }
    for (i=0; i<6; i++) {
        c->dq[i].mant = 1<<5;

		c->dq[i].sign = 0;
        c->dq[i].exp = 0;
		c->b[i] = 0;
    }
    c->yu = 544;
    c->yl = 34816;

    c->y = 544;


	c->sez = 0;
	c->ap = 0;
	c->dms = 0;
	c->dml = 0;
	c->td = 0;
	c->se = 0;

    return 0;
}

static int16_t g726_encode(G726Context* c, int16_t sig)
{
    uint8_t i;

    i = av_mod_uintp2(quant(c, sig/4 - c->se), c->code_size);
    g726_decode(c, i);

    return i;
}



/**
 * Initialize the PutBitContext s.
 *
 * @param buffer the buffer where to put bits
 * @param buffer_size the size in bytes of buffer
 */
static inline void init_put_bits(PutBitContext *s, uint8_t *buffer, int buffer_size)
{
    if (buffer_size < 0) {
        buffer_size = 0;
        buffer      = NULL;
    }

    s->buf          = buffer;
    s->buf_end      = s->buf + buffer_size;
    s->buf_ptr      = s->buf;
    s->bit_left     = BUF_BITS;
    s->bit_buf      = 0;
}

static inline void put_bits_le(PutBitContext *s, int n, uint32_t value)
{
    uint32_t bit_buf;
    int bit_left;

    bit_buf  = s->bit_buf;
    bit_left = s->bit_left;

    bit_buf |= value << (BUF_BITS - bit_left);
    if (n >= bit_left) {
        if (s->buf_end - s->buf_ptr >= sizeof(uint32_t)) {
            AV_WL32(s->buf_ptr, bit_buf);
			s->buf_ptr += sizeof(uint32_t);
        } else {
            printf("Internal error, put_bits buffer too small\n");
            return ;
        }
        bit_buf     = value >> bit_left;
        bit_left   += BUF_BITS;
    }
    bit_left -= n;

    s->bit_buf  = bit_buf;
    s->bit_left = bit_left;
}


static inline void flush_put_bits_le(PutBitContext *s)
{
    while (s->bit_left < BUF_BITS) {
        //av_assert0(s->buf_ptr < s->buf_end);
        *s->buf_ptr++ = s->bit_buf;
        s->bit_buf  >>= 8;
        s->bit_left  += 8;
    }
    s->bit_left = BUF_BITS;
    s->bit_buf  = 0;
}

/* Interfacing to the libavcodec */
int g726_encode_init(G726Context **context)
{
	G726Context *c = NULL;
	c = (G726Context *)malloc(sizeof(G726Context));
	if(!c){
		printf("param c is NULL\n");
		return -1;
	}

    c->little_endian = 1;
    c->code_size = 2;

    g726_reset(c);
	*context = c;

    return 0;
}

int g726_encode_frame(G726Context **context, uint8_t *out, int inS16Len, const int16_t *in)
{
	PutBitContext pb;
	int i = 0;

	G726Context *c = *context;
	if(!c){
		return -1;
	}

	if(!in || (inS16Len <= 0) || !out || (in && inS16Len <= 0)){
		printf(" in or inS16Len or out is NULL\n");
		return -1;
	}

	init_put_bits(&pb, out, inS16Len);
	for(i = 0; i < inS16Len; i++)
		put_bits_le(&pb, 2, g726_encode(c, *in++));

	flush_put_bits_le(&pb);

	return inS16Len/4;
}

int g726_encode_destroy(G726Context **context)
{
	G726Context *c = *context;
	if(!c){
		printf("param c is already NULL\n");
		return 0;
	}

	free(c);
    return 0;
}


/**
 * Initialize GetBitContext.
 * @param buffer bitstream buffer, must be AV_INPUT_BUFFER_PADDING_SIZE bytes
 *        larger than the actual read bits because some optimized bitstream
 *        readers read 32 or 64 bit at once and could read over the end
 * @param bit_size the size of the buffer in bits
 * @return 0 on success, AVERROR_INVALIDDATA if the buffer_size would overflow.
 */
static inline int init_get_bits(GetBitContext *s, const uint8_t *buffer,int bit_size)
{
    int buffer_size;
    int ret = 0;

    if (bit_size >= INT_MAX - FFMAX(7, 64*8) || bit_size < 0 || !buffer) {
        bit_size    = 0;
        buffer      = NULL;
        ret         = -1;
    }

    buffer_size = (bit_size + 7) >> 3;

    s->buffer             = buffer;
    s->size_in_bits       = bit_size;
    s->size_in_bits_plus8 = bit_size + 8;
    s->buffer_end         = buffer + buffer_size;
    s->index              = 0;
    return ret;
}

static inline unsigned int get_bits_le(GetBitContext *s, int n)
{
	int tmp;
	unsigned int index = s->index;
	unsigned int cache = AV_RL32(s->buffer + (index >> 3)) >> (index & 7);
	tmp = zero_extend(cache, n);
	index += n;
	s->index = index;

	return tmp;
}

int g726_decode_init(G726Context **context)
{
	G726Context *c = NULL;
	c = (G726Context *)malloc(sizeof(G726Context));
	if(!c){
		printf("param c is NULL\n");
		return -1;
	}

    c->little_endian = 1;
    c->code_size = 2;

    g726_reset(c);

	*context = c;
    return 0;
}

int g726_decode_frame(G726Context **context, const uint8_t *in, int inU8Len, int16_t *out)
{
	int out_samples;
	GetBitContext gb;
	int16_t *tmp = out;
	G726Context *c = *context;
	if(!c){
		printf(" c is NULL\n");
		return -1;
	}

	if(!in || (inU8Len <= 0) || !out || (in && inU8Len <= 0)){
		printf(" in or inU8Len or out is NULL\n");
		return -1;
	}

	if(-1 == init_get_bits(&gb, in, inU8Len * 8)){
		printf("init_get_bits failed!\n");
		return -1;
	}
	out_samples = inU8Len * 8 / c->code_size;
	while(out_samples--)
		*out++ = g726_decode(c, get_bits_le(&gb, c->code_size));

	out = tmp;

	return inU8Len;
}


void g726_decode_flush(G726Context **context)
{
	G726Context *c = *context;
	if(!c){
		printf(" c is NULL\n");
		return ;
	}

    g726_reset(c);
}

int g726_decode_destroy(G726Context **context)
{
	G726Context *c = *context;
	if(!c){
		printf("param c is already NULL\n");
		return 0;
	}

	free(c);
	c = NULL;

    return 0;
}

/******************  Next is the g726 codec test ********************/

//#define g726_test		//>>>Modify this comment to start or close the test

#ifdef g726_test
#define encodec_or_decodec 0	//>>>Modify this define to change type of test: \
												0:encodec;		1:decodec
int main(int argc, char **argv)
{
	const char *outfilename, *filename;
	if (argc <= 2) {
		fprintf(stderr, "Usage: %s <input file> <output file>\n", argv[0]);
		exit(0);
	}
	filename    = argv[1];
	outfilename = argv[2];

	FILE *f, *outfile;
	f = fopen(filename, "rb");
	if (!f) {
		fprintf(stderr, "Could not open %s\n", filename);
		exit(1);
	}
	outfile = fopen(outfilename, "wb");
	if (!outfile) {
		exit(1);
	}

	int   data_size;
	uint8_t *inbuf = (uint8_t *)malloc(8000 * sizeof(short) * 4096 * 4/ 1000/ 8);
	if(inbuf == NULL){
		fprintf(stderr, "Could not alloc inbuf \n");
		return -1;
	}

	int16_t *data = (int16_t *)malloc(8000 * sizeof(short) * 4096 *4/ 1000);
	if(data == NULL){
		fprintf(stderr, "Could not alloc inbuf \n");
		return -1;
	}

	int ret = 0;
	G726Context *handle = NULL;

#if encodec_or_decodec
//解码
	if(0 != g726_decode_init(&handle))
	{
		fprintf(stderr, "Could not init decodec \n");
		return -1;
	}

	while(1)
	{
		data_size = fread(inbuf, 1, 8000*sizeof(short) * 40 / 1000/8, f);
		if(data_size <= 0){
			fprintf(stderr, "Could not read :%s\n", strerror(errno));
			break;
		}
#if 0
		int i = 0;
		for(i = 0; i < data_size; i++){
			if(i % 8 == 0)
				printf("\n");
			printf("%02x ", inbuf[i]);
		}
		printf("\n");
#endif//END of printf

		ret = g726_decode_frame(&handle, (const uint8_t *)inbuf, data_size, (int16_t *)data);
		if(ret != data_size){
			printf("ret(%d) != datasize(%d)\n", ret, data_size);
			return -1;
		}

#if 0
		for(i = 0; i < ret*4; i++){
			if(i % 8 == 0)
				printf("\n");
			printf("%08x ", data[i]);
		}
		printf("\n");
#endif//END of printf

		ret = fwrite(data, 2, ret*4, outfile);
	}

	g726_decode_destroy(&handle);

#else//ELSE OF encodec_or_decodec 
//编码
	if(0 != g726_encode_init(&handle))
	{
		fprintf(stderr, "Could not init encodec \n");
		return -1;
	}

	while(1)
	{
		data_size = fread(data, 1, 640, f);
		if(data_size <= 0){
			fprintf(stderr, "Could not read :%s\n", strerror(errno));
			break;
		}

		int i = 0;
#if 0
		for(i = 0; i < data_size/2; i++){
			if(i % 8 == 0)
				printf("\n");
			printf("%04x ", (uint16_t)data[i]);
		}

		printf("\n");
#endif//END of printf

		ret = g726_encode_frame(&handle, (uint8_t *)inbuf, data_size/2, (const int16_t *)data);
		if(ret != data_size/8)
		{
			printf("ret(%d) != datasize(%d)\n", ret, data_size/8);
			return -1;
		}
#if 0
		for(i = 0; i < ret; i++){
			if(i % 8 == 0)
				printf("\n");
			printf("%02x ", inbuf[i]);
		}
		printf("\n");
#endif//END of printf

		ret = fwrite(inbuf, 1, data_size/8, outfile);
	}

	g726_encode_destroy(&handle);

#endif//END OF encodec_or_decodec

	fclose(f);
	fclose(outfile);

	free(data);
	free(inbuf);
	return 0;
}
#endif//END OF g726_test
