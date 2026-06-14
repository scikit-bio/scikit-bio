/*
 bit_array.c
 project: bit array C library
 url: https://github.com/noporpoise/BitArray/
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Aug 2014
*/

// 64 bit words
// Array length can be zero
// Unused top bits must be zero

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <limits.h> // ULONG_MAX
#include <errno.h>
#include <signal.h> // needed for abort()
#include <string.h> // memset()
#include <assert.h>
#include <time.h> // needed for seeding rand()
#include <unistd.h>  // need for getpid() for seeding rand number
#include <ctype.h>  // need for tolower()
#include <errno.h>  // perror()
#include <sys/time.h> // for seeding random

// Windows includes
#if defined(_WIN32)
#include <intrin.h>
#endif

#include "bit_array.h"
#include "bit_macros.h"

//
// Tables of constants
//

// byte reverse look up table
static const word_t reverse_table[256] =
{
  0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0,
  0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
  0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8,
  0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
  0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4,
  0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
  0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC,
  0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
  0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2,
  0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
  0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA,
  0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
  0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6,
  0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
  0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE,
  0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
  0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1,
  0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
  0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9,
  0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
  0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5,
  0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
  0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED,
  0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
  0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3,
  0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
  0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB,
  0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
  0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7,
  0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
  0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF,
  0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF,
};

// Morton table for interleaving bytes
static const word_t morton_table0[256] =
{
  0x0000, 0x0001, 0x0004, 0x0005, 0x0010, 0x0011, 0x0014, 0x0015,
  0x0040, 0x0041, 0x0044, 0x0045, 0x0050, 0x0051, 0x0054, 0x0055,
  0x0100, 0x0101, 0x0104, 0x0105, 0x0110, 0x0111, 0x0114, 0x0115,
  0x0140, 0x0141, 0x0144, 0x0145, 0x0150, 0x0151, 0x0154, 0x0155,
  0x0400, 0x0401, 0x0404, 0x0405, 0x0410, 0x0411, 0x0414, 0x0415,
  0x0440, 0x0441, 0x0444, 0x0445, 0x0450, 0x0451, 0x0454, 0x0455,
  0x0500, 0x0501, 0x0504, 0x0505, 0x0510, 0x0511, 0x0514, 0x0515,
  0x0540, 0x0541, 0x0544, 0x0545, 0x0550, 0x0551, 0x0554, 0x0555,
  0x1000, 0x1001, 0x1004, 0x1005, 0x1010, 0x1011, 0x1014, 0x1015,
  0x1040, 0x1041, 0x1044, 0x1045, 0x1050, 0x1051, 0x1054, 0x1055,
  0x1100, 0x1101, 0x1104, 0x1105, 0x1110, 0x1111, 0x1114, 0x1115,
  0x1140, 0x1141, 0x1144, 0x1145, 0x1150, 0x1151, 0x1154, 0x1155,
  0x1400, 0x1401, 0x1404, 0x1405, 0x1410, 0x1411, 0x1414, 0x1415,
  0x1440, 0x1441, 0x1444, 0x1445, 0x1450, 0x1451, 0x1454, 0x1455,
  0x1500, 0x1501, 0x1504, 0x1505, 0x1510, 0x1511, 0x1514, 0x1515,
  0x1540, 0x1541, 0x1544, 0x1545, 0x1550, 0x1551, 0x1554, 0x1555,
  0x4000, 0x4001, 0x4004, 0x4005, 0x4010, 0x4011, 0x4014, 0x4015,
  0x4040, 0x4041, 0x4044, 0x4045, 0x4050, 0x4051, 0x4054, 0x4055,
  0x4100, 0x4101, 0x4104, 0x4105, 0x4110, 0x4111, 0x4114, 0x4115,
  0x4140, 0x4141, 0x4144, 0x4145, 0x4150, 0x4151, 0x4154, 0x4155,
  0x4400, 0x4401, 0x4404, 0x4405, 0x4410, 0x4411, 0x4414, 0x4415,
  0x4440, 0x4441, 0x4444, 0x4445, 0x4450, 0x4451, 0x4454, 0x4455,
  0x4500, 0x4501, 0x4504, 0x4505, 0x4510, 0x4511, 0x4514, 0x4515,
  0x4540, 0x4541, 0x4544, 0x4545, 0x4550, 0x4551, 0x4554, 0x4555,
  0x5000, 0x5001, 0x5004, 0x5005, 0x5010, 0x5011, 0x5014, 0x5015,
  0x5040, 0x5041, 0x5044, 0x5045, 0x5050, 0x5051, 0x5054, 0x5055,
  0x5100, 0x5101, 0x5104, 0x5105, 0x5110, 0x5111, 0x5114, 0x5115,
  0x5140, 0x5141, 0x5144, 0x5145, 0x5150, 0x5151, 0x5154, 0x5155,
  0x5400, 0x5401, 0x5404, 0x5405, 0x5410, 0x5411, 0x5414, 0x5415,
  0x5440, 0x5441, 0x5444, 0x5445, 0x5450, 0x5451, 0x5454, 0x5455,
  0x5500, 0x5501, 0x5504, 0x5505, 0x5510, 0x5511, 0x5514, 0x5515,
  0x5540, 0x5541, 0x5544, 0x5545, 0x5550, 0x5551, 0x5554, 0x5555,
};

// Morton table for interleaving bytes, shifted left 1 bit
static const word_t morton_table1[256] =
{
  0x0000, 0x0002, 0x0008, 0x000A, 0x0020, 0x0022, 0x0028, 0x002A,
  0x0080, 0x0082, 0x0088, 0x008A, 0x00A0, 0x00A2, 0x00A8, 0x00AA,
  0x0200, 0x0202, 0x0208, 0x020A, 0x0220, 0x0222, 0x0228, 0x022A,
  0x0280, 0x0282, 0x0288, 0x028A, 0x02A0, 0x02A2, 0x02A8, 0x02AA,
  0x0800, 0x0802, 0x0808, 0x080A, 0x0820, 0x0822, 0x0828, 0x082A,
  0x0880, 0x0882, 0x0888, 0x088A, 0x08A0, 0x08A2, 0x08A8, 0x08AA,
  0x0A00, 0x0A02, 0x0A08, 0x0A0A, 0x0A20, 0x0A22, 0x0A28, 0x0A2A,
  0x0A80, 0x0A82, 0x0A88, 0x0A8A, 0x0AA0, 0x0AA2, 0x0AA8, 0x0AAA,
  0x2000, 0x2002, 0x2008, 0x200A, 0x2020, 0x2022, 0x2028, 0x202A,
  0x2080, 0x2082, 0x2088, 0x208A, 0x20A0, 0x20A2, 0x20A8, 0x20AA,
  0x2200, 0x2202, 0x2208, 0x220A, 0x2220, 0x2222, 0x2228, 0x222A,
  0x2280, 0x2282, 0x2288, 0x228A, 0x22A0, 0x22A2, 0x22A8, 0x22AA,
  0x2800, 0x2802, 0x2808, 0x280A, 0x2820, 0x2822, 0x2828, 0x282A,
  0x2880, 0x2882, 0x2888, 0x288A, 0x28A0, 0x28A2, 0x28A8, 0x28AA,
  0x2A00, 0x2A02, 0x2A08, 0x2A0A, 0x2A20, 0x2A22, 0x2A28, 0x2A2A,
  0x2A80, 0x2A82, 0x2A88, 0x2A8A, 0x2AA0, 0x2AA2, 0x2AA8, 0x2AAA,
  0x8000, 0x8002, 0x8008, 0x800A, 0x8020, 0x8022, 0x8028, 0x802A,
  0x8080, 0x8082, 0x8088, 0x808A, 0x80A0, 0x80A2, 0x80A8, 0x80AA,
  0x8200, 0x8202, 0x8208, 0x820A, 0x8220, 0x8222, 0x8228, 0x822A,
  0x8280, 0x8282, 0x8288, 0x828A, 0x82A0, 0x82A2, 0x82A8, 0x82AA,
  0x8800, 0x8802, 0x8808, 0x880A, 0x8820, 0x8822, 0x8828, 0x882A,
  0x8880, 0x8882, 0x8888, 0x888A, 0x88A0, 0x88A2, 0x88A8, 0x88AA,
  0x8A00, 0x8A02, 0x8A08, 0x8A0A, 0x8A20, 0x8A22, 0x8A28, 0x8A2A,
  0x8A80, 0x8A82, 0x8A88, 0x8A8A, 0x8AA0, 0x8AA2, 0x8AA8, 0x8AAA,
  0xA000, 0xA002, 0xA008, 0xA00A, 0xA020, 0xA022, 0xA028, 0xA02A,
  0xA080, 0xA082, 0xA088, 0xA08A, 0xA0A0, 0xA0A2, 0xA0A8, 0xA0AA,
  0xA200, 0xA202, 0xA208, 0xA20A, 0xA220, 0xA222, 0xA228, 0xA22A,
  0xA280, 0xA282, 0xA288, 0xA28A, 0xA2A0, 0xA2A2, 0xA2A8, 0xA2AA,
  0xA800, 0xA802, 0xA808, 0xA80A, 0xA820, 0xA822, 0xA828, 0xA82A,
  0xA880, 0xA882, 0xA888, 0xA88A, 0xA8A0, 0xA8A2, 0xA8A8, 0xA8AA,
  0xAA00, 0xAA02, 0xAA08, 0xAA0A, 0xAA20, 0xAA22, 0xAA28, 0xAA2A,
  0xAA80, 0xAA82, 0xAA88, 0xAA8A, 0xAAA0, 0xAAA2, 0xAAA8, 0xAAAA,
};

//
// Macros
//

// WORD_SIZE is the number of bits per word
// sizeof gives size in bytes (8 bits per byte)
#define WORD_SIZE 64
// #define WORD_SIZE (sizeof(word_t) * 8)

// POPCOUNT is number of bits set

#if defined(_WIN32)

// See http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
static word_t __inline windows_popcount(word_t w)
{
  w = w - ((w >> 1) & (word_t)~(word_t)0/3);
  w = (w & (word_t)~(word_t)0/15*3) + ((w >> 2) & (word_t)~(word_t)0/15*3);
  w = (w + (w >> 4)) & (word_t)~(word_t)0/255*15;
  c = (word_t)(w * ((word_t)~(word_t)0/255)) >> (sizeof(word_t) - 1) * 8;
}

static word_t __inline windows_parity(word_t w)
{
  w ^= w >> 1;
  w ^= w >> 2;
  w = (w & 0x1111111111111111UL) * 0x1111111111111111UL;
  return (w >> 60) & 1;
}

#define POPCOUNT(x) windows_popcountl(x)
#define PARITY(x) windows_parity(x)
#else
#define POPCOUNT(x) (unsigned)__builtin_popcountll(x)
#define PARITY(x) (unsigned)__builtin_parityll(x)
#endif

#define MIN(a, b)  (((a) <= (b)) ? (a) : (b))
#define MAX(a, b)  (((a) >= (b)) ? (a) : (b))

// Make this a power of two
#define INIT_CAPACITY_WORDS 2

// word of all 1s
#define WORD_MAX  (~(word_t)0)

#define SET_REGION(arr,start,len)    _set_region((arr),(start),(len),FILL_REGION)
#define CLEAR_REGION(arr,start,len)  _set_region((arr),(start),(len),ZERO_REGION)
#define TOGGLE_REGION(arr,start,len) _set_region((arr),(start),(len),SWAP_REGION)

// Have we initialised with srand() ?
static char rand_initiated = 0;

static void _seed_rand()
{
  if(!rand_initiated)
  {
    // Initialise random number generator
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((((time.tv_sec ^ getpid()) * 1000001) + time.tv_usec));
    rand_initiated = 1;
  }
}

//
// Common internal functions
//

#define bits_in_top_word(nbits) ((nbits) ? bitset64_idx((nbits) - 1) + 1 : 0)

// Mostly used for debugging
static inline void _print_word(word_t word, FILE* out)
{
  word_offset_t i;
  for(i = 0; i < WORD_SIZE; i++)
  {
    fprintf(out, "%c", ((word >> i) & (word_t)0x1) == 0 ? '0' : '1');
  }
}

// prints right to left
static inline char* _word_to_str(word_t word, char str[WORD_SIZE+1])
  __attribute__((unused));

static inline char* _word_to_str(word_t word, char str[WORD_SIZE+1])
{
  word_offset_t i;
  for(i = 0; i < WORD_SIZE; i++)
  {
    str[WORD_SIZE-i-1] = ((word >> i) & (word_t)0x1) == 0 ? '0' : '1';
  }
  str[WORD_SIZE] = '\0';
  return str;
}

// Used in debugging
#ifdef DEBUG
  #define DEBUG_PRINT(msg,...) printf("[%s:%i] "msg, __FILE__, __LINE__, ##__VA_ARGS__);
  #define DEBUG_VALIDATE(a) validate_bitarr((a), __FILE__, __LINE__)
#else
  #define DEBUG_PRINT(msg,...)
  #define DEBUG_VALIDATE(a)
#endif

void validate_bitarr(BIT_ARRAY *arr, const char *file, int lineno)
{
  // Check top word is masked
  word_addr_t tw = arr->num_of_words == 0 ? 0 : arr->num_of_words - 1;
  bit_index_t top_bits = bits_in_top_word(arr->num_of_bits);
  int err = 0;

  if(arr->words[tw] > bitmask64(top_bits))
  {
    _print_word(arr->words[tw], stderr);
    fprintf(stderr, "\n[%s:%i] Expected %i bits in top word[%i]\n",
            file, lineno, (int)top_bits, (int)tw);
    err = 1;
  }

  // Check num of words is correct
  word_addr_t num_words = roundup_bits2words64(arr->num_of_bits);
  if(num_words != arr->num_of_words)
  {
    fprintf(stderr, "\n[%s:%i] num of words wrong "
                    "[bits: %i, word: %i, actual words: %i]\n", file, lineno,
            (int)arr->num_of_bits, (int)num_words, (int)arr->num_of_words);
    err = 1;
  }

  if(err) abort();
}

// Reverse a word
static inline word_t _reverse_word(word_t word)
{
  word_t reverse = (reverse_table[(word)       & 0xff] << 56) |
                   (reverse_table[(word >>  8) & 0xff] << 48) |
                   (reverse_table[(word >> 16) & 0xff] << 40) |
                   (reverse_table[(word >> 24) & 0xff] << 32) |
                   (reverse_table[(word >> 32) & 0xff] << 24) |
                   (reverse_table[(word >> 40) & 0xff] << 16) |
                   (reverse_table[(word >> 48) & 0xff] << 8) |
                   (reverse_table[(word >> 56) & 0xff]);

  return reverse;
}

static inline void _mask_top_word(BIT_ARRAY* bitarr)
{
  // Mask top word
  word_addr_t num_of_words = MAX(1, bitarr->num_of_words);
  word_offset_t bits_active = bits_in_top_word(bitarr->num_of_bits);
  bitarr->words[num_of_words-1] &= bitmask64(bits_active);
}

//
// Get and set words (internal use only -- no bounds checking)
//

static inline word_t _get_word(const BIT_ARRAY* bitarr, bit_index_t start)
{
  word_addr_t word_index = bitset64_wrd(start);
  word_offset_t word_offset = bitset64_idx(start);

  word_t result = bitarr->words[word_index] >> word_offset;

  word_offset_t bits_taken = WORD_SIZE - word_offset;

  // word_offset is now the number of bits we need from the next word
  // Check the next word has at least some bits
  if(word_offset > 0 && start + bits_taken < bitarr->num_of_bits)
  {
    result |= bitarr->words[word_index+1] << (WORD_SIZE - word_offset);
  }

  return result;
}

// Set 64 bits from a particular start position
// Doesn't extend bit array
static inline void _set_word(BIT_ARRAY* bitarr, bit_index_t start, word_t word)
{
  word_addr_t word_index = bitset64_wrd(start);
  word_offset_t word_offset = bitset64_idx(start);

  if(word_offset == 0)
  {
    bitarr->words[word_index] = word;
  }
  else
  {
    bitarr->words[word_index]
      = (word << word_offset) |
        (bitarr->words[word_index] & bitmask64(word_offset));

    if(word_index+1 < bitarr->num_of_words)
    {
      bitarr->words[word_index+1]
        = (word >> (WORD_SIZE - word_offset)) |
          (bitarr->words[word_index+1] & (WORD_MAX << word_offset));
    }
  }

  // Mask top word
  _mask_top_word(bitarr);
  DEBUG_VALIDATE(bitarr);
}

static inline void _set_byte(BIT_ARRAY *bitarr, bit_index_t start, uint8_t byte)
{
  word_t w = _get_word(bitarr, start);
  _set_word(bitarr, start, (w & ~(word_t)0xff) | byte);
}

// 4 bits
static inline void _set_nibble(BIT_ARRAY *bitarr, bit_index_t start,
                               uint8_t nibble)
{
  word_t w = _get_word(bitarr, start);
  _set_word(bitarr, start, (w & ~(word_t)0xf) | nibble);
}

// Wrap around
static inline word_t _get_word_cyclic(const BIT_ARRAY* bitarr, bit_index_t start)
{
  word_t word = _get_word(bitarr, start);

  bit_index_t bits_taken = bitarr->num_of_bits - start;

  if(bits_taken < WORD_SIZE)
  {
    word |= (bitarr->words[0] << bits_taken);

    if(bitarr->num_of_bits < (bit_index_t)WORD_SIZE)
    {
      // Mask word to prevent repetition of the same bits
      word = word & bitmask64(bitarr->num_of_bits);
    }
  }

  return word;
}

// Wrap around
static inline void _set_word_cyclic(BIT_ARRAY* bitarr,
                                    bit_index_t start, word_t word)
{
  _set_word(bitarr, start, word);

  bit_index_t bits_set = bitarr->num_of_bits - start;

  if(bits_set < WORD_SIZE && start > 0)
  {
    word >>= bits_set;

    // Prevent overwriting the bits we've just set
    // by setting 'start' as the upper bound for the number of bits to write
    word_offset_t bits_remaining = MIN(WORD_SIZE - bits_set, start);
    word_t mask = bitmask64(bits_remaining);

    bitarr->words[0] = bitmask_merge(word, bitarr->words[0], mask);
  }
}

//
// Fill a region (internal use only)
//

// FillAction is fill with 0 or 1 or toggle
typedef enum {ZERO_REGION, FILL_REGION, SWAP_REGION} FillAction;

static inline void _set_region(BIT_ARRAY* bitarr, bit_index_t start,
                               bit_index_t length, FillAction action)
{
  if(length == 0) return;

  word_addr_t first_word = bitset64_wrd(start);
  word_addr_t last_word = bitset64_wrd(start+length-1);
  word_offset_t foffset = bitset64_idx(start);
  word_offset_t loffset = bitset64_idx(start+length-1);

  if(first_word == last_word)
  {
    word_t mask = bitmask64(length) << foffset;

    switch(action)
    {
      case ZERO_REGION: bitarr->words[first_word] &= ~mask; break;
      case FILL_REGION: bitarr->words[first_word] |=  mask; break;
      case SWAP_REGION: bitarr->words[first_word] ^=  mask; break;
    }
  }
  else
  {
    // Set first word
    switch(action)
    {
      case ZERO_REGION: bitarr->words[first_word] &=  bitmask64(foffset); break;
      case FILL_REGION: bitarr->words[first_word] |= ~bitmask64(foffset); break;
      case SWAP_REGION: bitarr->words[first_word] ^= ~bitmask64(foffset); break;
    }

    word_addr_t i;

    // Set whole words
    switch(action)
    {
      case ZERO_REGION:
        for(i = first_word + 1; i < last_word; i++)
          bitarr->words[i] = (word_t)0;
        break;
      case FILL_REGION:
        for(i = first_word + 1; i < last_word; i++)
          bitarr->words[i] = WORD_MAX;
        break;
      case SWAP_REGION:
        for(i = first_word + 1; i < last_word; i++)
          bitarr->words[i] ^= WORD_MAX;
        break;
    }

    // Set last word
    switch(action)
    {
      case ZERO_REGION: bitarr->words[last_word] &= ~bitmask64(loffset+1); break;
      case FILL_REGION: bitarr->words[last_word] |=  bitmask64(loffset+1); break;
      case SWAP_REGION: bitarr->words[last_word] ^=  bitmask64(loffset+1); break;
    }
  }
}



//
// Constructor
//

// If cannot allocate memory, set errno to ENOMEM, return NULL
BIT_ARRAY* bit_array_alloc(BIT_ARRAY* bitarr, bit_index_t nbits)
{
  bitarr->num_of_bits = nbits;
  bitarr->num_of_words = roundup_bits2words64(nbits);
  bitarr->capacity_in_words = MAX(8, roundup2pow(bitarr->num_of_words));
  bitarr->words = (word_t*)calloc(bitarr->capacity_in_words, sizeof(word_t));
  if(bitarr->words == NULL) {
    errno = ENOMEM;
    return NULL;
  }
  return bitarr;
}

void bit_array_dealloc(BIT_ARRAY* bitarr)
{
  free(bitarr->words);
  memset(bitarr, 0, sizeof(BIT_ARRAY));
}

// If cannot allocate memory, set errno to ENOMEM, return NULL
BIT_ARRAY* bit_array_create(bit_index_t nbits)
{
  BIT_ARRAY* bitarr = (BIT_ARRAY*)malloc(sizeof(BIT_ARRAY));

  // error if could not allocate enough memory
  if(bitarr == NULL || bit_array_alloc(bitarr, nbits) == NULL)
  {
    if(bitarr != NULL) free(bitarr);
    errno = ENOMEM;
    return NULL;
  }

  DEBUG_PRINT("Creating BIT_ARRAY (bits: %lu; allocated words: %lu; "
              "using words: %lu; WORD_SIZE: %i)\n",
              (unsigned long)nbits, (unsigned long)bitarr->capacity_in_words,
              (unsigned long)roundup_bits2words64(nbits), (int)WORD_SIZE);

  DEBUG_VALIDATE(bitarr);

  return bitarr;
}

//
// Destructor
//
void bit_array_free(BIT_ARRAY* bitarr)
{
  if(bitarr->words != NULL)
    free(bitarr->words);

  free(bitarr);
}

bit_index_t bit_array_length(const BIT_ARRAY* bit_arr)
{
  return bit_arr->num_of_bits;
}

// Change the size of a bit array. Enlarging an array will add zeros
// to the end of it. Returns 1 on success, 0 on failure (e.g. not enough memory)
char bit_array_resize(BIT_ARRAY* bitarr, bit_index_t new_num_of_bits)
{
  word_addr_t old_num_of_words = bitarr->num_of_words;
  word_addr_t new_num_of_words = roundup_bits2words64(new_num_of_bits);

  bitarr->num_of_bits = new_num_of_bits;
  bitarr->num_of_words = new_num_of_words;

  DEBUG_PRINT("Resize: old_num_of_words: %i; new_num_of_words: %i capacity: %i\n",
              (int)old_num_of_words, (int)new_num_of_words,
              (int)bitarr->capacity_in_words);

  if(new_num_of_words > bitarr->capacity_in_words)
  {
    // Need to change the amount of memory used
    word_addr_t old_capacity_in_words = bitarr->capacity_in_words;
    size_t old_capacity_in_bytes = old_capacity_in_words * sizeof(word_t);

    bitarr->capacity_in_words = roundup2pow(new_num_of_words);
    bitarr->capacity_in_words = MAX(8, bitarr->capacity_in_words);

    size_t new_capacity_in_bytes = bitarr->capacity_in_words * sizeof(word_t);
    bitarr->words = (word_t*)realloc(bitarr->words, new_capacity_in_bytes);

    if(bitarr->words == NULL)
    {
      // error - could not allocate enough memory
      perror("resize realloc");
      errno = ENOMEM;
      return 0;
    }

    // Need to zero new memory
    size_t num_bytes_to_zero = new_capacity_in_bytes - old_capacity_in_bytes;
    memset(bitarr->words + old_capacity_in_words, 0, num_bytes_to_zero);

    DEBUG_PRINT("zeroing from word %i for %i bytes\n", (int)old_capacity_in_words,
                (int)num_bytes_to_zero);
  }
  else if(new_num_of_words < old_num_of_words)
  {
    // Shrunk -- need to zero old memory
    size_t num_bytes_to_zero = (old_num_of_words - new_num_of_words)*sizeof(word_t);

    memset(bitarr->words + new_num_of_words, 0, num_bytes_to_zero);
  }

  // Mask top word
  _mask_top_word(bitarr);
  DEBUG_VALIDATE(bitarr);
  return 1;
}

void bit_array_resize_critical(BIT_ARRAY* bitarr, bit_index_t num_of_bits)
{
  bit_index_t old_num_of_bits = bitarr->num_of_bits;

  if(!bit_array_resize(bitarr, num_of_bits))
  {
    fprintf(stderr, "Ran out of memory resizing [%lu -> %lu]",
            (unsigned long)old_num_of_bits, (unsigned long)num_of_bits);
    abort();
  }
}

// If bitarr length < num_bits, resizes to num_bits
char bit_array_ensure_size(BIT_ARRAY* bitarr, bit_index_t ensure_num_of_bits)
{
  if(bitarr->num_of_bits < ensure_num_of_bits)
  {
    return bit_array_resize(bitarr, ensure_num_of_bits);
  }

  return 1;
}

void bit_array_ensure_size_critical(BIT_ARRAY* bitarr, bit_index_t num_of_bits)
{
  if(num_of_bits > bitarr->num_of_bits)
  {
    bit_array_resize_critical(bitarr, num_of_bits);
  }
}

static inline
void _bit_array_ensure_nwords(BIT_ARRAY* bitarr, word_addr_t nwords,
                              const char *file, int lineno, const char *func)
{
  size_t newmem, oldmem;
  if(bitarr->capacity_in_words < nwords) {
    oldmem = bitarr->capacity_in_words * sizeof(word_t);
    bitarr->capacity_in_words = roundup2pow(nwords);
    newmem = bitarr->capacity_in_words * sizeof(word_t);
    bitarr->words = (word_t*)realloc(bitarr->words, newmem);

    if(bitarr->words == NULL) {
      fprintf(stderr, "[%s:%i:%s()] Ran out of memory resizing [%zu -> %zu]",
              file, lineno, func, oldmem, newmem);
      abort();
    }

    DEBUG_PRINT("Ensure nwords realloc %zu -> %zu\n", oldmem, newmem);
  }
}


//
// Get, set, clear, assign and toggle individual bits
//

// Get the value of a bit (returns 0 or 1)
char bit_array_get_bit(const BIT_ARRAY* bitarr, bit_index_t b)
{
  assert(b < bitarr->num_of_bits);
  return bit_array_get(bitarr, b);
}

// set a bit (to 1) at position b
void bit_array_set_bit(BIT_ARRAY* bitarr, bit_index_t b)
{
  assert(b < bitarr->num_of_bits);
  bit_array_set(bitarr,b);
  DEBUG_VALIDATE(bitarr);
}

// clear a bit (to 0) at position b
void bit_array_clear_bit(BIT_ARRAY* bitarr, bit_index_t b)
{
  assert(b < bitarr->num_of_bits);
  bit_array_clear(bitarr, b);
  DEBUG_VALIDATE(bitarr);
}

// If bit is 0 -> 1, if bit is 1 -> 0.  AKA 'flip'
void bit_array_toggle_bit(BIT_ARRAY* bitarr, bit_index_t b)
{
  assert(b < bitarr->num_of_bits);
  bit_array_toggle(bitarr, b);
  DEBUG_VALIDATE(bitarr);
}

// If char c != 0, set bit; otherwise clear bit
void bit_array_assign_bit(BIT_ARRAY* bitarr, bit_index_t b, char c)
{
  assert(b < bitarr->num_of_bits);
  bit_array_assign(bitarr, b, c ? 1 : 0);
  DEBUG_VALIDATE(bitarr);
}

//
// Get, set etc with resize
//

// Get the value of a bit (returns 0 or 1)
char bit_array_rget(BIT_ARRAY* bitarr, bit_index_t b)
{
  bit_array_ensure_size_critical(bitarr, b+1);
  return bit_array_get(bitarr, b);
}

// set a bit (to 1) at position b
void bit_array_rset(BIT_ARRAY* bitarr, bit_index_t b)
{
  bit_array_ensure_size_critical(bitarr, b+1);
  bit_array_set(bitarr,b);
  DEBUG_VALIDATE(bitarr);
}

// clear a bit (to 0) at position b
void bit_array_rclear(BIT_ARRAY* bitarr, bit_index_t b)
{
  bit_array_ensure_size_critical(bitarr, b+1);
  bit_array_clear(bitarr, b);
  DEBUG_VALIDATE(bitarr);
}

// If bit is 0 -> 1, if bit is 1 -> 0.  AKA 'flip'
void bit_array_rtoggle(BIT_ARRAY* bitarr, bit_index_t b)
{
  bit_array_ensure_size_critical(bitarr, b+1);
  bit_array_toggle(bitarr, b);
  DEBUG_VALIDATE(bitarr);
}

// If char c != 0, set bit; otherwise clear bit
void bit_array_rassign(BIT_ARRAY* bitarr, bit_index_t b, char c)
{
  bit_array_ensure_size_critical(bitarr, b+1);
  bit_array_assign(bitarr, b, c ? 1 : 0);
  DEBUG_VALIDATE(bitarr);
}

//
// Get, set, clear and toggle several bits at once
//

// Get the offsets of the set bits (for offsets start<=offset<end)
// Returns the number of bits set
// It is assumed that dst is at least of length (end-start)
bit_index_t bit_array_get_bits(const BIT_ARRAY* bitarr,
                               bit_index_t start, bit_index_t end,
                               bit_index_t* dst)
{
  bit_index_t i, n = 0;
  assert(end <= bitarr->num_of_bits);
  for(i = start; i < end; i++) {
    if(bit_array_get(bitarr, i)) {
      dst[n++] = i;
    }
  }
  return n;
}

// Set multiple bits at once.
// e.g. set bits 1, 20 & 31: bit_array_set_bits(bitarr, 3, 1,20,31);
void bit_array_set_bits(BIT_ARRAY* bitarr, size_t n, ...)
{
  size_t i;
  va_list argptr;
  va_start(argptr, n);

  for(i = 0; i < n; i++)
  {
    unsigned int bit_index = va_arg(argptr, unsigned int);
    bit_array_set_bit(bitarr, bit_index);
  }

  va_end(argptr);
  DEBUG_VALIDATE(bitarr);
}

// Clear multiple bits at once.
// e.g. clear bits 1, 20 & 31: bit_array_clear_bits(bitarr, 3, 1,20,31);
void bit_array_clear_bits(BIT_ARRAY* bitarr, size_t n, ...)
{
  size_t i;
  va_list argptr;
  va_start(argptr, n);

  for(i = 0; i < n; i++)
  {
    unsigned int bit_index = va_arg(argptr, unsigned int);
    bit_array_clear_bit(bitarr, bit_index);
  }

  va_end(argptr);
  DEBUG_VALIDATE(bitarr);
}

// Toggle multiple bits at once
// e.g. toggle bits 1, 20 & 31: bit_array_toggle_bits(bitarr, 3, 1,20,31);
void bit_array_toggle_bits(BIT_ARRAY* bitarr, size_t n, ...)
{
  size_t i;
  va_list argptr;
  va_start(argptr, n);

  for(i = 0; i < n; i++)
  {
    unsigned int bit_index = va_arg(argptr, unsigned int);
    bit_array_toggle_bit(bitarr, bit_index);
  }

  va_end(argptr);
  DEBUG_VALIDATE(bitarr);
}


//
// Set, clear and toggle all bits in a region
//

// Set all the bits in a region
void bit_array_set_region(BIT_ARRAY* bitarr, bit_index_t start, bit_index_t len)
{
  assert(start + len <= bitarr->num_of_bits);
  SET_REGION(bitarr, start, len);
  DEBUG_VALIDATE(bitarr);
}


// Clear all the bits in a region
void bit_array_clear_region(BIT_ARRAY* bitarr, bit_index_t start, bit_index_t len)
{
  assert(start + len <= bitarr->num_of_bits);
  CLEAR_REGION(bitarr, start, len);
  DEBUG_VALIDATE(bitarr);
}

// Toggle all the bits in a region
void bit_array_toggle_region(BIT_ARRAY* bitarr, bit_index_t start, bit_index_t len)
{
  assert(start + len <= bitarr->num_of_bits);
  TOGGLE_REGION(bitarr, start, len);
  DEBUG_VALIDATE(bitarr);
}


//
// Set, clear and toggle all bits at once
//

// set all elements of data to one
void bit_array_set_all(BIT_ARRAY* bitarr)
{
  bit_index_t num_of_bytes = bitarr->num_of_words * sizeof(word_t);
  memset(bitarr->words, 0xFF, num_of_bytes);
  _mask_top_word(bitarr);
  DEBUG_VALIDATE(bitarr);
}

// set all elements of data to zero
void bit_array_clear_all(BIT_ARRAY* bitarr)
{
  memset(bitarr->words, 0, bitarr->num_of_words * sizeof(word_t));
  DEBUG_VALIDATE(bitarr);
}

// Set all 1 bits to 0, and all 0 bits to 1. AKA flip
void bit_array_toggle_all(BIT_ARRAY* bitarr)
{
  word_addr_t i;
  for(i = 0; i < bitarr->num_of_words; i++)
  {
    bitarr->words[i] ^= WORD_MAX;
  }

  _mask_top_word(bitarr);
  DEBUG_VALIDATE(bitarr);
}

//
// Get a word at a time
//

uint64_t bit_array_get_word64(const BIT_ARRAY* bitarr, bit_index_t start)
{
  assert(start < bitarr->num_of_bits);
  return (uint64_t)_get_word(bitarr, start);
}

uint32_t bit_array_get_word32(const BIT_ARRAY* bitarr, bit_index_t start)
{
  assert(start < bitarr->num_of_bits);
  return (uint32_t)_get_word(bitarr, start);
}

uint16_t bit_array_get_word16(const BIT_ARRAY* bitarr, bit_index_t start)
{
  assert(start < bitarr->num_of_bits);
  return (uint16_t)_get_word(bitarr, start);
}

uint8_t bit_array_get_word8(const BIT_ARRAY* bitarr, bit_index_t start)
{
  assert(start < bitarr->num_of_bits);
  return (uint8_t)_get_word(bitarr, start);
}

uint64_t bit_array_get_wordn(const BIT_ARRAY* bitarr, bit_index_t start, int n)
{
  assert(start < bitarr->num_of_bits);
  assert(n <= 64);
  return (uint64_t)(_get_word(bitarr, start) & bitmask64(n));
}

//
// Set a word at a time
//
// Doesn't extend bit array. However it is safe to TRY to set bits beyond the
// end of the array, as long as: `start` is < `bit_array_length(arr)`
//

void bit_array_set_word64(BIT_ARRAY* bitarr, bit_index_t start, uint64_t word)
{
  assert(start < bitarr->num_of_bits);
  _set_word(bitarr, start, (word_t)word);
}

void bit_array_set_word32(BIT_ARRAY* bitarr, bit_index_t start, uint32_t word)
{
  assert(start < bitarr->num_of_bits);
  word_t w = _get_word(bitarr, start);
  _set_word(bitarr, start, bitmask_merge(w, word, 0xffffffff00000000UL));
}

void bit_array_set_word16(BIT_ARRAY* bitarr, bit_index_t start, uint16_t word)
{
  assert(start < bitarr->num_of_bits);
  word_t w = _get_word(bitarr, start);
  _set_word(bitarr, start, bitmask_merge(w, word, 0xffffffffffff0000UL));
}

void bit_array_set_word8(BIT_ARRAY* bitarr, bit_index_t start, uint8_t byte)
{
  assert(start < bitarr->num_of_bits);
  _set_byte(bitarr, start, byte);
}

void bit_array_set_wordn(BIT_ARRAY* bitarr, bit_index_t start, uint64_t word, int n)
{
  assert(start < bitarr->num_of_bits);
  assert(n <= 64);
  word_t w = _get_word(bitarr, start), m = bitmask64(n);
  _set_word(bitarr, start, bitmask_merge(word,w,m));
}

//
// Number of bits set
//

// Get the number of bits set (hamming weight)
bit_index_t bit_array_num_bits_set(const BIT_ARRAY* bitarr)
{
  word_addr_t i;

  bit_index_t num_of_bits_set = 0;

  for(i = 0; i < bitarr->num_of_words; i++)
  {
    if(bitarr->words[i] > 0)
    {
      num_of_bits_set += POPCOUNT(bitarr->words[i]);
    }
  }

  return num_of_bits_set;
}

// Get the number of bits not set (1 - hamming weight)
bit_index_t bit_array_num_bits_cleared(const BIT_ARRAY* bitarr)
{
  return bitarr->num_of_bits - bit_array_num_bits_set(bitarr);
}


// Get the number of bits set in on array and not the other.  This is equivalent
// to hamming weight of the XOR when the two arrays are the same length.
// e.g. 10101 vs 00111 => hamming distance 2 (XOR is 10010)
bit_index_t bit_array_hamming_distance(const BIT_ARRAY* arr1,
                                       const BIT_ARRAY* arr2)
{
  word_addr_t min_words = MIN(arr1->num_of_words, arr2->num_of_words);
  word_addr_t max_words = MAX(arr1->num_of_words, arr2->num_of_words);

  bit_index_t hamming_distance = 0;
  word_addr_t i;

  for(i = 0; i < min_words; i++)
  {
    hamming_distance += POPCOUNT(arr1->words[i] ^ arr2->words[i]);
  }

  if(min_words != max_words)
  {
    const BIT_ARRAY* long_arr
      = (arr1->num_of_words > arr2->num_of_words ? arr1 : arr2);

    for(i = min_words; i < max_words; i++)
    {
      hamming_distance += POPCOUNT(long_arr->words[i]);
    }
  }

  return hamming_distance;
}

// Parity - returns 1 if odd number of bits set, 0 if even
char bit_array_parity(const BIT_ARRAY* bitarr)
{
  word_addr_t w;
  unsigned int parity = 0;

  for(w = 0; w < bitarr->num_of_words; w++)
  {
    parity ^= PARITY(bitarr->words[w]);
  }

  return (char)parity;
}

//
// Find indices of set/clear bits
//

// Find the index of the next bit that is set/clear, at or after `offset`
// Returns 1 if such a bit is found, otherwise 0
// Index is stored in the integer pointed to by `result`
// If no such bit is found, value at `result` is not changed
#define _next_bit_func_def(FUNC,GET) \
char FUNC(const BIT_ARRAY* bitarr, bit_index_t offset, bit_index_t* result) \
{ \
  assert(offset < bitarr->num_of_bits); \
  if(bitarr->num_of_bits == 0 || offset >= bitarr->num_of_bits) { return 0; } \
 \
  /* Find first word that is greater than zero */ \
  word_addr_t i = bitset64_wrd(offset); \
  word_t w = GET(bitarr->words[i]) & ~bitmask64(bitset64_idx(offset)); \
 \
  while(1) { \
    if(w > 0) { \
      bit_index_t pos = i * WORD_SIZE + trailing_zeros(w); \
      if(pos < bitarr->num_of_bits) { *result = pos; return 1; } \
      else { return 0; } \
    } \
    i++; \
    if(i >= bitarr->num_of_words) break; \
    w = GET(bitarr->words[i]); \
  } \
 \
  return 0; \
}

// Find the index of the previous bit that is set/clear, before `offset`.
// Returns 1 if such a bit is found, otherwise 0
// Index is stored in the integer pointed to by `result`
// If no such bit is found, value at `result` is not changed
#define _prev_bit_func_def(FUNC,GET) \
char FUNC(const BIT_ARRAY* bitarr, bit_index_t offset, bit_index_t* result) \
{ \
  assert(offset <= bitarr->num_of_bits); \
  if(bitarr->num_of_bits == 0 || offset == 0) { return 0; } \
 \
  /* Find prev word that is greater than zero */ \
  word_addr_t i = bitset64_wrd(offset-1); \
  word_t w = GET(bitarr->words[i]) & bitmask64(bitset64_idx(offset-1)+1); \
 \
  if(w > 0) { *result = (i+1) * WORD_SIZE - leading_zeros(w) - 1; return 1; } \
 \
  /* i is unsigned so have to use break when i == 0 */ \
  for(--i; i != BIT_INDEX_MAX; i--) { \
    w = GET(bitarr->words[i]); \
    if(w > 0) { \
      *result = (i+1) * WORD_SIZE - leading_zeros(w) - 1; \
      return 1; \
    } \
  } \
 \
  return 0; \
}

#define GET_WORD(x) (x)
#define NEG_WORD(x) (~(x))
_next_bit_func_def(bit_array_find_next_set_bit,  GET_WORD);
_next_bit_func_def(bit_array_find_next_clear_bit,NEG_WORD);
_prev_bit_func_def(bit_array_find_prev_set_bit,  GET_WORD);
_prev_bit_func_def(bit_array_find_prev_clear_bit,NEG_WORD);

// Find the index of the first bit that is set.
// Returns 1 if a bit is set, otherwise 0
// Index of first set bit is stored in the integer pointed to by result
// If no bits are set, value at `result` is not changed
char bit_array_find_first_set_bit(const BIT_ARRAY* bitarr, bit_index_t* result)
{
  return bit_array_find_next_set_bit(bitarr, 0, result);
}

// same same
char bit_array_find_first_clear_bit(const BIT_ARRAY* bitarr, bit_index_t* result)
{
  return bit_array_find_next_clear_bit(bitarr, 0, result);
}

// Find the index of the last bit that is set.
// Returns 1 if a bit is set, otherwise 0
// Index of last set bit is stored in the integer pointed to by `result`
// If no bits are set, value at `result` is not changed
char bit_array_find_last_set_bit(const BIT_ARRAY* bitarr, bit_index_t* result)
{
  return bit_array_find_prev_set_bit(bitarr, bitarr->num_of_bits, result);
}

// same same
char bit_array_find_last_clear_bit(const BIT_ARRAY* bitarr, bit_index_t* result)
{
  return bit_array_find_prev_clear_bit(bitarr, bitarr->num_of_bits, result);
}

//
// "Sorting" bits
//

// Put all the 0s before all the 1s
void bit_array_sort_bits(BIT_ARRAY* bitarr)
{
  bit_index_t num_of_bits_set = bit_array_num_bits_set(bitarr);
  bit_index_t num_of_bits_cleared = bitarr->num_of_bits - num_of_bits_set;
  bit_array_set_all(bitarr);
  CLEAR_REGION(bitarr, 0, num_of_bits_cleared);
  DEBUG_VALIDATE(bitarr);
}

// Put all the 1s before all the 0s
void bit_array_sort_bits_rev(BIT_ARRAY* bitarr)
{
  bit_index_t num_of_bits_set = bit_array_num_bits_set(bitarr);
  bit_array_clear_all(bitarr);
  SET_REGION(bitarr, 0, num_of_bits_set);
  DEBUG_VALIDATE(bitarr);
}


//
// Strings and printing
//

// Construct a BIT_ARRAY from a substring with given on and off characters.
void bit_array_from_substr(BIT_ARRAY* bitarr, bit_index_t offset,
                           const char *str, size_t len,
                           const char *on, const char *off,
                           char left_to_right)
{
  bit_array_ensure_size(bitarr, offset + len);
  bit_array_clear_region(bitarr, offset, len);

  // BitArray region is now all 0s -- just set the 1s
  size_t i;
  bit_index_t j;

  for(i = 0; i < len; i++)
  {
    if(strchr(on, str[i]) != NULL)
    {
      j = offset + (left_to_right ? i : len - i - 1);
      bit_array_set(bitarr, j);
    }
    else { assert(strchr(off, str[i]) != NULL); }
  }

  DEBUG_VALIDATE(bitarr);
}

// From string method
void bit_array_from_str(BIT_ARRAY* bitarr, const char* str)
{
  bit_array_from_substr(bitarr, 0, str, strlen(str), "1", "0", 1);
}

// Takes a char array to write to.  `str` must be bitarr->num_of_bits+1 in length
// Terminates string with '\0'
char* bit_array_to_str(const BIT_ARRAY* bitarr, char* str)
{
  bit_index_t i;

  for(i = 0; i < bitarr->num_of_bits; i++)
  {
    str[i] = bit_array_get(bitarr, i) ? '1' : '0';
  }

  str[bitarr->num_of_bits] = '\0';

  return str;
}

char* bit_array_to_str_rev(const BIT_ARRAY* bitarr, char* str)
{
  bit_index_t i;

  for(i = 0; i < bitarr->num_of_bits; i++)
  {
    str[i] = bit_array_get(bitarr, bitarr->num_of_bits-i-1) ? '1' : '0';
  }

  str[bitarr->num_of_bits] = '\0';

  return str;
}


// Get a string representations for a given region, using given on/off characters.
// Note: does not null-terminate
void bit_array_to_substr(const BIT_ARRAY* bitarr,
                         bit_index_t start, bit_index_t length,
                         char* str, char on, char off,
                         char left_to_right)
{
  assert(start + length <= bitarr->num_of_bits);

  bit_index_t i, j;
  bit_index_t end = start + length - 1;

  for(i = 0; i < length; i++)
  {
    j = (left_to_right ? start + i : end - i);
    str[i] = bit_array_get(bitarr, j) ? on : off;
  }

//  str[length] = '\0';
}

// Print this array to a file stream.  Prints '0's and '1'.  Doesn't print newline.
void bit_array_print(const BIT_ARRAY* bitarr, FILE* fout)
{
  bit_index_t i;

  for(i = 0; i < bitarr->num_of_bits; i++)
  {
    fprintf(fout, "%c", bit_array_get(bitarr, i) ? '1' : '0');
  }
}

// Print a string representations for a given region, using given on/off characters.
void bit_array_print_substr(const BIT_ARRAY* bitarr,
                            bit_index_t start, bit_index_t length,
                            FILE* fout, char on, char off,
                            char left_to_right)
{
  assert(start + length <= bitarr->num_of_bits);

  bit_index_t i, j;
  bit_index_t end = start + length - 1;

  for(i = 0; i < length; i++)
  {
    j = (left_to_right ? start + i : end - i);
    fprintf(fout, "%c", bit_array_get(bitarr, j) ? on : off);
  }
}

//
// Decimal
//

// Get bit array as decimal str (e.g. 0b1101 -> "13")
// len is the length of str char array -- will write at most len-1 chars
// returns the number of characters needed
// return is the same as strlen(str)
size_t bit_array_to_decimal(const BIT_ARRAY *bitarr, char *str, size_t len)
{
  size_t i = 0;

  if(bit_array_cmp_uint64(bitarr, 0) == 0)
  {
    if(len >= 2)
    {
      *str = '0';
      *(str+1) = '\0';
    }

    return 1;
  }

  BIT_ARRAY *tmp = bit_array_clone(bitarr);
  uint64_t rem;

  str[len-1] = '\0';

  while(bit_array_cmp_uint64(tmp, 0) != 0)
  {
    bit_array_div_uint64(tmp, 10, &rem);

    if(i < len-1)
    {
      str[len-2-i] = '0' + rem;
    }

    i++;
  }

  if(i < len-1)
  {
    // Moves null-terminator as well
    memmove(str, str+len-i-1, i+1);
  }

  bit_array_free(tmp);

  return i;
}

// Get bit array from decimal str (e.g. "13" -> 0b1101)
// Returns number of characters used
size_t bit_array_from_decimal(BIT_ARRAY *bitarr, const char* decimal)
{
  bit_array_clear_all(bitarr);
  size_t i = 0;

  if(decimal[0] == '\0' || decimal[0] < '0' || decimal[0] > '9')
  {
    return 0;
  }

  bit_array_add_uint64(bitarr, decimal[i] - '0');
  i++;

  while(decimal[i] != '\0' && decimal[i] >= '0' && decimal[i] <= '9')
  {
    bit_array_mul_uint64(bitarr, 10);
    bit_array_add_uint64(bitarr, decimal[i] - '0');
    i++;
  }

  return i;
}

//
// Hexidecimal
//

char bit_array_hex_to_nibble(char c, uint8_t *b)
{
  c = tolower(c);

  if(c >= '0' && c <= '9')
  {
    *b = c - '0';
    return 1;
  }
  else if(c >= 'a' && c <= 'f')
  {
    *b = 0xa + (c - 'a');
    return 1;
  }
  else
  {
    return 0;
  }
}

char bit_array_nibble_to_hex(uint8_t b, char uppercase)
{
  if(b <= 9)
  {
    return '0' + b;
  }
  else
  {
    return (uppercase ? 'A' : 'a') + (b - 0xa);
  }
}

// Loads array from hex string
// Returns the number of bits loaded (will be chars rounded up to multiple of 4)
// (0 on failure)
bit_index_t bit_array_from_hex(BIT_ARRAY* bitarr, bit_index_t offset,
                               const char* str, size_t len)
{
  if(str[0] == '0' && tolower(str[1]) == 'x')
  {
    str += 2;
    len -= 2;
  }

  size_t i;
  for(i = 0; i < len; i++, offset += 4)
  {
    uint8_t b;
    if(bit_array_hex_to_nibble(str[i], &b))
    {
      bit_array_ensure_size(bitarr, offset + 4);
      _set_nibble(bitarr, offset, b);
    }
    else
    {
      break;
    }
  }

  return 4 * i;
}

// Returns number of characters written
size_t bit_array_to_hex(const BIT_ARRAY* bitarr,
                        bit_index_t start, bit_index_t length,
                        char* str, char uppercase)
{
  assert(start + length <= bitarr->num_of_bits);

  size_t k = 0;
  bit_index_t offset, end = start + length;

  for(offset = start; offset + WORD_SIZE <= end; offset += WORD_SIZE)
  {
    word_t w = _get_word(bitarr, offset);

    word_offset_t j;
    for(j = 0; j < 64; j += 4)
    {
      str[k++] = bit_array_nibble_to_hex((w>>j) & 0xf, uppercase);
    }
  }

  if(offset < end)
  {
    // Remaining full nibbles (4 bits)
    word_t w = _get_word(bitarr, offset);

    for(; offset + 4 <= end; offset += 4)
    {
      str[k++] = bit_array_nibble_to_hex(w & 0xf, uppercase);
      w >>= 4;
    }

    if(offset < end)
    {
      // Remaining bits
      str[k++] = bit_array_nibble_to_hex(w & bitmask64(end - offset), uppercase);
    }
  }

  str[k] = '\0';

  // Return number of characters written
  return k;
}

// Print bit array as hex
size_t bit_array_print_hex(const BIT_ARRAY* bitarr,
                           bit_index_t start, bit_index_t length,
                           FILE* fout, char uppercase)
{
  assert(start + length <= bitarr->num_of_bits);

  size_t k = 0;
  bit_index_t offset, end = start + length;

  for(offset = start; offset + WORD_SIZE <= end; offset += WORD_SIZE)
  {
    word_t w = _get_word(bitarr, offset);

    word_offset_t j;
    for(j = 0; j < 64; j += 4)
    {
      fprintf(fout, "%c", bit_array_nibble_to_hex((w>>j) & 0xf, uppercase));
      k++;
    }
  }

  if(offset < end)
  {
    // Remaining full nibbles (4 bits)
    word_t w = _get_word(bitarr, offset);

    for(; offset + 4 <= end; offset += 4)
    {
      fprintf(fout, "%c", bit_array_nibble_to_hex(w & 0xf, uppercase));
      w >>= 4;
      k++;
    }

    if(offset < end)
    {
      // Remaining bits
      char hex = bit_array_nibble_to_hex(w & bitmask64(end - offset), uppercase);
      fprintf(fout, "%c", hex);
      k++;
    }
  }

  return k;
}

//
// Clone and copy
//

// Returns NULL if cannot malloc
BIT_ARRAY* bit_array_clone(const BIT_ARRAY* bitarr)
{
  BIT_ARRAY* cpy = bit_array_create(bitarr->num_of_bits);

  if(cpy == NULL)
  {
    return NULL;
  }

  // Copy across bits
  memcpy(cpy->words, bitarr->words, bitarr->num_of_words * sizeof(word_t));

  DEBUG_VALIDATE(cpy);
  return cpy;
}

// destination and source may be the same bit_array
// and src/dst regions may overlap
static void _array_copy(BIT_ARRAY* dst, bit_index_t dstindx,
                        const BIT_ARRAY* src, bit_index_t srcindx,
                        bit_index_t length)
{
  DEBUG_PRINT("bit_array_copy(dst: %zu, src: %zu, length: %zu)\n",
              (size_t)dstindx, (size_t)srcindx, (size_t)length);

  // Num of full words to copy
  word_addr_t num_of_full_words = length / WORD_SIZE;
  word_addr_t i;

  word_offset_t bits_in_last_word = bits_in_top_word(length);

  if(dst == src && srcindx > dstindx)
  {
    // Work left to right
    DEBUG_PRINT("work left to right\n");

    for(i = 0; i < num_of_full_words; i++)
    {
      word_t word = _get_word(src, srcindx+i*WORD_SIZE);
      _set_word(dst, dstindx+i*WORD_SIZE, word);
    }

    if(bits_in_last_word > 0)
    {
      word_t src_word = _get_word(src, srcindx+i*WORD_SIZE);
      word_t dst_word = _get_word(dst, dstindx+i*WORD_SIZE);

      word_t mask = bitmask64(bits_in_last_word);
      word_t word = bitmask_merge(src_word, dst_word, mask);

      _set_word(dst, dstindx+num_of_full_words*WORD_SIZE, word);
    }
  }
  else
  {
    // Work right to left
    DEBUG_PRINT("work right to left\n");

    for(i = 0; i < num_of_full_words; i++)
    {
      word_t word = _get_word(src, srcindx+length-(i+1)*WORD_SIZE);
      _set_word(dst, dstindx+length-(i+1)*WORD_SIZE, word);
    }

    DEBUG_PRINT("Copy %i,%i to %i\n", (int)srcindx, (int)bits_in_last_word,
                                      (int)dstindx);

    if(bits_in_last_word > 0)
    {
      word_t src_word = _get_word(src, srcindx);
      word_t dst_word = _get_word(dst, dstindx);

      word_t mask = bitmask64(bits_in_last_word);
      word_t word = bitmask_merge(src_word, dst_word, mask);
      _set_word(dst, dstindx, word);
    }
  }

  _mask_top_word(dst);
}

// destination and source may be the same bit_array
// and src/dst regions may overlap
void bit_array_copy(BIT_ARRAY* dst, bit_index_t dstindx,
                    const BIT_ARRAY* src, bit_index_t srcindx,
                    bit_index_t length)
{
  assert(srcindx + length <= src->num_of_bits);
  assert(dstindx <= dst->num_of_bits);
  _array_copy(dst, dstindx, src, srcindx, length);
  DEBUG_VALIDATE(dst);
}

// Clone `src` into `dst`. Resizes `dst`.
void bit_array_copy_all(BIT_ARRAY* dst, const BIT_ARRAY* src)
{
  bit_array_resize_critical(dst, src->num_of_bits);
  memmove(dst->words, src->words, src->num_of_words * sizeof(word_t));
  DEBUG_VALIDATE(dst);
}


//
// Logic operators
//

// Destination can be the same as one or both of the sources
void bit_array_and(BIT_ARRAY* dst, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
{
  // Ensure dst array is big enough
  word_addr_t max_bits = MAX(src1->num_of_bits, src2->num_of_bits);
  bit_array_ensure_size_critical(dst, max_bits);

  word_addr_t min_words = MIN(src1->num_of_words, src2->num_of_words);

  word_addr_t i;

  for(i = 0; i < min_words; i++)
  {
    dst->words[i] = src1->words[i] & src2->words[i];
  }

  // Set remaining bits to zero
  for(i = min_words; i < dst->num_of_words; i++)
  {
    dst->words[i] = (word_t)0;
  }

  DEBUG_VALIDATE(dst);
}

// Destination can be the same as one or both of the sources
static void _logical_or_xor(BIT_ARRAY* dst,
                            const BIT_ARRAY* src1,
                            const BIT_ARRAY* src2,
                            char use_xor)
{
  // Ensure dst array is big enough
  bit_array_ensure_size_critical(dst, MAX(src1->num_of_bits, src2->num_of_bits));

  word_addr_t min_words = MIN(src1->num_of_words, src2->num_of_words);
  word_addr_t max_words = MAX(src1->num_of_words, src2->num_of_words);

  word_addr_t i;

  if(use_xor)
  {
    for(i = 0; i < min_words; i++)
      dst->words[i] = src1->words[i] ^ src2->words[i];
  }
  else
  {
    for(i = 0; i < min_words; i++)
      dst->words[i] = src1->words[i] | src2->words[i];
  }

  // Copy remaining bits from longer src array
  if(min_words != max_words)
  {
    const BIT_ARRAY* longer = src1->num_of_words > src2->num_of_words ? src1 : src2;

    for(i = min_words; i < max_words; i++)
    {
      dst->words[i] = longer->words[i];
    }
  }

  // Set remaining bits to zero
  size_t size = (dst->num_of_words - max_words) * sizeof(word_t);
  memset(dst->words + max_words, 0, size);

  DEBUG_VALIDATE(dst);
}

void bit_array_or(BIT_ARRAY* dst, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
{
  _logical_or_xor(dst, src1, src2, 0);
}

// Destination can be the same as one or both of the sources
void bit_array_xor(BIT_ARRAY* dst, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
{
  _logical_or_xor(dst, src1, src2, 1);
}

// If dst is longer than src, top bits are set to 1
void bit_array_not(BIT_ARRAY* dst, const BIT_ARRAY* src)
{
  bit_array_ensure_size_critical(dst, src->num_of_bits);

  word_addr_t i;

  for(i = 0; i < src->num_of_words; i++)
  {
    dst->words[i] = ~(src->words[i]);
  }

  // Set remaining words to 1s
  for(i = src->num_of_words; i < dst->num_of_words; i++)
  {
    dst->words[i] = WORD_MAX;
  }

  _mask_top_word(dst);

  DEBUG_VALIDATE(dst);
}

//
// Comparisons
//

// Compare two bit arrays by value stored, with index 0 being the Least
// Significant Bit (LSB). Arrays do not have to be the same length.
// Example: ..0101 (5) > ...0011 (3) [index 0 is LSB at right hand side]
// Sorts on length if all zeros: (0,0) < (0,0,0)
// returns:
//  >0 iff bitarr1 > bitarr2
//   0 iff bitarr1 == bitarr2
//  <0 iff bitarr1 < bitarr2
int bit_array_cmp(const BIT_ARRAY* bitarr1, const BIT_ARRAY* bitarr2)
{
  word_addr_t i;
  word_t word1, word2;
  word_addr_t min_words = bitarr1->num_of_words;

  // i is unsigned so break when i == 0
  if(bitarr1->num_of_words > bitarr2->num_of_words) {
    min_words = bitarr2->num_of_words;
    for(i = bitarr1->num_of_words-1; ; i--) {
      if(bitarr1->words[i]) return 1;
      if(i == bitarr2->num_of_words) break;
    }
  }
  else if(bitarr1->num_of_words < bitarr2->num_of_words) {
    for(i = bitarr2->num_of_words-1; ; i--) {
      if(bitarr2->words[i]) return 1;
      if(i == bitarr1->num_of_words) break;
    }
  }

  if(min_words == 0) return 0;

  for(i = min_words-1; ; i--)
  {
    word1 = bitarr1->words[i];
    word2 = bitarr2->words[i];
    if(word1 != word2) return (word1 > word2 ? 1 : -1);
    if(i == 0) break;
  }

  if(bitarr1->num_of_bits == bitarr2->num_of_bits) return 0;
  return bitarr1->num_of_bits > bitarr2->num_of_bits ? 1 : -1;
}

// Compare two bit arrays by value stored, with index 0 being the Most
// Significant Bit (MSB). Arrays do not have to be the same length.
// Example: 10.. > 01.. [index 0 is MSB at left hand side]
// Sorts on length if all zeros: (0,0) < (0,0,0)
// returns:
//  >0 iff bitarr1 > bitarr2
//   0 iff bitarr1 == bitarr2
//  <0 iff bitarr1 < bitarr2
int bit_array_cmp_big_endian(const BIT_ARRAY* bitarr1, const BIT_ARRAY* bitarr2)
{
  word_addr_t min_words = MAX(bitarr1->num_of_words, bitarr2->num_of_words);

  word_addr_t i;
  word_t word1, word2;

  for(i = 0; i < min_words; i++) {
    word1 = _reverse_word(bitarr1->words[i]);
    word2 = _reverse_word(bitarr2->words[i]);
    if(word1 != word2) return (word1 > word2 ? 1 : -1);
  }

  // Check remaining words. Only one of these loops will execute
  for(; i < bitarr1->num_of_words; i++)
    if(bitarr1->words[i]) return 1;
  for(; i < bitarr2->num_of_words; i++)
    if(bitarr2->words[i]) return -1;

  if(bitarr1->num_of_bits == bitarr2->num_of_bits) return 0;
  return bitarr1->num_of_bits > bitarr2->num_of_bits ? 1 : -1;
}

// compare bitarr with (bitarr2 << pos)
// bit_array_cmp(bitarr1, bitarr2<<pos)
// returns:
//  >0 iff bitarr1 > bitarr2
//   0 iff bitarr1 == bitarr2
//  <0 iff bitarr1 < bitarr2
int bit_array_cmp_words(const BIT_ARRAY *arr1,
                        bit_index_t pos, const BIT_ARRAY *arr2)
{
  if(arr1->num_of_bits == 0 && arr2->num_of_bits == 0)
  {
    return 0;
  }

  bit_index_t top_bit1 = 0, top_bit2 = 0;

  char arr1_zero = !bit_array_find_last_set_bit(arr1, &top_bit1);
  char arr2_zero = !bit_array_find_last_set_bit(arr2, &top_bit2);

  if(arr1_zero && arr2_zero) return 0;
  if(arr1_zero) return -1;
  if(arr2_zero) return 1;

  bit_index_t top_bit2_offset = top_bit2 + pos;

  if(top_bit1 != top_bit2_offset) {
    return top_bit1 > top_bit2_offset ? 1 : -1;
  }

  word_addr_t i;
  word_t word1, word2;

  for(i = top_bit2 / WORD_SIZE; i > 0; i--)
  {
    word1 = _get_word(arr1, pos + i * WORD_SIZE);
    word2 = arr2->words[i];

    if(word1 > word2) return 1;
    if(word1 < word2) return -1;
  }

  word1 = _get_word(arr1, pos);
  word2 = arr2->words[0];

  if(word1 > word2) return 1;
  if(word1 < word2) return -1;

  // return 1 if arr1[0..pos] != 0, 0 otherwise

  // Whole words
  word_addr_t num_words = pos / WORD_SIZE;

  for(i = 0; i < num_words; i++)
  {
    if(arr1->words[i] > 0)
    {
      return 1;
    }
  }

  word_offset_t bits_remaining = pos - num_words * WORD_SIZE;

  if(arr1->words[num_words] & bitmask64(bits_remaining))
  {
    return 1;
  }

  return 0;
}


//
// Reverse -- coords may wrap around
//

// No bounds checking
// length cannot be zero
static void _reverse_region(BIT_ARRAY* bitarr,
                            bit_index_t start,
                            bit_index_t length)
{
  bit_index_t left = start;
  bit_index_t right = (start + length - WORD_SIZE) % bitarr->num_of_bits;

  while(length >= 2 * WORD_SIZE)
  {
    // Swap entire words
    word_t left_word = _get_word_cyclic(bitarr, left);
    word_t right_word = _get_word_cyclic(bitarr, right);

    // reverse words individually
    left_word = _reverse_word(left_word);
    right_word = _reverse_word(right_word);

    // Swap
    _set_word_cyclic(bitarr, left, right_word);
    _set_word_cyclic(bitarr, right, left_word);

    // Update
    left = (left + WORD_SIZE) % bitarr->num_of_bits;
    right = (right < WORD_SIZE ? right + bitarr->num_of_bits : right) - WORD_SIZE;
    length -= 2 * WORD_SIZE;
  }

  word_t word, rev;

  if(length == 0)
  {
    return;
  }
  else if(length > WORD_SIZE)
  {
    // Words overlap
    word_t left_word = _get_word_cyclic(bitarr, left);
    word_t right_word = _get_word_cyclic(bitarr, right);

    rev = _reverse_word(left_word);
    right_word = _reverse_word(right_word);

    // fill left 64 bits with right word rev
    _set_word_cyclic(bitarr, left, right_word);

    // Now do remaining bits (length is between 1 and 64 bits)
    left += WORD_SIZE;
    length -= WORD_SIZE;

    word = _get_word_cyclic(bitarr, left);
  }
  else
  {
    word = _get_word_cyclic(bitarr, left);
    rev = _reverse_word(word);
  }

  rev >>= WORD_SIZE - length;
  word_t mask = bitmask64(length);

  word = bitmask_merge(rev, word, mask);

  _set_word_cyclic(bitarr, left, word);
}

void bit_array_reverse_region(BIT_ARRAY* bitarr, bit_index_t start, bit_index_t len)
{
  assert(start + len <= bitarr->num_of_bits);
  if(len > 0) _reverse_region(bitarr, start, len);
  DEBUG_VALIDATE(bitarr);
}

void bit_array_reverse(BIT_ARRAY* bitarr)
{
  if(bitarr->num_of_bits > 0) _reverse_region(bitarr, 0, bitarr->num_of_bits);
  DEBUG_VALIDATE(bitarr);
}

//
// Shift left / right
//

// Shift towards MSB / higher index
void bit_array_shift_left(BIT_ARRAY* bitarr, bit_index_t shift_dist, char fill)
{
  if(shift_dist >= bitarr->num_of_bits)
  {
    fill ? bit_array_set_all(bitarr) : bit_array_clear_all(bitarr);
    return;
  }
  else if(shift_dist == 0)
  {
    return;
  }

  FillAction action = fill ? FILL_REGION : ZERO_REGION;

  bit_index_t cpy_length = bitarr->num_of_bits - shift_dist;
  _array_copy(bitarr, shift_dist, bitarr, 0, cpy_length);
  _set_region(bitarr, 0, shift_dist, action);
}

// shift left extend - don't truncate bits when shifting UP, instead
// make room for them.
void bit_array_shift_left_extend(BIT_ARRAY* bitarr, bit_index_t shift_dist,
                                 char fill)
{
   bit_index_t newlen = bitarr->num_of_bits + shift_dist;
   bit_index_t cpy_length = bitarr->num_of_bits;

  if(shift_dist == 0)
  {
    return;
  }

  bit_array_resize_critical(bitarr, newlen);

  FillAction action = fill ? FILL_REGION : ZERO_REGION;
  _array_copy(bitarr, shift_dist, bitarr, 0, cpy_length);
  _set_region(bitarr, 0, shift_dist, action);
}

// Shift towards LSB / lower index
void bit_array_shift_right(BIT_ARRAY* bitarr, bit_index_t shift_dist, char fill)
{
  if(shift_dist >= bitarr->num_of_bits)
  {
    fill ? bit_array_set_all(bitarr) : bit_array_clear_all(bitarr);
    return;
  }
  else if(shift_dist == 0)
  {
    return;
  }

  FillAction action = fill ? FILL_REGION : ZERO_REGION;

  bit_index_t cpy_length = bitarr->num_of_bits - shift_dist;
  bit_array_copy(bitarr, 0, bitarr, shift_dist, cpy_length);

  _set_region(bitarr, cpy_length, shift_dist, action);
}

//
// Cycle
//

// Cycle towards index 0
void bit_array_cycle_right(BIT_ARRAY* bitarr, bit_index_t cycle_dist)
{
  if(bitarr->num_of_bits == 0)
  {
    return;
  }

  cycle_dist = cycle_dist % bitarr->num_of_bits;

  if(cycle_dist == 0)
  {
    return;
  }

  bit_index_t len1 = cycle_dist;
  bit_index_t len2 = bitarr->num_of_bits - cycle_dist;

  _reverse_region(bitarr, 0, len1);
  _reverse_region(bitarr, len1, len2);
  bit_array_reverse(bitarr);
}

// Cycle away from index 0
void bit_array_cycle_left(BIT_ARRAY* bitarr, bit_index_t cycle_dist)
{
  if(bitarr->num_of_bits == 0)
  {
    return;
  }

  cycle_dist = cycle_dist % bitarr->num_of_bits;

  if(cycle_dist == 0)
  {
    return;
  }

  bit_index_t len1 = bitarr->num_of_bits - cycle_dist;
  bit_index_t len2 = cycle_dist;

  _reverse_region(bitarr, 0, len1);
  _reverse_region(bitarr, len1, len2);
  bit_array_reverse(bitarr);
}

//
// Next permutation
//

static word_t _next_permutation(word_t v)
{
  assert(v);
  // From http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
  word_t t = v | (v - 1); // t gets v's least significant 0 bits set to 1
  // Next set to 1 the most significant bit to change,
  // set to 0 the least significant ones, and add the necessary 1 bits.
  return (t+1) | (((~t & (t+1)) - 1) >> (trailing_zeros(v) + 1));
}

// Get the next permutation of an array with a fixed size and given number of
// bits set.  Also known as next lexicographic permutation.
// Given a bit array find the next lexicographic orginisation of the bits
// Number of possible combinations given by (size choose bits_set) i.e. nCk
// 00011 -> 00101 -> 00110 -> 01001 -> 01010 ->
// 01100 -> 10001 -> 10010 -> 10100 -> 11000 -> 00011 (back to start)
void bit_array_next_permutation(BIT_ARRAY* bitarr)
{
  if(bitarr->num_of_bits == 0)
  {
    return;
  }

  word_addr_t w;

  char carry = 0;
  word_offset_t top_bits = bitset64_idx(bitarr->num_of_bits);

  for(w = 0; w < bitarr->num_of_words; w++)
  {
    word_t mask
      = (w < bitarr->num_of_words - 1 || top_bits == 0) ? WORD_MAX
                                                        : bitmask64(top_bits);

    if(bitarr->words[w] > 0 &&
       (bitarr->words[w] | (bitarr->words[w]-1)) == mask)
    {
      // Bits in this word cannot be moved forward
      carry = 1;
    }
    else if(carry)
    {
      // 0111 -> 1000, 1000 -> 1001
      word_t tmp = bitarr->words[w] + 1;

      // Count bits previously set
      bit_index_t bits_previously_set = POPCOUNT(bitarr->words[w]);

      // set new word
      bitarr->words[w] = tmp;

      // note: w is unsigned
      // Zero words while counting bits set
      while(w > 0)
      {
        bits_previously_set += POPCOUNT(bitarr->words[w-1]);
        bitarr->words[w-1] = 0;
        w--;
      }

      // Set bits at the beginning
      SET_REGION(bitarr, 0, bits_previously_set - POPCOUNT(tmp));

      carry = 0;
      break;
    }
    else if(bitarr->words[w] > 0)
    {
      bitarr->words[w] = _next_permutation(bitarr->words[w]);
      break;
    }
  }

  if(carry)
  {
    // Loop around
    bit_index_t num_bits_set = bit_array_num_bits_set(bitarr);
    bit_array_clear_all(bitarr);
    SET_REGION(bitarr, 0, num_bits_set);
  }

  DEBUG_VALIDATE(bitarr);
}


//
// Interleave
//

// dst cannot point to the same bit array as src1 or src2
// src1, src2 may point to the same bit array
// abcd 1234 -> a1b2c3d4
// 0011 0000 -> 00001010
// 1111 0000 -> 10101010
// 0101 1010 -> 01100110
void bit_array_interleave(BIT_ARRAY* dst,
                          const BIT_ARRAY* src1,
                          const BIT_ARRAY* src2)
{
  // dst cannot be either src1 or src2
  assert(dst != src1 && dst != src2);
  // Behaviour undefined when src1 length != src2 length",
  assert(src1->num_of_bits == src2->num_of_bits);

  // Need at least src1->num_of_words + src2->num_of_words
  size_t nwords = MIN(src1->num_of_words + src2->num_of_words, 2);
  _bit_array_ensure_nwords(dst, nwords, __FILE__, __LINE__, __func__);
  dst->num_of_bits = src1->num_of_bits + src2->num_of_bits;
  dst->num_of_words = roundup_bits2words64(dst->num_of_bits);

  word_addr_t i, j;

  for(i = 0, j = 0; i < src1->num_of_words; i++)
  {
    word_t a = src1->words[i];
    word_t b = src2->words[i];

    dst->words[j++] =  morton_table0[(a      ) & 0xff] |
                       morton_table1[(b      ) & 0xff] |
                      (morton_table0[(a >>  8) & 0xff] << 16) |
                      (morton_table1[(b >>  8) & 0xff] << 16) |
                      (morton_table0[(a >> 16) & 0xff] << 32) |
                      (morton_table1[(b >> 16) & 0xff] << 32) |
                      (morton_table0[(a >> 24) & 0xff] << 48) |
                      (morton_table1[(b >> 24) & 0xff] << 48);

    dst->words[j++] =  morton_table0[(a >> 32) & 0xff] |
                       morton_table1[(b >> 32) & 0xff] |
                      (morton_table0[(a >> 40) & 0xff] << 16) |
                      (morton_table1[(b >> 40) & 0xff] << 16) |
                      (morton_table0[(a >> 48) & 0xff] << 32) |
                      (morton_table1[(b >> 48) & 0xff] << 32) |
                      (morton_table0[(a >> 56)       ] << 48) |
                      (morton_table1[(b >> 56)       ] << 48);
  }

  DEBUG_VALIDATE(dst);
}

//
// Random
//

// Set bits randomly with probability prob : 0 <= prob <= 1
void bit_array_random(BIT_ARRAY* bitarr, float prob)
{
  assert(prob >= 0 && prob <= 1);

  if(bitarr->num_of_bits == 0)
  {
    return;
  }
  else if(prob == 1)
  {
    bit_array_set_all(bitarr);
    return;
  }

  // rand() generates number between 0 and RAND_MAX inclusive
  // therefore we want to check if rand() <= p
  long p = RAND_MAX * prob;

  _seed_rand();

  word_addr_t w;
  word_offset_t o;

  // Initialise to zero
  memset(bitarr->words, 0, bitarr->num_of_words * sizeof(word_t));

  for(w = 0; w < bitarr->num_of_words - 1; w++)
  {
    for(o = 0; o < WORD_SIZE; o++)
    {
      if(rand() <= p)
      {
        bitarr->words[w] |= ((word_t)0x1 << o);
      }
    }
  }

  // Top word
  word_offset_t bits_in_last_word = bits_in_top_word(bitarr->num_of_bits);
  w = bitarr->num_of_words - 1;

  for(o = 0; o < bits_in_last_word; o++)
  {
    if(rand() <= p)
    {
      bitarr->words[w] |= ((word_t)0x1 << o);
    }
  }

  DEBUG_VALIDATE(bitarr);
}

// Shuffle the bits in an array randomly
void bit_array_shuffle(BIT_ARRAY* bitarr)
{
  if(bitarr->num_of_bits == 0)
    return;

  _seed_rand();

  bit_index_t i, j;

  for(i = bitarr->num_of_bits - 1; i > 0; i--)
  {
    j = (bit_index_t)rand() % i;

    // Swap i and j
    char x = (bitarr->words[bitset64_wrd(i)] >> bitset64_idx(i)) & 0x1;
    char y = (bitarr->words[bitset64_wrd(j)] >> bitset64_idx(j)) & 0x1;

    if(!y)
      bitarr->words[bitset64_wrd(i)] &= ~((word_t)0x1 << bitset64_idx(i));
    else
      bitarr->words[bitset64_wrd(i)] |= (word_t)0x1 << bitset64_idx(i);

    if(!x)
      bitarr->words[bitset64_wrd(j)] &= ~((word_t)0x1 << bitset64_idx(j));
    else
      bitarr->words[bitset64_wrd(j)] |= (word_t)0x1 << bitset64_idx(j);
  }

  DEBUG_VALIDATE(bitarr);
}

//
// Arithmetic
//

// Returns 1 on sucess, 0 if value in array is too big
char bit_array_as_num(const BIT_ARRAY* bitarr, uint64_t* result)
{
  if(bitarr->num_of_bits == 0)
  {
    *result = 0;
    return 1;
  }

  word_addr_t w;

  for(w = bitarr->num_of_words-1; w > 0; w--)
  {
    if(bitarr->words[w] > 0)
    {
      return 0;
    }
  }

  *result = bitarr->words[0];
  return 1;
}


// 1 iff bitarr > value
// 0 iff bitarr == value
// -1 iff bitarr < value
int bit_array_cmp_uint64(const BIT_ARRAY* bitarr, uint64_t value)
{
  uint64_t arr_num = 0;

  // If cannot put bitarr in uint64, it is > value
  if(!bit_array_as_num(bitarr, &arr_num)) return 1;

  if(arr_num > value)      return  1;
  else if(arr_num < value) return -1;
  else                     return  0;
}

// If value is zero, no change is made
void bit_array_add_uint64(BIT_ARRAY* bitarr, uint64_t value)
{
  if(value == 0)
  {
    return;
  }
  else if(bitarr->num_of_bits == 0)
  {
    bit_array_resize_critical(bitarr, WORD_SIZE - leading_zeros(value));
    bitarr->words[0] = (word_t)value;
    return;
  }

  char carry = 0;
  word_addr_t i;

  for(i = 0; i < bitarr->num_of_words; i++)
  {
    if(WORD_MAX - bitarr->words[i] < value)
    {
      carry = 1;
      bitarr->words[i] += value;
      value = 1;
    }
    else
    {
      // Carry is absorbed
      bitarr->words[i] += value;
      carry = 0;
      break;
    }
  }

  if(carry)
  {
    // Bit array full, need another bit after all words filled
    bit_array_resize_critical(bitarr, bitarr->num_of_words * WORD_SIZE + 1);

    // Set top word to 1
    bitarr->words[bitarr->num_of_words-1] = 1;
  }
  else
  {
    word_t final_word = bitarr->words[bitarr->num_of_words-1];
    word_offset_t expected_bits = bits_in_top_word(bitarr->num_of_bits);
    word_offset_t actual_bits = WORD_SIZE - leading_zeros(final_word);

    if(actual_bits > expected_bits)
    {
      // num_of_bits has increased -- num_of_words has not
      bitarr->num_of_bits += (actual_bits - expected_bits);
    }
  }
}

// If value is greater than bitarr, bitarr is not changed and 0 is returned
// Returns 1 on success, 0 if value > bitarr
char bit_array_sub_uint64(BIT_ARRAY* bitarr, uint64_t value)
{
  if(value == 0)
  {
    return 1;
  }
  else if(bitarr->words[0] >= value)
  {
    bitarr->words[0] -= value;
    return 1;
  }

  value -= bitarr->words[0];

  word_addr_t i;

  for(i = 1; i < bitarr->num_of_words; i++)
  {
    if(bitarr->words[i] > 0)
    {
      // deduct one
      bitarr->words[i]--;

      for(; i > 0; i--)
      {
        bitarr->words[i] = WORD_MAX;
      }

      // -1 since we've already deducted 1
      bitarr->words[0] = WORD_MAX - value - 1;

      return 1;
    }
  }

  // subtract value is greater than array
  return 0;
}

//
// Arithmetic between bit arrays
//

// src1, src2 and dst can all be the same BIT_ARRAY
static void _arithmetic(BIT_ARRAY* dst,
                        const BIT_ARRAY* src1,
                        const BIT_ARRAY* src2,
                        char subtract)
{
  word_addr_t max_words = MAX(src1->num_of_words, src2->num_of_words);

  // Adding: dst_words >= max(src1 words, src2 words)
  // Subtracting: dst_words is >= src1->num_of_words

  char carry = subtract ? 1 : 0;

  word_addr_t i;
  word_t word1, word2;

  for(i = 0; i < max_words; i++)
  {
    word1 = (i < src1->num_of_words ? src1->words[i] : 0);
    word2 = (i < src2->num_of_words ? src2->words[i] : 0);

    if(subtract)
      word2 = ~word2;

    dst->words[i] = word1 + word2 + carry;
    // Update carry
    carry = WORD_MAX - word1 < word2 || WORD_MAX - word1 - word2 < (word_t)carry;
  }

  if(subtract)
  {
    carry = 0;
  }
  else
  {
    // Check last word
    word_offset_t bits_on_last_word = bits_in_top_word(dst->num_of_bits);

    if(bits_on_last_word < WORD_SIZE)
    {
      word_t mask = bitmask64(bits_on_last_word);

      if(dst->words[max_words-1] > mask)
      {
        // Array has overflowed, increase size
        dst->num_of_bits++;
      }
    }
    else if(carry)
    {
      // Carry onto a new word
      if(dst->num_of_words == max_words)
      {
        // Need to resize for the carry bit
        bit_array_resize_critical(dst, dst->num_of_bits+1);
      }

      dst->words[max_words] = (word_t)1;
    }
  }

  // Zero the rest of dst array
  for(i = max_words+carry; i < dst->num_of_words; i++)
  {
    dst->words[i] = (word_t)0;
  }

  DEBUG_VALIDATE(dst);
}

// src1, src2 and dst can all be the same BIT_ARRAY
// If dst is shorter than either of src1, src2, it is enlarged
void bit_array_add(BIT_ARRAY* dst, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
{
  bit_array_ensure_size_critical(dst, MAX(src1->num_of_bits, src2->num_of_bits));
  _arithmetic(dst, src1, src2, 0);
}

// dst = src1 - src2
// src1, src2 and dst can all be the same BIT_ARRAY
// If dst is shorter than src1, it will be extended to be as long as src1
// src1 must be greater than or equal to src2 (src1 >= src2)
void bit_array_subtract(BIT_ARRAY* dst,
                          const BIT_ARRAY* src1, const BIT_ARRAY* src2)
{
  // subtraction by method of complements:
  // a - b = a + ~b + 1 = src1 + ~src2 +1

  assert(bit_array_cmp(src1, src2) >= 0); // Require src1 >= src2

  bit_array_ensure_size_critical(dst, src1->num_of_bits);
  _arithmetic(dst, src1, src2, 1);
}


// Add `add` to `bitarr` at `pos`
// Bounds checking not needed as out of bounds is valid
void bit_array_add_word(BIT_ARRAY *bitarr, bit_index_t pos, uint64_t add)
{
  DEBUG_VALIDATE(bitarr);

  if(add == 0)
  {
    return;
  }
  else if(pos >= bitarr->num_of_bits)
  {
    // Resize and add!
    bit_index_t num_bits_required = pos + (WORD_SIZE - leading_zeros(add));
    bit_array_resize_critical(bitarr, num_bits_required);
    _set_word(bitarr, pos, (word_t)add);
    return;
  }

  /*
  char str[1000];
  printf(" add_word: %s\n", bit_array_to_str_rev(bitarr, str));
  printf("     word: %s [pos: %i]\n", _word_to_str(add, str), (int)pos);
  */

  word_t w = _get_word(bitarr, pos);
  word_t sum = w + add;
  char carry = WORD_MAX - w < add;

  // Ensure array is big enough
  bit_index_t num_bits_required = pos + (carry ? WORD_SIZE + 1
                                               : (WORD_SIZE - leading_zeros(sum)));

  bit_array_ensure_size(bitarr, num_bits_required);

  _set_word(bitarr, pos, sum);
  pos += WORD_SIZE;

  if(carry)
  {
    word_offset_t offset = pos % WORD_SIZE;
    word_addr_t addr = bitset64_wrd(pos);

    add = (word_t)0x1 << offset;
    carry = (WORD_MAX - bitarr->words[addr] < add);
    sum = bitarr->words[addr] + add;

    num_bits_required = addr * WORD_SIZE +
                        (carry ? WORD_SIZE + 1 : (WORD_SIZE - leading_zeros(sum)));

    bit_array_ensure_size(bitarr, num_bits_required);

    bitarr->words[addr++] = sum;

    if(carry)
    {
      while(addr < bitarr->num_of_words && bitarr->words[addr] == WORD_MAX)
      {
        bitarr->words[addr++] = 0;
      }

      if(addr == bitarr->num_of_words)
      {
        bit_array_resize_critical(bitarr, addr * WORD_SIZE + 1);
      }
      else if(addr == bitarr->num_of_words-1 &&
              bitarr->words[addr] == bitmask64(bits_in_top_word(bitarr->num_of_bits)))
      {
        bit_array_resize_critical(bitarr, bitarr->num_of_bits + 1);
      }

      bitarr->words[addr]++;
    }
  }

  DEBUG_VALIDATE(bitarr);
}

// Add `add` to `bitarr` at `pos`
// Bounds checking not needed as out of bounds is valid
void bit_array_add_words(BIT_ARRAY *bitarr, bit_index_t pos, const BIT_ARRAY *add)
{
  assert(bitarr != add); // bitarr and add cannot point to the same bit array

  bit_index_t add_top_bit_set;

  if(!bit_array_find_last_set_bit(add, &add_top_bit_set))
  {
    // No bits set in add
    return;
  }
  else if(pos >= bitarr->num_of_bits)
  {
    // Just resize and copy!
    bit_index_t num_bits_required = pos + add_top_bit_set + 1;
    bit_array_resize_critical(bitarr, num_bits_required);
    _array_copy(bitarr, pos, add, 0, add->num_of_bits);
    return;
  }
  else if(pos == 0)
  {
    bit_array_add(bitarr, bitarr, add);
    return;
  }

  /*
  char str[1000];
  printf(" add_words1: %s\n", bit_array_to_str_rev(bitarr, str));
  printf(" add_words2: %s\n", bit_array_to_str_rev(add, str));
  printf(" [pos: %i]\n", (int)pos);
  */

  bit_index_t num_bits_required = pos + add_top_bit_set + 1;
  bit_array_ensure_size(bitarr, num_bits_required);

  word_addr_t first_word = bitset64_wrd(pos);
  word_offset_t first_offset = bitset64_idx(pos);

  word_t w = add->words[0] << first_offset;
  unsigned char carry = (WORD_MAX - bitarr->words[first_word] < w);

  bitarr->words[first_word] += w;

  word_addr_t i = first_word + 1;
  bit_index_t offset = WORD_SIZE - first_offset;

  for(; carry || offset <= add_top_bit_set; i++, offset += WORD_SIZE)
  {
    w = offset < add->num_of_bits ? _get_word(add, offset) : (word_t)0;

    if(i >= bitarr->num_of_words)
    {
      // Extend by a word
      bit_array_resize_critical(bitarr, (bit_index_t)(i+1)*WORD_SIZE+1);
    }

    word_t prev = bitarr->words[i];

    bitarr->words[i] += w + carry;

    carry = (WORD_MAX - prev < w || (carry && prev + w == WORD_MAX)) ? 1 : 0;
  }

  word_offset_t top_bits
    = WORD_SIZE - leading_zeros(bitarr->words[bitarr->num_of_words-1]);

  bit_index_t min_bits = (bitarr->num_of_words-1)*WORD_SIZE + top_bits;

  if(bitarr->num_of_bits < min_bits)
  {
    // Extend within the last word
    bitarr->num_of_bits = min_bits;
  }

  DEBUG_VALIDATE(bitarr);
}

char bit_array_sub_word(BIT_ARRAY* bitarr, bit_index_t pos, word_t minus)
{
  DEBUG_VALIDATE(bitarr);

  if(minus == 0)
  {
    return 1;
  }

  word_t w = _get_word(bitarr, pos);

  if(w >= minus)
  {
    _set_word(bitarr, pos, w - minus);
    DEBUG_VALIDATE(bitarr);
    return 1;
  }

  minus -= w;

  bit_index_t offset;
  for(offset = pos + WORD_SIZE; offset < bitarr->num_of_bits; offset += WORD_SIZE)
  {
    w = _get_word(bitarr, offset);

    if(w > 0)
    {
      // deduct one
      _set_word(bitarr, offset, w - 1);

      SET_REGION(bitarr, pos, offset-pos);

      // -1 since we've already deducted 1
      minus--;

      _set_word(bitarr, pos, WORD_MAX - minus);

      DEBUG_VALIDATE(bitarr);
      return 1;
    }
  }

  DEBUG_VALIDATE(bitarr);

  return 0;
}

char bit_array_sub_words(BIT_ARRAY* bitarr, bit_index_t pos, BIT_ARRAY* minus)
{
  assert(bitarr != minus); // bitarr and minus cannot point to the same bit array

  int cmp = bit_array_cmp_words(bitarr, pos, minus);

  if(cmp == 0)
  {
    bit_array_clear_all(bitarr);
    return 1;
  }
  else if(cmp < 0)
  {
    return 0;
  }

  bit_index_t bitarr_length = bitarr->num_of_bits;

  bit_index_t bitarr_top_bit_set;
  bit_array_find_last_set_bit(bitarr, &bitarr_top_bit_set);

  // subtraction by method of complements:
  // a - b = a + ~b + 1 = src1 + ~src2 +1

  bit_array_not(minus, minus);

  bit_array_add_words(bitarr, pos, minus);
  bit_array_add_word(bitarr, pos, (word_t)1);

  bit_array_sub_word(bitarr, pos+minus->num_of_bits, 1);
  bit_array_resize(bitarr, bitarr_length);

  bit_array_not(minus, minus);

  DEBUG_VALIDATE(bitarr);

  return 1;
}

void bit_array_mul_uint64(BIT_ARRAY *bitarr, uint64_t multiplier)
{
  if(bitarr->num_of_bits == 0 || multiplier == 1)
  {
    return;
  }
  else if(multiplier == 0)
  {
    bit_array_clear_all(bitarr);
    return;
  }

  bit_index_t i;

  for(i = bitarr->num_of_bits; i > 0; i--)
  {
    if(bit_array_get(bitarr, i-1))
    {
      bit_array_clear(bitarr, i-1);
      bit_array_add_word(bitarr, i-1, multiplier);
    }
  }

  DEBUG_VALIDATE(bitarr);
}

void bit_array_multiply(BIT_ARRAY *dst, BIT_ARRAY *src1, BIT_ARRAY *src2)
{
  if(src1->num_of_bits == 0 || src2->num_of_bits == 0)
  {
    bit_array_clear_all(dst);
    return;
  }

  // Cannot pass the same array as dst, src1 AND src2
  assert(dst != src1 || dst != src2);

  // Dev: multiplier == 1?

  BIT_ARRAY *read_arr, *add_arr;

  if(src1 == dst)
  {
    read_arr = src1;
    add_arr = src2;
  }
  else
  {
    read_arr = src2;
    add_arr = src1;
  }

  if(dst != src1 && dst != src2)
  {
    bit_array_clear_all(dst);
  }

  bit_index_t i;

  for(i = read_arr->num_of_bits; i > 0; i--)
  {
    if(bit_array_get(read_arr, i-1))
    {
      bit_array_clear(dst, i-1);
      bit_array_add_words(dst, i-1, add_arr);
    }
  }

  DEBUG_VALIDATE(dst);
}

// bitarr = round_down(bitarr / divisor)
// rem = bitarr % divisor
void bit_array_div_uint64(BIT_ARRAY *bitarr, uint64_t divisor, uint64_t *rem)
{
  assert(divisor != 0); // cannot divide by zero

  bit_index_t div_top_bit = 63 - leading_zeros(divisor);
  bit_index_t bitarr_top_bit;

  if(!bit_array_find_last_set_bit(bitarr, &bitarr_top_bit))
  {
    *rem = 0;
    return;
  }

  if(bitarr_top_bit < div_top_bit)
  {
    *rem = bitarr->words[0];
    bit_array_clear_all(bitarr);
    return;
  }

  // When div is shifted by offset, their top set bits are aligned
  bit_index_t offset = bitarr_top_bit - div_top_bit;

  uint64_t tmp = _get_word(bitarr, offset);
  _set_word(bitarr, offset, (word_t)0);

  // Carry if 1 if the top bit was set before left shift
  char carry = 0;

  // offset unsigned so break when offset == 0
  while(1)
  {
    if(carry)
    {
      // (carry:tmp) - divisor = (WORD_MAX+1+tmp)-divisor
      tmp = WORD_MAX - divisor + tmp + 1;
      bit_array_set(bitarr, offset);
    }
    else if(tmp >= divisor)
    {
      tmp -= divisor;
      bit_array_set(bitarr, offset);
    }
    else
    {
      bit_array_clear(bitarr, offset);
    }

    if(offset == 0)
      break;

    offset--;

    // Is the top bit set (that we're about to shift off)?
    carry = tmp & 0x8000000000000000;

    tmp <<= 1;
    tmp |= bit_array_get(bitarr, offset);
  }

  *rem = tmp;
}

// Results in:
//   quotient = dividend / divisor
//   dividend = dividend % divisor
// (dividend is used to return the remainder)
void bit_array_divide(BIT_ARRAY *dividend, BIT_ARRAY *quotient, BIT_ARRAY *divisor)
{
  assert(bit_array_cmp_uint64(divisor, 0) != 0); // Cannot divide by zero

  bit_array_clear_all(quotient);

  int cmp = bit_array_cmp(dividend, divisor);

  if(cmp == 0)
  {
    bit_array_ensure_size(quotient, 1);
    bit_array_set(quotient, 0);
    bit_array_clear_all(dividend);
    return;
  }
  else if(cmp < 0)
  {
    // dividend is < divisor, quotient is zero -- done
    return;
  }

  // now we know: dividend > divisor, quotient is zero'd,
  //              dividend != 0, divisor != 0
  bit_index_t dividend_top_bit = 0, div_top_bit = 0;

  bit_array_find_last_set_bit(dividend, &dividend_top_bit);
  bit_array_find_last_set_bit(divisor, &div_top_bit);

  // When divisor is shifted by offset, their top set bits are aligned
  bit_index_t offset = dividend_top_bit - div_top_bit;

  // offset unsigned so break when offset == 0
  for(; ; offset--)
  {
    if(bit_array_cmp_words(dividend, offset, divisor) >= 0)
    {
      bit_array_sub_words(dividend, offset, divisor);
      bit_array_ensure_size(quotient, offset+1);
      bit_array_set(quotient, offset);
    }

    if(offset == 0)
      break;
  }
}

//
// Read/Write from files
//
// file format is [8 bytes: for number of elements in array][data]
// data is written in little endian order (least sig byte first)
//

// Saves bit array to a file. Returns the number of bytes written
// number of bytes returned should be 8+(bitarr->num_of_bits+7)/8
bit_index_t bit_array_save(const BIT_ARRAY* bitarr, FILE* f)
{
  bit_index_t num_of_bytes = roundup_bits2bytes(bitarr->num_of_bits);
  bit_index_t bytes_written = 0;

  const int endian = 1;
  if(*(uint8_t*)&endian == 1)
  {
    // Little endian machine
    // Write 8 bytes to store the number of bits in the array
    bytes_written += fwrite(&bitarr->num_of_bits, 1, 8, f);

    // Write the array
    bytes_written += fwrite(bitarr->words, 1, num_of_bytes, f);
  }
  else
  {
    // Big endian machine
    uint64_t i, w, whole_words = num_of_bytes/sizeof(word_t);
    uint64_t rem_bytes = num_of_bytes - whole_words*sizeof(word_t);
    uint64_t n_bits = byteswap64(bitarr->num_of_bits);

    // Write 8 bytes to store the number of bits in the array
    bytes_written += fwrite(&n_bits, 1, 8, f);

    // Write the array
    for(i = 0; i < whole_words; i++) {
      w = byteswap64(bitarr->words[i]);
      bytes_written += fwrite(&w, 1, 8, f);
    }

    if(rem_bytes > 0) {
      w = byteswap64(bitarr->words[whole_words]);
      bytes_written += fwrite(&w, 1, rem_bytes, f);
    }
  }

  return bytes_written;
}

// Load a uint64 from little endian format.
// Works for both big and little endian architectures
static inline uint64_t le64_to_cpu(const uint8_t *x)
{
  return (((uint64_t)(x[0]))       | ((uint64_t)(x[1]) << 8)  |
          ((uint64_t)(x[2]) << 16) | ((uint64_t)(x[3]) << 24) |
          ((uint64_t)(x[4]) << 32) | ((uint64_t)(x[5]) << 40) |
          ((uint64_t)(x[6]) << 48) | ((uint64_t)(x[7]) << 56));
}

// Reads bit array from a file. bitarr is resized and filled.
// Returns 1 on success, 0 on failure
char bit_array_load(BIT_ARRAY* bitarr, FILE* f)
{
  // Read in number of bits, return 0 if we can't read in
  bit_index_t num_bits;
  if(fread(&num_bits, 1, 8, f) != 8) return 0;
  num_bits = le64_to_cpu((uint8_t*)&num_bits);

  // Resize
  bit_array_resize_critical(bitarr, num_bits);

  // Have to calculate how many bytes are needed for the file
  // (Note: this may be different from num_of_words * sizeof(word_t))
  bit_index_t num_of_bytes = roundup_bits2bytes(bitarr->num_of_bits);
  if(fread(bitarr->words, 1, num_of_bytes, f) != num_of_bytes) return 0;

  // Fix endianness
  word_addr_t i;
  for(i = 0; i < bitarr->num_of_words; i++)
    bitarr->words[i] = le64_to_cpu((uint8_t*)&bitarr->words[i]);

  // Mask top word
  _mask_top_word(bitarr);
  DEBUG_VALIDATE(bitarr);
  return 1;
}

//
// Hash function
//

/* From: lookup3.c, by Bob Jenkins, May 2006, Public Domain. */
#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

/* From: lookup3.c, by Bob Jenkins, May 2006, Public Domain. */
#define mix(a,b,c) \
{ \
  a -= c;  a ^= rot(c, 4);  c += b; \
  b -= a;  b ^= rot(a, 6);  a += c; \
  c -= b;  c ^= rot(b, 8);  b += a; \
  a -= c;  a ^= rot(c,16);  c += b; \
  b -= a;  b ^= rot(a,19);  a += c; \
  c -= b;  c ^= rot(b, 4);  b += a; \
}

/* From: lookup3.c, by Bob Jenkins, May 2006, Public Domain. */
#define final(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}

/*
From: lookup3.c, by Bob Jenkins, May 2006, Public Domain.
--------------------------------------------------------------------
hashword2() -- same as hashword(), but take two seeds and return two
32-bit values.  pc and pb must both be nonnull, and *pc and *pb must
both be initialized with seeds.  If you pass in (*pb)==0, the output
(*pc) will be the same as the return value from hashword().
--------------------------------------------------------------------
*/
static void hashword2 (
const uint32_t *k,                   /* the key, an array of uint32_t values */
size_t          length,               /* the length of the key, in uint32_ts */
uint32_t       *pc,                      /* IN: seed OUT: primary hash value */
uint32_t       *pb)               /* IN: more seed OUT: secondary hash value */
{
  uint32_t a,b,c;

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + ((uint32_t)(length<<2)) + *pc;
  c += *pb;

  /*------------------------------------------------- handle most of the key */
  while (length > 3)
  {
    a += k[0];
    b += k[1];
    c += k[2];
    mix(a,b,c);
    length -= 3;
    k += 3;
  }

  /*------------------------------------------- handle the last 3 uint32_t's */
  switch(length)                     /* all the case statements fall through */
  {
  case 3 : c+=k[2];
  case 2 : b+=k[1];
  case 1 : a+=k[0];
    final(a,b,c);
  case 0:     /* case 0: nothing left to add */
    break;
  }
  /*------------------------------------------------------ report the result */
  *pc=c; *pb=b;
}

// Pass seed as 0 on first call, pass previous hash value if rehashing due
// to a collision
// Using bob jenkins hash lookup3
uint64_t bit_array_hash(const BIT_ARRAY* bitarr, uint64_t seed)
{
  uint32_t seed32[2];
  memcpy(seed32, &seed, sizeof(uint32_t)*2);

  // Round up length to number 32bit words
  hashword2((uint32_t*)bitarr->words, (bitarr->num_of_bits + 31) / 32,
            &seed32[0], &seed32[1]);

  // XOR with array length. This ensures arrays with different length but same
  // contents have different hash values
  seed ^= bitarr->num_of_bits;

  return seed;
}


//
// Generally useful functions
//

// Generalised 'binary to string' function
// Adds bits to the string in order of lsb to msb
// e.g. 0b11010 (26 in decimal) would come out as "01011"
char* bit_array_word2str(const void *ptr, size_t num_of_bits, char *str)
{
  const uint8_t* d = (const uint8_t*)ptr;

  size_t i;
  for(i = 0; i < num_of_bits; i++)
  {
    uint8_t bit = (d[i/8] >> (i % 8)) & 0x1;
    str[i] = bit ? '1' : '0';
  }
  str[num_of_bits] = '\0';
  return str;
}

char* bit_array_word2str_rev(const void *ptr, size_t num_of_bits, char *str)
{
  const uint8_t* d = (const uint8_t*)ptr;

  size_t i;
  for(i = 0; i < num_of_bits; i++)
  {
    uint8_t bit = (d[i/8] >> (i % 8)) & 0x1;
    str[num_of_bits-1-i] = bit ? '1' : '0';
  }
  str[num_of_bits] = '\0';
  return str;
}
