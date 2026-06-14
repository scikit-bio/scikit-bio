/*
 bit_macros.h
 project: bit array C library
 url: https://github.com/noporpoise/BitArray/
 author: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Dec 2013
*/

#ifndef BITSET_H_
#define BITSET_H_

#include <inttypes.h>
#include <sched.h>

// trailing_zeros is number of least significant zeros
// leading_zeros is number of most significant zeros
#if defined(_WIN32)
  #define trailing_zeros(x) ({ __typeof(x) _r; _BitScanReverse64(&_r, x); _r; })
  #define leading_zeros(x) ({ __typeof(x) _r; _BitScanForward64(&_r, x); _r; })
#else
  #define trailing_zeros(x) ((x) ? (__typeof(x))__builtin_ctzll(x) : (__typeof(x))sizeof(x)*8)
  #define leading_zeros(x) ((x) ? (__typeof(x))__builtin_clzll(x) : (__typeof(x))sizeof(x)*8)
#endif

// Get index of top set bit. If x is 0 return nbits
#define top_set_bit(x) ((x) ? sizeof(x)*8-leading_zeros(x)-1 : sizeof(x)*8)

#define roundup_bits2bytes(bits)   (((bits)+7)/8)
#define roundup_bits2words32(bits) (((bits)+31)/32)
#define roundup_bits2words64(bits) (((bits)+63)/64)

// Round a number up to the nearest number that is a power of two
#define roundup2pow(x) (1UL << (64 - leading_zeros(x)))

#define rot32(x,r) (((x)<<(r)) | ((x)>>(32-(r))))
#define rot64(x,r) (((x)<<(r)) | ((x)>>(64-(r))))

// need to check for length == 0, undefined behaviour if uint64_t >> 64 etc
#define bitmask(nbits,type) ((nbits) ? ~(type)0 >> (sizeof(type)*8-(nbits)): (type)0)
#define bitmask32(nbits) bitmask(nbits,uint32_t)
#define bitmask64(nbits) bitmask(nbits,uint64_t)

// A possibly faster way to combine two words with a mask
//#define bitmask_merge(a,b,abits) ((a & abits) | (b & ~abits))
#define bitmask_merge(a,b,abits) (b ^ ((a ^ b) & abits))

// Swap lowest four bits. A nibble is 4 bits (i.e. half a byte)
#define rev_nibble(x) ((((x)&1)<<3)|(((x)&2)<<1)|(((x)&4)>>1)|(((x)&8)>>3))

//
// Bit array (bitset)
//
// bitsetX_wrd(): get word for a given position
// bitsetX_idx(): get index within word for a given position
#define _VOLPTR(x) ((volatile __typeof(x) *)(&(x)))
#define _VOLVALUE(x) (*_VOLPTR(x))

#define _TYPESHIFT(arr,word,shift) \
        ((__typeof(*(arr)))((__typeof(*(arr)))(word) << (shift)))

#define bitsetX_wrd(wrdbits,pos) ((pos) / (wrdbits))
#define bitsetX_idx(wrdbits,pos) ((pos) % (wrdbits))

#define bitset32_wrd(pos) ((pos) >> 5)
#define bitset32_idx(pos) ((pos) & 31)

#define bitset64_wrd(pos) ((pos) >> 6)
#define bitset64_idx(pos) ((pos) & 63)

//
// Bit functions on arrays
//
#define bitset2_get(arr,wrd,idx)     (((arr)[wrd] >> (idx)) & 0x1)
#define bitset2_set(arr,wrd,idx)     ((arr)[wrd] |=  _TYPESHIFT(arr,1,idx))
#define bitset2_del(arr,wrd,idx)     ((arr)[wrd] &=~ _TYPESHIFT(arr,1,idx))
#define bitset2_tgl(arr,wrd,idx)     ((arr)[wrd] ^=  _TYPESHIFT(arr,1,idx))
#define bitset2_or(arr,wrd,idx,bit)  ((arr)[wrd] |=  _TYPESHIFT(arr,bit,idx))
#define bitset2_xor(arr,wrd,idx,bit) ((arr)[wrd]  = ~((arr)[wrd] ^ (~_TYPESHIFT(arr,bit,idx))))
#define bitset2_and(arr,wrd,idx,bit) ((arr)[wrd] &= (_TYPESHIFT(arr,bit,idx) | ~_TYPESHIFT(arr,1,idx)))
#define bitset2_cpy(arr,wrd,idx,bit) ((arr)[wrd]  = ((arr)[wrd] &~ _TYPESHIFT(arr,1,idx)) | _TYPESHIFT(arr,bit,idx))

//
// Thread safe versions
//
// They return the value of the bit (0 or 1) before it was updated
#define bitset2_get_mt(arr,wrd,idx)     bitset2_get(_VOLPTR(*(arr)),wrd,idx)
#define bitset2_set_mt(arr,wrd,idx)     ((__sync_fetch_and_or (_VOLPTR((arr)[wrd]),  _TYPESHIFT(arr,1,idx)) >> (idx))&1)
#define bitset2_del_mt(arr,wrd,idx)     ((__sync_fetch_and_and(_VOLPTR((arr)[wrd]), ~_TYPESHIFT(arr,1,idx)) >> (idx))&1)
#define bitset2_tgl_mt(arr,wrd,idx)     ((__sync_fetch_and_xor(_VOLPTR((arr)[wrd]),  _TYPESHIFT(arr,1,idx)) >> (idx))&1)
#define bitset2_or_mt(arr,wrd,idx,bit)  ((__sync_fetch_and_or (_VOLPTR((arr)[wrd]),  _TYPESHIFT(arr,bit,idx)) >> (idx))&1)
#define bitset2_xor_mt(arr,wrd,idx,bit) ((__sync_fetch_and_xor(_VOLPTR((arr)[wrd]),  _TYPESHIFT(arr,bit,idx)) >> (idx))&1)
#define bitset2_and_mt(arr,wrd,idx,bit) ((__sync_fetch_and_and(_VOLPTR((arr)[wrd]), (_TYPESHIFT(arr,bit,idx) | ~_TYPESHIFT(arr,1,idx))) >> (idx))&1)
#define bitset2_cpy_mt(arr,wrd,idx,bit) ((bit) ? bitset2_set_mt(arr,wrd,idx) : bitset2_del_mt(arr,wrd,idx))

//
// Auto detect size of type from pointer
//
#define bitset_wrd(arr,pos) bitsetX_wrd(sizeof(*(arr))*8,pos)
#define bitset_idx(arr,pos) bitsetX_idx(sizeof(*(arr))*8,pos)
#define bitset_op(func,arr,pos)      func(arr, bitset_wrd(arr,pos), bitset_idx(arr,pos))
#define bitset_op2(func,arr,pos,bit) func(arr, bitset_wrd(arr,pos), bitset_idx(arr,pos), bit)

// Auto-detect type size: bit functions
#define bitset_get(arr,pos)     bitset_op(bitset2_get, arr, pos)
#define bitset_set(arr,pos)     bitset_op(bitset2_set, arr, pos)
#define bitset_del(arr,pos)     bitset_op(bitset2_del, arr, pos)
#define bitset_tgl(arr,pos)     bitset_op(bitset2_tgl, arr, pos)
#define bitset_or(arr,pos,bit)  bitset_op2(bitset2_or, arr, pos, bit)
#define bitset_xor(arr,pos,bit) bitset_op2(bitset2_xor, arr, pos, bit)
#define bitset_and(arr,pos,bit) bitset_op2(bitset2_and, arr, pos, bit)
#define bitset_cpy(arr,pos,bit) bitset_op2(bitset2_cpy, arr, pos, bit)

// Auto-detect type size: thread safe bit functions
// They return the value of the bit (0 or 1) before it was updated
#define bitset_get_mt(arr,pos)     bitset_op(bitset2_get_mt,  arr, pos)
#define bitset_set_mt(arr,pos)     bitset_op(bitset2_set_mt,  arr, pos)
#define bitset_del_mt(arr,pos)     bitset_op(bitset2_del_mt,  arr, pos)
#define bitset_tgl_mt(arr,pos)     bitset_op(bitset2_tgl_mt,  arr, pos)
#define bitset_or_mt(arr,pos,bit)  bitset_op2(bitset2_or_mt,  arr, pos, bit)
#define bitset_xor_mt(arr,pos,bit) bitset_op2(bitset2_xor_mt, arr, pos, bit)
#define bitset_and_mt(arr,pos,bit) bitset_op2(bitset2_and_mt, arr, pos, bit)
#define bitset_cpy_mt(arr,pos,bit) bitset_op2(bitset2_cpy_mt, arr, pos, bit)

// Clearing a word does not return a meaningful value
#define bitset_clear_word(arr,pos) ((arr)[bitset_wrd(arr,pos)] = 0)
#define bitset_clear_word_mt(arr,pos) (_VOLVALUE((arr)[bitset_wrd(arr,pos)]) = 0)

//
// Compact bit array of spin locks
// These are most effecient when arr is of type: volatile char*
//
// Acquire a lock
#define bitlock_acquire_block(arr,pos,wait,abandon) do {                       \
  size_t _w = bitset_wrd(arr,pos);                                             \
  __typeof(*(arr)) _o, _n, _b = _TYPESHIFT(arr, 1, bitset_idx(arr,pos));       \
  do {                                                                         \
    while((_o = _VOLVALUE((arr)[_w])) & _b) { wait }                           \
    abandon                                                                    \
    _n = _o | _b;                                                              \
  } while(!__sync_bool_compare_and_swap(_VOLPTR((arr)[_w]), _o, _n));          \
  __sync_synchronize(); /* Must not move commands to before acquiring lock */  \
} while(0)

// Undefined behaviour if you do not already hold the lock
#define bitlock_release(arr,pos) do {                                          \
  size_t _w = bitset_wrd(arr,pos);                                             \
  __typeof(*(arr)) _mask = ~_TYPESHIFT(arr, 1, bitset_idx(arr,pos));           \
  __sync_synchronize(); /* Must get the lock before releasing it */            \
  __sync_and_and_fetch(_VOLPTR((arr)[_w]), _mask);                             \
} while(0)

#define bitlock_acquire(arr,pos) bitlock_acquire_block(arr,pos,{},{})

// calls yield if cannot acquire the lock
#define bitlock_yield_acquire(arr,pos) bitlock_acquire_block(arr,pos,sched_yield();,{})

// Block until we get the lock or someone else does
// sets the memory pointed to by retptr to 1 if we got the lock, 0 otherwise
#define bitlock_try_acquire(arr,pos,retptr) do {                               \
  *(retptr) = 1; /* default to success, set to zero if locked */               \
  bitlock_acquire_block(arr,pos,{*(retptr)=0;break;},if(!*(retptr)){break;});  \
} while(0)

/*
 * Byteswapping
 */

/* clang uses these to check for features */
#ifndef __has_feature
#define __has_feature(x) 0
#endif

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

/* GCC versions < 4.3 do not have __builtin_bswapX() */
#if ( defined(__clang__) && !__has_builtin(__builtin_bswap64) ) ||             \
    ( !defined(__clang__) && defined(__GNUC__) && defined(__GNUC_MINOR__) &&   \
      ( (__GNUC__ < 4)  || (__GNUC__ == 4 && __GNUC_MINOR__ < 3)) )
  #define byteswap64(x) ( (((uint64_t)(x) << 56))                       | \
                          (((uint64_t)(x) << 40) & 0xff000000000000ULL) | \
                          (((uint64_t)(x) << 24) & 0xff0000000000ULL)   | \
                          (((uint64_t)(x) <<  8) & 0xff00000000ULL)     | \
                          (((uint64_t)(x) >>  8) & 0xff000000ULL)       | \
                          (((uint64_t)(x) >> 24) & 0xff0000ULL)         | \
                          (((uint64_t)(x) >> 40) & 0xff00ULL)           | \
                          (((uint64_t)(x) >> 56)) )

  #define byteswap32(x) ( (((uint32_t)(x) << 24))                       | \
                          (((uint32_t)(x) <<  8) & 0xff0000U)           | \
                          (((uint32_t)(x) >>  8) & 0xff00U)             | \
                          (((uint32_t)(x) >> 24)) )

  /* uint16_t type might be bigger than 2 bytes, so need to mask */
  #define byteswap16(x) ( (((uint16_t)(x) & 0xff) << 8) | \
                          (((uint16_t)(x) >> 8) & 0xff) )
#else
  #define byteswap64(x) __builtin_bswap64(x)
  #define byteswap32(x) __builtin_bswap64(x)
  #define byteswap16(x) __builtin_bswap64(x)
#endif

#endif /* BITLOCK_H_ */
