#ifndef _KALLOC_H_
#define _KALLOC_H_

#include <stddef.h> /* for size_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	size_t capacity, available, n_blocks, n_cores, largest;
} ri_km_stat_t;

void *ri_kmalloc(void *km, size_t size);
void *ri_krealloc(void *km, void *ptr, size_t size);
void *ri_kcalloc(void *km, size_t count, size_t size);
void ri_kfree(void *km, void *ptr);

void *ri_km_init(void);
void *ri_km_init2(void *km_par, size_t min_core_size);
void ri_km_destroy(void *km);
void ri_km_stat(const void *_km, ri_km_stat_t *s);

#ifdef __cplusplus
}
#endif

#define KMALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))ri_kmalloc((km), (len) * sizeof(*(ptr))))
#define RI_KCALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))ri_kcalloc((km), (len), sizeof(*(ptr))))
#define RI_KREALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))ri_krealloc((km), (ptr), (len) * sizeof(*(ptr))))

#define RI_KEXPAND(km, a, m) do { \
		(m) = (m) >= 4? (m) + ((m)>>1) : 16; \
		RI_KREALLOC((km), (a), (m)); \
	} while (0)

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

#define RI_KALLOC_POOL_INIT2(SCOPE, name, kmptype_t) \
	typedef struct { \
		size_t cnt, n, max; \
		kmptype_t **buf; \
		void *km; \
	} kmp_##name##_t; \
	SCOPE kmp_##name##_t *kmp_init_##name(void *km) { \
		kmp_##name##_t *mp; \
		RI_KCALLOC(km, mp, 1); \
		mp->km = km; \
		return mp; \
	} \
	SCOPE void kmp_destroy_##name(kmp_##name##_t *mp) { \
		size_t k; \
		for (k = 0; k < mp->n; ++k) ri_kfree(mp->km, mp->buf[k]); \
		ri_kfree(mp->km, mp->buf); ri_kfree(mp->km, mp); \
	} \
	SCOPE kmptype_t *kmp_alloc_##name(kmp_##name##_t *mp) { \
		++mp->cnt; \
		if (mp->n == 0) return (kmptype_t*)ri_kcalloc(mp->km, 1, sizeof(kmptype_t)); \
		return mp->buf[--mp->n]; \
	} \
	SCOPE void kmp_free_##name(kmp_##name##_t *mp, kmptype_t *p) { \
		--mp->cnt; \
		if (mp->n == mp->max) RI_KEXPAND(mp->km, mp->buf, mp->max); \
		mp->buf[mp->n++] = p; \
	}

#define RI_KALLOC_POOL_INIT(name, kmptype_t) \
	RI_KALLOC_POOL_INIT2(static inline klib_unused, name, kmptype_t)

#endif
