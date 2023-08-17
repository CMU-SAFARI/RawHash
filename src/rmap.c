#include "rmap.h"
#include <assert.h>
#include "kthread.h"
#include "rh_kvec.h"
#include "rutils.h"
#include "rsketch.h"
#include "revent.h"
#include "rseed.h"
#include "sequence_until.h"
#include "chain.h"

#ifdef PROFILERH
double ri_filereadtime = 0.0;
double ri_signaltime = 0.0;
double ri_sketchtime = 0.0;
double ri_seedtime = 0.0;
double ri_chaintime = 0.0;
double ri_maptime = 0.0;
double ri_maptime_multithread = 0.0;
#endif

static ri_tbuf_t *ri_tbuf_init(void)
{
	ri_tbuf_t *b;
	b = (ri_tbuf_t*)calloc(1, sizeof(ri_tbuf_t));
	b->km = ri_km_init();
	return b;
}

static void ri_tbuf_destroy(ri_tbuf_t *b)
{
	if (b == 0) return;
	ri_km_destroy(b->km);
	free(b);
}

/**
 * Find seed matches between a chunk of a raw signal and a reference genome
 *
 * @param km     		thread-local memory pool; using NULL falls back to malloc()
 * @return seed_hits    seed hits
 *               		a[i].x = strand | ref_id | ref_pos
 *               		a[i].y = flags  | q_span | q_pos
 */
static mm128_t *collect_seed_hits(void *km,
								//   const ri_mapopt_t *opt,
								  int max_occ,
								  int max_max_occ,
                              	  int dist,
								  const ri_idx_t *ri,
								//   const char *qname,
								  ri_reg1_t* reg,
								  const mm128_v *riv,
								  int qlen,
								  int64_t *n_seed_pos,
								  int *rep_len)
								//   int *n_seed_mini,
								//   uint64_t **seed_mini)
{
	int i,j, n_seed_m;
	ri_seed_t *seed_hits0;
	mm128_t *seed_hits;
	uint32_t mask_pos = (1ULL<<31)-1;
	uint64_t mask_id_shift = (((1ULL<<RI_ID_SHIFT) - 1)<<RI_ID_SHIFT)>>RI_POS_SHIFT;

	seed_hits0 = ri_collect_matches(km, &n_seed_m, qlen, max_occ, max_max_occ, dist, ri, riv, n_seed_pos, rep_len);
	seed_hits = (mm128_t*)ri_kmalloc(km, (*n_seed_pos + reg->n_prev_anchors) * sizeof(mm128_t));
	for (i = 0, *n_seed_pos = 0; i < n_seed_m; ++i) {
		ri_seed_t *s_match = &seed_hits0[i];
		const uint64_t *hits = s_match->cr;
		uint32_t k;
		for (k = 0; k < s_match->n; ++k) {
			uint32_t is_self = 0, ref_pos = (uint32_t)(hits[k] >> RI_POS_SHIFT)&mask_pos;
			mm128_t *p;
			// if (skip_seed(opt->flag, hits[k], s_match, qname, qlen, ri, &is_self)) continue;
			p = &seed_hits[(*n_seed_pos)++];

			if (!(hits[k]&1)) {
				p->x = hits[k]&mask_id_shift | ref_pos; //leftmost bit is 0 to indicate forward strand match
			}else{
				uint32_t len = ri->seq[hits[k]>>RI_ID_SHIFT].len;
				// p->x = 1ULL<<63 | hits[k]&mask_id_shift | (len - (ref_pos + 1 - s_match->q_span) - 1); // reverse strand
				// p->x = 1ULL<<63 | hits[k]&mask_id_shift | (len - ref_pos + 1); // reverse strand
				p->x = 1ULL<<63 | hits[k]&mask_id_shift | ref_pos; // reverse strand
			}

			// // forward strand because s_match->q_pos&1 is always 0 as we use single strand (0) for reads
			// if ((hits[k]&1) == (s_match->q_pos&1)) {
			// 	p->x = hits[k]&mask_id_shift | ref_pos; //leftmost bit is 0 to indicate forward strand match
			// 	p->y = (uint64_t)s_match->q_span << RI_ID_SHIFT | ((s_match->q_pos>>RI_POS_SHIFT)+reg->offset);
			// }
			// // else if (!(opt->flag & RI_F_QSTRAND)) { // reverse strand and not in the query-strand mode
			// // 	p->x = 1ULL<<63 | hits[k]&mask_id_shift | ref_pos;
			// // 	p->y = (uint64_t)s_match->q_span << RI_ID_SHIFT | (qlen - ((s_match->q_pos>>1) + 1 - s_match->q_span) - 1);
			// // }
			// else { // reverse strand; query-strand
			// 	// int32_t len = ri->seq[hits[k]>>RI_ID_SHIFT].len;
			// 	// // coordinate only accurate for non-HPC seeds
				// p->x = 1ULL<<63 | hits[k]&mask_id_shift | (len - (ref_pos + 1 - s_match->q_span) - 1);
			// 	p->x = 1ULL<<63 | hits[k]&mask_id_shift | ref_pos; // TODO: see what happens when you enable the above line
			// 	p->y = (uint64_t)s_match->q_span << RI_ID_SHIFT | ((s_match->q_pos>>RI_POS_SHIFT)+reg->offset);
			// }

			p->y = (uint64_t)s_match->q_span << RI_ID_SHIFT | ((s_match->q_pos>>RI_POS_SHIFT)+reg->offset);
			p->y |= (uint64_t)s_match->seg_id << RI_SEED_SEG_SHIFT;
			if (s_match->is_tandem) p->y |= RI_SEED_TANDEM;
			if (is_self) p->y |= RI_SEED_SELF;
		}
	}
	ri_kfree(km, seed_hits0);

	//memcpy reg->prev_anchors to seed_hits starting from index *n_seed_pos
	if(reg->n_prev_anchors > 0){
		memcpy(seed_hits + *n_seed_pos, reg->prev_anchors, reg->n_prev_anchors * sizeof(mm128_t));
		*n_seed_pos += reg->n_prev_anchors;
		if(reg->prev_anchors) ri_kfree(km, reg->prev_anchors); reg->prev_anchors = NULL; reg->n_prev_anchors = 0;
	}

	radix_sort_128x(seed_hits, seed_hits + *n_seed_pos);
	return seed_hits;
}

//returns n_regs // currently we report one mapping
void ri_map_frag(const ri_idx_t *ri,
				const uint32_t s_len,
				const float *sig,
				ri_reg1_t* reg,
				ri_tbuf_t *b,
				const ri_mapopt_t *opt,
				const char *qname)
{	
	uint32_t n_events = 0;

	#ifdef PROFILERH
	double signal_t = ri_realtime();
	#endif
	float* events = detect_events(b->km, s_len, sig, opt, &n_events);
	#ifdef PROFILERH
	ri_signaltime += ri_realtime() - signal_t;
	#endif

	if(n_events < opt->min_events) {
		if(events)ri_kfree(b->km, events);
		return;
	}

	// double new_chain_st = ri_realtime();
	// gen_chains(b->km, ri, events, n_events, reg, opt);

	#ifdef PROFILERH
	double sketch_t = ri_realtime();
	#endif

	//Sketching
	mm128_v riv = {0,0,0};
	ri_sketch(b->km, events, 0, 0, n_events, ri->w, ri->e, ri->n, ri->q, ri->lq, ri->k, &riv);
	if(events)ri_kfree(b->km, events);
	if (opt->q_occ_frac > 0.0f) ri_seed_mz_flt(b->km, &riv, opt->mid_occ, opt->q_occ_frac);

	#ifdef PROFILERH
	ri_sketchtime += ri_realtime() - sketch_t;
	#endif

	#ifdef PROFILERH
	double seed_t = ri_realtime();
	#endif

	// int n_seed_m, n_seed_mini;
	int rep_len;
	int64_t n_seed_pos;
	// ri_seed_t *seed_matches;
	mm128_t *seed_hits;
	uint64_t *u;
	uint32_t hash;
	// uint64_t *seed_mini;

	//Seeding
	seed_hits = collect_seed_hits(b->km, opt->mid_occ, opt->max_max_occ, opt->occ_dist, ri, reg, &riv, n_events, &n_seed_pos, &rep_len);
	ri_kfree(b->km, riv.a);
	// ri_kfree(b->km, seed_mini);

	#ifdef PROFILERH
	ri_seedtime += ri_realtime() - seed_t;
	#endif

	#ifdef PROFILERH
	double chain_t = ri_realtime();
	#endif

	//Chaining (DP or RMQ)
	float chn_pen_gap = opt->chain_gap_scale * 0.01 * ri->e, chn_pen_skip = opt->chain_skip_scale * 0.01 * ri->e;
	if(!(opt->flag & RI_M_RMQ))
		seed_hits = mg_lchain_dp(opt->max_target_gap_length, opt->max_query_gap_length, opt->bw, opt->max_num_skips, opt->max_chain_iter,
							 	opt->min_num_anchors,opt->min_chaining_score,chn_pen_gap,chn_pen_skip,&n_seed_pos,seed_hits, 
							 	&(reg->prev_anchors), &(reg->n_cregs), &u, b->km);
	else
		seed_hits = mg_lchain_rmq(opt->max_target_gap_length, opt->rmq_inner_dist, opt->bw, opt->max_num_skips, opt->rmq_size_cap,
							 	  opt->min_num_anchors,opt->min_chaining_score,chn_pen_gap,chn_pen_skip,&n_seed_pos,seed_hits, 
							 	  &(reg->prev_anchors), &(reg->n_cregs), &u, b->km);

	reg->n_prev_anchors = 0;
	if(n_seed_pos>0)reg->n_prev_anchors = n_seed_pos;
	else if(reg->prev_anchors) ri_kfree(b->km, reg->prev_anchors);

	hash = 0;
	hash ^= __ac_Wang_hash(reg->offset+n_events) + __ac_Wang_hash(11);
	hash  = __ac_Wang_hash(hash);

	//Find primary chains
	reg->creg = mm_gen_regs(b->km, hash, reg->offset+n_events, reg->n_cregs, u, seed_hits);
	mm_set_parent(b->km, opt->mask_level, opt->mask_len, reg->n_cregs, reg->creg, opt->a * 2 + opt->b, opt->flag&RI_M_HARD_MLEVEL, opt->alt_drop);
	mm_select_sub(b->km, opt->pri_ratio, ri->k*2, opt->best_n, 1, opt->max_target_gap_length * 0.8, &(reg->n_cregs), reg->creg);
	//Set MAPQ
	mm_set_mapq(b->km, reg->n_cregs, reg->creg, opt->min_chaining_score, opt->a, rep_len);

	ri_kfree(b->km, seed_hits);
	if(u)ri_kfree(b->km, u);

	#ifdef PROFILERH
	ri_chaintime += ri_realtime() - chain_t;
	#endif

	reg->offset += n_events;
}

static void map_worker_for(void *_data,
						   long i,
						   int tid) // kt_for() callback
{
    step_mt *s = (step_mt*)_data; //s->sig and s->n_sig (signals read in this step and num of them)
	const ri_mapopt_t *opt = s->p->opt;
	ri_tbuf_t* b = s->buf[tid];
	ri_reg1_t* reg0 = s->reg[i];
	reg0->prev_anchors = NULL;
	reg0->n_prev_anchors = 0;
	reg0->creg = NULL;
	ri_sig_t* sig = s->sig[i];

	uint32_t l_chunk = opt->chunk_size;
	uint32_t max_chunk =  opt->max_num_chunk;
	uint32_t qlen = sig->l_sig;
	uint32_t s_qs, s_qe = l_chunk;

	uint32_t c_count = 0;
	int is_mapped = 0;
	
	double t = ri_realtime();
	for (s_qs = c_count = 0; s_qs < qlen && c_count < max_chunk; s_qs += l_chunk, ++c_count) {
		s_qe = s_qs + l_chunk;
		if(s_qe > qlen) s_qe = qlen;

		if(reg0->creg)free(reg0->creg);
		reg0->creg = NULL;
		reg0->n_cregs = 0;

		ri_map_frag(s->p->ri, (const uint32_t)s_qe-s_qs, (const float*)&(sig->sig[s_qs]), reg0, b, opt, sig->name);
		
		//TODO: There should be much better automatic ways to make decisions here.
		if (reg0->n_cregs >= 2) {
			if (reg0->creg[0].mapq > opt->min_bestmapq && reg0->creg[1].mapq == 0) {is_mapped = 1; break;}
			
			float mean_chain_score = 0;
			float mean_mapq = 0;
			int non_zero_mapq = 0;
			for (uint32_t c_ind = 0; c_ind < reg0->n_cregs; ++c_ind){
				mean_chain_score += reg0->creg[c_ind].score;
				mean_mapq += reg0->creg[c_ind].mapq;
				if(reg0->creg[c_ind].mapq > 0) ++non_zero_mapq;
			}

			mean_chain_score /= reg0->n_cregs;
			mean_mapq /= reg0->n_cregs;

			if(non_zero_mapq < 3 && reg0->creg[1].mapq > 0 && 
			   reg0->creg[0].mapq / reg0->creg[1].mapq >= opt->min_bestmapq_ratio) {is_mapped = 1; break;}
			if(non_zero_mapq < 3 && reg0->creg[1].score > 0 &&
			   reg0->creg[0].score / reg0->creg[1].score >= opt->min_bestchain_ratio) {is_mapped = 1; break;}

			if (non_zero_mapq > 2 && mean_mapq > 0 && reg0->creg[0].mapq >= opt->min_meanmapq_ratio * mean_mapq) {is_mapped = 1; break;}
			if (reg0->creg[0].mapq > 0 && mean_chain_score > 0 && reg0->creg[0].score >= opt->min_meanchain_ratio * mean_chain_score) {is_mapped = 1; break;}

			//TODO increase score for those with query ending towards the end (i.e., meaning they benefited from more recent sequencing)
		} else if (reg0->n_cregs == 1 && reg0->creg[0].mapq >= opt->min_mapq) {is_mapped = 1; break;}
	} double mapping_time = ri_realtime() - t;

	#ifdef PROFILERH
	ri_maptime += mapping_time;
	#endif

	if (c_count > 0 && (s_qs >= qlen || c_count == max_chunk)) --c_count;

	float read_position_scale = ((float)(c_count+1)*l_chunk/reg0->offset) / opt->sample_per_base;
	// Save results in vector and output PAF
	
	mm_reg1_t* chains = reg0->creg;

	if(!chains) {reg0->n_cregs = 0;}

	// uint32_t n_anchors0 = (reg0->n_cregs)?chains[0].n_anchors:0;
	uint32_t n_anchors0 = (reg0->n_cregs)?chains[0].cnt:0;
	float mean_chain_score = 0;
	for (uint32_t c_ind = 0; c_ind < reg0->n_cregs; ++c_ind)
		mean_chain_score += chains[c_ind].score;

	mean_chain_score /= reg0->n_cregs;
	if (is_mapped || (reg0->creg && reg0->creg[0].mapq > opt->min_mapq)) {
							
		// float anchor_ref_gap_avg_length = 0;
		// float anchor_read_gap_avg_length = 0;
		// for (size_t ai = 0; ai < n_anchors0; ++ai) {
		// 	if (ai < n_anchors0 - 1) {
		// 		// anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
		// 		// anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
		// 		anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
		// 		anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
		// 	}
		// }
		// anchor_ref_gap_avg_length /= n_anchors0;
		// anchor_read_gap_avg_length /= n_anchors0;
		// std::string tags;
		char *tags = (char *)malloc(1024 * sizeof(char));
		tags[0] = '\0'; // make it an empty string
		char buffer[256]; // temporary buffer
		sprintf(buffer, "mt:f:%.6f", mapping_time * 1000); strcat(tags, buffer);
		sprintf(buffer, "\tci:i:%d", c_count + 1); strcat(tags, buffer);
		sprintf(buffer, "\tsl:i:%d", qlen); strcat(tags, buffer);
		sprintf(buffer, "\tcm:i:%d", n_anchors0); strcat(tags, buffer);
		sprintf(buffer, "\tnc:i:%d", reg0->n_cregs); strcat(tags, buffer);
		sprintf(buffer, "\ts1:i:%d", chains[0].score); strcat(tags, buffer);
		sprintf(buffer, "\ts2:i:%d", reg0->n_cregs > 1 ? chains[1].score : 0); strcat(tags, buffer);
		sprintf(buffer, "\tsm:f:%.2f", mean_chain_score); strcat(tags, buffer);
		
		// tags.append("mt:f:" + std::to_string(mapping_time * 1000));
		// tags.append("\tci:i:" + std::to_string(c_count + 1));
		// tags.append("\tsl:i:" + std::to_string(qlen));
		// tags.append("\tcm:i:" + std::to_string(n_anchors0));
		// tags.append("\tnc:i:" + std::to_string(reg0->n_cregs));
		// tags.append("\ts1:i:" + std::to_string(chains[0].score));
		// tags.append("\ts2:f:" + std::to_string(reg0->n_cregs > 1 ? chains[1].score : 0));
		// tags.append("\tsm:f:" + std::to_string(mean_chain_score));
		// tags.append("\tat:f:" + std::to_string((float)anchor_ref_gap_avg_length));
		// tags.append("\taq:f:" + std::to_string((float)anchor_read_gap_avg_length));

		reg0->read_id = sig->rid;
		reg0->ref_id = chains[0].rid;
		reg0->read_name = sig->name;
		reg0->read_length = (uint32_t)(read_position_scale*chains[0].qe);
		reg0->read_start_position = (uint32_t)(read_position_scale*chains[0].qs);
		reg0->read_end_position = (uint32_t)(read_position_scale*chains[0].qe);
		reg0->fragment_start_position = chains[0].rev?(uint32_t)(s->p->ri->seq[chains[0].rid].len+1-chains[0].re):chains[0].rs;
		reg0->fragment_length = (uint32_t)(chains[0].re - chains[0].rs + 1);
		reg0->mapq = chains[0].mapq;
		reg0->rev = (chains[0].rev == 1)?1:0;
		reg0->mapped = 1; 
		// reg0->tags = strdup(tags.c_str());
		reg0->tags = tags;
	} else {
		// std::string tags;
		char *tags = (char *)malloc(1024 * sizeof(char));
		tags[0] = '\0'; // make it an empty string
		char buffer[256]; // temporary buffer

		sprintf(buffer, "mt:f:%.6f", mapping_time * 1000); strcat(tags, buffer);
		sprintf(buffer, "\tci:i:%d", c_count + 1); strcat(tags, buffer);
		sprintf(buffer, "\tsl:i:%d", qlen); strcat(tags, buffer);
		// tags.append("mt:f:" + std::to_string(mapping_time * 1000));
		// tags.append("\tci:i:" + std::to_string(c_count + 1));
		// tags.append("\tsl:i:" + std::to_string(qlen));
		if (reg0->n_cregs >= 1) {
			// float anchor_ref_gap_avg_length = 0;
			// float anchor_read_gap_avg_length = 0;
			// for (size_t ai = 0; ai < n_anchors0; ++ai) {
			// 	if (ai < n_anchors0 - 1) {
			// 		anchor_ref_gap_avg_length += chains[0].anchors[ai].target_position - chains[0].anchors[ai + 1].target_position;
			// 		anchor_read_gap_avg_length += chains[0].anchors[ai].query_position - chains[0].anchors[ai + 1].query_position;
			// 	}
			// }
			// if(n_anchors0)anchor_ref_gap_avg_length /= n_anchors0;
			// if(n_anchors0)anchor_read_gap_avg_length /= n_anchors0;

			sprintf(buffer, "\tcm:i:%d", n_anchors0); strcat(tags, buffer);
			sprintf(buffer, "\tnc:i:%d", reg0->n_cregs); strcat(tags, buffer);
			sprintf(buffer, "\ts1:i:%d", chains[0].score); strcat(tags, buffer);
			sprintf(buffer, "\ts2:i:%d", reg0->n_cregs > 1 ? chains[1].score : 0); strcat(tags, buffer);
			sprintf(buffer, "\tsm:f:%.2f", mean_chain_score); strcat(tags, buffer);
			// tags.append("\tcm:i:" + std::to_string(n_anchors0));
			// tags.append("\tnc:i:" + std::to_string(reg0->n_cregs));
			// tags.append("\ts1:i:" + std::to_string(chains[0].score));
			// tags.append("\ts2:i:" + std::to_string(reg0->n_cregs > 1 ? chains[1].score: 0));
			// tags.append("\tsm:f:" + std::to_string(mean_chain_score));
			// tags.append("\tat:f:" + std::to_string((float)anchor_ref_gap_avg_length));
			// tags.append("\taq:f:" + std::to_string((float)anchor_read_gap_avg_length));
		}else {
			sprintf(buffer, "\tcm:i:0"); strcat(tags, buffer);
			sprintf(buffer, "\tnc:i:0"); strcat(tags, buffer);
			sprintf(buffer, "\ts1:i:0"); strcat(tags, buffer);
			sprintf(buffer, "\ts2:i:0"); strcat(tags, buffer);
			sprintf(buffer, "\tsm:f:0"); strcat(tags, buffer);
			// tags.append("\tcm:i:0");
			// tags.append("\tnc:i:0");
			// tags.append("\ts1:i:0");
			// tags.append("\ts2:i:0");
			// tags.append("\tsm:f:0");
			// tags.append("\tat:f:0");
			// tags.append("\taq:f:0");
		}

		reg0->read_id = sig->rid;
		reg0->ref_id = 0;
		reg0->read_name = sig->name;
		reg0->read_length = (uint32_t)(read_position_scale * reg0->offset);
		reg0->read_start_position = 0;
		reg0->read_end_position = 0;
		reg0->fragment_start_position = 0;
		reg0->fragment_length = 0;
		reg0->mapq = 0;
		reg0->rev = 0;
		reg0->mapped = 0;
		// reg0->tags = strdup(tags.c_str());
		reg0->tags = tags;
	}

	if(reg0->prev_anchors) {ri_kfree(b->km, reg0->prev_anchors); reg0->prev_anchors = NULL; reg0->n_prev_anchors = 0;}
	if(reg0->creg)free(reg0->creg);
	reg0->creg = NULL;
	reg0->n_cregs = 0;

	if (b->km) {
		ri_km_stat_t kmst;
		ri_km_stat(b->km, &kmst);
		// assert(kmst.n_blocks == kmst.n_cores);
		ri_km_destroy(b->km);
		b->km = ri_km_init();
	}
}

typedef struct { size_t n, m; ri_sig_t **a; } rhsig_v;

ri_sig_t** ri_sig_read_frag(pipeline_mt *pl,
							int64_t chunk_size,
							int *n_)
{
	int64_t size = 0;
	rhsig_v rsigv = {0,0,0};
	rh_kv_resize(ri_sig_t*, 0, rsigv, 4000);
	*n_ = 0;
	if (pl->n_fp < 1) return 0;
	while (pl->fp) {
		//Reading data in bulk if the buffer is emptied
		while(pl->fp && pl->fp->cur_read == pl->fp->num_read){
			ri_sig_close(pl->fp);
			if(pl->cur_f < pl->n_f){
				if((pl->fp = open_sig(pl->f[pl->cur_f++])) == 0) break;
			}else if(pl->cur_fp < pl->n_fp){
				if(pl->f){
					for(int i = 0; i < pl->n_f; ++i) if(pl->f[i])free(pl->f[i]);
					free(pl->f);
				}
				pl->n_f = 0; pl->cur_f = 0;

				ri_char_v fnames = {0,0,0};
				rh_kv_resize(char*, 0, fnames, 256);
				find_sfiles(pl->fn[pl->cur_fp++], &fnames);
				pl->f =  fnames.a;
				if(!fnames.n || ((pl->fp = open_sig(pl->f[pl->cur_f++])) == 0)) break;
				pl->n_f = fnames.n;
				// ++n_read;
			}else {pl->fp = 0; break;}
		}

		if(!pl->fp || pl->fp->cur_read == pl->fp->num_read) break;
		
		ri_sig_t *s = (ri_sig_t*)calloc(1, sizeof(ri_sig_t));
		rh_kv_push(ri_sig_t*, 0, rsigv, s);
		ri_read_sig(pl->fp, s);
		size += s->l_sig;

		if(size >= chunk_size) break;
	}

	ri_sig_t** a = 0;
	if(rsigv.n) a = rsigv.a;
	*n_ = rsigv.n;

	return a;
}

static void *map_worker_pipeline(void *shared,
								int step,
								void *in)
{
	int i, k;
    pipeline_mt *p = (pipeline_mt*)shared;
    if (step == 0) { // step 0: read sequences
		#ifdef PROFILERH
		double file_t = ri_realtime();
		#endif
        step_mt *s;
        s = (step_mt*)calloc(1, sizeof(step_mt));
		s->sig = ri_sig_read_frag(p, p->mini_batch_size, &s->n_sig);
		if (s->n_sig && !p->su_stop) {
			s->p = p;
			for (i = 0; i < s->n_sig; ++i)
				(*(s->sig[i])).rid = p->n_processed++;
			s->buf = (ri_tbuf_t**)calloc(p->n_threads, sizeof(ri_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = ri_tbuf_init();
			s->reg = (ri_reg1_t**)calloc(s->n_sig, sizeof(ri_reg1_t*));
			for(i = 0; i < s->n_sig; ++i)
				s->reg[i] = (ri_reg1_t*)calloc(1, sizeof(ri_reg1_t));
			return s;
		} else if(s){
			if(s->sig) free(s->sig);
			free(s);
		}
		#ifdef PROFILERH
		ri_filereadtime += ri_realtime() - file_t;
		#endif
    } else if (step == 1) { // step 1: detect events
		step_mt *s = (step_mt*)in;
		#ifdef PROFILERH
		double map_multit = ri_realtime();
		#endif
		if(!p->su_stop) kt_for(p->n_threads, map_worker_for, in, s->n_sig);
		#ifdef PROFILERH
		ri_maptime_multithread += ri_realtime() - map_multit;
		#endif

		if(p->opt->flag & RI_M_SEQUENCEUNTIL && !p->su_stop){
			const ri_idx_t *ri = p->ri;
			for (k = 0; k < s->n_sig; ++k) {
				if(s->reg[k] && s->reg[k]->ref_id < ri->n_seq && s->reg[k]->read_name && s->reg[k]->mapped){
					ri_reg1_t* reg0 = s->reg[k];
					p->su_c_estimations[reg0->ref_id] += reg0->fragment_length;
					p->ab_count += reg0->fragment_length;
					p->su_nreads++;
					if(p->su_nreads > p->opt->tmin_reads && !(p->su_nreads%p->opt->ttest_freq)){
						//calculate abundance
						for(uint32_t ce = 0; ce < ri->n_seq; ++ce){
							p->su_estimations[p->su_cur][ce] =  (float)p->su_c_estimations[ce]/p->ab_count;
						}
						
						if(++p->su_cur >= p->opt->tn_samples) p->su_cur = 0;
						if(p->su_nestimations++ >= p->opt->tn_samples){
							if(find_outlier((const float**)p->su_estimations, ri->n_seq, p->opt->tn_samples) <= p->opt->t_threshold){
								//sending the stop signal.
								p->su_stop = k+1;
								fprintf(stderr, "[M::%s] Sequence Until is activated, stopping sequencing after processing %d mapped reads\n", __func__, p->su_nreads);
								break;
							}
						}
					}
				}
			}
		}
		return in;
    } else if (step == 2) { // step 2: output
		// void *km = 0;
        step_mt *s = (step_mt*)in;
		const ri_idx_t *ri = p->ri;
		for (k = 0; k < s->n_sig; ++k) {
			if(s->reg[k]){
				ri_reg1_t* reg0 = s->reg[k];
				char strand = reg0->rev?'-':'+';
				
				if(reg0->ref_id < ri->n_seq && reg0->read_name){
					if(reg0->mapped && (!p->su_stop || k < p->su_stop) ){
						fprintf(stdout, "%s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%u\t%u\t%d\t%s\n", 
										reg0->read_name, reg0->read_length, reg0->read_start_position, reg0->read_end_position, 
										strand, ri->seq[reg0->ref_id].name, ri->seq[reg0->ref_id].len, reg0->fragment_start_position, reg0->fragment_start_position + reg0->fragment_length, reg0->read_end_position-reg0->read_start_position-1, reg0->fragment_length, reg0->mapq, reg0->tags);
					}else{
						fprintf(stdout, "%s\t%u\t*\t*\t*\t*\t*\t*\t*\t*\t*\t%d\t%s\n", reg0->read_name, reg0->read_length, reg0->mapq, reg0->tags);
					}
				}

				if(reg0->tags) free(reg0->tags);
				free(reg0);
			}
		}

		for (i = 0; i < p->n_threads; ++i) ri_tbuf_destroy(s->buf[i]);
		free(s->buf);

		if(s->reg)free(s->reg);

		for(int i = 0; i < s->n_sig; ++i){
			ri_sig_t *curS = s->sig[i];
			free(curS->sig);
			free(curS->name);
			free(curS);
		} if(s->sig)free(s->sig); free(s);
	}
    return 0;
}

int ri_map_file(const ri_idx_t *idx,
				const char *fn,
				const ri_mapopt_t *opt,
				int n_threads)
{
	return ri_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int ri_map_file_frag(const ri_idx_t *idx,
					int n_segs,
					const char **fn,
					const ri_mapopt_t *opt,
					int n_threads)
{
	int pl_threads;
	pipeline_mt pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_mt));
	pl.n_fp = n_segs;
	pl.n_f = 0; pl.cur_f = 0;
	ri_char_v fnames = {0,0,0};
	rh_kv_resize(char*, 0, fnames, 256);
	find_sfiles(fn[0], &fnames);
	pl.f =  fnames.a;
	if(!fnames.n || ((pl.fp = open_sig(pl.f[0])) == 0)) return -1;
	if (pl.fp == 0) return -1;
	pl.fn = fn;
	pl.n_f = fnames.n;
	pl.cur_fp = 1;
	pl.cur_f = 1;
	pl.opt = opt, pl.ri = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl_threads = pl.n_threads == 1?1:2;
	pl.su_stop = 0;

	if(opt->flag & RI_M_SEQUENCEUNTIL){
		pl.su_nreads = 0;
		pl.su_nestimations = 0;
		pl.ab_count = 0;
		pl.su_cur = 0;
		pl.su_estimations = (float**)calloc(opt->tn_samples, sizeof(float*));
		for(uint32_t i = 0; i < opt->tn_samples; ++i) pl.su_estimations[i] = (float*)calloc(idx->n_seq, sizeof(float));

		pl.su_c_estimations = (uint32_t*)calloc(idx->n_seq, sizeof(uint32_t));
	}
	
	kt_pipeline(pl_threads, map_worker_pipeline, &pl, 3);

	if(opt->flag & RI_M_SEQUENCEUNTIL){
		// pl.su_nreads = 0;
		// pl.su_nestimations = 0;
		// pl.su_stop = 0;
		if(pl.su_estimations){
			for(uint32_t i = 0; i < opt->tn_samples; ++i) if(pl.su_estimations[i]) free(pl.su_estimations[i]);
			free(pl.su_estimations);
		}

		if(pl.su_c_estimations)free(pl.su_c_estimations);
	}

	#ifdef PROFILERH
	fprintf(stderr, "\n[M::%s] File read: %.6f sec; Signal-to-event: %.6f sec; Sketching: %.6f sec; Seeding: %.6f sec; Chaining: %.6f sec; Mapping: %.6f sec; Mapping (multi-threaded): %.6f sec\n", __func__, ri_filereadtime, ri_signaltime, ri_sketchtime, ri_seedtime, ri_chaintime, ri_maptime, ri_maptime_multithread);
	#endif

	return 0;
}