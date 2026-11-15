/* Streamed FASTQ filter for the FLT1 isoform cohort scan.
 *
 * Reads FASTQ on stdin and, in one pass:
 *   1. counts reads containing any pre-specified exact-unique marker 50-mer (same semantics as
 *      scan_streamed_fastq_markers.py: a read is counted once per marker);
 *   2. writes every read that shares at least one exact BAIT_K-mer (either strand) with the FLT1
 *      cDNA bait to the candidate FASTQ, for splice-aware genome alignment.
 *
 * Usage: flt1_stream_filter marker-kmers.tsv flt1-bait.fa RUN candidates.fastq counts.json < reads.fastq
 * Build: cc -O3 -o flt1_stream_filter flt1_stream_filter.c
 */
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MARKER_K 50
#define BAIT_K 20
#define MAX_MARKERS 64
#define LINE_MAX_LEN (1 << 16)

typedef unsigned __int128 u128;

static int code_of(int c) {
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return -1;
    }
}

/* ---- 128-bit marker table (open addressing) ---- */
#define MT_BITS 16
#define MT_SIZE (1u << MT_BITS)
static u128 mt_key[MT_SIZE];
static uint64_t mt_mask[MT_SIZE]; /* bitmask of marker ids; 0 = empty */

static uint32_t hash128(u128 x) {
    uint64_t h = (uint64_t)x ^ ((uint64_t)(x >> 64) * 0x9E3779B97F4A7C15ULL);
    h ^= h >> 31;
    h *= 0xBF58476D1CE4E5B9ULL;
    h ^= h >> 29;
    return (uint32_t)h & (MT_SIZE - 1);
}

static void mt_add(u128 key, int marker) {
    uint32_t i = hash128(key);
    while (mt_mask[i] && mt_key[i] != key) i = (i + 1) & (MT_SIZE - 1);
    mt_key[i] = key;
    mt_mask[i] |= 1ULL << marker;
}

static inline uint64_t mt_get(u128 key) {
    uint32_t i = hash128(key);
    while (mt_mask[i]) {
        if (mt_key[i] == key) return mt_mask[i];
        i = (i + 1) & (MT_SIZE - 1);
    }
    return 0;
}

/* ---- 64-bit bait table ---- */
#define BT_BITS 20
#define BT_SIZE (1u << BT_BITS)
static uint64_t bt_key[BT_SIZE]; /* stored as key+1 so 0 = empty */

static uint32_t hash64(uint64_t h) {
    h ^= h >> 33;
    h *= 0xFF51AFD7ED558CCDULL;
    h ^= h >> 33;
    return (uint32_t)h & (BT_SIZE - 1);
}

static void bt_add(uint64_t key) {
    uint32_t i = hash64(key);
    while (bt_key[i] && bt_key[i] != key + 1) i = (i + 1) & (BT_SIZE - 1);
    bt_key[i] = key + 1;
}

static inline int bt_has(uint64_t key) {
    uint32_t i = hash64(key);
    while (bt_key[i]) {
        if (bt_key[i] == key + 1) return 1;
        i = (i + 1) & (BT_SIZE - 1);
    }
    return 0;
}

static void add_bait_sequence(const char *s, size_t n) {
    const uint64_t mask = (1ULL << (2 * BAIT_K)) - 1;
    uint64_t fw = 0, rv = 0;
    int run = 0;
    for (size_t i = 0; i < n; i++) {
        int c = code_of(toupper((unsigned char)s[i]));
        if (c < 0) { run = 0; fw = rv = 0; continue; }
        fw = ((fw << 2) | (uint64_t)c) & mask;
        rv = (rv >> 2) | ((uint64_t)(3 - c) << (2 * (BAIT_K - 1)));
        if (++run >= BAIT_K) { bt_add(fw); bt_add(rv); }
    }
}

static char marker_names[MAX_MARKERS][128];
static int n_markers = 0;

static int marker_id(const char *name) {
    for (int i = 0; i < n_markers; i++)
        if (strcmp(marker_names[i], name) == 0) return i;
    if (n_markers >= MAX_MARKERS) { fprintf(stderr, "too many markers\n"); exit(2); }
    strncpy(marker_names[n_markers], name, 127);
    return n_markers++;
}

static void chomp(char *s) {
    size_t n = strlen(s);
    while (n && (s[n - 1] == '\n' || s[n - 1] == '\r')) s[--n] = 0;
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: %s marker-kmers.tsv flt1-bait.fa RUN candidates.fastq counts.json\n", argv[0]);
        return 2;
    }
    static char line[LINE_MAX_LEN], seq[LINE_MAX_LEN], plus[LINE_MAX_LEN], qual[LINE_MAX_LEN];

    /* markers */
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }
    long n_marker_kmers = 0;
    while (fgets(line, sizeof line, f)) {
        chomp(line);
        char *tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = 0;
        const char *kmer = tab + 1;
        if (strlen(kmer) != MARKER_K) { fprintf(stderr, "bad marker k-mer length\n"); return 2; }
        u128 key = 0;
        for (int i = 0; i < MARKER_K; i++) {
            int c = code_of(kmer[i]);
            if (c < 0) { fprintf(stderr, "bad marker base\n"); return 2; }
            key = (key << 2) | (u128)c;
        }
        mt_add(key, marker_id(line));
        n_marker_kmers++;
    }
    fclose(f);

    /* bait */
    f = fopen(argv[2], "r");
    if (!f) { perror(argv[2]); return 2; }
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '>') continue;
        chomp(line);
        add_bait_sequence(line, strlen(line));
    }
    fclose(f);

    FILE *cand = fopen(argv[4], "w");
    if (!cand) { perror(argv[4]); return 2; }
    static char inbuf[1 << 22];
    setvbuf(stdin, inbuf, _IOFBF, sizeof inbuf);

    const u128 mmask = (((u128)1) << (2 * MARKER_K)) - 1;
    const uint64_t bmask = (1ULL << (2 * BAIT_K)) - 1;
    uint64_t total = 0, invalid = 0, n_cand = 0, total_bases = 0;
    uint64_t counts[MAX_MARKERS] = {0};
    uint64_t len_min = UINT64_MAX, len_max = 0;

    while (fgets(line, sizeof line, stdin)) {
        if (!fgets(seq, sizeof seq, stdin)) { invalid++; break; }
        if (!fgets(plus, sizeof plus, stdin)) { invalid++; break; }
        if (!fgets(qual, sizeof qual, stdin)) { invalid++; break; }
        chomp(seq);
        chomp(qual);
        size_t n = strlen(seq);
        if (n != strlen(qual)) { invalid++; continue; }
        total++;
        total_bases += n;
        if (n < len_min) len_min = n;
        if (n > len_max) len_max = n;

        u128 mk = 0;
        uint64_t bk = 0, hits = 0;
        int run = 0, is_cand = 0;
        for (size_t i = 0; i < n; i++) {
            int c = code_of(toupper((unsigned char)seq[i]));
            if (c < 0) { run = 0; mk = 0; bk = 0; continue; }
            mk = ((mk << 2) | (u128)c) & mmask;
            bk = ((bk << 2) | (uint64_t)c) & bmask;
            run++;
            if (run >= BAIT_K && !is_cand && bt_has(bk)) is_cand = 1;
            if (run >= MARKER_K) hits |= mt_get(mk);
        }
        if (hits)
            for (int m = 0; m < n_markers; m++)
                if (hits & (1ULL << m)) counts[m]++;
        if (is_cand || hits) {
            n_cand++;
            fputs(line, cand);
            fputs(seq, cand);
            fputs("\n+\n", cand);
            fputs(qual, cand);
            fputc('\n', cand);
        }
    }
    if (fclose(cand) != 0) { perror("close candidates"); return 3; }

    FILE *out = fopen(argv[5], "w");
    if (!out) { perror(argv[5]); return 2; }
    fprintf(out, "{\n  \"run\": \"%s\",\n", argv[3]);
    fprintf(out, "  \"method\": \"compiled exact marker 50-mer count plus FLT1 cDNA %d-mer candidate extraction\",\n", BAIT_K);
    fprintf(out, "  \"total_reads\": %llu,\n  \"total_bases\": %llu,\n", (unsigned long long)total, (unsigned long long)total_bases);
    fprintf(out, "  \"read_length_min\": %llu,\n  \"read_length_max\": %llu,\n",
            (unsigned long long)(total ? len_min : 0), (unsigned long long)len_max);
    fprintf(out, "  \"invalid_fastq_records\": %llu,\n  \"candidate_reads\": %llu,\n",
            (unsigned long long)invalid, (unsigned long long)n_cand);
    fprintf(out, "  \"marker_kmers_loaded\": %ld,\n  \"exact_marker_counts\": {\n", n_marker_kmers);
    for (int m = 0; m < n_markers; m++)
        fprintf(out, "    \"%s\": %llu%s\n", marker_names[m], (unsigned long long)counts[m], m + 1 < n_markers ? "," : "");
    fprintf(out, "  }\n}\n");
    fclose(out);
    return 0;
}
