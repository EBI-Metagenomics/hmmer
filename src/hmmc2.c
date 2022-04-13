#include "p7_config.h"

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/socket.h>
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h> /* On FreeBSD, you need netinet/in.h for struct sockaddr_in            */
#endif /* On OpenBSD, netinet/in.h is required for (must precede) arpa/inet.h  \
        */
#include <arpa/inet.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "hmmpgmd.h"

#include "hmmc2_opts.h"

#define SERVER_PORT 51371
#define MAX_READ_LEN 4096

static void sig_int(int signo)
{
    fprintf(stderr, "Exiting due to ctrl-c\n");
    exit(1);
}

static int usage(char *pgm)
{
    fprintf(stderr, "Usage: %s [-i addr] [-p port] [-A] [-S]\n", pgm);
    fprintf(stderr, "    -S      : print sequence scores\n");
    fprintf(stderr, "    -A      : print sequence alignments\n");
    fprintf(stderr,
            "    -i addr : ip address running daemon (default: 127.0.0.1)\n");
    fprintf(stderr, "    -p port : port daemon listens to clients on "
                    "(default: 51371)\n");
    exit(1);
}

static char buffer[MAX_READ_LEN] = {0};
static char seq[1024 * 1024] = {0};
static uint8_t *buf = 0;
static uint32_t buf_offset = 0;
static uint32_t hits_start = 0;
static char serv_ip[64] = "127.0.0.1";
static unsigned short serv_port = SERVER_PORT;
static int const ali = 1;
static int const scores = 1;

struct h3conn
{
    struct sockaddr_in addr;
    int sockfd;
    ESL_GETOPTS *go;
} conn = {0};

static FILE *fp_targets = 0;
static FILE *fp_domains = 0;
static FILE *tblfp = 0;
static FILE *domtblfp = 0;
static FILE *pfamtblfp = 0;

// hmmdb_idx is in [0 to n-1], but the output files index in [1, n]
static char phony_args[] = "--hmmdb 1 --hmmdb_idx 12492 --acc --cut_ga";
static char phony_cmd[] = "X --hmmdb 1 --hmmdb_idx 12492 --acc --cut_ga";

static int read_input(FILE *fp)
{
    seq[0] = 0;
    strcpy(seq, "@");
    strcpy(seq, phony_args);
    strcpy(seq, "\n");

    while (fgets(buffer, MAX_READ_LEN, fp))
    {
        int n = strlen(buffer);
        strcat(seq, buffer);
    }
    strcat(seq, "//");
    return !!feof(fp);
}

static int parse_input_options(char *ptr)
{
    char t;
    char *s = ++ptr;

    /* skip to the end of the line */
    while (*ptr && (*ptr != '\n' && *ptr != '\r'))
        ++ptr;
    t = *ptr;
    *ptr = 0;

    int rc = esl_getopts_Reuse(conn.go);
    assert(rc == eslOK);

    if (esl_opt_ProcessSpoof(conn.go, phony_cmd) != eslOK) return 1;
    assert(esl_opt_VerifyConfig(conn.go) == eslOK);

    /* skip remaining white spaces */
    *ptr = t;
    return 0;
}

static int perform(void)
{
    P7_HMM *hmm = NULL;
    ESL_SQ *sq = NULL;
    ESL_ALPHABET *abc = esl_alphabet_Create(eslAMINO);
    HMMD_SEARCH_STATS *stats = 0;
    HMMD_SEARCH_STATUS sstatus = {0};
    P7_PIPELINE *pli = NULL;
    P7_TOPHITS *th = NULL;
    int size = 0;

    int total = 0;

    /* Send the string to the server */
    uint64_t n = strlen(seq);
    printf("[%d] Sending data %" PRIu64 ": [%.*s]\n", __LINE__, n, (int)n, seq);
    if (writen(conn.sockfd, seq, n) != n)
    {
        fprintf(stderr, "[%s:%d] write (size %" PRIu64 ") error %d - %s\n",
                __FILE__, __LINE__, n, errno, strerror(errno));
        exit(1);
    }

    // Get the status structure back from the server
    buf = malloc(HMMD_SEARCH_STATUS_SERIAL_SIZE);
    buf_offset = 0;
    n = HMMD_SEARCH_STATUS_SERIAL_SIZE;
    if (buf == NULL)
    {
        printf("Unable to allocate memory for search status "
               "structure\n");
        exit(1);
    }

    if ((size = readn(conn.sockfd, buf, n)) == -1)
    {
        fprintf(stderr, "[%s:%d] read error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }

    if (hmmd_search_status_Deserialize(buf, &buf_offset, &sstatus) != eslOK)
    {
        printf("Unable to deserialize search status object \n");
        exit(1);
    }

    if (sstatus.status != eslOK)
    {
        char *ebuf;
        n = sstatus.msg_size;
        total += n;
        ebuf = malloc(n);
        if ((size = readn(conn.sockfd, ebuf, n)) == -1)
        {
            fprintf(stderr, "[%s:%d] read error %d - %s\n", __FILE__, __LINE__,
                    errno, strerror(errno));
            exit(1);
        }
        fprintf(stderr, "ERROR (%d): %s\n", sstatus.status, ebuf);
        free(ebuf);
        goto COMPLETE;
    }

    free(buf);      // clear this out
    buf_offset = 0; // reset to beginning for next serialized object
    n = sstatus.msg_size;

    if ((buf = malloc(n)) == NULL)
    {
        fprintf(stderr, "[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }
    // Grab the serialized search results
    if ((size = readn(conn.sockfd, buf, n)) == -1)
    {
        fprintf(stderr, "[%s:%d] read error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }

    if ((stats = malloc(sizeof(HMMD_SEARCH_STATS))) == NULL)
    {
        fprintf(stderr, "[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }
    stats->hit_offsets =
        NULL; // force allocation of memory for this in _Deserialize
    if (p7_hmmd_search_stats_Deserialize(buf, &buf_offset, stats) != eslOK)
    {
        printf("Unable to deserialize search stats object \n");
        exit(1);
    }

    // Create the structures we'll deserialize the hits into
    pli = p7_pipeline_Create(conn.go, 100, 100, FALSE, p7_SCAN_MODELS);
    pli->nmodels = stats->nmodels;
    pli->nseqs = stats->nseqs;
    pli->n_past_msv = stats->n_past_msv;
    pli->n_past_bias = stats->n_past_bias;
    pli->n_past_vit = stats->n_past_vit;
    pli->n_past_fwd = stats->n_past_fwd;

    pli->Z = stats->Z;
    pli->domZ = stats->domZ;
    pli->Z_setby = stats->Z_setby;
    pli->domZ_setby = stats->domZ_setby;

    th = p7_tophits_Create();

    free(th->unsrt);
    free(th->hit);

    th->N = stats->nhits;
    if ((th->unsrt = malloc(stats->nhits * sizeof(P7_HIT))) == NULL)
    {
        fprintf(stderr, "[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }
    th->nreported = stats->nreported;
    th->nincluded = stats->nincluded;
    th->is_sorted_by_seqidx = FALSE;
    th->is_sorted_by_sortkey = TRUE;

    if ((th->hit = malloc(sizeof(void *) * stats->nhits)) == NULL)
    {
        fprintf(stderr, "[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }
    hits_start = buf_offset;
    // deserialize the hits
    for (int i = 0; i < stats->nhits; ++i)
    {
        // set all internal pointers of the hit to NULL before
        // deserializing into it
        th->unsrt[i].name = NULL;
        th->unsrt[i].acc = NULL;
        th->unsrt[i].desc = NULL;
        th->unsrt[i].dcl = NULL;
        if ((buf_offset - hits_start) != stats->hit_offsets[i])
        {
            printf("Hit offset %d did not match expected.  Found "
                   "%d, expected %" PRIu64 "\n",
                   i, (buf_offset - hits_start), stats->hit_offsets[i]);
        }
        if (p7_hit_Deserialize(buf, &buf_offset, &(th->unsrt[i])) != eslOK)
        {
            printf("Unable to deserialize hit %d\n", i);
            exit(0);
        }
        th->hit[i] = &(th->unsrt[i]);
    }

    /* adjust the reported and included hits */
    // th->is_sorted = FALSE;
    // p7_tophits_Sort(th);

    /* Print the results.  */
    if (scores)
    {
        p7_tophits_Targets(fp_targets, th, pli, 120);
        fprintf(fp_targets, "\n\n");
    }
    if (ali)
    {
        p7_tophits_Domains(fp_domains, th, pli, 120);
        fprintf(fp_domains, "\n\n");
    }
    p7_tophits_TabularTargets(tblfp, "QNAME", "QACC", th, pli, 1);
    p7_tophits_TabularDomains(domtblfp, "QNAME", "QACC", th, pli, 1);
    p7_tophits_TabularXfam(pfamtblfp, "QNAME", "QACC", th, pli);

    p7_pipeline_Destroy(pli);
    p7_tophits_Destroy(th);
    free(buf);

    fprintf(stdout, "//\n");
    fflush(stdout);

    fprintf(stdout, "Total bytes received %" PRId64 "\n", sstatus.msg_size);

COMPLETE:
    if (abc) esl_alphabet_Destroy(abc);
    if (hmm) p7_hmm_Destroy(hmm);
    if (sq) esl_sq_Destroy(sq);

    return 0;
}

enum ipv
{
    IPV4,
    IPV6
};

static void setup_address(enum ipv ipv, char const *ip, uint16_t port)
{
    conn.addr.sin_family = ipv == IPV4 ? AF_INET : AF_INET6;
    conn.addr.sin_port = htons(port);
    if ((inet_pton(conn.addr.sin_family, ip, &conn.addr.sin_addr)) != 1)
        exit(1);
}

void h3conn_open(enum ipv ipv, char const *ip, uint16_t port)
{
    setup_address(ipv, ip, port);

    if ((conn.sockfd = socket(conn.addr.sin_family, SOCK_STREAM, 0)) == -1)
        exit(1);

    if (connect(conn.sockfd, (struct sockaddr *)&conn.addr,
                sizeof(conn.addr)) == -1)
        exit(1);

    conn.go = esl_getopts_Create(searchOpts);
}

void h3conn_close(void)
{
    esl_getopts_Destroy(conn.go);
    if (close(conn.sockfd)) exit(1);
}

int main(void)
{
    fp_targets = fopen("targets.txt", "wb");
    fp_domains = fopen("domains.txt", "wb");
    tblfp = fopen("tbl.txt", "w");
    domtblfp = fopen("domtbl.txt", "w");
    pfamtblfp = fopen("pfamtbl.txt", "w");

    P7_PIPELINE *pli = NULL;
    P7_TOPHITS *th = NULL;

    HMMD_SEARCH_STATS *stats;
    HMMD_SEARCH_STATUS sstatus;
    char *ptr;

    buf = NULL;
    buf_offset = 0;

    // /* set up a signal handler for broken pipes */
    if (signal(SIGINT, sig_int) == SIG_ERR) exit(1);

    h3conn_open(IPV4, serv_ip, serv_port);

    seq[0] = 0;
    read_input(stdin);

    ptr = seq;
    while (*ptr && isspace(*ptr))
        ++ptr;

    parse_input_options(ptr);
    while (*ptr && isspace(*ptr))
        ++ptr;

    perform();

    h3conn_close();

    fclose(fp_targets);
    fclose(fp_domains);
    fclose(tblfp);
    fclose(domtblfp);
    fclose(pfamtblfp);
    return 0;
}
