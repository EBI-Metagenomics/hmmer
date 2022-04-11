#include "p7_config.h"

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
static int sock = 0;
static uint8_t *buf = 0;
static uint32_t buf_offset = 0;
static uint32_t hits_start = 0;
static ESL_STOPWATCH *w = 0;
static char serv_ip[64] = "127.0.0.1";
static unsigned short serv_port = SERVER_PORT;
static int ali = 0;
static int scores = 0;
static char opts[MAX_READ_LEN] = {0};

static FILE *fp_targets = 0;
static FILE *fp_domains = 0;
static FILE *fp_stats = 0;

static int read_input(void)
{
    int eod = 0;
    seq[0] = 0;

    while (!eod)
    {

        if (fgets(buffer, MAX_READ_LEN, stdin) != NULL)
        {
            int n = strlen(buffer);
            strcat(seq, buffer);

            eod = (strncmp(buffer, "//", 2) == 0);
        }
        else
        {
            eod = 1;
        }
    }
    return 0;
}

static int parse_input_options(char *ptr, ESL_GETOPTS *go)
{
    char t;
    char *s = ++ptr;

    /* skip to the end of the line */
    while (*ptr && (*ptr != '\n' && *ptr != '\r'))
        ++ptr;
    t = *ptr;
    *ptr = 0;

    /* create a commandline string with dummy program name for
     * the esl_opt_ProcessSpoof() function to parse.
     */
    snprintf(opts, MAX_READ_LEN, "X %s\n", s);

    if (esl_getopts_Reuse(go) != eslOK)
        p7_Die("Internal failure reusing options object");
    if (esl_opt_ProcessSpoof(go, opts) != eslOK)
    {
        printf("Failed to parse options string: %s\n", go->errbuf);
        return 1;
    }
    if (esl_opt_VerifyConfig(go) != eslOK)
    {
        printf("Failed to parse options string: %s\n", go->errbuf);
        return 1;
    }

    /* the options string can handle an optional database */
    if (esl_opt_ArgNumber(go) != 0)
    {
        printf("Incorrect number of command line arguments.");
        return 1;
    }

    /* skip remaining white spaces */
    *ptr = t;
    return 0;
}

static int perform(ESL_GETOPTS *go)
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
    if (writen(sock, seq, n) != n)
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

    if ((size = readn(sock, buf, n)) == -1)
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
        if ((size = readn(sock, ebuf, n)) == -1)
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
    if ((size = readn(sock, buf, n)) == -1)
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
    printf("Ponto 1\n");
    fflush(stdout);

    // Create the structures we'll deserialize the hits into
    pli = p7_pipeline_Create(go, 100, 100, FALSE,
                             (esl_opt_IsUsed(go, "--seqdb")) ? p7_SEARCH_SEQS
                                                             : p7_SCAN_MODELS);
    printf("Ponto 2\n");
    fflush(stdout);

    /* copy the search stats */
    w->elapsed = stats->elapsed;
    w->user = stats->user;
    w->sys = stats->sys;

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
    printf("Ponto 3\n");
    fflush(stdout);
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
    p7_pli_Statistics(fp_stats, pli, w);

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

static void parse_cmd_args(int argc, char *argv[])
{
    int i = 1;
    while (i < argc)
    {
        if (argv[i][0] != '-') usage(argv[0]);
        if (argv[i][1] == 0 || argv[i][2] != 0) usage(argv[0]);
        switch (argv[i][1])
        {
        case 'p':
            if (i + 1 >= argc)
            {
                printf("Missing port number\n");
                usage(argv[0]);
            }
            serv_port = atoi(argv[i + 1]);
            ++i;
            break;
        case 'i':
            if (i + 1 >= argc)
            {
                printf("Missing ip address\n");
                usage(argv[0]);
            }
            strcpy(serv_ip, argv[i + 1]);
            ++i;
            break;
        case 'A':
            ali = 1;
            scores = 1;
            break;
        case 'S':
            scores = 1;
            break;
        default:
            usage(argv[0]);
        }
        ++i;
    }
}

int main(int argc, char *argv[])
{
    fp_targets = fopen("targets.txt", "wb");
    fp_domains = fopen("domains.txt", "wb");
    fp_stats = fopen("stats.txt", "wb");

    ESL_GETOPTS *go = NULL;
    P7_PIPELINE *pli = NULL;
    P7_TOPHITS *th = NULL;

    HMMD_SEARCH_STATS *stats;
    HMMD_SEARCH_STATUS sstatus;
    char *ptr;

    struct sockaddr_in serv_addr;

    buf = NULL;
    buf_offset = 0;

    parse_cmd_args(argc, argv);

    /* set up a signal handler for broken pipes */
    if (signal(SIGINT, sig_int) == SIG_ERR) exit(1);

    /* Create a reliable, stream socket using TCP */
    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        fprintf(stderr, "[%s:%d] socket error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }

    /* Construct the server address structure */
    memset(&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = inet_addr(serv_ip);
    serv_addr.sin_port = htons(serv_port);
    if ((inet_pton(AF_INET, serv_ip, &serv_addr.sin_addr)) < 0)
    {
        fprintf(stderr, "[%s:%d] inet_pton error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }

    /* Establish the connection to the echo server */
    if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
        fprintf(stderr, "[%s:%d] connect error %d - %s\n", __FILE__, __LINE__,
                errno, strerror(errno));
        exit(1);
    }

    w = esl_stopwatch_Create();
    go = esl_getopts_Create(searchOpts);

    seq[0] = 0;
    while (strncmp(seq, "//", 2) != 0)
    {
        // int rem;
        int total = 0;
        read_input();

        /* skip all leading white spaces */
        ptr = seq;
        while (*ptr && isspace(*ptr))
            ++ptr;

        /* process search specific options */
        if (*ptr == '@')
        {
            if (parse_input_options(ptr, go)) continue;
            while (*ptr && isspace(*ptr))
                ++ptr;
        }

        if (*ptr && strncmp(ptr, "//", 2) != 0)
        {
            perform(go);
        }
    }

    esl_getopts_Destroy(go);
    esl_stopwatch_Destroy(w);

    close(sock);

    fclose(fp_targets);
    fclose(fp_domains);
    fclose(fp_stats);
    return 0;
}
