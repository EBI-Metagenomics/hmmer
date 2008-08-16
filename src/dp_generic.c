/* HMMER's standard implementation of the DP algorithms.
 * 
 * This implementation is modified from an optimized implementation
 * contributed by Jeremy D. Buhler (Washington University in
 * St. Louis).
 * 
 * Relative to the implementation in HMMER2, Jeremy rearranged data
 * structures to reduce the number of registers needed in the inner
 * loop; eliminated branches from the inner loop by unrolling the Mth
 * iteration in Viterbi and by replacing a bunch of "if" tests with
 * MAX; and exposed opportunities for hoisting and strength reduction
 * to the compiler. (The preceding sentence is nearly verbatim from
 * Jeremy's notes.) I then uplifted the JB code to H3, most notably 
 * involving conversion from H2's scaled integers to H3's floating
 * point calculations.
 * 
 * Contents:
 *     1. HMM alignment algorithms.
 *     2. Traceback algorithms. 
 *     3. Benchmark driver.
 *     4. Unit tests.
 *     5. Test driver.
 *     6. Copyright and license information.
 * 
 * SRE, Tue Jan 30 10:49:43 2007 [at Einstein's in St. Louis]
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"		/* when StochasticTrace() moves, don't need random or vectorops here */
#include "esl_vectorops.h"

#include "hmmer.h"





/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7DP_GENERIC_BENCHMARK
/*
   gcc -o benchmark-generic -g -O2 -I. -L. -I../easel -L../easel -Dp7DP_GENERIC_BENCHMARK dp_generic.c -lhmmer -leasel -lm
   icc -O3 -static -o benchmark-generic -I. -L. -I../easel -L../easel -Dp7DP_GENERIC_BENCHMARK dp_generic.c -lhmmer -leasel -lm
   ./benchmark-generic <hmmfile>
 */
/* As of Fri Dec 28 14:48:39 2007
 *    Viterbi  = 61.8 Mc/s
 *    Forward  =  8.6 Mc/s
 *   Backward  =  7.1 Mc/s
 *        MSV  = 55.9 Mc/s
 * (gcc -g -O2, 3.2GHz Xeon, N=50K, L=400, M=72 RRM_1 model)
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing: don't do DP",                   0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use the Backward algorithm",                     0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use the Forward algorithm",                      0 },
  { "-M",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use the MSV algorithm",                          0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",             0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  float           nullsc;
  float           bitscore;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);

  gx = p7_gmx_Create(gm->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      if (esl_opt_GetBoolean(go, "-b")) continue;

      if      (esl_opt_GetBoolean(go, "-F"))           p7_GForward     (dsq, L, gm, gx, &sc);
      else if (esl_opt_GetBoolean(go, "-B"))           p7_GBackward    (dsq, L, gm, gx, &sc);
      else if (esl_opt_GetBoolean(go, "-M"))           p7_GMSV         (dsq, L, gm, gx, &sc);
      else                                             p7_GViterbi     (dsq, L, gm, gx, &sc);

      p7_bg_NullOne(bg, dsq, L, &nullsc);
      bitscore = (sc - nullsc) / eslCONST_LOG2;
      
      if (esl_opt_GetBoolean(go, "-v")) printf("%.4f bits  (%.4f raw)\n", bitscore, sc); 
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  free(dsq);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7DP_GENERIC_BENCHMARK*/


/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7DP_GENERIC_TESTDRIVE
#include <string.h>

#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_msa.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"


/* The "basic" utest is a minimal driver for making a small DNA profile and a small DNA sequence,
 * then running Viterbi and Forward. It's useful for dumping DP matrices and profiles for debugging.
 */
static void
utest_basic(ESL_GETOPTS *go)
{
  char           *query= "# STOCKHOLM 1.0\n\nseq1 GAATTC\nseq2 GAATTC\n//\n";
  int             fmt  = eslMSAFILE_STOCKHOLM;
  char           *targ = "GAATTC";
  ESL_ALPHABET   *abc  = NULL;
  ESL_MSA        *msa  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  P7_DPRIOR      *pri  = NULL;	
  ESL_DSQ        *dsq  = NULL;
  P7_GMX         *gx   = NULL;
  P7_TRACE        *tr  = NULL;
  int             L    = strlen(targ);
  float           vsc, vsc2, fsc;

  if ((abc = esl_alphabet_Create(eslDNA))          == NULL)  esl_fatal("failed to create alphabet");
  if ((pri = p7_dprior_CreateNucleic())            == NULL)  esl_fatal("failed to create prior");
  if ((msa = esl_msa_CreateFromString(query, fmt)) == NULL)  esl_fatal("failed to create MSA");
  if (esl_msa_Digitize(abc, msa)                   != eslOK) esl_fatal("failed to digitize MSA");
  if (p7_Fastmodelmaker(msa, 0.5, &hmm, NULL)      != eslOK) esl_fatal("failed to create GAATTC model");
  if (p7_ParameterEstimation(hmm, pri)             != eslOK) esl_fatal("failed to parameterize GAATTC model");
  if ((bg = p7_bg_Create(abc))                     == NULL)  esl_fatal("failed to create DNA null model");
  if ((gm = p7_profile_Create(hmm->M, abc))        == NULL)  esl_fatal("failed to create GAATTC profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL)!= eslOK) esl_fatal("failed to config profile");
  if (p7_profile_Validate(gm, NULL, 0.0001)        != eslOK) esl_fatal("whoops, profile is bad!");
  if (esl_abc_CreateDsq(abc, targ, &dsq)           != eslOK) esl_fatal("failed to create GAATTC digital sequence");
  if ((gx = p7_gmx_Create(gm->M, L))               == NULL)  esl_fatal("failed to create DP matrix");
  if ((tr = p7_trace_Create())                     == NULL)  esl_fatal("trace creation failed");

  p7_GViterbi   (dsq, L, gm, gx, &vsc);
  if (esl_opt_GetBoolean(go, "-v")) printf("Viterbi score: %.4f\n", vsc);
  if (esl_opt_GetBoolean(go, "-v")) p7_gmx_Dump(stdout, gx);

  p7_GTrace     (dsq, L, gm, gx, tr);
  p7_trace_Score(tr, dsq, gm, &vsc2);
  if (esl_opt_GetBoolean(go, "-v")) p7_trace_Dump(stdout, tr, gm, dsq);
  
  if (esl_FCompare(vsc, vsc2, 1e-5) != eslOK)  esl_fatal("trace score and Viterbi score don't agree.");

  p7_GForward   (dsq, L, gm, gx, &fsc);
  if (esl_opt_GetBoolean(go, "-v")) printf("Forward score: %.4f\n", fsc);
  if (esl_opt_GetBoolean(go, "-v")) p7_gmx_Dump(stdout, gx);

  p7_trace_Destroy(tr);
  p7_gmx_Destroy(gx);
  free(dsq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_msa_Destroy(msa);
  p7_dprior_Destroy(pri);
  esl_alphabet_Destroy(abc);
  return;
}

/* Viterbi validation is done by comparing the returned score
 * to the score of the optimal trace. Not foolproof, but catches
 * many kinds of errors.
 * 
 * Another check is that the average score should be <= 0,
 * since the random sequences are drawn from the null model.
 */ 
static void
utest_viterbi(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq, int L)
{
  float     avg_sc = 0.;
  char      errbuf[eslERRBUFSIZE];
  ESL_DSQ  *dsq = NULL;
  P7_GMX   *gx  = NULL;
  P7_TRACE *tr  = NULL;
  int       idx;
  float     sc1, sc2;

  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((tr     = p7_trace_Create())              == NULL)  esl_fatal("trace creation failed");
  if ((gx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");

  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
      if (p7_GViterbi(dsq, L, gm, gx, &sc1)       != eslOK) esl_fatal("viterbi failed");
      if (p7_GTrace  (dsq, L, gm, gx, tr)         != eslOK) esl_fatal("trace failed");
      if (p7_trace_Validate(tr, abc, dsq, errbuf) != eslOK) esl_fatal("trace invalid:\n%s", errbuf);
      if (p7_trace_Score(tr, dsq, gm, &sc2)       != eslOK) esl_fatal("trace score failed");
      if (esl_FCompare(sc1, sc2, 1e-6)            != eslOK) esl_fatal("Trace score != Viterbi score"); 
      if (p7_bg_NullOne(bg, dsq, L, &sc2)         != eslOK) esl_fatal("null score failed");

      avg_sc += (sc1 - sc2);

      if (esl_opt_GetBoolean(go, "--vv"))       
	printf("utest_viterbi: Viterbi score: %.4f (null %.4f) (total so far: %.4f)\n", sc1, sc2, avg_sc);
    }

  avg_sc /= (float) nseq;
  if (avg_sc > 0.) esl_fatal("Viterbi scores have positive expectation (%f nats)", avg_sc);

  p7_gmx_Destroy(gx);
  p7_trace_Destroy(tr);
  free(dsq);
  return;
}


/* Forward is harder to validate. 
 * We do know that the Forward score is >= Viterbi.
 * We also know that the expected score on random seqs is <= 0 (not
 * exactly - we'd have to sample the random length from the background
 * model too, not just use a fixed L - but it's close enough to
 * being true to be a useful test.)
 */
static void
utest_forward(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq, int L)
{
  float     avg_sc;
  ESL_DSQ  *dsq = NULL;
  P7_GMX   *gx  = NULL;
  P7_TRACE *tr  = NULL;
  int       idx;
  float     vsc, fsc, nullsc;

  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((gx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");

  avg_sc = 0.;
  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
      if (p7_GViterbi(dsq, L, gm, gx, &vsc)       != eslOK) esl_fatal("viterbi failed");
      if (p7_GForward(dsq, L, gm, gx, &fsc)       != eslOK) esl_fatal("forward failed");
      if (fsc < vsc) esl_fatal("Forward score can't be less than Viterbi score");
      if (p7_bg_NullOne(bg, dsq, L, &nullsc)      != eslOK) esl_fatal("null score failed");

      avg_sc += fsc - nullsc;

      if (esl_opt_GetBoolean(go, "--vv")) 
	printf("utest_forward: Forward score: %.4f (total so far: %.4f)\n", fsc, avg_sc);
    }

  avg_sc /= (float) nseq;
  if (avg_sc > 0.) esl_fatal("Forward scores have positive expectation (%f nats)", avg_sc);

  p7_gmx_Destroy(gx);
  p7_trace_Destroy(tr);
  free(dsq);
  return;
}


/* The MSV score can be validated against Viterbi (provided we trust
 * Viterbi), by creating a multihit local profile in which:
 *   1. All t_MM scores = 0
 *   2. All other core transitions = -inf
 *   3. All t_BMk entries uniformly log 2/(M(M+1))
 */
static void
utest_msv(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq, int L)
{
  P7_PROFILE *g2 = NULL;
  ESL_DSQ   *dsq = NULL;
  P7_GMX    *gx  = NULL;
  float     sc1, sc2;
  int       k, idx;

  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((gx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");
  if ((g2     = p7_profile_Clone(gm))           == NULL)  esl_fatal("profile clone failed");

  /* Make g2's scores appropriate for simulating the MSV algorithm in Viterbi */
  esl_vec_FSet(g2->tsc, p7P_NTRANS * g2->M, -eslINFINITY);
  for (k = 1; k <  g2->M; k++) p7P_TSC(g2, k, p7P_MM) = 0.0f;
  for (k = 0; k <  g2->M; k++) p7P_TSC(g2, k, p7P_BM) = log(2.0f / ((float) g2->M * (float) (g2->M+1)));

  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");

      if (p7_GMSV    (dsq, L, gm, gx, &sc1)       != eslOK) esl_fatal("MSV failed");
      if (p7_GViterbi(dsq, L, g2, gx, &sc2)       != eslOK) esl_fatal("viterbi failed");
      if (fabs(sc1-sc2) > 0.0001) esl_fatal("MSV score not equal to Viterbi score");
    }

  p7_gmx_Destroy(gx);
  p7_profile_Destroy(g2);
  free(dsq);
  return;
}


/* The "generation" test scores sequences generated by the same profile.
 * Each Viterbi and Forward score should be >= the trace score of the emitted seq.
 * The expectation of Forward scores should be positive.
 */
static void
utest_generation(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc,
		 P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, int nseq)
{
  ESL_SQ   *sq = esl_sq_CreateDigital(abc);
  P7_GMX   *gx = p7_gmx_Create(gm->M, 100);
  P7_TRACE *tr = p7_trace_Create();
  float     vsc, fsc, nullsc, tracesc;
  float     avg_fsc;
  int       idx;

  avg_fsc = 0.0;
  for (idx = 0; idx < nseq; idx++)
    {
      if (p7_ProfileEmit(r, hmm, gm, bg, sq, tr)     != eslOK) esl_fatal("profile emission failed");

      if (p7_gmx_GrowTo(gx, gm->M, sq->n)            != eslOK) esl_fatal("failed to reallocate gmx");
      if (p7_GViterbi(sq->dsq, sq->n, gm, gx, &vsc)  != eslOK) esl_fatal("viterbi failed");
      if (p7_GForward(sq->dsq, sq->n, gm, gx, &fsc)  != eslOK) esl_fatal("forward failed");
      if (p7_trace_Score(tr, sq->dsq, gm, &tracesc)  != eslOK) esl_fatal("trace score failed");
      if (p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc) != eslOK) esl_fatal("null score failed");

      if (vsc < tracesc) esl_fatal("viterbi score is less than trace");
      if (fsc < tracesc) esl_fatal("forward score is less than trace");

      if (esl_opt_GetBoolean(go, "--vv")) 
	printf("generated:  len=%d v=%8.4f  f=%8.4f  t=%8.4f\n", (int) sq->n, vsc, fsc, tracesc);
      
      avg_fsc += (fsc - nullsc);
    }
  
  avg_fsc /= (float) nseq;
  if (avg_fsc < 0.) esl_fatal("generation: Forward scores have negative expectation (%f nats)", avg_fsc);

  p7_gmx_Destroy(gx);
  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
}

/* The "enumeration" test samples a random enumerable HMM (transitions to insert are 0,
 * so the generated seq space only includes seqs of L<=M). 
 *
 * The test scores all seqs of length <=M by both Viterbi and Forward, verifies that 
 * the two scores are identical, and verifies that the sum of all the probabilities is
 * 1.0. It also verifies that the score of a sequence of length M+1 is indeed -infinity.
 * 
 * Because this function is going to work in unscaled probabilities, adding them up,
 * all P(seq) terms must be >> DBL_EPSILON.  That means M must be small; on the order 
 * of <= 10. 
 */
static void
utest_enumeration(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, int M)
{
  char            errbuf[eslERRBUFSIZE];
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  ESL_DSQ        *dsq  = NULL;
  P7_GMX         *gx   = NULL;
  float  vsc, fsc;
  float  bg_ll;   		/* log P(seq | bg) */
  double vp, fp;		/* P(seq,\pi | model) and P(seq | model) */
  int L;
  int i;
  double total_p;
  char   *seq;
    
  /* Sample an enumerable HMM & profile of length M.
   */
  if (p7_hmm_SampleEnumerable(r, M, abc, &hmm)      != eslOK) esl_fatal("failed to sample an enumerable HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, 0, p7_UNILOCAL) != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm, errbuf, 0.0001)      != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);
  if (p7_profile_Validate(gm, errbuf, 0.0001)       != eslOK) esl_fatal("whoops, profile is bad!: %s", errbuf);

  if (  (dsq = malloc(sizeof(ESL_DSQ) * (M+3)))     == NULL)  esl_fatal("allocation failed");
  if (  (seq = malloc(sizeof(char)    * (M+2)))     == NULL)  esl_fatal("allocation failed");
  if ((gx     = p7_gmx_Create(hmm->M, M+3))         == NULL)  esl_fatal("matrix creation failed");

  /* Enumerate all sequences of length L <= M
   */
  total_p = 0;
  for (L = 0; L <= M; L++)
    {
      /* Initialize dsq of length L at 0000... */
      dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
      for (i = 1; i <= L; i++) dsq[i] = 0;

      while (1) 		/* enumeration of seqs of length L*/
	{
	  if (p7_GViterbi(dsq, L, gm, gx, &vsc)  != eslOK) esl_fatal("viterbi failed");
	  if (p7_GForward(dsq, L, gm, gx, &fsc)  != eslOK) esl_fatal("forward failed");
 
	  /* calculate bg log likelihood component of the scores */
	  for (bg_ll = 0., i = 1; i <= L; i++)  bg_ll += log(bg->f[dsq[i]]);
	  
	  /* convert to probabilities, adding the bg LL back to the LLR */
	  vp =  exp(vsc + bg_ll);
	  fp =  exp(fsc + bg_ll);

	  if (esl_opt_GetBoolean(go, "--vv")) {
	    esl_abc_Textize(abc, dsq, L, seq);
	    printf("probability of sequence: %10s   %16g  (lod v=%8.4f f=%8.4f)\n", seq, fp, vsc, fsc);
	  }
	  total_p += fp;

	  /* Increment dsq like a reversed odometer */
	  for (i = 1; i <= L; i++) 
	    if (dsq[i] < abc->K-1) { dsq[i]++; break; } else { dsq[i] = 0; }
	  if (i > L) break;	/* we're done enumerating sequences */
	}
    }

  /* That sum is subject to significant numerical error because of
   * discretization error in FLogsum(); don't expect it to be too close.
   */
  if (total_p < 0.999 || total_p > 1.001) esl_fatal("Enumeration unit test failed: total Forward p isn't near 1.0 (%g)", total_p);
  if (esl_opt_GetBoolean(go, "-v")) {
    printf("enumeration test: total p is %g\n", total_p);
  }
  
  p7_gmx_Destroy(gx);
  p7_bg_Destroy(bg);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  free(dsq);
  free(seq);
}
#endif /*p7DP_GENERIC_TESTDRIVE*/




/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7DP_GENERIC_TESTDRIVE -I. -I../easel -L. -L../easel -o dp_generic_utest dp_generic.c -lhmmer -leasel -lm
 */
#ifdef p7DP_GENERIC_TESTDRIVE
#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  char            errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS *r    = NULL;
  ESL_ALPHABET   *abc  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  int             M    = 100;
  int             L    = 200;
  int             nseq = 20;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  utest_basic(go);
  
  if ((abc = esl_alphabet_Create(eslAMINO))         == NULL)  esl_fatal("failed to create alphabet");
  if (p7_hmm_Sample(r, M, abc, &hmm)                != eslOK) esl_fatal("failed to sample an HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)    != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm, errbuf, 0.0001)      != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);
  if (p7_profile_Validate(gm,  errbuf, 0.0001)      != eslOK) esl_fatal("whoops, profile is bad!: %s", errbuf);

  utest_viterbi    (go, r, abc, bg, gm, nseq, L);
  utest_forward    (go, r, abc, bg, gm, nseq, L);
  utest_msv        (go, r, abc, bg, gm, nseq, L);
  utest_generation (go, r, abc, gm, hmm, bg, nseq);
  utest_enumeration(go, r, abc, 4);	/* can't go much higher than 5; enumeration test is cpu-intensive. */
  
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7DP_GENERIC_TESTDRIVE*/

/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7DP_GENERIC_EXAMPLE
/* 
   gcc -g -O2 -Dp7DP_GENERIC_EXAMPLE -I. -I../easel -L. -L../easel -o example dp_generic.c -lhmmer -leasel -lm
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_dmatrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of a forward/backward posterior probability heat map";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *r       = esl_randomness_Create(42);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMX         *bck     = NULL;
  ESL_DMATRIX    *pp      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  P7_TRACE       *tr      = NULL;
  P7_ALIDISPLAY  *ad      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           sc;
  int             i,d,k;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");
  esl_sqfile_Close(sqfp);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  
  /* allocate DP matrices for forward and backward */
  fwd = p7_gmx_Create(gm->M, sq->n);
  bck = p7_gmx_Create(gm->M, sq->n);

  /* run Forward, Backward */
  tr = p7_trace_Create();

  p7_GViterbi (sq->dsq, sq->n, gm, fwd, &sc);
  p7_GTrace   (sq->dsq, sq->n, gm, fwd, tr);
  p7_trace_Index(tr);
  printf("# Viterbi: %d domains : ", tr->ndom);
  for (d = 0; d < tr->ndom; d++) printf("%6d %6d %6d %6d  ", tr->sqfrom[d], tr->sqto[d], tr->hmmfrom[d], tr->hmmto[d]);
  printf("\n");
  p7_trace_Reuse(tr);

  p7_GForward (sq->dsq, sq->n, gm, fwd, &sc);
  p7_GBackward(sq->dsq, sq->n, gm, bck, &sc);

  for (i = 0; i < 1000; i++)
    {
      p7_GStochasticTrace(r, sq->dsq, sq->n, gm, fwd, tr);
      p7_trace_Index(tr);
      /* printf("%3d  ", tr->ndom); */

      for (d = 0; d < tr->ndom; d++) printf("%6d %6d %6d %6d %6d %6d\n", 
					    tr->sqfrom[d], tr->sqto[d], tr->hmmfrom[d], tr->hmmto[d], 
					    tr->sqfrom[d]-tr->hmmfrom[d], tr->sqto[d]-tr->hmmto[d]);
	
#if 0
      for (d = 0; d < tr->ndom; d++) {
	printf("domain %d of %d\n", d+1, tr->ndom);
	ad = p7_alidisplay_Create(tr, d, gm, sq);
	p7_alidisplay_Print(stdout, ad, 40, 80);
	p7_alidisplay_Destroy(ad);
      }
#endif
      p7_trace_Reuse(tr);
    }
  p7_trace_Destroy(tr);
  
#if 0
  /* construct a LxM matrix holding posterior probs for each Match state emitting residue i */
  pp = esl_dmatrix_Create(sq->n, gm->M);
  for (i = 1; i <= sq->n; i++)
    for (k = 1; k <= gm->M; k++)
      pp->mx[i-1][k-1] = fwd->dp[i][k*3] + bck->dp[i][k*3] - sc;

  /* output a heatmap */
  dmx_Visualize(stdout, pp, -8.0, 0.0);
#endif

#if 0
  printf("min = %g\n",       esl_dmx_Min(pp));
  printf("max = %g\n",       esl_dmx_Max(pp));
  printf("sc  = %g nats\n",  sc);
#endif

  esl_dmatrix_Destroy(pp);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7DP_GENERIC_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
