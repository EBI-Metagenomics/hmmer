/* Building profile HMMs from alignments.
 *
 * Contents: 
 *    x. h4_Build():      build new profile from alignment
 *    x. Internal routines for profile construction
 *    x. H4_BUILD_CONFIG: customization of h4_Build()
 *    x. _experiment:     save counts files for training priors
 *    x. _experiment2:    compare old vs. new fragment marking
 * 
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_matrixops.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_vectorops.h"


#include "h4_path.h"
#include "h4_prior.h"
#include "h4_profile.h"

#include "build.h"
#include "eweight.h"
#include "parameterize.h"


static int consensus_by_symfrac(const ESL_MSA *msa, float symfrac, const int8_t *fragassign, int8_t *matassign);
static int consensus_by_hand   (const ESL_MSA *msa, int8_t *matassign, char *errbuf);
static int mark_fragments      (const ESL_MSA *msa, float fragthresh, int8_t *fragassign);
static int collect_counts      (const ESL_MSA *msa, const int8_t *fragassign, const int8_t *matassign, H4_PROFILE *hmm);

/*****************************************************************
 * 1. h4_Build(): build new profile from alignment
 *****************************************************************/

/* Function:  h4_Build()
 * Synopsis:  Build new profile from alignment
 * Incept:    SRE, Sun 15 Jul 2018 [Benasque]
 *
 * Purpose:   Given alignment <msa>, build a new profile and return it
 *            thru <*ret_hmm>. Caller can optionally provide a custom
 *            configuration for build options in <cfg>, or pass <NULL>
 *            for defaults. If anything goes wrong that's the user's
 *            fault, an informative error message is left in optional
 *            <errbuf>, if caller provides one.
 *            
 * Args:      cfg     : OPTIONAL: custom options, or <NULL> to use defaults
 *            msa     : alignment to build profile from
 *            ret_hmm : RETURN  : ptr to new profile
 *            errbuf  : OPTIONAL: allocated for <eslERRBUFSIZE>; or <NULL>
 *            
 * Returns:   <eslOK> on success. <*ret_hmm> is a new profile; 
 *            <errbuf> is '\0'.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_Build(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE **ret_hmm, char *errbuf)
{
  int         arch_strategy  = (cfg ? cfg->arch_strategy : h4BUILD_ARCH_RULES);
  float       symfrac        = (cfg ? cfg->symfrac       : h4BUILD_SYMFRAC);
  float       fragthresh     = (cfg ? cfg->fragthresh    : h4BUILD_FRAGTHRESH);
  int         wgt_strategy   = (cfg ? cfg->wgt_strategy  : h4BUILD_WGT_PB);
  float       wid            = (cfg ? cfg->wid           : h4BUILD_WID);
  int         effn_strategy  = (cfg ? cfg->effn_strategy : h4BUILD_EFFN_EWGT);
  float       re_target;
  float       re_sigma       = (cfg ? cfg->re_sigma      : h4BUILD_ESIGMA);
  float       n_eff          = (cfg ? cfg->effn_set      : (float) msa->nseq);    // default eweighting will reset this.
  H4_PRIOR   *pri            = ((cfg && cfg->pri) ? cfg->pri : h4_prior_Create(msa->abc));
  int         stop_early     = (cfg ? cfg->stop_early    : FALSE);
  ESL_MSAWEIGHT_CFG *wgt_cfg = NULL;
  H4_PROFILE *hmm            = NULL;
  int8_t     *fragassign     = NULL;
  int8_t     *matassign      = NULL;
  int         M              = 0;

  int         apos;
  int         status = eslOK;

  ESL_DASSERT1(( msa->flags && eslMSA_DIGITAL ));
  ESL_DASSERT1(( !cfg || !cfg->abc || cfg->abc->type == msa->abc->type ));

  if      (cfg)                        re_target = cfg->re_target;
  else if (msa->abc->type == eslAMINO) re_target = h4BUILD_ETARG_PRT;
  else if (msa->abc->type == eslDNA)   re_target = h4BUILD_ETARG_NUC;
  else if (msa->abc->type == eslRNA)   re_target = h4BUILD_ETARG_NUC;
  else                                 re_target = h4BUILD_ETARG_OTH;

  ESL_ALLOC(fragassign, sizeof(int8_t) * msa->nseq);
  ESL_ALLOC(matassign,  sizeof(int8_t) * (msa->alen + 1));
  if (( wgt_cfg = esl_msaweight_cfg_Create() ) == NULL) { status = eslEMEM; goto ERROR; }
  matassign[0] = 0; // unused; convention.c
  if (errbuf) errbuf[0] = '\0';


  // validate_msa()
  // esl_msa_Checksum()

  /* 1. Define which sequences are considered to be fragments (local alignments).
   */
  if (( status = mark_fragments(msa, fragthresh, fragassign)) != eslOK) goto ERROR;

  /* 2. Set relative weights, in msa->wgt.
   *    PB weighting determines consensus columns, so share the fragthresh/symfrac params.
   */
  wgt_cfg->fragthresh = fragthresh;
  wgt_cfg->symfrac    = symfrac;

  switch (wgt_strategy) {
  case h4BUILD_WGT_NONE:   esl_vec_DSet(msa->wgt, msa->nseq, 1.);                                         break;
  case h4BUILD_WGT_GIVEN:                                                                                 break;
  case h4BUILD_WGT_PB:     if ((status = esl_msaweight_PB_adv(wgt_cfg, msa, NULL))  != eslOK) goto ERROR; break;
  case h4BUILD_WGT_GSC:    if ((status = esl_msaweight_GSC(msa))                    != eslOK) goto ERROR; break;
  case h4BUILD_WGT_BLOSUM: if ((status = esl_msaweight_BLOSUM(msa, wid))            != eslOK) goto ERROR; break;
  default: esl_fatal("no such weighting strategy");
  }
  
  /* 3. Define which columns are considered to be consensus.
   */
  switch (arch_strategy) {
  case h4BUILD_ARCH_RULES: if ((status = consensus_by_symfrac(msa, symfrac, fragassign, matassign)) != eslOK) goto ERROR; break;
  case h4BUILD_ARCH_GIVEN: if ((status = consensus_by_hand(msa, matassign, errbuf))                 != eslOK) goto ERROR; break;
  }

  /* Allocate the new profile.
   */
  for (apos = 1; apos <= msa->alen; apos++) if (matassign[apos]) M++;
  hmm = h4_profile_Create(msa->abc, M);

  /* Collect observed (relative-weighted) counts from alignment in hmm->t[] and ->e[]
   */
  if ((status = collect_counts(msa, fragassign, matassign, hmm)) != eslOK) goto ERROR;

  if (stop_early) goto DONE;

  // mask_columns();

  /* Determine and apply effective sequence number
   */
  if (effn_strategy == h4BUILD_EFFN_EWGT)
    {
      re_target = ESL_MAX(re_target, (re_sigma - log2f( 2.0 / ((float) M * (float) (M+1)))) / (float) M); // assure a minimum total expected score, for short models. [J5/36]

      if (( status = h4_EntropyWeight(hmm, pri, msa->nseq, re_target, &n_eff)) != eslOK) goto ERROR;
  
      esl_mat_FScale(hmm->e, hmm->M+1, hmm->abc->K,  (n_eff / (float) msa->nseq));
      esl_mat_FScale(hmm->t, hmm->M+1, h4_NT,        (n_eff / (float) msa->nseq));
    }
  else if (effn_strategy == h4BUILD_EFFN_GIVEN)
    {
      esl_mat_FScale(hmm->e, hmm->M+1, hmm->abc->K,  (n_eff / (float) msa->nseq));
      esl_mat_FScale(hmm->t, hmm->M+1, h4_NT,        (n_eff / (float) msa->nseq));
    }

  /* Convert counts to mean posterior probability parameters
   */
  h4_Parameterize(hmm, pri);

  // annotate();
  // calibrate();
  // make_post_msa()


 DONE:
 ERROR:
  if (status == eslOK) *ret_hmm = hmm; else { *ret_hmm = NULL; h4_profile_Destroy(hmm); }
  esl_msaweight_cfg_Destroy(wgt_cfg);
  if (! cfg || ! cfg->pri) h4_prior_Destroy(pri); // if <cfg> provided the prior, <cfg> is managing that memory, not us.
  free(fragassign);
  free(matassign);
  return status;
}


/*****************************************************************
 * 2. Internal routines for profile construction
 *****************************************************************/

/* consensus_by_symfrac()
 * Sets <matassign[1..alen]> to 1/0 flags, defining consensus columns.
 * 
 * Args:    msa        : multiple sequence alignment
 *          symfrac    : define col as consensus if weighted residue fraction >= symfrac
 *          fragassign : [0..nseq-1] 1/0 flags marking fragments (local alignments)
 *          matassign  : RETURN: [1..alen] 1/0 flags marking consensus columns
 * 
 * Returns: <eslOK> on success;
 *          and matassign[(0)1..alen] has 1/0 for consensus vs. non.
 *          
 * Throws:  <eslEMEM> on allocation error.
 *          Now state of <matassign> is undefined.         
 */
static int
consensus_by_symfrac(const ESL_MSA *msa, float symfrac, const int8_t *fragassign, int8_t *matassign)
{
  float *r      = NULL;  // weighted residue count for each column 1..alen
  float *totwgt = NULL;  // weighted residue+gap count for each column 1..alen (not constant across cols, because of fragments)
  int    apos, idx;
  int    lpos, rpos;
  int    status;

  ESL_ALLOC(r,      sizeof(float) * (msa->alen+1));
  ESL_ALLOC(totwgt, sizeof(float) * (msa->alen+1));
  esl_vec_FSet(r,      msa->alen+1, 0.);
  esl_vec_FSet(totwgt, msa->alen+1, 0.);

  for (idx = 0; idx < msa->nseq; idx++)
    {
      lpos = 1;
      rpos = msa->alen;
      if (fragassign[idx])
	{
	  for (lpos = 1;         lpos <= msa->alen; lpos++) if (! esl_abc_XIsGap(msa->abc, msa->ax[idx][lpos])) break;  // find first residue of this seq
	  for (rpos = msa->alen; rpos >= 1;         rpos--) if (! esl_abc_XIsGap(msa->abc, msa->ax[idx][rpos])) break;  //  ... and last.
	  // if sequence is empty: lpos = msa->alen+1, rpos = 0, and the loop body below doesn't execute.
	}
      for (apos = lpos; apos <= rpos; apos++)
	{
	  if      (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) { r[apos] += msa->wgt[idx]; totwgt[apos] += msa->wgt[idx]; }
	  else if (esl_abc_XIsGap(msa->abc,     msa->ax[idx][apos])) {                           totwgt[apos] += msa->wgt[idx];  }
	  // missing ~, nonresidue * don't count either way toward the calculation.
	}
    }

  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (totwgt[apos] > 0. && r[apos] / totwgt[apos] >= symfrac) ? 1 : 0;

  free(r);
  free(totwgt);
  return eslOK;

 ERROR:
  free(r);
  free(totwgt);
  return status;
}


/* consensus_by_hand()
 * Define consensus columns using provided alignment annotation (#=GC RF or seq_cons)
 */
static int
consensus_by_hand(const ESL_MSA *msa, int8_t *matassign, char *errbuf)
{
  int apos;
  
  if (! msa->rf)
    ESL_FAIL(eslEFORMAT, errbuf, "no consensus column (#=GC RF, #=GC seq_cons) annotation on MSA");
  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1])? 0 : 1);  // watch off-by-one. rf is 0..alen-1, matassign is 1..alen
  return eslOK;
}



/* mark_fragments()
 * Set fragassign[i] TRUE | FALSE to mark local alignment frags
 * SRE, Fri 06 Jul 2018 [JB1685 BOS-PIT]
 * 
 * Heuristically define sequence fragments (as opposed to "full
 * length" sequences) in <msa>. Set <fragassign[i]> to TRUE if seq <i>
 * is a fragment, else FALSE.
 * 
 * Args:    msa        : msa with msa->nseq seqs
 *          fragthresh : if alen[i]/alen <= fragthresh, seq is a fragment
 *          fragassign : RESULT: fragassign[i]=1|0 if seq i is|isn't a fragment.
 *                       Caller provides allocation for <msa->nseq> flags.
 * 
 * Xref:  See h4_build.md notes for why we use this ad hoc rule, as
 *        opposed to others you might imagine.
 */
static int
mark_fragments(const ESL_MSA *msa, float fragthresh, int8_t *fragassign)
{
  int   idx,lpos,rpos;
  float alispan;

  for (idx = 0; idx < msa->nseq; idx++)
    {
      for (lpos = 1;         lpos <= msa->alen; lpos++) if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][lpos])) break;
      for (rpos = msa->alen; rpos >= 1;         rpos--) if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][rpos])) break;
      /* L=0 sequence? lpos == msa->alen+1, rpos == 0; lpos > rpos. 
       * alen=0 alignment? lpos == 1, rpos == 0; lpos > rpos.
       */
      alispan = (lpos > rpos ? 0. : (float) (rpos - lpos + 1));
      alispan /= (float) msa->alen;
      fragassign[idx] = ( alispan < fragthresh ? TRUE : FALSE );
    }
  return eslOK;
}



static int
collect_counts(const ESL_MSA *msa, const int8_t *fragassign, const int8_t *matassign, H4_PROFILE *hmm)
{
  H4_PATH *pi = h4_path_Create();
  int      idx;
  int      status;
  
  for (idx = 0; idx < msa->nseq; idx++)
    {
      if (fragassign[idx]) { if ((status = h4_path_InferLocal (msa->abc, msa->ax[idx], msa->alen, matassign, pi)) != eslOK) goto ERROR; }
      else                 { if ((status = h4_path_InferGlocal(msa->abc, msa->ax[idx], msa->alen, matassign, pi)) != eslOK) goto ERROR; }

      // SRE DEBUGGING
      //status = h4_path_Validate(pi, msa->abc, hmm->M, esl_abc_dsqrlen(msa->abc, msa->ax[idx]), errbuf);
      //if (status != eslOK) esl_fatal(errbuf);

      //h4_path_Dump(stdout, pi);

      if ((status = h4_path_Count(pi, msa->ax[idx], msa->wgt[idx], hmm)) != eslOK) goto ERROR;

      h4_path_Reuse(pi);
    }

  h4_path_Destroy(pi);
  return eslOK;

 ERROR:
  h4_path_Destroy(pi);
  return status;
}


#if 0
static int
annotate(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE *hmm, char *errbuf, int8_t *matassign)
{
  /* transfer from MSA to new HMM: */
  /* msa->rf   */
  /* msa->mm   */
  /* msa->ss_cons */
  /* msa->sa_cons */
  /* map */
  

  return eslOK;
}
#endif


/*****************************************************************
 * x. H4_BUILD_CONFIG
 *****************************************************************/

H4_BUILD_CONFIG *
h4_build_config_Create(const ESL_ALPHABET *abc)
{
  H4_BUILD_CONFIG *cfg = NULL;
  int              status;

  ESL_ALLOC(cfg, sizeof(H4_BUILD_CONFIG));

  cfg->arch_strategy       = h4BUILD_ARCH_RULES;
  cfg->symfrac             = h4BUILD_SYMFRAC;
  cfg->fragthresh          = h4BUILD_FRAGTHRESH;

  cfg->wgt_strategy        = h4BUILD_WGT_PB;
  cfg->wid                 = h4BUILD_WID;

  cfg->effn_strategy       = h4BUILD_EFFN_EWGT;
  switch (abc->type) {
  case eslAMINO: cfg->re_target = h4BUILD_ETARG_PRT;  break;
  case eslDNA:   cfg->re_target = h4BUILD_ETARG_NUC;  break;
  case eslRNA:   cfg->re_target = h4BUILD_ETARG_NUC;  break;
  default:       cfg->re_target = h4BUILD_ETARG_OTH;  break;
  }
  cfg->re_sigma            = h4BUILD_ESIGMA;
  cfg->effn_set            = -1.;
  cfg->pri                 = NULL;
  cfg->stop_early          = FALSE;
  cfg->abc                 = abc;
  return cfg;

 ERROR:
  h4_build_config_Destroy(cfg);
  return NULL;
}

void
h4_build_config_Destroy(H4_BUILD_CONFIG *cfg)
{
  if (cfg)
    {
      h4_prior_Destroy(cfg->pri);
      free(cfg);
    }
      
}
/*----------------- end, H4_BUILD_CONFIG ------------------------*/




/*****************************************************************
 * x. _experiment: save counts files for training priors
 *****************************************************************/
#ifdef h4BUILD_EXPERIMENT
#include "h4_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"

#include "h4_profile.h"

#include "build.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--dna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                            0 },
  { "--rna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                            0 },
  { "--amino",   eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile> <outpfx>";
static char banner[] = "utility for saving counts files for training priors";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *msafile = esl_opt_GetArg(go, 1);
  char           *outpfx  = esl_opt_GetArg(go, 2);
  int             infmt   = eslMSAFILE_UNKNOWN;
  ESL_ALPHABET   *abc     = NULL;
  ESL_MSAFILE    *afp     = NULL;
  ESL_MSA        *msa     = NULL;
  H4_BUILD_CONFIG *cfg    = NULL;
  H4_PROFILE     *hmm     = NULL;
  char           *efile   = NULL;                      // counts file for match emissions
  char           *mtfile  = NULL;                      //             ... match transitions
  char           *itfile  = NULL;                      //             ... insert transitions 
  char           *dtfile  = NULL;                      //             ... delete transitions
  FILE           *efp     = NULL;
  FILE           *mtfp    = NULL;
  FILE           *itfp    = NULL;
  FILE           *dtfp    = NULL;
  int             nali    = 0;
  int             k,a,z;
  char            errbuf[eslERRBUFSIZE];
  int             status;

  

  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  esl_sprintf(&efile,  "%s.ect", outpfx);
  esl_sprintf(&mtfile, "%s.mct", outpfx);
  esl_sprintf(&itfile, "%s.ict", outpfx);
  esl_sprintf(&dtfile, "%s.dct", outpfx);

  if (( efp  = fopen(efile,  "w")) == NULL) esl_fatal("failed to open %s", efile);
  if (( mtfp = fopen(mtfile, "w")) == NULL) esl_fatal("failed to open %s", mtfile);
  if (( itfp = fopen(itfile, "w")) == NULL) esl_fatal("failed to open %s", itfile);
  if (( dtfp = fopen(dtfile, "w")) == NULL) esl_fatal("failed to open %s", dtfile);

  status = esl_msafile_Open(&abc, msafile, NULL, infmt, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  cfg = h4_build_config_Create(abc);
  cfg->stop_early = TRUE;

  while ((status = esl_msafile_Read(afp, &msa)) == eslOK)
    {
      status = h4_Build(cfg, msa, &hmm, errbuf);

      for (k = 1; k <= hmm->M; k++)
	if (esl_vec_FSum(hmm->e[k], hmm->abc->K) > 0.)
	  for (a = 0; a < abc->K; a++)  fprintf(efp,  "%10.2f%c", hmm->e[k][a], a == abc->K-1 ? '\n' : ' ');

      for (k = 1; k < hmm->M; k++)
	{
	  if (esl_vec_FSum(hmm->t[k]+h4_TMM, 3) > 0.)
	    for (z = 0; z < 3; z++) fprintf(mtfp, "%10.2f%c", hmm->t[k][h4_TMM+z], z == 2 ? '\n' : ' ');
	  if (esl_vec_FSum(hmm->t[k]+h4_TIM, 3) > 0.)
	    for (z = 0; z < 3; z++) fprintf(itfp, "%10.2f%c", hmm->t[k][h4_TIM+z], z == 2 ? '\n' : ' ');
	  if (esl_vec_FSum(hmm->t[k]+h4_TDM, 3) > 0.)
	    for (z = 0; z < 3; z++) fprintf(dtfp, "%10.2f%c", hmm->t[k][h4_TDM+z], z == 2 ? '\n' : ' ');
	}
      
      esl_msa_Destroy(msa);    msa = NULL;
      h4_profile_Destroy(hmm); hmm = NULL;
      nali++;
    }
  if (nali == 0 || status != eslEOF)
    esl_msafile_ReadFailure(afp, status);

  if (efp)  fclose(efp);
  if (mtfp) fclose(mtfp);
  if (itfp) fclose(itfp);
  if (dtfp) fclose(dtfp);
  free(efile);
  free(mtfile);
  free(itfile);
  free(dtfile);
  esl_msafile_Close(afp);    
  esl_alphabet_Destroy(abc);
  h4_build_config_Destroy(cfg);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /* h4BUILD_EXPERIMENT */



/*****************************************************************
 * x. _experiment2: compare old vs. new fragment marking 
 *****************************************************************/

#ifdef h4BUILD_EXPERIMENT2
#include "h4_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",        0 },
  { "--dna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                            0 },
  { "--rna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                            0 },
  { "--amino",   eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "testing old v. new fragment-marking strategy";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char        *msafile = esl_opt_GetArg(go, 1);
  int          infmt   = eslMSAFILE_UNKNOWN;
  ESL_ALPHABET *abc    = NULL;
  ESL_MSAFILE  *afp    = NULL;
  ESL_MSA      *msa    = NULL;
  int           nali   = 0;
  float         fragthresh = 0.5;
  int           nold, nnew;
  int           idx;
  int           status;
  
  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  status = esl_msafile_Open(&abc, msafile, NULL, infmt, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  esl_dataheader(stdout,
		 20, "name",  10,  "nseq",  10, "alen",     
		 10, "n_old", 10, "n_new",  10, "frac_old",
		 10, "frac_new", 0);

  while ((status = esl_msafile_Read(afp, &msa)) == eslOK)
    {
      int8_t *old_fragassign = NULL;
      int8_t *new_fragassign = NULL;

      ESL_ALLOC(old_fragassign, sizeof(int8_t) * msa->nseq);
      ESL_ALLOC(new_fragassign, sizeof(int8_t) * msa->nseq);
      
      if ((status = mark_fragments(msa, h4BUILD_FRAGTHRESH, new_fragassign)) != eslOK) goto ERROR;

      /* Reproduce the H3 calculation, while setting flag instead of marking ~ in the msa */
      for (idx = 0; idx < msa->nseq; idx++)
	old_fragassign[idx] = 
	  ((float) esl_abc_dsqrlen(msa->abc, msa->ax[idx]) / (float) msa->alen <= fragthresh) ?  
	  TRUE : FALSE;

      for (idx = 0, nold = 0, nnew = 0; idx < msa->nseq; idx++)
	{
	  if (new_fragassign[idx]) nnew++;
	  if (old_fragassign[idx]) nold++;
	}
      
      printf("%20s %10d %10d %10d %10d %10.4f %10.4f\n",
	     msa->name, msa->nseq, (int) msa->alen, nold, nnew,
	     (float) nold / (float) msa->nseq,
	     (float) nnew / (float) msa->nseq);
      
      nali++;
      esl_msa_Destroy(msa);  msa            = NULL;
      free(old_fragassign);  old_fragassign = NULL;
      free(new_fragassign);  new_fragassign = NULL;
    }
  if (nali == 0 || status != eslEOF)
    esl_msafile_ReadFailure(afp, status);

  esl_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR:
  return status;
}
#endif /* h4BUILD_EXPERIMENT2 */
/*--------------- end, experiment driver ------------------------*/
