/* Manipulate AFNI design matrices using the shared thd_xmat API. */
#include "mrilib.h"
#include "thd_xmat.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

static void usage(void)
{
   puts(
"1dXmatTool: manipulate AFNI X.xmat.1D design matrices\n"
"\n"
"Split a multi-run 3dDeconvolve design into run-local matrices:\n"
"  1dXmatTool -infile X.xmat.1D -uncensored X.nocensor.xmat.1D \\\n"
"             -split_runs -prefix split/X\n"
"\n"
"Operations (one at a time; -describe and -validate may be combined):\n"
"  -describe / -show  Full-matrix and per-run rank, DF, censor/stimulus summary.\n"
"  -validate         As above; exit 2 for failed checks, 0 for pass.\n"
"  -for lss|reml     With validation: check per-run LSS trial-model ranks, or\n"
"                    full/per-run GLM rank and residual degrees of freedom.\n"
"  -rank_tol T       Relative singular-value cutoff on unit-norm columns\n"
"                    (default 1e-7). Duplicate pairs use 1-abs(cosine)<=T.\n"
"  -json FILE        Optional JSON report with describe/validate.\n"
"  -select_columns LIST  Zero-based column indices, e.g. 0..3,8,10.\n"
"  -select_labels LIST   Exact comma-separated ColumnLabels.\n"
"  -select_stimuli LIST  Exact comma-separated StimLabels.\n"
"  -select_groups LIST   baseline (-1), nuisance (0), stimuli (>0), or group IDs.\n"
"                    Selections combine by union; original column order stays.\n"
"                    Writes PREFIX.xmat.1D, .rows.1D and .columns.1D.\n"
"  -drop_glts        Required for selection when input contains GLTs.\n"
"  -compare FILE     Match unique column labels/full TR indices; distinguish\n"
"                    reorder, censor/run/stimulus changes and changed values.\n"
"                    Exit 2 for differences, 0 for equivalence (reordering OK).\n"
"  -atol A -rtol R   Comparison tolerance A+R*max(abs(a),abs(b)); default 1e-6.\n"
"  -export_manifest  Write PREFIX.trials.1D from an IM design plus explicit events.\n"
"  -events FILE      Header/table: Run Volume Condition [Item]; complete unique\n"
"                    keys; Run one-based, Volume zero-based within each run.\n"
"  -subject ID       Subject label for the manifest.\n"
"  -beta FILE        Single-run beta bucket, or use -beta_pattern for many runs.\n"
"  -beta_pattern P   One {run} placeholder, replaced by 01,02,... (quote it).\n"
"  -manifest_format manifest|selected  Default manifest: whole InputFile plus\n"
"                    Volume for 3dMVPA -trial_manifest. selected: InputFile[n]\n"
"                    for direct 3dMVPA/3dRSA/3dSRM tables. Files need not yet exist.\n"
"                    Conditions are supplied by you, never guessed from IM labels.\n"
"Inputs/splitting:\n"
"  -infile FILE       NIML-comment matrix from 3dDeconvolve.\n"
"  -split_runs        Extract runs, preserving values and column order.\n"
"  -prefix PREFIX    Output prefix; its directory must already exist.\n"
"Optional:\n"
"  -uncensored FILE  Matching full-row design, recommended with censoring.\n"
"                    Determines column support before censoring. Dimensions,\n"
"                    labels, stimuli, run starts and retained values must match.\n"
"  -runs LIST        One-based run numbers, e.g. '1,3..5' (default: all).\n"
"  -zero_tol T       Relative support threshold, 0 <= T < 1 (default: 1e-6).\n"
"                    A column's peak within a run must exceed T times its\n"
"                    global peak. Use 0 for exact nonzero support.\n"
"  -help             This help.\n"
"\n"
"Outputs (r01 denotes the original, one-based run number):\n"
"  PREFIX.r01.xmat.1D    Run matrix with local GoodList, NRowFull, RunStart=0,\n"
"                       ColumnGroups and remapped stimulus/Basis metadata.\n"
"  PREFIX.r01.rows.1D   Row, SourceRow, RunTR, SourceTR (all zero-based).\n"
"  PREFIX.r01.columns.1D  Every input column, including dropped columns;\n"
"                       numeric columns only; labels in # ColumnLabel comments.\n"
"                       local Column/Stim/IMVolume are -1 when unavailable.\n"
"                       IMVolume is within the corresponding IM stimulus,\n"
"                       and is the LSS volume when there is one IM stimulus.\n"
"  PREFIX.runs.1D       Original run number, start TR, full length, retained\n"
"                       rows and columns; output paths are in # RunFiles comments.\n"
"Comment text escapes backslash, tab, newline and CR as \\\\, \\t, \\n and \\r.\n"
"All mappings refer to the immediate input file. Matrix headers also record\n"
"XmatSource* indices and the support file/tolerance. No files are overwritten.\n"
"All selected runs are checked before writing. Existing outputs cause failure.\n"
"\n"
"Policy:\n"
"* Requires explicit RunStart, NRowFull and GoodList. No timing is guessed.\n"
"* Drops columns inactive in the selected run; retains shared nuisance columns.\n"
"  Does not refit, normalize or demean. Values retain thd_xmat float precision;\n"
"  original double-precision text is not preserved beyond that precision.\n"
"* IM columns must belong to exactly one run. Zero/cross-run IM support, or\n"
"  loss of an IM column's support after censoring, is an error. Fully censored\n"
"  runs also fail. Without -uncensored, support is known only in retained rows.\n"
"* Rebuilds stimulus ranges; preserves column labels and per-stimulus Basis\n"
"  attributes. Drops GLTs and ColumnMasks (and records their removal).\n"
"* Extraction does not establish full rank or satisfy every downstream tool.\n"
"  3dLSS still requires exactly one IM stimulus and sufficient observations.\n"
"\n"
"Example downstream use (run 1, full-length per-run scaled EPI):\n"
"  3dLSS -matrix split/X.r01.xmat.1D \\\n"
"        -input pb03.SUBJ.r01.scale+orig -prefix LSS.r01\n"
"Use the matching run's FULL time series: 3dLSS applies the local GoodList.\n"
"Do not concatenate censored rows or feed the entire multi-run dataset.\n"
"This tool splits an existing design; it does not create new regressors.\n"
   );
}
/* Fill caller-owned arrays with unique one-based runs; reject trailing garbage. */
static int parse_run_values(const char *p,int limit,NI_int_array *a,int *seen)
{
   while(*p){
      char *end; long first,last,v; int step;
      while(isspace((unsigned char)*p)) p++;
      errno=0; first=strtol(p,&end,10);
      if(errno || p==end || first<1 || first>limit) return 0;
      p=end; last=first;
      if(p[0]=='.' && p[1]=='.'){
         p+=2; errno=0; last=strtol(p,&end,10);
         if(errno || p==end || last<1 || last>limit) return 0;
         p=end;
      }
      while(isspace((unsigned char)*p)) p++;
      if(*p && *p!=',') return 0;
      step=last>=first?1:-1;
      for(v=first;;v+=step){
         if(seen[v-1]) return 0;
         seen[v-1]=1; a->ar[a->num++]=(int)v;
         if(v==last) break;
      }
      if(*p==','){p++;if(!*p) return 0;}
   }
   return 1;
}

/* Own parsing storage; return a complete list or release it on invalid input. */
static NI_int_array *parse_runs(const char *p,int limit)
{
   NI_int_array *a=calloc(1,sizeof(*a));
   int *seen=calloc((size_t)limit,sizeof(int)),ok=0;
   if(a && seen){
      a->ar=malloc(sizeof(int)*(size_t)limit);
      if(a->ar) ok=parse_run_values(p,limit,a,seen);
   }
   free(seen);
   if(!ok){if(a)NI_delete_int_array(a);return NULL;}
   return a;
}
static char *outname(const char *prefix,int run,const char *suffix)
{
   size_t n=strlen(prefix)+strlen(suffix)+40; char *p=malloc(n);
   if(p) snprintf(p,n,"%s.r%02d.%s",prefix,run,suffix);
   return p;
}
static void comment_text(FILE *f,const char *p)
{
   for(;*p;p++) switch(*p){
      case '\\':fputs("\\\\",f);break; case '\t':fputs("\\t",f);break;
      case '\n':fputs("\\n",f);break; case '\r':fputs("\\r",f);break;
      default:fputc(*p,f);
   }
}
/* Additional operations: all tabular outputs use # Columns: and whitespace. */
typedef struct {
   int describe,validate,select,manifest,drop_glts;
   const char *compare,*json,*for_tool,*columns,*labels,*stimuli,*groups;
   const char *events,*subject,*beta,*beta_pattern,*manifest_format;
   double rank_tol,atol,rtol;
} XT_options;

static int xt_stim(const THD_xmat *a,int c)
{
   int s;for(s=0;s<a->nstim;s++)if(c>=a->stim[s].bot && c<=a->stim[s].top)return s;return -1;
}
static int xt_is_im(const THD_xmat *a,int s)
{return s>=0 && a->stim[s].option && !strcmp(a->stim[s].option,"-stim_times_IM");}
static FILE *xt_open(const char *path)
{
   int fd=open(path,O_WRONLY|O_CREAT|O_EXCL,0666);FILE *f;
   if(fd<0){fprintf(stderr,"** cannot create new output %s: %s\n",path,strerror(errno));return NULL;}
   f=fdopen(fd,"w");if(!f){close(fd);unlink(path);}return f;
}
static int xt_close(FILE *f,const char *path)
{
   int bad=ferror(f);if(fclose(f))bad=1;
   if(bad){unlink(path);fprintf(stderr,"** failed writing %s\n",path);return 0;}return 1;
}
static void xt_json_text(FILE *f,const char *s)
{
   const unsigned char *p=(const unsigned char *)(s?s:"");fputc('"',f);
   for(;*p;p++){
      if(*p=='"' || *p=='\\')fprintf(f,"\\%c",*p);
      else if(*p<32)fprintf(f,"\\u%04x",*p);else fputc(*p,f);
   }fputc('"',f);
}
static int xt_metadata(const THD_xmat *a,char *err,int nerr)
{
   const char *v=THD_xmat_attr(a,"Nstim");int s;
   if(v && atoi(v)!=a->nstim){snprintf(err,nerr,"Nstim does not match stimulus ranges");return 0;}
   for(s=0;s<a->nattr;s++){
      const char *k=a->attr_name[s],*p=strchr(k,'_');
      if(!strncasecmp(k,"Basis",5) && p && (atoi(p+1)<1 || atoi(p+1)>a->nstim)){
         snprintf(err,nerr,"Basis attribute refers to an unknown stimulus");return 0;
      }
   }
   v=THD_xmat_attr(a,"ColumnGroups");
   if(v){
      NI_int_array *g=NI_decode_int_list((char *)v,",");int valid=g && g->num==a->ncol;
      if(valid)for(s=0;s<a->ncol;s++){
         int st=xt_stim(a,s);
         if(g->ar[s]>0 && g->ar[s]!=st+1)valid=0;
      }
      if(g)NI_delete_int_array(g);
      if(!valid){snprintf(err,nerr,"ColumnGroups disagrees with columns/stimulus ranges");return 0;}
   }
   for(s=0;s<a->nstim;s++){
      int b,t,n=0;char k[64];snprintf(k,sizeof(k),"BasisColumns_%06d",s+1);v=THD_xmat_attr(a,k);
      if(v && (sscanf(v,"%d:%d%n",&b,&t,&n)!=2 || v[n] || b!=a->stim[s].bot || t!=a->stim[s].top)){
         snprintf(err,nerr,"%s does not match stimulus range",k);return 0;
      }
   }return 1;
}
static int xt_lss(const THD_xmat *a,double tol,char *err,int nerr)
{
   THD_xmat z=*a;THD_xmat_diagnostic d;int s,im=-1,c,r,j,k,nt,nn,bad=0;
   for(s=0;s<a->nstim;s++)if(xt_is_im(a,s)){if(im>=0){snprintf(err,nerr,"LSS requires exactly one IM family");return 0;}im=s;}
   if(im<0){snprintf(err,nerr,"LSS requires an IM family");return 0;}
   nt=a->stim[im].top-a->stim[im].bot+1;nn=a->ncol-nt;
   if(nt<2 || a->nrow<a->ncol+2){snprintf(err,nerr,"LSS needs >=2 IM columns and nrow>=ncol+2");return 0;}
   z.ncol=nn+2;z.nstim=0;z.stim=NULL;z.mat=calloc((size_t)z.nrow*z.ncol,sizeof(float));
   z.collab=calloc((size_t)z.ncol,sizeof(char *));
   if(!z.mat || !z.collab){free(z.mat);free(z.collab);return 0;}
   for(c=0;c<z.ncol;c++)z.collab[c]="LSS";
   for(c=0,k=0;c<a->ncol;c++)if(c<a->stim[im].bot || c>a->stim[im].top)
      memcpy(z.mat+(size_t)k++*z.nrow,a->mat+(size_t)c*a->nrow,sizeof(float)*a->nrow);
   for(j=0;j<nt;j++){
      for(r=0;r<a->nrow;r++){
         double sum=0;for(c=0;c<nt;c++)if(c!=j)sum+=a->mat[(size_t)(a->stim[im].bot+c)*a->nrow+r];
         z.mat[(size_t)nn*z.nrow+r]=(float)sum;
         z.mat[(size_t)(nn+1)*z.nrow+r]=a->mat[(size_t)(a->stim[im].bot+j)*a->nrow+r];
      }
      if(!THD_xmat_diagnostics(&z,tol,&d,err,nerr) || d.rank<z.ncol)bad++;
   }
   free(z.mat);free(z.collab);
   if(bad)snprintf(err,nerr,"%d LSS trial model(s) rank deficient at requested tolerance",bad);
   return bad==0;
}
static int xt_report(const THD_xmat *a,const THD_xmat *u,const XT_options *o,double zero_tol)
{
   THD_xmat_diagnostic d;FILE *json=NULL;char err[1024]="";int run,s,ok,failed=0;
   if(!xt_metadata(a,err,sizeof(err)) || !THD_xmat_diagnostics(a,o->rank_tol,&d,err,sizeof(err))){fprintf(stderr,"** %s\n",err);return 1;}
   if(o->json){json=xt_open(o->json);if(!json)return 1;}
   printf("# 1dXmatTool %s; rank uses unit-norm columns, tolerance %.9g\n",o->validate?"validation":"description",o->rank_tol);
   printf("# Matrix: %s\n# FullModel Rows=%d Columns=%d FullRows=%d Runs=%d Stimuli=%d Rank=%d DF=%d Condition=%.9g ZeroColumns=%d DuplicatePairs=%d MaxAbsCosine=%.9g\n",
      a->fname,a->nrow,a->ncol,a->nfull,a->nrun,a->nstim,d.rank,d.residual_df,d.condition,d.zero_columns,d.duplicate_pairs,d.max_abs_cosine);
   for(s=0;s<a->nstim;s++)printf("# Stimulus %d %s columns=%d..%d option=%s\n",s,a->stim[s].label,a->stim[s].bot,a->stim[s].top,a->stim[s].option?a->stim[s].option:"none");
   if(!u && a->nrow<a->nfull)printf("# Support is inferred only from retained rows; supply -uncensored to detect support hidden by censoring.\n");
   if(!o->for_tool || strcmp(o->for_tool,"lss"))if(d.rank<a->ncol || d.residual_df<1)failed=1;
   if(o->for_tool && !strcmp(o->for_tool,"reml") && a->nrow<9)failed=1;
   if(json){
      fprintf(json,"{\n  \"input\": ");xt_json_text(json,a->fname);
      fprintf(json,",\n  \"rows\": %d, \"columns\": %d, \"full_rows\": %d, \"rank\": %d, \"residual_df\": %d, \"zero_columns\": %d, \"duplicate_pairs\": %d, \"condition\": ",a->nrow,a->ncol,a->nfull,d.rank,d.residual_df,d.zero_columns,d.duplicate_pairs);
      if(isfinite(d.condition))fprintf(json,"%.17g",d.condition);else fputs("null",json);
      fprintf(json,",\n  \"rank_tolerance\": %.17g, \"runs\": [\n",o->rank_tol);
   }
   puts("# Columns: Run Start Length NKept NColumns Rank ResidualDF Condition ZeroColumns DuplicatePairs Valid");
   for(run=0;run<a->nrun;run++){
      THD_xmat_run_map *map=NULL;THD_xmat *x=THD_xmat_extract_run(a,u,run,zero_tol,&map,err,sizeof(err));
      THD_xmat_diagnostic rd={0};int nr=0,nc=0;ok=x!=NULL;
      if(ok){nr=x->nrow;nc=x->ncol;ok=THD_xmat_diagnostics(x,o->rank_tol,&rd,err,sizeof(err));}
      if(ok && o->for_tool && !strcmp(o->for_tool,"lss"))ok=xt_lss(x,o->rank_tol,err,sizeof(err));
      else if(ok && (rd.rank<nc || rd.residual_df<1)){ok=0;snprintf(err,sizeof(err),"rank deficient or no residual degrees of freedom");}
      if(ok && o->for_tool && !strcmp(o->for_tool,"reml") && x->nrow<9){ok=0;snprintf(err,sizeof(err),"3dREMLfit requires at least 9 rows");}
      if(!ok){failed=1;printf("# RunIssue %d %s\n",run+1,err);}
      printf("%d %d %d %d %d %d %d %.9g %d %d %d\n",run+1,a->runstart[run],THD_xmat_run_length(a,run),nr,nc,x?rd.rank:-1,x?rd.residual_df:-1,x?rd.condition:NAN,rd.zero_columns,rd.duplicate_pairs,ok);
      if(json){
         fprintf(json,"    %s{\"run\": %d, \"start\": %d, \"length\": %d, \"retained\": %d, \"columns\": %d, \"rank\": %d, \"residual_df\": %d, \"valid\": %s, \"issue\": ",run?",":"",run+1,a->runstart[run],THD_xmat_run_length(a,run),nr,nc,x?rd.rank:-1,x?rd.residual_df:-1,ok?"true":"false");
         xt_json_text(json,ok?"":err);fputs("}\n",json);
      }
      THD_xmat_free(x);THD_xmat_run_map_free(map);
   }
   if(o->for_tool)printf("# %s preflight checks matrix structure/estimability; dataset, spatial and runtime checks remain with the fitting program.\n",o->for_tool);
   if(json){fprintf(json,"  ], \"valid\": %s\n}\n",failed?"false":"true");if(!xt_close(json,o->json))return 1;}
   if(o->validate)printf("# Validation %s\n",failed?"FAIL":"PASS");
   return o->validate && failed?2:0;
}
static int xt_compare(const THD_xmat *a,const XT_options *o)
{
   THD_xmat *b;NI_int_array *ga=NULL,*gb=NULL;const char *v;
   char err[1024];int c,j,r,t,missing=0,reorder=0,changed=0,runs=1,stim=1,censor=1,labels=1,groups=1;double maximum=0;
   b=THD_xmat_read(o->compare,err,sizeof(err));if(!b){fprintf(stderr,"** %s\n",err);return 1;}
   if(!xt_metadata(a,err,sizeof(err)) || !xt_metadata(b,err,sizeof(err))){fprintf(stderr,"** %s\n",err);THD_xmat_free(b);return 1;}
   if(a->nfull!=b->nfull || a->tr!=b->tr || a->nrun!=b->nrun)runs=0;
   if(runs)for(r=0;r<a->nrun;r++)if(a->runstart[r]!=b->runstart[r])runs=0;
   if(a->nrow!=b->nrow)censor=0;
   if(censor)for(r=0;r<a->nrow;r++)if(a->good[r]!=b->good[r])censor=0;
   /* Unique labels are necessary for unambiguous identity matching. */
   for(c=0;c<a->ncol;c++)for(j=0;j<c;j++)if(!strcmp(a->collab[c],a->collab[j]))labels=0;
   for(c=0;c<b->ncol;c++)for(j=0;j<c;j++)if(!strcmp(b->collab[c],b->collab[j]))labels=0;
   if(!labels){fprintf(stderr,"** compare requires unique ColumnLabels in both matrices\n");THD_xmat_free(b);return 1;}
   v=THD_xmat_attr(a,"ColumnGroups");if(v)ga=NI_decode_int_list((char *)v,",");
   v=THD_xmat_attr(b,"ColumnGroups");if(v)gb=NI_decode_int_list((char *)v,",");
   if((ga==NULL)!=(gb==NULL))groups=0;
   puts("# Match columns by label and rows by full-timeline GoodList index.");
   puts("# Columns: Column OtherColumn CommonRows ChangedValues MaxAbsDifference");
   for(c=0;c<a->ncol;c++){
      int diff=0,common=0,sa,sb;double md=0;
      for(j=0;j<b->ncol;j++)if(!strcmp(a->collab[c],b->collab[j]))break;
      if(j==b->ncol){missing++;printf("%d -1 0 0 nan\n",c);continue;}
      if(c!=j)reorder=1;sa=xt_stim(a,c);sb=xt_stim(b,j);
      if(ga && gb && (ga->ar[c]<=0 || gb->ar[j]<=0) && ga->ar[c]!=gb->ar[j])groups=0;
      if((sa<0)!=(sb<0))stim=0;
      else if(sa>=0 && (strcmp(a->stim[sa].label,b->stim[sb].label) || strcmp(a->stim[sa].option?a->stim[sa].option:"",b->stim[sb].option?b->stim[sb].option:"") || strcmp(a->stim[sa].formula?a->stim[sa].formula:"",b->stim[sb].formula?b->stim[sb].formula:"")))stim=0;
      for(r=0,t=0;r<a->nrow && t<b->nrow;){
         if(a->good[r]<b->good[t]){r++;continue;}if(a->good[r]>b->good[t]){t++;continue;}
         {double av=a->mat[(size_t)c*a->nrow+r],bv=b->mat[(size_t)j*b->nrow+t],d=fabs(av-bv);
          if(!isfinite(av) || !isfinite(bv) || d>o->atol+o->rtol*MAX(fabs(av),fabs(bv)))diff++;
          if(d>md)md=d;common++;}r++;t++;
      }
      changed+=diff;if(md>maximum)maximum=md;printf("%d %d %d %d %.17g\n",c,j,common,diff,md);
   }
   if(a->ncol!=b->ncol)missing++;
   printf("# SameColumns %d\n# ColumnOrderChanged %d\n# SameCensor %d\n# SameRunsAndTR %d\n# SameStimulusMetadata %d\n# ChangedValues %d\n# MaxAbsDifference %.17g\n",missing==0,reorder,censor,runs,stim,changed,maximum);
   printf("# SameColumnGroups %d\n",groups);
   if(ga)NI_delete_int_array(ga);if(gb)NI_delete_int_array(gb);
   THD_xmat_free(b);return missing || !censor || !runs || !stim || !groups || changed?2:0;
}
static int xt_mark_names(const THD_xmat *a,const char *list,int stimuli,unsigned char *keep)
{
   NI_str_array *v=NI_decode_string_list((char *)list,",");int i,c,s,found;
   if(!v || !v->num){if(v)NI_delete_str_array(v);return 0;}
   for(i=0;i<v->num;i++){
      found=0;for(c=0;c<a->ncol;c++){
         s=xt_stim(a,c);
         if((!stimuli && !strcmp(v->str[i],a->collab[c])) || (stimuli && s>=0 && !strcmp(v->str[i],a->stim[s].label))){keep[c]=1;found=1;}
      }
      if(!found){fprintf(stderr,"** unknown %s label: %s\n",stimuli?"stimulus":"column",v->str[i]);NI_delete_str_array(v);return 0;}
   }NI_delete_str_array(v);return 1;
}
/* Translate zero-based column ranges to the strict one-based run parser.
   On success, add columns to the caller's mask; always release parsing storage. */
static int xt_mark_columns(const THD_xmat *a,const char *columns,
                           unsigned char *keep,char *err,int nerr)
{
   const char *q=columns;
   NI_int_array *v=NULL;char *one=NULL;int i,ok=1;
   while(*q && ok){
      char *end;long n;errno=0;n=strtol(q,&end,10);
      if(errno || end==q || n<0 || n>=a->ncol){
         snprintf(err,nerr,"invalid zero-based -select_columns");ok=0;
      }else{
         one=THD_zzprintf(one,"%ld",n+1);q=end;
         if(q[0]=='.' && q[1]=='.'){one=THD_zzprintf(one,"..");q+=2;}
         else if(*q==','){one=THD_zzprintf(one,",");q++;}
         else if(*q){snprintf(err,nerr,"invalid -select_columns");ok=0;}
      }
   }
   if(ok){
      v=one?parse_runs(one,a->ncol):NULL;
      if(!v || !v->num){snprintf(err,nerr,"invalid -select_columns");ok=0;}
      else for(i=0;i<v->num;i++)keep[v->ar[i]-1]=1;
   }
   if(v)NI_delete_int_array(v);
   free(one);return ok;
}

/* Resolve group names/numbers into the caller's mask, releasing decoded lists. */
static int xt_mark_groups(const THD_xmat *a,const char *list,
                          unsigned char *keep,char *err,int nerr)
{
   const char *v=THD_xmat_attr(a,"ColumnGroups");
   NI_int_array *g=v?NI_decode_int_list((char *)v,","):NULL;
   NI_str_array *groups=NI_decode_string_list((char *)list,",");
   int i,c,ok=g && g->num==a->ncol && groups;
   if(!ok)snprintf(err,nerr,"group selection needs valid ColumnGroups");
   for(i=0;ok && i<groups->num;i++){
      char *end;long n=strtol(groups->str[i],&end,10);
      int all=!strcmp(groups->str[i],"stimuli"),found=0;
      if(!strcmp(groups->str[i],"baseline"))n=-1;
      else if(!strcmp(groups->str[i],"nuisance"))n=0;
      else if(!all && (end==groups->str[i] || *end)){
         snprintf(err,nerr,"unknown group: %s",groups->str[i]);ok=0;
      }
      if(ok){
         for(c=0;c<a->ncol;c++)if(all?g->ar[c]>0:g->ar[c]==n){keep[c]=1;found=1;}
         if(!found){snprintf(err,nerr,"group selects no columns: %s",groups->str[i]);ok=0;}
      }
   }
   if(g)NI_delete_int_array(g);
   if(groups)NI_delete_str_array(groups);
   return ok;
}

/* Select and write a matrix plus provenance; roll back the matrix if maps fail. */
static int xt_select(const THD_xmat *a,const XT_options *o,const char *prefix)
{
   unsigned char *keep=calloc((size_t)a->ncol,1);THD_xmat *x=NULL;THD_xmat_run_map *m=NULL;
   char err[1024]="",*paths[3]={NULL,NULL,NULL};int i,ok,wrote=0;
   if(!keep || !prefix){fprintf(stderr,"** selection requires -prefix\n");free(keep);return 1;}
   ok=xt_metadata(a,err,sizeof(err));
   if(ok && o->columns)ok=xt_mark_columns(a,o->columns,keep,err,sizeof(err));
   if(ok && o->labels)ok=xt_mark_names(a,o->labels,0,keep);
   if(ok && o->stimuli)ok=xt_mark_names(a,o->stimuli,1,keep);
   if(ok && o->groups)ok=xt_mark_groups(a,o->groups,keep,err,sizeof(err));
   if(ok){x=THD_xmat_select_columns(a,keep,o->drop_glts,&m,err,sizeof(err));ok=x!=NULL;}
   if(ok){
      paths[0]=THD_zzprintf(NULL,"%s.xmat.1D",prefix);paths[1]=THD_zzprintf(NULL,"%s.rows.1D",prefix);paths[2]=THD_zzprintf(NULL,"%s.columns.1D",prefix);
      wrote=THD_xmat_write(x,paths[0],err,sizeof(err));
      ok=wrote && THD_xmat_write_run_map(a,x,m,paths[1],paths[2],err,sizeof(err));
   }
   if(ok)printf("# Selected %d/%d columns -> %s\n",x->ncol,a->ncol,paths[0]);
   else{if(wrote)unlink(paths[0]);fprintf(stderr,"** selection: %s\n",*err?err:"invalid selection");}
   for(i=0;i<3;i++)free(paths[i]);
   free(keep);THD_xmat_free(x);THD_xmat_run_map_free(m);return !ok;
}
static char *xt_beta_name(const XT_options *o,int run)
{
   const char *p=o->beta_pattern,*at;
   if(o->beta)return strdup(o->beta);
   if(!p || !(at=strstr(p,"{run}")) || strstr(at+5,"{run}"))return NULL;
   return THD_zzprintf(NULL,"%.*s%02d%s",(int)(at-p),p,run,at+5);
}
static int xt_integer(const char *s,int *n)
{
   char *e;long v;errno=0;v=strtol(s,&e,10);if(errno || e==s || *e || v<0 || v>INT_MAX)return 0;*n=(int)v;return 1;
}
/* Validate all event keys before opening a manifest output file. */
static int xt_events_valid(const THD_xmat *a,const THD_datatable *events,
                           int ri,int vi,int ci,char *err,int nerr)
{
   int r;
   if(ri<0 || vi<0 || ci<0){
      snprintf(err,nerr,"events need Run Volume Condition (optional Item)");return 0;
   }
   for(r=0;r<events->nrow;r++){
      int rr,vv,j,rr2,vv2;
      if(!xt_integer(DT_CELL(events,r,ri),&rr) || rr<1 || rr>a->nrun || !xt_integer(DT_CELL(events,r,vi),&vv) || !*DT_CELL(events,r,ci)){
         snprintf(err,nerr,"invalid events Run/Volume/Condition on row %d",r+1);return 0;}
      for(j=0;j<r;j++){xt_integer(DT_CELL(events,j,ri),&rr2);xt_integer(DT_CELL(events,j,vi),&vv2);
         if(rr==rr2 && vv==vv2){snprintf(err,nerr,"duplicate Run/Volume in events");return 0;}}
   }
   return 1;
}

/* Write one run using borrowed resources, marking matched rows in seen. */
static int xt_manifest_run(const THD_xmat *a,const THD_xmat *x,
   const THD_xmat_run_map *m,const XT_options *o,const THD_datatable *events,
   int ri,int vi,int ci,int ii,int im,int run,const char *beta,FILE *fp,
   int *seen,int *total,char *err,int nerr)
{
   int c,r,local=0;
   for(c=0;c<x->ncol;c++)if(m->columns[c]>=a->stim[im].bot && m->columns[c]<=a->stim[im].top){
      for(r=0;r<events->nrow;r++){int rr,vv;xt_integer(DT_CELL(events,r,ri),&rr);xt_integer(DT_CELL(events,r,vi),&vv);if(rr==run+1 && vv==local)break;}
      if(r==events->nrow){snprintf(err,nerr,"missing event for run %d volume %d",run+1,local);return 0;}
      seen[r]=1;(*total)++;
      {char *sample=THD_zzprintf(NULL,"%s_r%02d_v%04d",o->subject,run+1,local),*sel=NULL;
       THD_1D_print_token(fp,sample);fputc(' ',fp);THD_1D_print_token(fp,o->subject);
       fprintf(fp," %d ",run+1);THD_1D_print_token(fp,ii>=0?DT_CELL(events,r,ii):a->collab[m->columns[c]]);
       fputc(' ',fp);THD_1D_print_token(fp,DT_CELL(events,r,ci));fprintf(fp," %d ",local);THD_1D_print_token(fp,a->fname);
       fprintf(fp," %d ",m->columns[c]);
       if(o->manifest_format && !strcmp(o->manifest_format,"selected"))sel=THD_zzprintf(NULL,"%s[%d]",beta,local);
       THD_1D_print_token(fp,sel?sel:beta);fputc('\n',fp);free(sel);free(sample);}
      local++;
   }
   return 1;
}

/* Own event, run and output resources; remove an incomplete manifest on failure. */
static int xt_manifest(const THD_xmat *a,const THD_xmat *u,const XT_options *o,const char *prefix,double tol)
{
   THD_datatable *events=NULL;FILE *fp=NULL;int ri,vi,ci,ii,run,r,im=-1,ok=0,total=0;
   int *seen=NULL;char *path=NULL,err[1024]="";
   if(!prefix || !o->events || !o->subject || (!o->beta && !o->beta_pattern) || (o->beta && o->beta_pattern)){
      fprintf(stderr,"** export needs -prefix -events -subject and exactly one of -beta/-beta_pattern\n");return 1;
   }
   if(o->beta && a->nrun!=1){fprintf(stderr,"** -beta requires a single-run matrix; use -beta_pattern for per-run buckets\n");return 1;}
   for(r=0;r<a->nstim;r++)if(xt_is_im(a,r)){if(im>=0){fprintf(stderr,"** export requires exactly one IM family\n");return 1;}im=r;}
   if(im<0){fprintf(stderr,"** export requires an IM family\n");return 1;}
   events=THD_read_datatable_file((char *)o->events);
   ri=THD_datatable_column(events,"Run");vi=THD_datatable_column(events,"Volume");ci=THD_datatable_column(events,"Condition");ii=THD_datatable_column(events,"Item");
   ok=xt_events_valid(a,events,ri,vi,ci,err,sizeof(err));
   if(ok){
      seen=calloc((size_t)events->nrow,sizeof(int));
      if(!seen){snprintf(err,sizeof(err),"out of memory");ok=0;}
   }
   if(ok){path=THD_zzprintf(NULL,"%s.trials.1D",prefix);fp=xt_open(path);ok=fp!=NULL;}
   if(ok){
      fprintf(fp,"# 1dXmatTool trial manifest; Run is one-based, Volume/XColumn are zero-based.\n");
      fprintf(fp,"# Columns: Sample Subj Run Item Condition Volume %s InputFile\n",
         o->manifest_format && !strcmp(o->manifest_format,"selected")?"SourceXmat SourceColumn":"Xmat XColumn");
   }
   for(run=0;ok && run<a->nrun;run++){
      THD_xmat_run_map *m=NULL;
      THD_xmat *x=THD_xmat_extract_run(a,u,run,tol,&m,err,sizeof(err));
      char *beta=xt_beta_name(o,run+1);
      if(!x || !beta){
         if(!*err)snprintf(err,sizeof(err),"beta pattern must contain exactly one {run}");
         ok=0;
      }else ok=xt_manifest_run(a,x,m,o,events,ri,vi,ci,ii,im,run,beta,fp,
                               seen,&total,err,sizeof(err));
      THD_xmat_free(x);THD_xmat_run_map_free(m);free(beta);
   }
   for(r=0;ok && r<events->nrow;r++)if(!seen[r]){snprintf(err,sizeof(err),"events contain unmatched rows");ok=0;}
   if(ok){ok=xt_close(fp,path);fp=NULL;}
   if(ok)printf("# Exported %d trials -> %s\n",total,path);
   if(fp){fclose(fp);unlink(path);}if(!ok)fprintf(stderr,"** manifest: %s\n",*err?err:"cannot write output");
   free(path);free(seen);THD_free_datatable(events);return !ok;
}

/* Validate operation combinations before acquiring input or output resources. */
static const char *xt_options_error(const XT_options *ops,const char *in,
   const char *unc,const char *prefix,const char *runlist,int split)
{
   if(!in) return "requires -infile (see -help)";
   if(split+ops->select+ops->manifest+(ops->compare!=NULL)+((ops->describe||ops->validate)!=0)!=1)
      return "choose one operation: split, describe/validate, select, compare, or export_manifest";
   if(ops->for_tool && strcmp(ops->for_tool,"lss") && strcmp(ops->for_tool,"reml"))return "-for must be lss or reml";
   if(ops->manifest_format && strcmp(ops->manifest_format,"selected") && strcmp(ops->manifest_format,"manifest"))return "-manifest_format must be manifest or selected";
   if(ops->json && !(ops->describe||ops->validate))return "-json is available with describe/validate";
   if(runlist && !split)return "-runs applies to -split_runs";
   if(prefix && !*prefix)return "-prefix cannot be empty";
   if(ops->for_tool && !(ops->describe||ops->validate))return "-for applies to describe/validate";
   if(ops->drop_glts && !ops->select)return "-drop_glts applies to column selection";
   if((ops->events || ops->subject || ops->beta || ops->beta_pattern || ops->manifest_format) && !ops->manifest)
      return "event/beta options apply to -export_manifest";
   if(unc && (ops->select || ops->compare))return "-uncensored applies to split, describe/validate and export_manifest";
   return NULL;
}

/* Resources owned by one split operation, including only outputs we created. */
typedef struct {
   THD_xmat **x;
   THD_xmat_run_map **maps;
   NI_int_array *runs;
   char **paths,*summary;
   FILE *fp;
   int n,done,summary_owned;
} XT_split;

static int xt_error(char *err,int nerr,const char *message)
{
   snprintf(err,nerr,"%s",message);return 1;
}

/* Prepare and write run subsets. The caller releases resources and rolls back
   completed files on failure; done tracks precisely which paths we own. */
static int xt_split_write(const THD_xmat *a,const THD_xmat *u,
   const char *prefix,const char *runlist,double tol,XT_split *out,
   char *err,int nerr)
{
   int i,j,fd;struct stat st;
   if(!prefix || !*prefix)return xt_error(err,nerr,"split requires -prefix");
   if(runlist){
      out->runs=parse_runs(runlist,a->nrun);
      if(!out->runs || out->runs->num<1) return xt_error(err,nerr,"invalid -runs list: use unique run numbers in range, e.g. 1,3..5");
   }else{
      out->runs=calloc(1,sizeof(*out->runs));
      if(!out->runs)return xt_error(err,nerr,"out of memory");
      out->runs->num=a->nrun;
      out->runs->ar=malloc(sizeof(int)*(size_t)a->nrun);
      if(!out->runs->ar)return xt_error(err,nerr,"out of memory");
      for(i=0;i<a->nrun;i++) out->runs->ar[i]=i+1;
   }
   out->n=out->runs->num;
   out->x=calloc((size_t)out->n,sizeof(*out->x)); out->maps=calloc((size_t)out->n,sizeof(*out->maps));
   out->paths=calloc((size_t)out->n*3,sizeof(*out->paths)); out->summary=malloc(strlen(prefix)+16);
   if(!out->x || !out->maps || !out->paths || !out->summary) return xt_error(err,nerr,"out of memory");
   sprintf(out->summary,"%s.runs.1D",prefix);
   if(!u && a->nrow<a->nfull)
      fprintf(stderr,"++ Support inferred from retained rows only; use -uncensored for pre-censor support.\n");
   for(i=0;i<out->n;i++){
      if(out->runs->ar[i]<1 || out->runs->ar[i]>a->nrun) return xt_error(err,nerr,"-runs index outside the input run range");
      for(j=0;j<i;j++) if(out->runs->ar[i]==out->runs->ar[j]) return xt_error(err,nerr,"duplicate run in -runs");
      out->x[i]=THD_xmat_extract_run(a,u,out->runs->ar[i]-1,tol,out->maps+i,err,nerr);
      if(!out->x[i]) return 1;
      out->paths[3*i]=outname(prefix,out->runs->ar[i],"xmat.1D");
      out->paths[3*i+1]=outname(prefix,out->runs->ar[i],"rows.1D");
      out->paths[3*i+2]=outname(prefix,out->runs->ar[i],"columns.1D");
      if(!out->paths[3*i] || !out->paths[3*i+1] || !out->paths[3*i+2]) return xt_error(err,nerr,"out of memory");
   }
   for(i=0;i<=3*out->n;i++){
      const char *p=i==3*out->n?out->summary:out->paths[i];
      if(!lstat(p,&st)){snprintf(err,nerr,"output already exists: %s",p);return 1;}
      if(errno!=ENOENT){snprintf(err,nerr,"cannot access output path: %s",p);return 1;}
   }
   fd=open(out->summary,O_WRONLY|O_CREAT|O_EXCL,0666);
   if(fd<0) return xt_error(err,nerr,"cannot create summary (prefix directory must exist)");
   out->summary_owned=1; out->fp=fdopen(fd,"w");
   if(!out->fp){close(fd);return xt_error(err,nerr,"cannot open summary stream");}
   fputs("# Columns: Run Start Length NKept NColumns\n",out->fp);
   for(i=0;i<out->n;i++){
      if(!THD_xmat_write(out->x[i],out->paths[3*i],err,nerr)) return 1;
      out->done=3*i+1;
      if(!THD_xmat_write_run_map(a,out->x[i],out->maps[i],out->paths[3*i+1],out->paths[3*i+2],err,nerr)) return 1;
      out->done=3*i+3;
      fprintf(out->fp,"# RunFiles %d ",out->runs->ar[i]);
      comment_text(out->fp,out->paths[3*i]); fputc(' ',out->fp);comment_text(out->fp,out->paths[3*i+1]);fputc(' ',out->fp);comment_text(out->fp,out->paths[3*i+2]);fputc('\n',out->fp);
      fprintf(out->fp,"%d %d %d %d %d\n",out->runs->ar[i],out->maps[i]->start,out->x[i]->nfull,out->x[i]->nrow,out->x[i]->ncol);
   }
   j=ferror(out->fp); if(fclose(out->fp)) j=1; out->fp=NULL;
   if(j) return xt_error(err,nerr,"failed writing summary");
   for(i=0;i<out->n;i++) fprintf(stderr,"++ run %d: %d/%d rows, %d/%d columns -> %s\n",
      out->runs->ar[i],out->x[i]->nrow,out->x[i]->nfull,out->x[i]->ncol,a->ncol,out->paths[3*i]);
   fprintf(stderr,"++ Wrote %d run(s); provenance: %s\n",out->n,out->summary); return 0;
}

/* Release split storage once, removing only newly created outputs on failure. */
static int xt_split_runs(const THD_xmat *a,const THD_xmat *u,
   const char *prefix,const char *runlist,double tol)
{
   XT_split out={0};char err[1024]="";int i;
   int rc=xt_split_write(a,u,prefix,runlist,tol,&out,err,sizeof(err));
   if(out.fp)fclose(out.fp);
   if(rc){
      fprintf(stderr,"** 1dXmatTool: %s\n",err);
      for(i=0;i<out.done;i++)unlink(out.paths[i]);
      if(out.summary_owned)unlink(out.summary);
   }
   for(i=0;i<out.n;i++){
      if(out.x)THD_xmat_free(out.x[i]);
      if(out.maps)THD_xmat_run_map_free(out.maps[i]);
   }
   if(out.paths){for(i=0;i<3*out.n;i++)free(out.paths[i]);free(out.paths);}
   free(out.x);free(out.maps);free(out.summary);
   if(out.runs)NI_delete_int_array(out.runs);
   return rc;
}

int main(int argc,char **argv)
{
   XT_options ops={0};
   const char *in=NULL,*unc=NULL,*prefix=NULL,*runlist=NULL;
   THD_xmat *a=NULL,*u=NULL;char err[1024]="";const char *problem;
   double tol=THD_XMAT_SPLIT_TOL;int split=0,i,rc=1;
   ops.rank_tol=1.e-7;ops.atol=1.e-6;ops.rtol=1.e-6;
   if(argc==1){usage();return 0;}
   for(i=1;i<argc;i++){
      if(!strcmp(argv[i],"-help") || !strcmp(argv[i],"--help")){usage();return 0;}
      if(!strcmp(argv[i],"-describe") || !strcmp(argv[i],"-show")){ops.describe=1;continue;}
      if(!strcmp(argv[i],"-validate")){ops.validate=1;continue;}
      if(!strcmp(argv[i],"-export_manifest")){ops.manifest=1;continue;}
      if(!strcmp(argv[i],"-drop_glts")){ops.drop_glts=1;continue;}
      if(!strcmp(argv[i],"-split_runs")){split=1;continue;}
      if(i+1>=argc){fprintf(stderr,"** missing argument after %s\n",argv[i]);return 1;}
      if(!strcmp(argv[i],"-infile")) in=argv[++i];
      else if(!strcmp(argv[i],"-uncensored")) unc=argv[++i];
      else if(!strcmp(argv[i],"-prefix")) prefix=argv[++i];
      else if(!strcmp(argv[i],"-runs")) runlist=argv[++i];
      else if(!strcmp(argv[i],"-compare")) ops.compare=argv[++i];
      else if(!strcmp(argv[i],"-json")) ops.json=argv[++i];
      else if(!strcmp(argv[i],"-for")) ops.for_tool=argv[++i];
      else if(!strcmp(argv[i],"-select_columns")) ops.columns=argv[++i];
      else if(!strcmp(argv[i],"-select_labels")) ops.labels=argv[++i];
      else if(!strcmp(argv[i],"-select_stimuli")) ops.stimuli=argv[++i];
      else if(!strcmp(argv[i],"-select_groups")) ops.groups=argv[++i];
      else if(!strcmp(argv[i],"-events")) ops.events=argv[++i];
      else if(!strcmp(argv[i],"-subject")) ops.subject=argv[++i];
      else if(!strcmp(argv[i],"-beta")) ops.beta=argv[++i];
      else if(!strcmp(argv[i],"-beta_pattern")) ops.beta_pattern=argv[++i];
      else if(!strcmp(argv[i],"-manifest_format")) ops.manifest_format=argv[++i];
      else if(!strcmp(argv[i],"-rank_tol") || !strcmp(argv[i],"-atol") || !strcmp(argv[i],"-rtol")){
         const char *opt=argv[i],*v=argv[++i];char *end;double d;errno=0;d=strtod(v,&end);
         if(errno || end==v || *end || !isfinite(d) || d<0 || (!strcmp(opt,"-rank_tol") && (d<=0 || d>=1))){fprintf(stderr,"** invalid tolerance\n");return 1;}
         if(!strcmp(opt,"-rank_tol"))ops.rank_tol=d;else if(!strcmp(opt,"-atol"))ops.atol=d;else ops.rtol=d;
      }
      else if(!strcmp(argv[i],"-zero_tol")){
         char *end; const char *v=argv[++i]; errno=0; tol=strtod(v,&end);
         if(errno || end==v || *end || !isfinite(tol) || tol<0 || tol>=1){
            fprintf(stderr,"** -zero_tol must be a number in [0,1)\n");return 1;
         }
      }else{fprintf(stderr,"** unknown option: %s\n",argv[i]);return 1;}
   }
   ops.select=ops.columns || ops.labels || ops.stimuli || ops.groups;
   problem=xt_options_error(&ops,in,unc,prefix,runlist,split);
   if(problem){fprintf(stderr,"** 1dXmatTool: %s\n",problem);return 1;}
   a=THD_xmat_read(in,err,sizeof(err));
   if(a && unc)u=THD_xmat_read(unc,err,sizeof(err));
   if(!a || (unc && !u))fprintf(stderr,"** 1dXmatTool: %s\n",err);
   else if(split)rc=xt_split_runs(a,u,prefix,runlist,tol);
   else if(ops.describe||ops.validate)rc=xt_report(a,u,&ops,tol);
   else if(ops.compare)rc=xt_compare(a,&ops);
   else if(ops.select)rc=xt_select(a,&ops,prefix);
   else if(ops.manifest)rc=xt_manifest(a,u,&ops,prefix,tol);
   THD_xmat_free(a);THD_xmat_free(u);
   return rc;
}
