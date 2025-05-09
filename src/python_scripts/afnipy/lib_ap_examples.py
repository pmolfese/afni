#!/usr/bin/env python

from afnipy import afni_util          as UTIL
from afnipy import lib_format_cmd_str as lfcs
import copy

# ----------------------------------------------------------------------
# This is a library for storing basic information regarding options and
# examples for afni_proc.py.
#
#    - examples are stored as individual dictionaries
#    - examples from afni_proc.py -help are created by:
#       - set noglob
#       - add '-optlist_show_argv_array dict' to command
#
#    - for a given example, want:
#       - name (e.g. Example 11) - want case insensitive checking
#       - source (from help, from AFNI_data6, etc.)
#       - description
#       - header  (before the example in verbose -help output)
#       - trailer (after  the example in verbose -help output)
#       - keywords = []; e.g. 'obsolete', 'complete', 'registration'. 'ME',
#                             'surface', 'rest', 'task', 'physio'
#                             ('noshow' : default to not showing with help)
#
# --------------------------
# preferred option ordering: 
#
#   - id and script specs
#   - EPI data and specs, including blip
#   - anat data and specs (followers)
#   - non-block data inputs
#   - block order (-blocks)
#   - general block options (e.g. -radial_correlate_blocks)
#   - block ordered options (e.g. -blur_size)
#      -surf dests go in surf block section
#      -blip dsets do not go in blip block section (as not user specified)
#   -regress order:
#       ROI (regs of interest), opts_3dD, polort, bandpass,
#       motion, AICOR/RONI (regs of no interest), censoring,
#       non-reg opts (3dD_stop, reml_exec, then post-reg options)

# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# main array of APExample instances
# (will rename later to _old and examples)
ap_examples = []
ap_examples_new = []

# ----------------------------------------------------------------------
# class definition for instances in ap_examples array
class APExample:
   def __init__(self, name, olist, source='', descrip='', moddate='UNKNOWN',
                header='', trailer='', keywords=[]):
      self.name     = name          # used to reference example
      self.source   = source        # from AP help, AD6, etc.
      self.descrip  = descrip       # very short description
      self.moddate  = moddate       # most recent modification date
      self.header   = header        # shown before example (in -help)
      self.trailer  = trailer       # shown after example (in -help)
      self.keywords = keywords      # list:

      self.keys     = []            # convenience: olist[][0] entries
      self.olist    = olist         # list of options [opt, [params]]
      self.aphelp   = (source == 'afni_proc.py -help')

      if not self.valid_olist():
         return None

      self.keys = [o[0] for o in olist]

   def valid_olist(self):
      valid = 1
      for oind, opt in enumerate(self.olist):
         if len(opt) != 2:
            print("** APExample '%s' olist #%d needs opt and list, have:\n" \
                  "   %s" % (self.name, oind, opt))
            valid = 0
            continue
         if type(opt[0]) != str:
            print("** APExample '%s' olist #%d should start with string:\n" \
                  "   %s" % (self.name, oind, opt))
            valid = 0
         if type(opt[1]) != list:
            print("** APExample '%s' olist #%d should have a list:\n" \
                  "   %s" % (self.name, oind, opt))
            valid = 0
      return valid

   def copy(self, quotize=0):
      """return a deep copy
            quotize - if set, quotize any needed option parameters
      """
      cc = copy.deepcopy(self)

      if quotize:
         for entry in cc.olist:
            entry[1] = UTIL.quotize_list(entry[1])

      return cc

   def compare(self, target, eskip=[], verb=1):
      """compare against some target, which can be a 'name' string
         or another APExample instance

         for any key in eskip, do not print all details of element differences
      """
      itarg = self.find_instance(target, verb=verb)
      if itarg is None:
         return

      return self.compare_v_instance(itarg, eskip=eskip, verb=verb)

   def find_instance(self, target, verb):
      """target might already be an instance, but allow for a string,
         possibly with spaces, to be found as an example name

         on error, return None
      """
      if isinstance(target, APExample):
         return target

      # target must be either an instance or a string
      if type(target) != str:
         print("** find_instance: target not instance or string, %s" % target)
         return None

      # otherwise, try to find an instance in the ap_examples list
      eg = find_eg(target)
      # if not found, try replacing any '_' with ' '
      if eg is None and target.find('_') >= 0:
          eg = find_eg(target.replace('_', ' '))

      # success
      if eg is not None:
         if verb > 2:
            print("-- found instance for target string '%s'" % target)
         return eg

      # failure
      print("** find_instance : failed for target '%s'" % target)
      return None

   def compare_v_instance(self, target, eskip=[], verb=1):
      """compare against another APExample instance

            -show missing, extra and different categories

            target  can be name (string) or APExample instance
            eskip list of keys for which elements should not be compared
            verb    0 - use short list notation
                    1 - show missing, extra and diffs
                    2 - show diffs, but after trying to remove paths
                    3 - include all parameter lists
      """
      print("="*75)
      print("==== comparing (current) '%s' vs (target) '%s'" \
            % (self.name, target.name))
      if verb == 1:
         print('     (for more detail, consider "-verb 2")')
      print("")

      # use more generic name
      source = self

      # get global max key len
      maxs = max([len(key) for key in source.keys])
      maxt = max([len(key) for key in target.keys])
      maxk = max(maxs, maxt)

      # first look for missing and extra
      emiss =  [e for e in target.olist if e[0] not in source.keys]
      eextra = [e for e in source.olist if e[0] not in target.keys]

      # common keys, might be fewer, more, and/or have diffs
      efewer = []   # fewer of keys than target
      emore = []    # more than keys than target
      pdiff = []    # index pairs, where keys are the same but lists differ
      ncommon = 0

      # loop over unique keys, counting number in both key lists
      uniq_target_keys = UTIL.get_unique_sublist(target.keys)

      # rcr - right now, options of the same name need to be in the same order
      #       to be equal, we could try to remove the requirement...

      for key in uniq_target_keys:

         if key not in source.keys: continue

         ksource = [ind for ind, e in enumerate(source.olist) if e[0]==key]
         ktarget = [ind for ind, e in enumerate(target.olist) if e[0]==key]

         # how many of each?  both must be positive
         nsource = len(ksource)
         ntarget = len(ktarget)

         # get any efewer or emore entries
         if nsource < ntarget:
            efewer.extend([target.olist[ktarget[ind]] \
                          for ind in range(nsource,ntarget)])
         elif ntarget < nsource:
            emore.extend([source.olist[ksource[ind]] \
                         for ind in range(ntarget,nsource)])

         # get matching e[0]'s where e[1]'s differ
         ncomp = min(nsource, ntarget)
         ncommon += ncomp
         for ind in range(ncomp):
            # rcr - should we separate path diffs from true diffs here?
            if source.olist[ksource[ind]] != target.olist[ktarget[ind]]:
               # just store source and target indices here
               pdiff.append([ksource[ind], ktarget[ind]])

      nindent = 3
      ind1 = ' '*nindent
      ind2 = ' '*(2*nindent)

      # possibly print compact output
      if verb == 0:
         print("== missing (%d): %s\n" % (len(emiss),  [e[0] for e in emiss ]))
         print("== extra   (%d): %s\n" % (len(eextra), [e[0] for e in eextra]))
         print("== fewer   (%d): %s\n" % (len(efewer), [e[0] for e in efewer]))
         print("== more    (%d): %s\n" % (len(emore),  [e[0] for e in emore ]))
         print("== common  (%d - with %d differing)\n" % (ncommon, len(pdiff)))
         print("== diffs   (%d): %s\n" % (len(pdiff),
                                       [source.olist[p[0]][0] for p in pdiff ]))
         return

      # more verbose output allows for full printing
      if verb > 1: lmax = 0
      else:        lmax = 45

      print("==========  missing option(s) : %d" % len(emiss))
      if len(emiss) > 0:
         for e in emiss:
             estr = ' '.join(e[1])
             self._print_opt_lin(ind1, e[0], maxk, estr, lmax=lmax)
      print("")

      print("==========  extra option(s) : %d" % len(eextra))
      if len(eextra) > 0:
         for e in eextra:
             estr = ' '.join(e[1])
             self._print_opt_lin(ind1, e[0], maxk, estr, lmax=lmax)
      print("")

      print("==========  differing option(s) : %d" % len(pdiff))
      if len(pdiff) > 0:
         samestr = '  : SAME, but path diff (skipping details)'
         skips = []
         for pair in pdiff:
             oname = source.olist[pair[0]][0]
             if oname in eskip and verb < 2:
                # only show skipping details once
                if oname in skips:
                   continue
                skips.append(oname)

                # if only path diff, say so
                if self.depath_lists_equal(source.olist[pair[0]][1],
                                           target.olist[pair[1]][1]):
                   print("%s%-*s%s" % (ind1, maxk, oname, samestr))
                else:
                   print("%s%-*s  (verb=1; skipping details)" \
                         % (ind1, maxk, oname))

             elif verb <= 2:
                # try to remove any paths
                kslist = self.depath_list(source.olist[pair[0]][1])
                ktlist = self.depath_list(target.olist[pair[1]][1])
                ksstr = ' '.join(kslist)
                ktstr = ' '.join(ktlist)
                # if removing paths makes them the same, say so
                if ksstr == ktstr:
                    print("%s%-*s%s" % (ind1, maxk, oname, samestr))
                else:
                    print("%s%-*s" % (ind1, maxk, oname))
                    self._print_diff_line(ind2, 'current', ksstr, lmax=lmax)
                    self._print_diff_line(ind2, 'target',  ktstr, lmax=lmax)
             else:
                # full details
                print("%s%-*s" % (ind1, maxk, oname))

                ksstr = ' '.join(source.olist[pair[0]][1])
                ktstr = ' '.join(target.olist[pair[1]][1])
                self._print_diff_line(ind2, 'current', ksstr, lmax=lmax)
                self._print_diff_line(ind2, 'target',  ktstr, lmax=lmax)
             print("")

      print("==========  fewer applied option(s) : %d" % len(efewer))
      if len(efewer) > 0:
         for e in efewer:
             estr = ' '.join(e[1])
             self._print_opt_lin(ind1, e[0], maxk, estr, lmax=lmax)
      print("")

      print("==========  more applied option(s) : %d" % len(emore))
      if len(emore) > 0:
         for e in emore:
             estr = ' '.join(e[1])
             self._print_opt_lin(ind1, e[0], maxk, estr, lmax=lmax)
      print("")

      print("")

   def merge_w_instance(self, target, eskip=[], verb=1):
      """merge with another APExample instance

         Expand the current instance with options in the target.

         For each option in the target:

            if it exists in current
              ignore (for repeated options, do we require all?)
            else
              add (if in eskip, prepend -CHECK)

         Do not need much verb, since we can compare opts afterwards.

         return 0 on success
      """

      # use more generic name
      source = self

      if verb > 2: print("-- merge_w_instance() ..., target %s" % target)

      # be sure we have an instance, and not a string
      target = self.find_instance(target, verb=verb)
      if target is None:
         return 1

      # loop over unique keys, counting number in both key lists
      uniq_target_keys = UTIL.get_unique_sublist(target.keys)

      for key in uniq_target_keys:

         # ignore if already in current (ponder repeated options?)
         if key in source.keys:
            if verb > 2:
               print("-- skip already applied, %s" % key)
            continue

         # key is missing in source, so add all such target options to source

         # make note of all target options
         # ksource = [ind for ind, e in enumerate(source.olist) if e[0]==key]
         ktarget = [ind for ind, e in enumerate(target.olist) if e[0]==key]
         ntarget = len(ktarget)

         # if the key is in eskip, modify the option name with prefix -CHECK
         prefix = ''
         if key in eskip:
            if verb > 2:
               print("-- add with prefix, %s" % key)
            prefix = '-CHECK'
         else:
            if verb > 2:
               print("-- add directly, %s" % key)

         # for each such option, duplicate, possibly modify, and insert
         for ind in range(ntarget):
            newopt = copy.deepcopy(target.olist[ktarget[ind]])
            if prefix:
               newopt[0] = prefix+newopt[0]
            source.olist.append(newopt)

         # and add the new key, without prefix
         source.keys.append(key)

      return 0

   def depath_lists_equal(self, list0, list1):
       """if we depath_list() both lists, are they equal"""
       kslist = self.depath_list(list0)
       ktlist = self.depath_list(list1)
       ksstr = ' '.join(kslist)
       ktstr = ' '.join(ktlist)

       return (ksstr == ktstr)

   def depath_list(self, arglist):
       """return a copy of arglist, removing anything up to each final '/'
       """
       # rfind returns the highest index of '/', else -1
       # --> so always start at the return index+1
       newlist = [s[s.rfind('/')+1:] for s in arglist]

       return newlist

   def _print_opt_lin(self, indent, oname, maxlen, parstr, lmax=50):
       """if lmax > 0: restrict parstr to given length, plus ellipsis
       """
       estr = ''
       kprint = parstr
       if lmax > 0:
          if len(parstr) > lmax:
             kprint = parstr[0:lmax]
             estr = ' ...'

       print("%s%-*s  %s%s" % (indent, maxlen, oname, kprint, estr))

   def _print_diff_line(self, indent, tstr, parstr, lmax=50):
       """if lmax > 0: restrict parstr to given length, plus ellipsis
       """
       estr = ''
       kprint = parstr
       if lmax > 0:
          if len(parstr) > lmax:
             kprint = parstr[0:lmax]
             estr = ' ...'

       print("%s%-8s : %s%s" % (indent, tstr, kprint, estr))

   def command_string(self, wrap=0, windent=8):
      """return a string that is an afni_proc.py command
            if wrap: - return indented command
                     - indent by windent
      """

      # minor format conversion to a list of lists per line
      newolist = [['afni_proc.py']] + [([o[0]] + o[1]) for o in self.olist]

      # set default left-padding and line length
      if wrap:
         leftstr = ' '*windent
         linewid = 78
         is_diff, str_nice = lfcs.afni_niceify_cmd_str('',
                                                       comment_start=leftstr,
                                                       max_lw=linewid,
                                                       big_list=newolist)
      else:
         # join nested list into single string
         str_nice = ' '.join([' '.join(entry) for entry in newolist])

      return str_nice

   def display(self, verb=0, sphinx=1):
      """display a single example - use a copy for quoting
         verb: verbosity level
           0: only show example, with no wrap or formatting
           1: include name, source, descrip
           2: full verbosity: as ap -help: include descrip, header, trailer
              - terminate descrip with ~2~, for sphinxificaiton

         ponder indentation
      """
      cc = self.copy(quotize=1)

      # handle verb 0 right away, as the command string will differ
      if verb == 0:
         print("%s" % cc.command_string())
         return

      cmd = cc.command_string(wrap=1)

      indent = ' '*8

      if sphinx: sstr = ' ~2~'
      else:      sstr = ''

      # possibly show the name/description line,
      # (with trailing ~2~ for sphinx help)
      if verb > 0:
         print("%s%s. %s%s" % (indent, cc.name, cc.descrip, sstr))

      # header
      if verb > 1:
         print("%s" % cc.header)
         print("%s--------------------------" % indent)
         print("%slast mod date : %s" % (indent, cc.moddate))
         print("%skeywords      : %s" % (indent, ', '.join(cc.keywords)))
         print("%s--------------------------" % indent)

      print("")

      # print the actual example
      print("%s" % cmd)

      # any trailer
      if verb > 1 and cc.trailer != '':
         print("%s" % cc.trailer)

      print("")

      del(cc)


def get_all_examples():
   """return all known examples"""

   examples = []

   examples.extend(egs_example())
   examples.extend(egs_class())
   examples.extend(egs_demo())
   examples.extend(egs_short())
   examples.extend(egs_publish())

   return examples


def populate_examples(keys_keep=[], keys_rm=[], verb=1):
   """only populate the examples array if someone wants it

        keys_keep : if not empty, keep only entries with any of these keywords
        keys_rm   : if not empty, remove entries with any of these keywords
                    - default to 'noshow', those examples not shown by default
                      with the -help output

      populate the global ap_examples array from some of the egs_* functions
   """

   global ap_examples

   if len(ap_examples) > 0:
      return

   # use a local name to populate the global list
   examples = get_all_examples()

   # keys_keep = ['task', 'surface']
   # keys_rm   = ['obsolete', 'surface']

   if verb > 2:
      show_example_keywords(examples, mesg='initial examples')

   # --------------------------------------------------
   # keep only those entries with at least one of the given keep keys
   if len(keys_keep) > 0:
      newex = []
      for key in keys_keep:
         for ex in examples:
            if key in ex.keywords:
               if ex not in newex:
                  newex.append(ex)

      examples = newex

      if verb > 2:
         show_example_keywords(examples, mesg='post-keep examples')

   # --------------------------------------------------
   # remove any entry with any given rm key
   poplist = []
      # counting from the end, pop unwanted indices
   for eind, eg in enumerate(examples):
      # if there is any bad key, add to poplist and break
      for key in keys_rm:
         if key in eg.keywords:
            poplist.append(eg)
            break
   # remove everything in poplist
   if len(poplist) > 0:
      examples = [e for e in examples if e not in poplist]

   ap_examples = examples

def show_example_keywords(elist, mesg='', verb=1):
   """
   """
   if mesg: mstr = '(%s) ' % mesg
   else:    mstr = ''

   if len(elist) == 0:
      print("-- %d examples %s:" % (len(elist), mstr))
      print()
      return

   # special cases, get all keywords, and show a unique list
   if 'ALL' in elist:
      elist = get_all_examples()

   # if being quiet just show a compact key list
   if verb <= 1:
      # now python has a strange left to right reading to nest 2 lists
      klist = [e for ex in elist for e in ex.keywords]
      klist = UTIL.get_unique_sublist(klist)
      klist.sort()
      print('   ' + '\n   '.join(klist))
      return

   # else verbose, so show all examples and their keys
   print("-- %d examples %s:" % (len(elist), mstr))
   # align ':' using max name length
   nlen = max([len(ex.name) for ex in elist])
   for ex in elist:
      print("    %-*s : %s" % (nlen, ex.name, ', '.join(ex.keywords)))
   print()

def egs_ap_run(keys_keep=[], keys_rm=[]):
   """
            ***** this is probably garbage *****

      only populate the examples array if someone wants it

      return an array of APExample objects

       # sample example generation
       sample = APExample( 'my help example name',
         source='place to find this example in use',
         descrip='short description w/name in help output',
         moddate='1946.02.30',
         header='''
                  (recommended?  no, not intended for a complete analysis)
                  (              merely shows how simple a command can be)

               Purpose of this command, shown before example command.
               ''',
         trailer='''Any parts to show after the example command.''',
         olist = [
            ['-dsets',                 ['epiRT*.HEAD']],
            ['-regress_stim_files',    ['stims.1D']],
            # any other needed options
           ],
         )

   """

   examples = []

   examples.append( APExample( 'run_ap FT',
     source='ap_run_simple_rest.tcsh on class subject FT',
     descrip='Basic example to generate APQC for quick review.',
     moddate='2023.12.18',
     keywords=['rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (            * might be done immediately after acquisition)

           This example was generated by running ap_run_simple_rest.tcsh,
           providing a single subject anat and (3 runs of) EPI.  It could
           be generated (and run) using the following:

              cd AFNI_data6/FT_analysis/FT
              ap_run_simple_rest.tcsh -subjid FT -run_proc \\
                -anat FT_anat+orig -epi FT_epi_r*.HEAD
           """,
     trailer="""
           This quick processing will generate results in a new FT.results
           directory, which contains the APQC HTML report under QC_FT (where
           FT is the subject ID).  The report can be viewed using:

              afni_open -b FT.results/QC_FT/index.html

           or using interactive afni and niivue:

              open_apqc.py -infiles FT.results/QC_FT/index.html

           """,
     olist = [
        ['-subj_id',                 ['FT.run_ap']],
        ['-script',                  ['proc.FT.run_ap']],
        ['-out_dir',                 [ 'FT.1.results']],
        ['-blocks',                  [ 'tshift', 'align', 'tlrc', 'volreg',
                                       'mask', 'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', [ 'tcat volreg']],
        ['-copy_anat',               [ 'FT_anat+orig']],
        ['-dsets',                   [ 'FT_epi_r1+orig.HEAD',
                'FT_epi_r2+orig.HEAD', 'FT_epi_r3+orig.HEAD']],
        ['-tcat_remove_first_trs',   [ '2']],
        ['-align_unifize_epi',       [ 'local']],
        ['-align_opts_aea',          [ '-cost lpc+ZZ -giant_move -check_flip']],
        ['-tlrc_base',               [ 'MNI152_2009_template_SSW.nii.gz']],
        ['-volreg_align_to',         [ 'MIN_OUTLIER']],
        ['-volreg_align_e2a',        []],
        ['-volreg_tlrc_warp',        []],
        ['-volreg_compute_tsnr',     [ 'yes']],
        ['-mask_epi_anat',           [ 'yes']],
        ['-blur_size',               [ '6']],
        ['-regress_apply_mot_types', [ 'demean deriv']],
        ['-regress_motion_per_run',  []],
        ['-regress_censor_motion',   [ '0.25']],
        ['-regress_censor_outliers', [ '0.05']],
        ['-regress_est_blur_epits',  []],
        ['-regress_est_blur_errts',  []],
        ['-regress_make_ideal_sum',  [ 'sum_ideal.1D']],
        ['-html_review_style',       [ 'pythonic']],
       ],
     ))

   examples.append( APExample('simple align orig EPI',
     source='basic example',
     descrip='Perform orig space alignment to EPI base.',
     moddate='2023.12.18',
     keywords=['registration'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (            - only to do basic registration to EPI base)

           Align the anat and EPI time series to the EPI MIN_OUTLIER volume.
           The regress block is included merely to provide extra QC.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',                 ['FT.align.orig']],
        ['-blocks',                  [ 'align', 'volreg', 'regress']],
        ['-dsets',                   [ 'FT_epi_r1+orig.HEAD',
                'FT_epi_r2+orig.HEAD', 'FT_epi_r3+orig.HEAD']],
        ['-tcat_remove_first_trs',   ['2']],
        ['-align_unifize_epi',       [ 'local']],
        ['-align_opts_aea',          [ '-cost lpc+ZZ -giant_move -check_flip']],
        ['-volreg_align_to',         [ 'MIN_OUTLIER']],
        ['-html_review_style',       [ 'pythonic']],
       ],
     ))

   return examples

def egs_example():
   """general "Example" list, corresponding to the 2023 list"""

   examples = []

   examples.append( APExample( 'Example 1',
     source='afni_proc.py -help',
     descrip='Minimum use.',
     moddate='2008.12.10',
     keywords=['obsolete', 'task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              merely shows how simple a command can be)

           Provide datasets and stim files (or stim_times files).  Note that a
           dataset suffix (e.g. HEAD) must be used with wildcards, so that
           datasets are not applied twice.  In this case, a stim_file with many
           columns is given, where the script to changes it to stim_times files.
           """,
     trailer=""" """,
     olist = [
        ['-dsets',                 ['epiRT*.HEAD']],
        ['-regress_stim_files',    ['stims.1D']],
       ],
     ))

   examples.append( APExample('Example 2',
     source='afni_proc.py -help',
     descrip='Very simple.',
     moddate='2009.05.28',
     keywords=['obsolete', 'task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              many missing preferences, e.g. @SSwarper)

           Use all defaults, except remove 3 TRs and use basis
           function BLOCK(30,1).  The default basis function is GAM.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['sb23.e2.simple']],
        ['-dsets',                 ['sb23/epi_r??+orig.HEAD']],
        ['-tcat_remove_first_trs', ['3']],
        ['-regress_stim_times',    ['sb23/stim_files/blk_times.*.1D']],
        ['-regress_basis',         ['BLOCK(30,1)']],
       ],
     ))

   examples.append( APExample('Example 3',
     source='afni_proc.py -help',
     descrip='Formerly a simple class example.',
     moddate='2009.05.28',
     keywords=['obsolete', 'task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              many missing preferences, e.g. @SSwarper)

           Copy the anatomy into the results directory, register EPI data to
           the last TR, specify stimulus labels, compute blur estimates, and
           provide GLT options directly to 3dDeconvolve.  The GLTs will be
           ignored after this, as they take up too many lines.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['sb23.blk']],
        ['-dsets',                 ['sb23/epi_r??+orig.HEAD']],
        ['-copy_anat',             ['sb23/sb23_mpra+orig']],
        ['-tcat_remove_first_trs', ['3']],
        ['-volreg_align_to',       ['last']],
        ['-regress_stim_times',    ['sb23/stim_files/blk_times.*.1D']],
        ['-regress_stim_labels',   ['tneg', 'tpos', 'tneu', 'eneg', 'epos',
                                    'eneu', 'fneg', 'fpos', 'fneu']],
        ['-regress_basis',         ['BLOCK(30,1)']],
        ['-regress_opts_3dD',      ['-gltsym', 'SYM: +eneg -fneg',
                '-glt_label', '1', 'eneg_vs_fneg', '-gltsym',
                'SYM: 0.5*fneg 0.5*fpos -1.0*fneu', '-glt_label', '2',
                'face_contrast', '-gltsym',
                'SYM: tpos epos fpos -tneg -eneg -fneg',
                '-glt_label', '3', 'pos_vs_neg']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample( 'Example 4',
     source='afni_proc.py -help',
     descrip='Similar to 3, but specify the processing blocks.',
     moddate='2009.05.28',
     keywords=['obsolete', 'task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              many missing preferences, e.g. @SSwarper)

           Adding despike and tlrc, and removing tshift.  Note that
           the tlrc block is to run @auto_tlrc on the anat.  Ignore the GLTs.
           """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['sb23.e4.blocks']],
        ['-dsets',                 ['sb23/epi_r??+orig.HEAD']],
        ['-blocks',                ['despike', 'volreg', 'blur', 'mask',
                                    'scale', 'regress', 'tlrc']],
        ['-copy_anat',             ['sb23/sb23_mpra+orig']],
        ['-tcat_remove_first_trs', ['3']],
        ['-regress_stim_times',    ['sb23/stim_files/blk_times.*.1D']],
        ['-regress_stim_labels',   ['tneg', 'tpos', 'tneu', 'eneg', 'epos',
                                    'eneu', 'fneg', 'fpos', 'fneu']],
        ['-regress_basis',         ['BLOCK(30,1)']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample( 'Example 5a',
     source='afni_proc.py -help',
     descrip='RETROICOR, resting state data.',
     moddate='2009.05.28',
     keywords=['obsolete', 'physio', 'rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              just a terribly simple example using ricor)

           Assuming the class data is for resting-state and that we have the
           appropriate slice-based regressors from RetroTS.py, apply the
           despike and ricor processing blocks.  Note that '-do_block' is used
           to add non-default blocks into their default positions.  Here the
           'despike' and 'ricor' processing blocks would come before 'tshift'.

           Remove 3 TRs from the ricor regressors to match the EPI data.  Also,
           since degrees of freedom are not such a worry, regress the motion
           parameters per-run (each run gets a separate set of 6 regressors).

           The regression will use 81 basic regressors (all of "no interest"),
           with 13 retroicor regressors being removed during preprocessing:

                 27 baseline  regressors ( 3 per run * 9 runs)
                 54 motion    regressors ( 6 per run * 9 runs)

           To example #3, add -do_block, -ricor_* and -regress_motion_per_run.
           """,
     trailer="""
           If tshift, blurring and masking are not desired, consider replacing
           the -do_block option with an explicit list of blocks:

                -blocks despike ricor volreg regress
           """,
     olist = [
        ['-subj_id',               ['sb23.e5a.ricor']],
        ['-dsets',                 ['sb23/epi_r??+orig.HEAD']],
        ['-do_block',              ['despike', 'ricor']],
        ['-tcat_remove_first_trs', ['3']],
        ['-ricor_regs_nfirst',     ['3']],
        ['-ricor_regs',            ['sb23/RICOR/r*.slibase.1D']],
        ['-regress_motion_per_run', []],
       ],
     ))

   examples.append( APExample( 'Example 5b',
     source='afni_proc.py -help',
     descrip='RETROICOR, while running a normal regression.',
     moddate='2009.05.28',
     keywords=['obsolete', 'physio', 'rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              another overly simple example using ricor)

           Add the ricor regressors to a normal regression-based processing
           stream.  Apply the RETROICOR regressors across runs (so using 13
           concatenated regressors, not 13*9).  Note that concatenation is
           normally done with the motion regressors too.

           To example #3, add -do_block and three -ricor options.
           """,
     trailer="""
           Also consider adding -regress_bandpass.
           """,
     olist = [
        ['-subj_id',               ['sb23.e5b.ricor']],
        ['-dsets',                 ['sb23/epi_r??+orig.HEAD']],
        ['-do_block',              ['despike', 'ricor']],
        ['-copy_anat',             ['sb23/sb23_mpra+orig']],
        ['-tcat_remove_first_trs', ['3']],
        ['-ricor_regs_nfirst',     ['3']],
        ['-ricor_regs',            ['sb23/RICOR/r*.slibase.1D']],
        ['-ricor_regress_method',  ['across-runs']],
        ['-volreg_align_to',       ['last']],
        ['-regress_stim_times',    ['sb23/stim_files/blk_times.*.1D']],
        ['-regress_stim_labels',   ['tneg', 'tpos', 'tneu', 'eneg', 'epos',
                                    'eneu', 'fneg', 'fpos', 'fneu']],
        ['-regress_basis',         ['BLOCK(30,1)']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample( 'Example 5c',
     source='afni_proc.py -help',
     descrip='RETROICOR: censor and band pass.',
     moddate='2016.05.03',
     keywords=['obsolete', 'physio', 'task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              many missing preferences, e.g. @SSwarper, no BP)

           This is an example of how we might currently suggest analyzing
           resting state data.  If no RICOR regressors exist, see example 9
           (or just remove any ricor options).

           Censoring due to motion has long been considered appropriate in
           BOLD FMRI analysis, but is less common for those doing bandpass
           filtering in RS FMRI because the FFT requires one to either break
           the time axis (evil) or to replace the censored data with something
           probably inappropriate.

           Instead, it is slow (no FFT, but maybe SFT :) but effective to
           regress frequencies within the regression model, where censoring
           is simple.

           Note: band passing in the face of RETROICOR is questionable.  It may
                 be questionable in general.  To skip bandpassing, remove the
                 -regress_bandpass option line.

           Also, align EPI to anat and warp to standard space.
           """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['sb23.e5a.ricor']],
        ['-dsets',                 ['sb23/epi_r??+orig.HEAD']],
        ['-blocks',                ['despike', 'ricor', 'tshift', 'align',
                            'tlrc', 'volreg', 'blur', 'mask', 'regress']],
        ['-copy_anat',             ['sb23/sb23_mpra+orig']],
        ['-tcat_remove_first_trs', ['3']],
        ['-ricor_regs_nfirst',     ['3']],
        ['-ricor_regs',            ['sb23/RICOR/r*.slibase.1D']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-blur_size',             ['6']],
        ['-regress_bandpass',      ['0.01', '0.1']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_run_clustsim',  ['no']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample( 'Example 6',
     source='afni_proc.py -help',
     descrip='A simple task example, based on AFNI_data6.',
     moddate='2020.02.15',
     keywords=['task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              meant to be fast, but not complete, e.g. NL warp)
              (              prefer: see Example 6b)

           This example has changed to more closely correspond with the
           the class analysis example, AFNI_data6/FT_analysis/s05.ap.uber.

           The tshift block will interpolate each voxel time series to adjust
           for differing slice times, where the result is more as if each
           entire volume were acquired at the beginning of the TR.

           The 'align' block implies using align_epi_anat.py to align the
           anatomy with the EPI.  Here, the EPI base is first unifized locally.
           Additional epi/anat alignment options specify using lpc+ZZ for the
           cost function (more robust than simply lpc), -giant_move (in case
           the anat and EPI start a bit far apart), and -check_flip, to try to
           verify whether EPI left and right agree with the anatomy.
           This block computes the anat to EPI transformation matrix, which
           will be inverted in the volreg block, based on -volreg_align_e2a.

           Also, compute the transformation of the anatomy to MNI space, using
           affine registration (for speed in this simple example) to align to
           the 2009c template.

           In the volreg block, align the EPI to the MIN_OUTLIER volume (a
           low-motion volume, determined based on the data).  Then concatenate
           all EPI transformations, warping the EPI to standard space in one
           step (without multiple resampling operations), combining:

              EPI  ->  EPI base  ->  anat  ->  MNI 2009c template

           The standard space transformation is included by specifying option
           -volreg_tlrc_warp.

           A 4 mm blur is applied, to keep it very light (about 1.5 times the
           voxel size).

           The regression model is based on 2 conditions, each lasting 20 s
           per event, modeled by convolving a 20 s boxcar function with the
           BLOCK basis function, specified as BLOCK(20,1) to make the regressor
           unit height (height 1).

           One extra general linear test (GLT) is included, contrasting the
           visual reliable condition (vis) with auditory reliable (aud).

           Motion regression will be per run (using one set of 6 regressors for
           each run, i.e. 18 regressors in this example).

           The regression includes censoring of large motion (>= 0.3 ~mm
           between successive time points, based on the motion parameters),
           as well as censoring of outlier time points, where at least 5% of
           the brain voxels are computed as outliers.

           The regression model starts as a full time series, for time
           continuity, before censored time points are removed.  The output
           errts will be zero at censored time points (no error there), and so
           the output fit times series (fitts) will match the original data.

           The model fit time series (fitts) will be computed AFTER the linear
           regression, to save RAM on class laptops.

           Create sum_ideal.1D, as the sum of all non-baseline regressors, for
           quality control.

           Estimate the blur in the residual time series.  The resulting 3 ACF
           parameters can be averaged across subjects for cluster correction at
           the group level.

           Skip running the Monte Carlo cluster simulation example (which would
           specify minimum cluster sizes for cluster significance, based on the
           ACF parameters and mask), for speed.

           Once the proc script is created, execute it.
           """,
     trailer="""
         * One could also use ANATICOR with task (e.g. -regress_anaticor_fast)
           in the case of -regress_reml_exec.  3dREMLfit supports voxelwise
           regression, but 3dDeconvolve does not.
             """,
     olist = [
        ['-subj_id',               ['FT.e6']],
        ['-copy_anat',             ['FT/FT_anat+orig']],
        ['-dsets',                 ['FT/FT_epi_r?+orig.HEAD']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-volreg_compute_tsnr',   ['yes']],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['4.0']],
        ['-regress_stim_times',    ['FT/AV1_vis.txt', 'FT/AV2_aud.txt']],
        ['-regress_stim_labels',   ['vis', 'aud']],
        ['-regress_basis',         ['BLOCK(20,1)']],
        ['-regress_opts_3dD',      ['-jobs', '2', '-gltsym', 'SYM: vis -aud',
                                    '-glt_label', '1', 'V-A']],
        ['-regress_motion_per_run',[]],
        ['-regress_censor_motion', ['0.3']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_compute_fitts', []],
        ['-regress_make_ideal_sum',['sum_ideal.1D']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-regress_run_clustsim',  ['no']],
        ['-html_review_style',     ['pythonic']],
        ['-execute',               []],
       ],
     ))

   examples.append( APExample( 'Example 6b',
     source='afni_proc.py -help',
     descrip='A modern task example, with preferable options.',
     moddate='2020.02.15',
     keywords=['complete', 'task'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)

           GOOD TO CONSIDER

           This is based on Example 6, but is more complete.
           Example 6 is meant to run quickly, as in an AFNI bootcamp setting.
           Example 6b is meant to process more as we might suggest.

              - apply -check_flip in align_epi_anat.py, to monitor consistency
              - apply non-linear registration to MNI template, using output
                from @SSwarper:
                  o apply skull-stripped anat in -copy_anat
                  o apply original anat as -anat_follower (QC, for comparison)
                  o pass warped anat and transforms via -tlrc_NL_warped_dsets,
                    to apply those already computed transformations
              - use -mask_epi_anat to tighten the EPI mask (for QC),
                intersecting it (full_mask) with the anat mask (mask_anat)
              - use 3dREMLfit for the regression, to account for temporal
                autocorrelation in the noise
                (-regress_3dD_stop, -regress_reml_exec)
              - generate the HTML QC report using the nicer pythonic functions
                (requires matplotlib)
           """,
     trailer="""
           To compare one's own command against this one, consider adding
                -compare_opts 'example 6b'
           to the end of (or anywhere in) the current command, as in:

                afni_proc.py ... my options ...   -compare_opts 'example 6b'
           """,
     olist = [
        ['-subj_id',               ['FT.e6b']],
        ['-copy_anat',             ['Qwarp/anat_warped/anatSS.FT.nii']],
        ['-anat_has_skull',        ['no']],
        ['-anat_follower',         ['anat_w_skull', 'anat', 'FT/FT_anat+orig']],
        ['-dsets',                 ['FT/FT_epi_r?+orig.HEAD']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template_SSW.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-tlrc_NL_warped_dsets',  ['Qwarp/anat_warped/anatQQ.FT.nii',
                                    'Qwarp/anat_warped/anatQQ.FT.aff12.1D',
                                    'Qwarp/anat_warped/anatQQ.FT_WARP.nii']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-volreg_compute_tsnr',   ['yes']],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['4.0']],
        ['-regress_stim_times',    ['FT/AV1_vis.txt', 'FT/AV2_aud.txt']],
        ['-regress_stim_labels',   ['vis', 'aud']],
        ['-regress_basis',         ['BLOCK(20,1)']],
        ['-regress_opts_3dD',      ['-jobs', '2', '-gltsym', 'SYM: vis -aud',
                                    '-glt_label', '1', 'V-A']],
        ['-regress_motion_per_run',[]],
        ['-regress_censor_motion', ['0.3']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_3dD_stop',      []],
        ['-regress_reml_exec',     []],
        ['-regress_compute_fitts', []],
        ['-regress_make_ideal_sum',['sum_ideal.1D']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-regress_run_clustsim',  ['no']],
        ['-html_review_style',     ['pythonic']],
        ['-execute',               []],
       ],
     ))

   examples.append( APExample( 'Example 7',
     source='afni_proc.py -help',
     descrip='Apply some esoteric options.',
     moddate='2020.01.08',
     keywords=['task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              e.g. NL warp without @SSwarper)
              (              prefer: see Example 6b)

           a. Blur only within the brain, as far as an automask can tell.  So
              add -blur_in_automask to blur only within an automatic mask
              created internally by 3dBlurInMask (akin to 3dAutomask).

           b. Let the basis functions vary.  For some reason, we expect the
              BOLD responses to the telephone classes to vary across the brain.
              So we have decided to use TENT functions there.  Since the TR is
              3.0s and we might expect up to a 45 second BOLD response curve,
              use 'TENT(0,45,16)' for those first 3 out of 9 basis functions.

              This means using -regress_basis_multi instead of -regress_basis,
              and specifying all 9 basis functions appropriately.

           c. Use amplitude modulation.

              We expect responses to email stimuli to vary proportionally with
              the number of punctuation characters used in the message (in
              certain brain regions).  So we will use those values as auxiliary
              parameters 3dDeconvolve by marrying the parameters to the stim
              times (using 1dMarry).

              Use -regress_stim_types to specify that the epos/eneg/eneu stim
              classes should be passed to 3dDeconvolve using -stim_times_AM2.

           d. Not only censor motion, but censor TRs when more than 10% of the
              automasked brain are outliers.  So add -regress_censor_outliers.

           e. Include both de-meaned and derivatives of motion parameters in
              the regression.  So add '-regress_apply_mot_types demean deriv'.

           f. Output baseline parameters so we can see the effect of motion.
              So add -bout under option -regress_opts_3dD.

           g. Save on RAM by computing the fitts only after 3dDeconvolve.
              So add -regress_compute_fitts.

           h. Speed things up.  Have 3dDeconvolve use 4 CPUs and skip the
              single subject 3dClustSim execution.  So add '-jobs 4' to the
              -regress_opts_3dD option and add '-regress_run_clustsim no'.
           """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['sb23.e7.esoteric']],
        ['-dsets',                 ['sb23/epi_r??+orig.HEAD']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'blur', 'mask', 'scale', 'regress']],
        ['-copy_anat',             ['sb23/sb23_mpra+orig']],
        ['-tcat_remove_first_trs', ['3']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['4']],
        ['-blur_in_automask',      []],
        ['-regress_stim_times',    ['sb23/stim_files/blk_times.*.1D']],
        ['-regress_stim_types',    ['times', 'times', 'times',
                                    'AM2', 'AM2', 'AM2',
                                    'times', 'times', 'times']],
        ['-regress_stim_labels',   ['tneg', 'tpos', 'tneu',
                                    'eneg', 'epos', 'eneu',
                                    'fneg', 'fpos', 'fneu']],
        ['-regress_basis_multi',['BLOCK(30,1)', 'TENT(0,45,16)','BLOCK(30,1)',
                                 'BLOCK(30,1)', 'TENT(0,45,16)','BLOCK(30,1)',
                                 'BLOCK(30,1)', 'TENT(0,45,16)','BLOCK(30,1)']],
        ['-regress_opts_3dD',      ['-bout', '-gltsym', 'SYM: +eneg -fneg',
                            '-glt_label', '1', 'eneg_vs_fneg', '-jobs', '4']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.3']],
        ['-regress_censor_outliers', ['0.1']],
        ['-regress_compute_fitts', []],
        ['-regress_run_clustsim',  ['no']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample('Example 8',
     source='afni_proc.py -help',
     descrip='Surface-based analysis.',
     moddate='2017.09.12',
     keywords=['complete', 'surface', 'task'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)

           This example is intended to be run from AFNI_data6/FT_analysis.
           It is provided with the class data in file s03.ap.surface.

           Add -surf_spec and -surf_anat to provide the required spec and
           surface volume datasets.  The surface volume will be aligned to
           the current anatomy in the processing script.  Two spec files
           (lh and rh) are provided, one for each hemisphere (via wildcard).

           Also, specify a (resulting) 6 mm FWHM blur via -blur_size.  This
           does not add a blur, but specifies a resulting blur level.  So
           6 mm can be given directly for correction for multiple comparisons
           on the surface.

           Censor per-TR motion above 0.3 mm.

           Note that no -regress_est_blur_errts option is given, since that
           applies to the volume only (and since the 6 mm blur is a resulting
           blur level, so the estimates are not needed).

           The -blocks option is provided, but it is the same as the default
           for surface-based analysis, so is not really needed here.  Note that
           the 'surf' block is added and the 'mask' block is removed from the
           volume-based defaults.

           important options:

                -blocks         : includes surf, but no mask
                                  (default blocks for surf, so not needed)
                -surf_anat      : volume aligned with surface
                -surf_spec      : spec file(s) for surface

           Note: one would probably want to use standard mesh surfaces here.
                 This example will be updated with them in the future.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT.surf']],
        ['-blocks',                ['tshift', 'align', 'volreg', 'surf',
                                    'blur', 'scale', 'regress']],
        ['-copy_anat',             ['FT/FT_anat+orig']],
        ['-dsets',                 ['FT/FT_epi_r?+orig.HEAD']],
        ['-surf_anat',             ['FT/SUMA/FTmb_SurfVol+orig']],
        ['-surf_spec',             ['FT/SUMA/FTmb_?h.spec']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-blur_size',             ['6']],
        ['-regress_stim_times',    ['FT/AV1_vis.txt', 'FT/AV2_aud.txt']],
        ['-regress_stim_labels',   ['vis', 'aud']],
        ['-regress_basis',         ['BLOCK(20,1)']],
        ['-regress_opts_3dD',      ['-jobs', '2', '-gltsym', 'SYM: vis -aud',
                                    '-glt_label', '1', 'V-A']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.3']],
       ],
     ))

   examples.append( APExample('Example 9',
     source='afni_proc.py -help',
     descrip='Resting state analysis with censoring and band passing.',
     moddate='2019.02.26',
     keywords=['rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              e.g. has band pass, no @SSwarper)
              (              prefer: see Example 11)

           With censoring and bandpass filtering.

           This is our suggested way to do preprocessing for resting state
           analysis, under the assumption that no cardio/physio recordings
           were made (see example 5 for cardio files).

           Censoring due to motion has long been considered appropriate in
           BOLD FMRI analysis, but is less common for those doing bandpass
           filtering in RS FMRI because the FFT requires one to either break
           the time axis (evil) or to replace the censored data with something
           probably inappropriate.

           Instead, it is slow (no FFT, but maybe SFT :) but effective to
           regress frequencies within the regression model, where censoring
           is simple.

           inputs: anat, EPI
           output: errts dataset (to be used for correlation)

           special processing:
              - despike, as another way to reduce motion effect
                 (see block despike)
              - censor motion TRs at the same time as bandpassing data
                 (see -regress_censor_motion, -regress_bandpass)
              - regress motion parameters AND derivatives
                 (see -regress_apply_mot_types)

           Note: for resting state data, a more strict threshold may be a good
                 idea, since motion artifacts should play a bigger role than in
                 a task-based analysis.

                 So the typical suggestion of motion censoring at 0.3 for task
                 based analysis has been changed to 0.2 for this resting state
                 example, and censoring of outliers has also been added, at a
                 value of 5% of the brain mask.

                 Outliers are typically due to motion, and may capture motion
                 in some cases where the motion parameters do not, because
                 motion is not generally a whole-brain-between-TRs event.

           Note: if regressing out regions of interest, either create the ROI
                 time series before the blur step, or remove blur from the list
                 of blocks (and apply any desired blur after the regression).

           Note: it might be reasonable to estimate the blur using epits rather
                 than errts in the case of bandpassing.  Both options are
                 included here.

           Note: scaling is optional here.  While scaling has no direct effect
                 on voxel correlations, it does have an effect on ROI averages
                 used for correlations.

           Other options to consider: -tlrc_NL_warp, -anat_uniform_method
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['subj123']],
        ['-dsets',                 ['epi_run1+orig.HEAD']],
        ['-copy_anat',             ['anat+orig']],
        ['-blocks',                ['despike', 'tshift', 'align', 'tlrc',
                          'volreg', 'blur', 'mask', 'scale', 'regress']],
        ['-tcat_remove_first_trs', ['3']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['4']],
        ['-regress_bandpass',      ['0.01', '0.1']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample('Example 9b',
     source='afni_proc.py -help',
     descrip='Resting state analysis with ANATICOR.',
     moddate='2020.01.08',
     keywords=['rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              e.g. has band pass, no @SSwarper)
              (              prefer: see Example 11)

           Like example #9, but also regress out the signal from locally
           averaged white matter.  The only change is adding the option
           -regress_anaticor.

           Note that -regress_anaticor implies options -mask_segment_anat and
           -mask_segment_erode.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['subj123']],
        ['-dsets',                 ['epi_run1+orig.HEAD']],
        ['-copy_anat',             ['anat+orig']],
        ['-blocks',                ['despike', 'tshift', 'align', 'tlrc',
                          'volreg', 'blur', 'mask', 'scale', 'regress']],
        ['-tcat_remove_first_trs', ['3']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['4']],
        ['-regress_bandpass',      ['0.01', '0.1']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_anaticor',      []],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample('Example 10',
     source='afni_proc.py -help',
     descrip='Resting state analysis, with tissue-based regressors.',
     moddate='2020.01.08',
     keywords=['rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              e.g. missing @SSwarper)
              (              prefer: see Example 11)

           Like example #9, but also regress the eroded white matter averages.
           The WMe mask come from the Classes dataset, created by 3dSeg via the
           -mask_segment_anat and -mask_segment_erode options.

        ** While -mask_segment_anat also creates a CSF mask, that mask is ALL
           CSF, not just restricted to the ventricles, for example.  So it is
           probably not appropriate for use in tissue-based regression.

           CSFe was previously used as an example of what one could do, but as
           it is not advised, it has been removed.

           Also, align to minimum outlier volume, and align to the anatomy
           using cost function lpc+ZZ.

           Note: it might be reasonable to estimate the blur using epits rather
                 than errts in the case of bandpassing.  Both options are
                 included here.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['subj123']],
        ['-dsets',                 ['epi_run1+orig.HEAD']],
        ['-copy_anat',             ['anat+orig']],
        ['-blocks',                ['despike', 'tshift', 'align', 'tlrc',
                                'volreg', 'blur', 'mask', 'scale', 'regress']],
        ['-tcat_remove_first_trs', ['3']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-blur_size',             ['4']],
        ['-mask_epi_anat',         ['yes']],
        ['-mask_segment_anat',     ['yes']],
        ['-mask_segment_erode',    ['yes']],
        ['-regress_bandpass',      ['0.01', '0.1']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_ROI',           ['WMe']],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample('Example 10b',
     source='afni_proc.py -help',
     descrip='Resting state analysis, as 10a with 3dRSFC.',
     moddate='2019.02.13',
     keywords=['rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              prefer: see Example 11)
              (         ***        : use censoring and 3dLombScargle)

            This is for band passing and computation of ALFF, etc.

          * This will soon use a modified 3dRSFC.

            Like example #10, but add -regress_RSFC to bandpass via 3dRSFC.
            Skip censoring and regression band passing because of the bandpass
            operation in 3dRSFC.

            To correspond to common tractography, this example stays in orig
            space (no 'tlrc' block, no -volreg_tlrc_warp option).  Of course,
            going to standard space is an option.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['subj123']],
        ['-dsets',                 ['epi_run1+orig.HEAD']],
        ['-copy_anat',             ['anat+orig']],
        ['-blocks',                ['despike', 'tshift', 'align', 'volreg',
                                    'blur', 'mask', 'scale', 'regress']],
        ['-tcat_remove_first_trs', ['3']],
        ['-volreg_align_e2a',      []],
        ['-blur_size',             ['6.0']],
        ['-mask_apply',            ['epi']],
        ['-mask_segment_anat',     ['yes']],
        ['-mask_segment_erode',    ['yes']],
        ['-regress_bandpass',      ['0.01', '0.1']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_ROI',           ['WMe']],
        ['-regress_RSFC',          []],
        ['-regress_run_clustsim',  ['no']],
        ['-regress_est_blur_errts', []],
       ],
     ))

   examples.append( APExample( 'Example 11',
     source='afni_proc.py -help',
     descrip='Resting state analysis (now even more modern :).',
     moddate='2022.10.06',
     keywords=['complete', 'rest'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)

         o Yes, censor (outliers and motion) and despike.
         o Align the anatomy and EPI using the lpc+ZZ cost function, rather
           than the default lpc one.  Apply -giant_move, in case the datasets
           do not start off well-aligned.  Include -check_flip for good measure.
           A locally unifized EPI base is used for anatomical registration.
         o Register EPI volumes to the one which has the minimum outlier
              fraction (so hopefully the least motion).
         o Use non-linear registration to MNI template (non-linear 2009c).
           * This adds a lot of processing time.
           * Let @SSwarper align to template MNI152_2009_template_SSW.nii.gz.
             Then use the resulting datasets in the afni_proc.py command below
             via -tlrc_NL_warped_dsets.
                  @SSwarper -input FT_anat+orig        \\
                            -subid FT                  \\
                            -odir  FT_anat_warped      \\
                            -base  MNI152_2009_template_SSW.nii.gz

            - The SS (skull-stripped) can be given via -copy_anat, and the
              with-skull unifized anatU can be given as a follower.
         o No bandpassing.
         o Use fast ANATICOR method (slightly different from default ANATICOR).
         o Use FreeSurfer segmentation for:
             - regression of first 3 principal components of lateral ventricles
             - ANATICOR white matter mask (for local white matter regression)
           * For details on how these masks were created, see "FREESURFER NOTE"
             in the help, as it refers to this "Example 11".
         o Erode FS white matter and ventricle masks before application.
         o Bring along FreeSurfer parcellation datasets:
             - aaseg : NN interpolated onto the anatomical grid
             - aeseg : NN interpolated onto the EPI        grid
           * These 'aseg' follower datasets are just for visualization,
             they are not actually required for the analysis.
         o Compute average correlation volumes of the errts against the
           the gray matter (aeseg) and ventricle (FSVent) masks.
         o Run @radial_correlate at the ends of the tcat, volreg and regress
           blocks.  If ANATICOR is being used to remove a scanner artifact,
           the errts radcor images might show the effect of this.

           Note: it might be reasonable to use either set of blur estimates
                 here (from epits or errts).  The epits (uncleaned) dataset
                 has all of the noise (though what should be considered noise
                 in this context is not clear), while the errts is motion
                 censored.  For consistency in resting state, it would be
                 reasonable to stick with epits.  They will likely be almost
                 identical.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT.11.rest']],
        ['-blocks',                ['despike', 'tshift', 'align', 'tlrc',
                          'volreg', 'blur', 'mask', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
        ['-copy_anat',             ['anatSS.FT.nii']],
        ['-anat_has_skull',        ['no']],
        ['-anat_follower',         ['anat_w_skull', 'anat', 'anatU.FT.nii']],
        ['-anat_follower_ROI',     ['aaseg', 'anat', 'aparc.a2009s+aseg_REN_all.nii.gz']],
        ['-anat_follower_ROI',     ['aeseg', 'epi', 'aparc.a2009s+aseg_REN_all.nii.gz']],
        ['-anat_follower_ROI',     ['FSvent', 'epi', 'fs_ap_latvent.nii.gz']],
        ['-anat_follower_ROI',     ['FSWe', 'epi', 'fs_ap_wm.nii.gz']],
        ['-anat_follower_erode',   ['FSvent', 'FSWe']],
        ['-dsets',                 ['FT_epi_r?+orig.HEAD']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template_SSW.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-tlrc_NL_warped_dsets',  ['anatQQ.FT.nii', 'anatQQ.FT.aff12.1D',
                                    'anatQQ.FT_WARP.nii']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['4']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_anaticor_fast', []],
        ['-regress_anaticor_label', ['FSWe']],
        ['-regress_ROI_PC',        ['FSvent', '3']],
        ['-regress_ROI_PC_per_run', ['FSvent']],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_make_corr_vols', ['aeseg', 'FSvent']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-html_review_style',     ['pythonic']],
       ],
     ))

   examples.append( APExample('Example 11b',
     source='afni_proc.py -help',
     descrip='Similar to 11, but without FreeSurfer.',
     moddate='2020.01.17',
     keywords=['complete', 'rest'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)
              (              if this ventricle extraction method seems okay)

         AFNI currently does not have a good program to extract ventricles.
         But it can make a CSF mask that includes them.  So without FreeSurfer,
         one could import a ventricle mask from the template (e.g. for TT space,
         using TT_desai_dd_mpm+tlrc).  For example, assuming Talairach space
         (and a 2.5 mm^3 final voxel grid) for the analysis, one could create a
         ventricle mask as follows:

                3dcalc -a ~/abin/TT_desai_dd_mpm+tlrc                       \\
                       -expr 'amongst(a,152,170)' -prefix template_ventricle
                3dresample -dxyz 2.5 2.5 2.5 -inset template_ventricle+tlrc \\
                       -prefix template_ventricle_2.5mm

         o Be explicit with 2.5mm, using '-volreg_warp_dxyz 2.5'.
         o Use template TT_N27+tlrc, to be aligned with the desai atlas.
         o No -anat_follower options, but use -mask_import to import the
           template_ventricle_2.5mm dataset (and call it Tvent).
         o Use -mask_intersect to intersect ventricle mask with the subject's
           CSFe mask, making a more reliable subject ventricle mask (Svent).
         o Ventricle principle components are created as per-run regressors.
         o Make WMe and Svent correlation volumes, which are just for
           entertainment purposes anyway.
         o Run the cluster simulation.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT.11b.rest']],
        ['-blocks',                ['despike', 'tshift', 'align', 'tlrc',
                          'volreg', 'blur', 'mask', 'scale', 'regress']],
        ['-copy_anat',             ['FT_anat+orig']],
        ['-dsets',                 ['FT_epi_r?+orig.HEAD']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['TT_N27+tlrc']],
        ['-tlrc_NL_warp',          []],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-volreg_warp_dxyz',      ['2.5']],
        ['-blur_size',             ['4']],
        ['-mask_segment_anat',     ['yes']],
        ['-mask_segment_erode',    ['yes']],
        ['-mask_import',           ['Tvent', 'template_ventricle_2.5mm+tlrc']],
        ['-mask_intersect',        ['Svent', 'CSFe', 'Tvent']],
        ['-mask_epi_anat',         ['yes']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_anaticor_fast', []],
        ['-regress_ROI_PC',        ['Svent', '3']],
        ['-regress_ROI_PC_per_run', ['Svent']],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_make_corr_vols', ['WMe', 'Svent']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-regress_run_clustsim',  ['yes']],
       ],
     ))

   examples.append( APExample('Example 12',
     source='afni_proc.py -help',
     descrip='background: Multi-echo data processing.',
     moddate='2018.02.27',
     keywords=['ME', 'rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              incomplete - just shows basic ME options)
              (              prefer: see Example 13)

         Processing multi-echo data should be similar to single echo data,
         except for perhaps:

            combine         : the addition of a 'combine' block
            -dsets_me_echo  : specify ME data, per echo
            -dsets_me_run   : specify ME data, per run (alternative to _echo)
            -echo_times     : specify echo times (if needed)
            -combine_method : specify method to combine echoes (if any)

         An afni_proc.py command might be updated to include something like:
            """,
     trailer=""" """,
     olist = [
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                            'mask', 'combine', 'blur', 'scale', 'regress']],
        ['-dsets_me_echo',         ['epi_run*_echo_01.nii']],
        ['-dsets_me_echo',         ['epi_run*_echo_02.nii']],
        ['-dsets_me_echo',         ['epi_run*_echo_03.nii']],
        ['-echo_times',            ['15', '30.5', '41']],
        ['-mask_epi_anat',         ['yes']],
        ['-combine_method',        ['OC']],
       ],
     ))

   examples.append( APExample('Example 12a',
     source='afni_proc.py -help',
     descrip='Multi-echo data processing - very simple.',
     moddate='2018.02.27',
     keywords=['ME', 'rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              many missing preferences, e.g. @SSwarper)
              (              prefer: see Example 13)

         Keep it simple and just focus on the basic ME options, plus a few
         for controlling registration.

         o This example uses 3 echoes of data across just 1 run.
            - so use a single -dsets_me_run option to input EPI datasets
         o Echo 2 is used to drive registration for all echoes.
            - That is the default, but it is good to be explicit.
         o The echo times are not needed, as the echoes are never combined.
         o The echo are never combined (in this example), so that there
           are always 3 echoes, even until the end.
            - Note that the 'regress' block is not valid for multiple echoes.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT.12a.ME']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'blur']],
        ['-copy_anat',             ['FT_anat+orig']],
        ['-dsets_me_run',          ['epi_run1_echo*.nii']],
        ['-reg_echo',              ['2']],
        ['-tcat_remove_first_trs', ['2']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
       ],
     ))

   examples.append( APExample('Example 12b',
     source='afni_proc.py -help',
     descrip='Multi-echo data processing - OC resting state.',
     moddate='2020.01.08',
     keywords=['ME', 'rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              many missing preferences, e.g. @SSwarper)
              (              prefer: see Example 13)

         Still keep this simple, mostly focusing on ME options, plus standard
         ones for resting state.

         o This example uses 3 echoes of data across just 1 run.
            - so use a single -dsets_me_run option to input EPI datasets
         o Echo 2 is used to drive registration for all echoes.
            - That is the default, but it is good to be explicit.
         o The echoes are combined via the 'combine' block.
         o So -echo_times is used to provided them.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT.12a.ME']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                            'mask', 'combine', 'blur', 'scale', 'regress']],
        ['-copy_anat',             ['FT_anat+orig']],
        ['-dsets_me_run',          ['epi_run1_echo*.nii']],
        ['-echo_times',            ['15', '30.5', '41']],
        ['-reg_echo',              ['2']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-combine_method',        ['OC']],
        ['-blur_size',             ['4']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits', []],
       ]
     ))

   examples.append( APExample('Example 12c',
     source='afni_proc.py -help',
     descrip='Multi-echo data processing - ME-ICA resting state.',
     moddate='2020.01.08',
     keywords=['ME', 'rest'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              many missing preferences, e.g. @SSwarper)
              (              prefer: see Example 13)

         As above, but run tedana.py for MEICA denoising.

         o Since tedana.py will mask the data, it may be preferable to
           blur only within that mask (-blur_in_mask yes).
         o A task analysis using tedana might look much the same,
           but with the extra -regress options for the tasks.
            """,
     trailer="""
         Consider an alternative combine method, 'tedana_OC_tedort'.
           """,
     olist = [
        ['-subj_id',               ['FT.12a.ME']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                            'mask', 'combine', 'blur', 'scale', 'regress']],
        ['-copy_anat',             ['FT_anat+orig']],
        ['-dsets_me_run',          ['epi_run1_echo*.nii']],
        ['-echo_times',            ['15', '30.5', '41']],
        ['-reg_echo',              ['2']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-combine_method',        ['tedana']],
        ['-blur_size',             ['4']],
        ['-blur_in_mask',          ['yes']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run',[]],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits',[]],
       ],
     ))

   examples.append( APExample( 'Example 13',
     source='afni_proc.py -help',
     descrip='Complicated ME, surface-based resting state example.',
     moddate='2019.09.06',
     keywords=['complete', 'ME', 'physio', 'rest', 'surface'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)

         Example 'publish 3d' might be preferable.

         Key aspects of this example:

            - multi-echo data, using "optimally combined" echoes
            - resting state analysis (without band passing)
            - surface analysis
            - blip up/blip down distortion correction
            - slice-wise regression of physiological parameters (RETROICOR)
            - ventricle principal component regression (3 PCs)
            - EPI volreg to per-run MIN_OUTLIER, with across-runs allineate
            - QC: @radial_correlate on tcat and volreg block results
            - QC: pythonic html report

            * since this is a surface-based example, the are no tlrc options

         Minor aspects:

            - a FWHM=6mm blur is applied, since blur on surface is TO is size

         Note: lacking good sample data for this example, it is simply faked
               for demonstration (echoes are identical, fake ricor parameters
               are not part of this data tree).
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT.complicated']],
        ['-dsets_me_echo',         ['FT/FT_epi_r?+orig.HEAD']],
        ['-dsets_me_echo',         ['FT/FT_epi_r?+orig.HEAD']],
        ['-dsets_me_echo',         ['FT/FT_epi_r?+orig.HEAD']],
        ['-echo_times',            ['11', '22.72', '34.44']],
        ['-blip_forward_dset',     ['FT/FT_epi_r1+orig.HEAD[0]']],
        ['-blip_reverse_dset',     ['FT/FT_epi_r1+orig.HEAD[0]']],
        ['-copy_anat',             ['FT/FT_anat+orig']],
        ['-anat_follower_ROI',     ['FSvent', 'epi', 'FT/SUMA/FT_vent.nii']],
        ['-anat_follower_erode',   ['FSvent']],
        ['-blocks',                ['despike', 'ricor', 'tshift', 'align',
                                    'volreg', 'mask', 'combine', 'surf',
                                    'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg']],
        ['-tcat_remove_first_trs', ['2']],
        ['-ricor_regs_nfirst',     ['2']],
        ['-ricor_regs',            ['FT/fake.slibase.FT.r?.1D']],
        ['-ricor_regress_method',  ['per-run']],
        ['-tshift_interp',         ['-wsinc9']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_post_vr_allin',  ['yes']],
        ['-volreg_pvra_base_index', ['MIN_OUTLIER']],
        ['-volreg_warp_final_interp', ['wsinc5']],
        ['-mask_epi_anat',         ['yes']],
        ['-combine_method',        ['OC']],
        ['-surf_anat',             ['FT/SUMA/FT_SurfVol.nii']],
        ['-surf_spec',             ['FT/SUMA/std.141.FT_?h.spec']],
        ['-blur_size',             ['6']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_ROI_PC',        ['FSvent', '3']],
        ['-regress_ROI_PC_per_run', ['FSvent']],
        ['-regress_make_corr_vols', ['FSvent']],
        ['-regress_censor_motion', ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-html_review_style',     ['pythonic']],
       ]
     ))

   return examples

def egs_class():
   """AP class examples
   """

   examples =  []

   examples.append( APExample('AP class 3',
     source='FT_analysis',
     descrip='s03.ap.surface - basic surface analysis',
     moddate='2025.03.11',
     keywords=['complete', 'surface', 'task'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)
                             (though it is a very simple example)

           This is the surface analysis run during an AFNI bootcamp.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT.surf']],
        ['-dsets',                 ['FT/FT_epi_r?+orig.HEAD']],
        ['-copy_anat',             ['FT/FT_anat+orig']],
        ['-blocks',                ['tshift', 'align', 'volreg', 'surf',
                                    'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_compute_tsnr',   ['yes']],
        ['-surf_anat',             ['FT/SUMA/FT_SurfVol.nii']],
        ['-surf_spec',             ['FT/SUMA/std.60.FT_?h.spec']],
        ['-blur_size',             ['6']],
        ['-regress_stim_times',    ['FT/AV1_vis.txt', 'FT/AV2_aud.txt']],
        ['-regress_stim_labels',   ['vis', 'aud']],
        ['-regress_basis',         ['BLOCK(20,1)']],
        ['-regress_opts_3dD',      ['-jobs', '2', '-gltsym', 'SYM: vis -aud',
                                    '-glt_label', '1', 'V-A']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.3']],
        ['-regress_censor_outliers', ['0.05']],
       ]
     ))

   examples.append( APExample('AP class 5',
     source='FT_analysis',
     descrip='s05.ap.uber - basic task analysis',
     moddate='2025.03.11',
     keywords=['task'],
     header="""
              (recommended?  no, not intended for a complete analysis)
              (              prefer: see Example publish 3b for NL warp)

           A basic task analysis with a pair of visual and auditory tasks.

           notable options include :
                - affine registration to MNI152_2009_template.nii.gz template
                - censoring based on both motion params and outliers
                - '-regress_compute_fitts' to reduce RAM needs in 3dD
                - mask_epi_anat - intersect full_mask (epi) with mask_anat
                - QC: computing radial correlation volumes at the end
                      of the tcat, volreg and regress processing blocks
                - QC: include -check_flip left/right consistency check
                - QC: compute sum of ideals, for evaluation
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT']],
        ['-dsets',                 ['FT/FT_epi_r1+orig.HEAD',
                                    'FT/FT_epi_r2+orig.HEAD',
                                    'FT/FT_epi_r3+orig.HEAD']],
        ['-copy_anat',             ['FT/FT_anat+orig']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template.nii.gz']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-volreg_compute_tsnr',   ['yes']],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['4.0']],
        ['-regress_stim_times',    ['FT/AV1_vis.txt', 'FT/AV2_aud.txt']],
        ['-regress_stim_labels',   ['vis', 'aud']],
        ['-regress_basis',         ['BLOCK(20,1)']],
        ['-regress_opts_3dD',      ['-jobs', '2', '-gltsym', 'SYM: vis -aud',
                                    '-glt_label', '1', 'V-A',
                                    '-gltsym', 'SYM: 0.5*vis +0.5*aud',
                                    '-glt_label', '2', 'mean.VA']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.3']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_compute_fitts', []],
        ['-regress_make_ideal_sum', ['sum_ideal.1D']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-regress_run_clustsim',  ['no']],
        ['-html_review_style',     ['pythonic']],
        ['-execute',               []],
       ]
     ))

   return examples

def egs_publish():
   """AP publish examples
   """

   examples =  []

   examples.append( APExample('AP publish 1',
     source='AFNI_demos',
     descrip='pamenc, ds000030.v16 parametric encoding task analysis.',
     moddate='2024.08.26',
     keywords=['complete', 'publish', 'task'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)

           While this example is reasonable, 'publish 3b' has more QC options,
           as well as updates for anat/EPI alignment and grid size.

           Events are modeled using duration modulation, so AM1 is applied.

           original analysis was from:
               Gorgolewski KJ, Durnez J and Poldrack RA.
               Preprocessed Consortium for Neuropsychiatric Phenomics dataset.
               F1000Research 2017, 6:1262
               https://doi.org/10.12688/f1000research.11964.2

           downloadable from https://legacy.openfmri.org/dataset/ds000030
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['SID']],
        ['-script',                ['proc.SID']],
        ['-scr_overwrite',         []],
        ['-dsets',                 ['func/SID_task-pamenc_bold.nii.gz']],
        ['-copy_anat',             ['anatSS.SID.nii']],
        ['-anat_has_skull',        ['no']],
        ['-anat_follower',         ['anat_w_skull', 'anat', 'anatU.SID.nii']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'blur', 'scale', 'regress']],
        ['-radial_correlate',      ['yes']],
        ['-tcat_remove_first_trs', ['0']],
        ['-tshift_opts_ts',        ['-tpattern', 'alt+z2']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template_SSW.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-tlrc_NL_warped_dsets',  ['anatQQ.SID.nii', 'anatQQ.SID.aff12.1D',
                                    'anatQQ.SID_WARP.nii']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['6']],
        ['-blur_in_mask',          ['yes']],
        ['-regress_stim_times',    ['timing/times.CONTROL.txt',
                                    'timing/times.TASK.txt']],
        ['-regress_stim_labels',   ['CONTROL', 'TASK']],
        ['-regress_stim_types',    ['AM1']],
        ['-regress_basis_multi',   ['dmBLOCK']],
        ['-regress_opts_3dD',      ['-jobs', '8']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.3']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_compute_fitts', []],
        ['-regress_fout',          ['no']],
        ['-regress_3dD_stop',      []],
        ['-regress_reml_exec',     []],
        ['-regress_make_ideal_sum', ['sum_ideal.1D']],
        ['-regress_est_blur_errts', []],
        ['-regress_run_clustsim',  ['no']],
        ['-html_review_style',     ['pythonic']],
       ],
     ))

   examples.append( APExample('AP publish 2',
     source='eventually mention paper reference?',
     descrip='NARPS analysis from AFNI.',
     moddate='2020.02.10',
     keywords=['complete', 'publish', 'task'],
     header="""
              (recommended?  yes, reasonable for a complete analysis)

           An amplitude modulation task analysis.  AM1 is used for NoResp
           merely to consistently apply duration modulation.
            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['sid']],
        ['-script',                ['proc.sid']],
        ['-scr_overwrite',         []],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'blur', 'scale', 'regress']],
        ['-copy_anat',             ['anatSS.sid.nii']],
        ['-anat_has_skull',        ['no']],
        ['-anat_follower',         ['anat_w_skull', 'anat', 'anatU.sid.nii']],
        ['-anat_follower_ROI',     ['FS_wm_e', 'epi',
                                    'SUMA/mask.aseg.wm.e1.nii.gz']],
        ['-anat_follower_ROI',     ['FS_REN_epi', 'epi',
                                    'SUMA/aparc+aseg_REN_all.nii.gz']],
        ['-anat_follower_ROI',     ['FS_REN_anat', 'anat',
                                    'SUMA/aparc+aseg_REN_all.nii.gz']],
        ['-anat_follower_erode',   ['FS_wm_e']],
        ['-dsets',                 ['func/sid_task-MGT_run-01_bold.nii.gz',
                                    'func/sid_task-MGT_run-02_bold.nii.gz',
                                    'func/sid_task-MGT_run-03_bold.nii.gz',
                                    'func/sid_task-MGT_run-04_bold.nii.gz']],
        ['-tcat_remove_first_trs', ['0']],
        ['-tshift_opts_ts',        ['-tpattern', 'alt+z2']],
        ['-radial_correlate',      ['yes']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template_SSW.nii.gz']],
        ['-tlrc_NL_warp',          []],
        ['-tlrc_NL_warped_dsets',  ['anatQQ.sid.nii', 'anatQQ.sid.aff12.1D',
                                    'anatQQ.sid_WARP.nii']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['5']],
        ['-test_stim_files',       ['no']],
        ['-regress_stim_times',    ['timing/times.Resp.txt',
                                    'timing/times.NoResp.txt']],
        ['-regress_stim_labels',   ['Resp', 'NoResp']],
        ['-regress_stim_types',    ['AM2', 'AM1']],
        ['-regress_basis_multi',   ['dmBLOCK']],
        ['-regress_anaticor_fast', []],
        ['-regress_anaticor_fwhm', ['20']],
        ['-regress_anaticor_label', ['FS_wm_e']],
        ['-regress_censor_motion', ['0.3']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_motion_per_run', []],
        ['-regress_compute_fitts', []],
        ['-regress_opts_3dD',      ['-jobs', '8',
                        '-gltsym', 'SYM: Resp[1] -Resp[2]',
                        '-glt_label', '1', 'gain-loss', '-GOFORIT', '10']],
        ['-regress_opts_reml',     ['-GOFORIT']],
        ['-regress_3dD_stop',      []],
        ['-regress_reml_exec',     []],
        ['-regress_make_ideal_sum', ['sum_ideal.1D']],
        ['-regress_make_corr_vols', ['FS_wm_e']],
        ['-regress_est_blur_errts', []],
        ['-regress_run_clustsim',  ['no']],
        ['-html_review_style',     ['pythonic']],
       ],
     ))

   examples.append( APExample('AP publish 3a',
     source='AP_paper/scripts_rest/do_21_ap_ex1_align.tcsh',
     descrip='do_21_ap_ex1_align.tcsh - only perform alignment steps.',
     moddate='2024.01.26',
     keywords=['partial', 'publish', 'align'],
     header="""
              (recommended?  somewhat, for alignment only)

         This example is based on the APMULTI_Demo1_rest tree, but will come
         with a new demo package.  Probably.  Maybe.

         This is a full analysis, including:

            - reverse phase encoding (blip) distortion correction
              (-blip_forward_dset, -blip_reverse_dset)
            - EPI motion registration (to MIN_OUTLIER)
            - EPI to anatomical registration
            - non-linear anatomical to MNI template registration
              (precomputed affine+non-linear warp is provided)
            * the regress block is included only for QC

            - QC options:
                -anat_follower (with skull), (-align_opts_aea) -check_flip

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',                 ['sub-005.ex1']],
        ['-dsets',                   ['func/sub-005_rest_echo-2_bold.nii.gz']],
        ['-blip_forward_dset',       ['func/sub-005_blip-match.nii.gz[0]']],
        ['-blip_reverse_dset',       ['func/sub-005_blip-opp.nii.gz[0]']],
        ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
        ['-anat_has_skull',          ['no']],
        ['-anat_follower',           ['anat_w_skull', 'anat',
                                      'ssw/anatU.sub-005.nii']],
        ['-blocks',                  ['align', 'tlrc', 'volreg', 'regress']],
        ['-tcat_remove_first_trs',   ['4']],
        ['-align_unifize_epi',       ['local']],
        ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                      '-check_flip']],
        ['-tlrc_base',               ['MNI152_2009_template_SSW.nii.gz']],
        ['-tlrc_NL_warp',            []],
        ['-tlrc_NL_warped_dsets',    ['ssw/anatQQ.sub-005.nii',
                                      'ssw/anatQQ.sub-005.aff12.1D',
                                      'ssw/anatQQ.sub-005_WARP.nii']],
        ['-volreg_align_to',         ['MIN_OUTLIER']],
        ['-volreg_align_e2a',        []],
        ['-volreg_tlrc_warp',        []],
        ['-volreg_warp_dxyz',        ['3']],
       ],
     ))

   examples.append( APExample('AP publish 3b',
     source='AP_paper/scripts_task/do_22_ap_ex2_task.tcsh',
     descrip='do_22_ap_ex2_task.tcsh - pamenc task analysis.',
     moddate='2024.02.20',
     keywords=['complete', 'publish', 'task'],
     header="""
              (recommended?  yes, for a volumetric task analysis)

         This example is based on the AFNI_demos/AFNI_pamenc data.

         This is a full analysis, including:
            - slice time correction (alt+z2 timing pattern)
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - NL warp to MNI152_2009 template, as computed by @SSwarper
            - all registration transformations are concatenated
            - computing an EPI mask intersected with the anatomical mask
              for blurring and QC (-mask_epi_anat)
            - applying a 6 mm FWHM Gaussian blur, restricted to the EPI mask
            - voxelwise scaling to percent signal change
            - linear regression of task events using duration modulation with
              the BLOCK basis function (dmUBLOCK(-1)), where the ideal response
              height is unit for a 1 s event; stim_type AM1 is required here
            - censoring time points where motion exceeds 0.3 mm or the outlier
              fraction exceeds 5%
            - regression is performed by 3dREMLfit, accounting for voxelwise
              temporal autocorrelation in the noise
            - estimate data blur from the regression residuals using
              the mixed-model ACF function

            - QC options:
                -anat_follower (with skull), (-align_opts_aea) -check_flip,
                -radial_correlate_blocks, -volreg_compute_tsnr,
                -regress_make_ideal_sum, -html_review_style

         * input dataset names have been shortened

            """,
     trailer=""" """,
     olist = [
      ['-subj_id',                      ['sub-10506.ex2']],
      ['-dsets',                        ['func/sub-10506_pamenc_bold.nii.gz']],
      ['-copy_anat',                    ['ssw/anatSS.sub-10506.nii']],
      ['-anat_has_skull',               ['no']],
      ['-anat_follower',                ['anat_w_skull', 'anat',
                                         'ssw/anatU.sub-10506.nii']],
      ['-blocks',                       ['tshift', 'align', 'tlrc', 'volreg',
                                         'mask', 'blur', 'scale', 'regress']],
      ['-radial_correlate_blocks',      ['tcat', 'volreg', 'regress']],
      ['-tcat_remove_first_trs',        ['0']],
      ['-tshift_opts_ts',               ['-tpattern', 'alt+z2']],
      ['-align_unifize_epi',            ['local']],
      ['-align_opts_aea',               ['-giant_move', '-cost', 'lpc+ZZ',
                                         '-check_flip']],
      ['-tlrc_base',                    ['MNI152_2009_template_SSW.nii.gz']],
      ['-tlrc_NL_warp',                 []],
      ['-tlrc_NL_warped_dsets',         ['ssw/anatQQ.sub-10506.nii',
                                         'ssw/anatQQ.sub-10506.aff12.1D',
                                         'ssw/anatQQ.sub-10506_WARP.nii']],
      ['-volreg_align_to',              ['MIN_OUTLIER']],
      ['-volreg_align_e2a',             []],
      ['-volreg_tlrc_warp',             []],
      ['-volreg_warp_dxyz',             ['3.0']],
      ['-volreg_compute_tsnr',          ['yes']],
      ['-mask_epi_anat',                ['yes']],
      ['-blur_size',                    ['6']],
      ['-blur_in_mask',                 ['yes']],
      ['-regress_stim_times',           ['timing/times.CONTROL.txt',
                                         'timing/times.TASK.txt']],
      ['-regress_stim_labels',          ['CONTROL', 'TASK']],
      ['-regress_stim_types',           ['AM1']],
      ['-regress_basis_multi',          ['dmUBLOCK(-1)']],
      ['-regress_opts_3dD',             ['-jobs', '8',
                                         '-gltsym', 'SYM: TASK -CONTROL',
                                         '-glt_label', '1', 'T-C',
                                         '-gltsym',
                                         'SYM: 0.5*TASK +0.5*CONTROL',
                                         '-glt_label', '2', 'meanTC']],
      ['-regress_motion_per_run',       []],
      ['-regress_censor_motion',        ['0.3']],
      ['-regress_censor_outliers',      ['0.05']],
      ['-regress_compute_fitts',        []],
      ['-regress_fout',                 ['no']],
      ['-regress_3dD_stop',             []],
      ['-regress_reml_exec',            []],
      ['-regress_make_ideal_sum',       ['sum_ideal.1D']],
      ['-regress_est_blur_errts',       []],
      ['-regress_run_clustsim',         ['no']],
      ['-html_review_style',            ['pythonic']],
       ],
     ))

   examples.append( APExample('AP publish 3c',
     source='AP_paper/scripts_rest/do_23_ap_ex3_vol.tcsh',
     descrip='do_23_ap_ex3_vol.tcsh - rest analysis.',
     moddate='2024.08.09',
     keywords=['complete', 'physio', 'publish', 'rest'],
     header="""
              (recommended?  yes, an example of resting state analysis)

         This example is based on the APMULTI_Demo1_rest tree, to perform a
         resting state analysis with a single echo time series.

         This is a resting state processing command, including:
            - physio regression, slicewise, before any temporal or volumetric
              alterations (and per-run, though there is only 1 run here)
            - slice timing correction (notably after physio regression)
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - NL warp to MNI152_2009 template, as computed by @SSwarper
            - apply 5 mm FWHM Gaussian blur, approx 1.5*voxel size
            - all registration transformations are concatenated
            - voxelwise scaling to percent signal change
            - regression (projection) of:
                - per run motion and first differences
                - censor motion exceeding 0.2 ~mm from enorm time series,
                  or outliers exceeding 5% of brain 
            - estimate data blur from the regression residuals and the
              regression input (separately) using the mixed-model ACF function

            - QC options:
                -anat_follower (with skull), -anat_follower_ROI (FS GM mask),
                -radial_correlate_blocks, (-align_opts_aea) -check_flip,
                -volreg_compute_tsnr, -regress_make_corr_vols,
                -html_review_style

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
      ['-subj_id',                 ['sub-005.ex3']],
      ['-dsets',                   ['func/sub-005_rest_echo-2_bold.nii.gz']],
      ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
      ['-anat_has_skull',          ['no']],
      ['-anat_follower',           ['anat_w_skull', 'anat',
                                   'ssw/anatU.sub-005.nii']],
      ['-anat_follower_ROI',       ['aagm09', 'anat',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-anat_follower_ROI',       ['aegm09', 'epi',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-ROI_import',              ['BrodPijn', 'Brodmann_pijn_afni.nii.gz']],
      ['-ROI_import',              ['SchYeo7N', 'Schaefer_7N_400.nii.gz']],
      ['-blocks',                  ['ricor', 'tshift', 'align', 'tlrc',
                                    'volreg', 'mask', 'blur', 'scale',
                                    'regress']],
      ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
      ['-tcat_remove_first_trs',   ['4']],
      ['-ricor_regs',              ['physio/sub-005_rest_physio.slibase.1D']],
      ['-ricor_regs_nfirst',       ['4']],
      ['-ricor_regress_method',    ['per-run']],
      ['-align_unifize_epi',       ['local']],
      ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
      ['-tlrc_base',               ['MNI152_2009_template_SSW.nii.gz']],
      ['-tlrc_NL_warp',            []],
      ['-tlrc_NL_warped_dsets',    ['ssw/anatQQ.sub-005.nii',
                                    'ssw/anatQQ.sub-005.aff12.1D',
                                    'ssw/anatQQ.sub-005_WARP.nii']],
      ['-volreg_align_to',         ['MIN_OUTLIER']],
      ['-volreg_align_e2a',        []],
      ['-volreg_tlrc_warp',        []],
      ['-volreg_warp_dxyz',        ['3']],
      ['-volreg_compute_tsnr',     ['yes']],
      ['-mask_epi_anat',           ['yes']],
      ['-blur_size',               ['5']],
      ['-regress_apply_mot_types', ['demean', 'deriv']],
      ['-regress_motion_per_run',  []],
      ['-regress_censor_motion',   ['0.2']],
      ['-regress_censor_outliers', ['0.05']],
      ['-regress_make_corr_vols',  ['aegm09']],
      ['-regress_est_blur_epits',  []],
      ['-regress_est_blur_errts',  []],
      ['-regress_compute_tsnr_stats', ['BrodPijn', '7', '10', '12', '39',
                                       '107', '110', '112', '139']],
      ['-regress_compute_tsnr_stats', ['SchYeo7N', '161', '149', '7', '364',
                                       '367', '207']],
      ['-html_review_style',       ['pythonic']],
       ],
     ))

   examples.append( APExample('AP publish 3d',
     source='AP_paper/scripts_rest/do_24_ap_ex4_mesurf.tcsh',
     descrip='do_24_ap_ex4_mesurf.tcsh - multi-echo surface-based analysis.',
     moddate='2024.05.30',
     keywords=['complete', 'blip', 'ME', 'publish', 'rest',
               'surface', 'tedana'],
     header="""
              (recommended?  yes)

         This example is based on the APMULTI_Demo1_rest tree, to perform a
         resting state analysis on the surface with multi-echo data.

         This is a surface-based resting state processing command, including:
            - slice timing correction (using wsinc9 interpolation)
            - distortion correction using reverse blip phase encoding
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - all registration transformations are concatenated, and
              based on echo 2 (as we did not specify), but applied to all
              echoes, and resampled using a wsinc9 interpolant
            - compute a mask dataset to give to tedana (-mask_epi_anat)
              (having tedana do the projection results in masked EPI data)
            - echos are combined and then "cleaned" by tedana
            - the EPI time series are then  projected onto the surface
              (a previously computed set of surfaces, registered to the
              current anat, making a new SurfVol_Alnd_Exp anat dset)
                - might have surf data gaps, due to coverage or tedana masking
            - (light) blurring _to_ of FWHM of 4 mm is applied on the surface
            - nodewise scaling to percent signal change
            - (light, since tedana) regression (projection) of:
                - per run motion and first differences
                - censor motion exceeding 0.2 ~mm from enorm time series,
                  or outliers exceeding 5% of brain 

            - QC options:
                -anat_follower (with skull), -radial_correlate_blocks,
                (-align_opts_aea) -check_flip, -volreg_compute_tsnr,
                -html_review_style

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',                  ['sub-005.ex4']],
        ['-dsets_me_run',             ['func/sub-005_rest_echo-1_bold.nii.gz',
                                       'func/sub-005_rest_echo-2_bold.nii.gz',
                                       'func/sub-005_rest_echo-3_bold.nii.gz']],
        ['-echo_times',               ['12.5', '27.6', '42.7']],
        ['-blip_forward_dset',        ['func/sub-005_blip-match.nii.gz[0]']],
        ['-blip_reverse_dset',        ['func/sub-005_blip-opp.nii.gz[0]']],
        ['-copy_anat',                ['ssw/anatSS.sub-005.nii']],
        ['-anat_has_skull',           ['no']],
        ['-anat_follower',            ['anat_w_skull', 'anat',
                                       'ssw/anatU.sub-005.nii']],
        ['-blocks',                   ['tshift', 'align', 'volreg', 'mask',
                                       'combine', 'surf', 'blur', 'scale',
                                       'regress']],
        ['-radial_correlate_blocks',  ['tcat', 'volreg']],
        ['-tcat_remove_first_trs',    ['4']],
        ['-tshift_interp',            ['-wsinc9']],
        ['-align_unifize_epi',        ['local']],
        ['-align_opts_aea',           ['-cost', 'lpc+ZZ', '-giant_move',
                                       '-check_flip']],
        ['-volreg_align_to',          ['MIN_OUTLIER']],
        ['-volreg_align_e2a',         []],
        ['-volreg_warp_final_interp', ['wsinc5']],
        ['-volreg_compute_tsnr',      ['yes']],
        ['-mask_epi_anat',            ['yes']],
        ['-combine_method',           ['m_tedana']],
        ['-surf_anat',                ['SUMA/sub-005_SurfVol.nii']],
        ['-surf_spec',                ['SUMA/std.141.sub-005_lh.spec',
                                       'SUMA/std.141.sub-005_rh.spec']],
        ['-blur_size',                ['4']],
        ['-regress_apply_mot_types',  ['demean', 'deriv']],
        ['-regress_motion_per_run',   []],
        ['-regress_censor_motion',    ['0.2']],
        ['-regress_censor_outliers',  ['0.05']],
        ['-html_review_style',        ['pythonic']],

       ],
     ))

   examples.append( APExample('AP publish 3e',
     source='AP_paper/scripts_rest/do_35_ap_ex5_vol.tcsh',
     descrip='do_35_ap_ex5_vol.tcsh - rest analysis, with bandpassing.',
     moddate='2024.08.27',
     keywords=['complete', 'noshow', 'physio', 'publish', 'rest'],
     header="""
              (recommended?  almost, bandpassing is not preferred)

         This example is based on the APMULTI_Demo1_rest tree, to perform a
         resting state analysis with a single echo time series.

         This is the same as AP publish 3c, but with bandpassing added.

         This is a resting state processing command, including:
            - physio regression, slicewise, before any temporal or volumetric
              alterations (and per-run, though there is only 1 run here)
            - slice timing correction (notably after physio regression)
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - NL warp to MNI152_2009 template, as computed by @SSwarper
            - apply 5 mm FWHM Gaussian blur, approx 1.5*voxel size
            - all registration transformations are concatenated
            - voxelwise scaling to percent signal change
            - regression (projection) of:
                - per run motion and first differences
                - bandpassing to include frequencies in [0.01, 0.1]
                - censor motion exceeding 0.2 ~mm from enorm time series,
                  or outliers exceeding 5% of brain 
            - estimate data blur from the regression residuals and the
              regression input (separately) using the mixed-model ACF function

            - QC options:
                -anat_follower (with skull), -anat_follower_ROI (FS GM mask),
                -radial_correlate_blocks, (-align_opts_aea) -check_flip,
                -volreg_compute_tsnr, -regress_make_corr_vols,
                -html_review_style

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
      ['-subj_id',                 ['sub-005.ex5']],
      ['-dsets',                   ['func/sub-005_rest_echo-2_bold.nii.gz']],
      ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
      ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
      ['-anat_has_skull',          ['no']],
      ['-anat_follower',           ['anat_w_skull', 'anat',
                                   'ssw/anatU.sub-005.nii']],
      ['-anat_follower_ROI',       ['aagm09', 'anat',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-anat_follower_ROI',       ['aegm09', 'epi',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-ROI_import',              ['BrodPijn', 'Brodmann_pijn_afni.nii.gz']],
      ['-ROI_import',              ['SchYeo7N', 'Schaefer_7N_400.nii.gz']],
      ['-blocks',                  ['ricor', 'tshift', 'align', 'tlrc',
                                    'volreg', 'mask', 'blur', 'scale',
                                    'regress']],
      ['-tcat_remove_first_trs',   ['4']],
      ['-ricor_regs',              ['physio/sub-005_rest_physio.slibase.1D']],
      ['-ricor_regs_nfirst',       ['4']],
      ['-ricor_regress_method',    ['per-run']],
      ['-align_unifize_epi',       ['local']],
      ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
      ['-tlrc_base',               ['MNI152_2009_template_SSW.nii.gz']],
      ['-tlrc_NL_warp',            []],
      ['-tlrc_NL_warped_dsets',    ['ssw/anatQQ.sub-005.nii',
                                    'ssw/anatQQ.sub-005.aff12.1D',
                                    'ssw/anatQQ.sub-005_WARP.nii']],
      ['-volreg_align_to',         ['MIN_OUTLIER']],
      ['-volreg_align_e2a',        []],
      ['-volreg_tlrc_warp',        []],
      ['-volreg_warp_dxyz',        ['3']],
      ['-volreg_compute_tsnr',     ['yes']],
      ['-mask_epi_anat',           ['yes']],
      ['-blur_size',               ['5']],
      ['-regress_bandpass',        ['0.01', '0.1']],
      ['-regress_apply_mot_types', ['demean', 'deriv']],
      ['-regress_motion_per_run',  []],
      ['-regress_censor_motion',   ['0.2']],
      ['-regress_censor_outliers', ['0.05']],
      ['-regress_make_corr_vols',  ['aegm09']],
      ['-regress_est_blur_epits',  []],
      ['-regress_est_blur_errts',  []],
      ['-regress_compute_tsnr_stats', ['BrodPijn', '7', '10', '12', '39',
                                       '107', '110', '112', '139']],
      ['-regress_compute_tsnr_stats', ['SchYeo7N', '161', '149', '7', '364',
                                       '367', '207']],
      ['-html_review_style',       ['pythonic']],
       ],
     ))

   examples.append( APExample('AP publish 3f',
     source='AP_paper/scripts_rest/do_36_ap_ex6_vol.tcsh',
     descrip='do_36_ap_ex6_vol.tcsh - rest analysis.',
     moddate='2024.08.30',
     keywords=['complete', 'noshow', 'physio', 'publish', 'rest'],
     header="""
              (recommended?  almost, this has extra regressors)

         This example is based on the APMULTI_Demo1_rest tree, to perform a
         resting state analysis with a single echo time series.

         This is the same as AP publish 3c, but with:
            - anat_follower_erode, ROI_PC, ANATICOR, extra followers

         This is a resting state processing command, including:
            - physio regression, slicewise, before any temporal or volumetric
              alterations (and per-run, though there is only 1 run here)
            - slice timing correction (notably after physio regression)
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - NL warp to MNI152_2009 template, as computed by @SSwarper
            - apply 8 mm FWHM Gaussian blur using -blur_to_fwhm to account
              for data from multiple sites
              * this results in the EPI being masked
            - all registration transformations are concatenated
            - voxelwise scaling to percent signal change
            - regression (projection) of:
                - per run motion and first differences
                - first 3 principle components from volreg ventricles
                  (per run, though only 1 run here)
                - fast ANATICOR: weighted local white matter (voxelwise regs)
                - censor motion exceeding 0.2 ~mm from enorm time series,
                  or outliers exceeding 5% of brain 
            - estimate data blur from the regression residuals and the
              regression input (separately) using the mixed-model ACF function

            - QC options:
                -anat_follower (with skull), -anat_follower_ROI (FS GM mask),
                -radial_correlate_blocks, (-align_opts_aea) -check_flip,
                -volreg_compute_tsnr, -regress_make_corr_vols,
                -html_review_style

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
      ['-subj_id',                 ['sub-005.ex6']],
      ['-dsets',                   ['func/sub-005_rest_echo-2_bold.nii.gz']],
      ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
      ['-anat_has_skull',          ['no']],
      ['-anat_follower',           ['anat_w_skull', 'anat',
                                   'ssw/anatU.sub-005.nii']],
      ['-anat_follower_ROI',       ['aaseg', 'anat',
                                    'SUMA/aparc.a2009s+aseg_REN_all.nii.gz']],
      ['-anat_follower_ROI',       ['aeseg', 'epi',
                                    'SUMA/aparc.a2009s+aseg_REN_all.nii.gz']],
      ['-anat_follower_ROI',       ['aagm09', 'anat',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-anat_follower_ROI',       ['aegm09', 'epi',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-anat_follower_ROI',       ['aagm00', 'anat',
                                   'SUMA/aparc+aseg_REN_gmrois.nii.gz']],
      ['-anat_follower_ROI',       ['aegm00', 'epi',
                                   'SUMA/aparc+aseg_REN_gmrois.nii.gz']],
      ['-anat_follower_ROI',       ['FSvent', 'epi',
                                    'SUMA/fs_ap_latvent.nii.gz']],
      ['-anat_follower_ROI',       ['FSWe', 'epi', 'SUMA/fs_ap_wm.nii.gz']],
      ['-anat_follower_erode',     ['FSvent', 'FSWe']],
      ['-ROI_import',              ['BrodPijn', 'Brodmann_pijn_afni.nii.gz']],
      ['-ROI_import',              ['SchYeo7N', 'Schaefer_7N_400.nii.gz']],
      ['-blocks',                  ['ricor', 'tshift', 'align', 'tlrc',
                                    'volreg', 'mask', 'blur', 'scale',
                                    'regress']],
      ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
      ['-tcat_remove_first_trs',   ['4']],
      ['-ricor_regs',              ['physio/sub-005_rest_physio.slibase.1D']],
      ['-ricor_regs_nfirst',       ['4']],
      ['-ricor_regress_method',    ['per-run']],
      ['-align_unifize_epi',       ['local']],
      ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
      ['-tlrc_base',               ['MNI152_2009_template_SSW.nii.gz']],
      ['-tlrc_NL_warp',            []],
      ['-tlrc_NL_warped_dsets',    ['ssw/anatQQ.sub-005.nii',
                                    'ssw/anatQQ.sub-005.aff12.1D',
                                    'ssw/anatQQ.sub-005_WARP.nii']],
      ['-volreg_align_to',         ['MIN_OUTLIER']],
      ['-volreg_align_e2a',        []],
      ['-volreg_tlrc_warp',        []],
      ['-volreg_warp_dxyz',        ['3']],
      ['-volreg_compute_tsnr',     ['yes']],
      ['-mask_epi_anat',           ['yes']],
      ['-blur_size',               ['8']],
      ['-blur_to_fwhm',            []],
      ['-regress_apply_mot_types', ['demean', 'deriv']],
      ['-regress_motion_per_run',  []],
      ['-regress_anaticor_fast',   []],
      ['-regress_anaticor_label',  ['FSWe']],
      ['-regress_ROI_PC',          ['FSvent', '3']],
      ['-regress_ROI_PC_per_run',  ['FSvent']],
      ['-regress_censor_motion',   ['0.2']],
      ['-regress_censor_outliers', ['0.05']],
      ['-regress_make_corr_vols',  ['aeseg', 'FSvent']],
      ['-regress_est_blur_epits',  []],
      ['-regress_est_blur_errts',  []],
      ['-regress_compute_tsnr_stats', ['BrodPijn', '7', '10', '12', '39',
                                       '107', '110', '112', '139']],
      ['-regress_compute_tsnr_stats', ['SchYeo7N', '161', '149', '7', '364',
                                       '367', '207']],
      ['-html_review_style',       ['pythonic']],
       ],
     ))

   examples.append( APExample('AP publish 3g',
     source='AP_paper/scripts_rest/do_37_ap_ex7_roi.tcsh',
     descrip='do_37_ap_ex7_roi.tcsh - rest analysis.',
     moddate='2024.08.29',
     keywords=['complete', 'noshow', 'physio', 'publish', 'rest'],
     header="""
              (recommended?  yes, an example of resting state analysis)

         This example is based on the APMULTI_Demo1_rest tree, to perform a
         resting state analysis with a single echo time series.

         This is the same as AP publish 3c, but it is for an ROI analysis,
         and so does not have any 'blur' block or blur_size.

         This is a resting state processing command, including:
            - physio regression, slicewise, before any temporal or volumetric
              alterations (and per-run, though there is only 1 run here)
            - slice timing correction (notably after physio regression)
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - NL warp to MNI152_2009 template, as computed by @SSwarper
            - do not apply any blur, to avoid corrupting regions of interest
            - all registration transformations are concatenated
            - voxelwise scaling to percent signal change
            - regression (projection) of:
                - per run motion and first differences
                - censor motion exceeding 0.2 ~mm from enorm time series,
                  or outliers exceeding 5% of brain 
            - estimate data blur from the regression residuals and the
              regression input (separately) using the mixed-model ACF function

            - QC options:
                -anat_follower (with skull), -anat_follower_ROI (FS GM mask),
                -radial_correlate_blocks, (-align_opts_aea) -check_flip,
                -volreg_compute_tsnr, -regress_make_corr_vols,
                -html_review_style

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
      ['-subj_id',                 ['sub-005.ex7']],
      ['-dsets',                   ['func/sub-005_rest_echo-2_bold.nii.gz']],
      ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
      ['-anat_has_skull',          ['no']],
      ['-anat_follower',           ['anat_w_skull', 'anat',
                                   'ssw/anatU.sub-005.nii']],
      ['-anat_follower_ROI',       ['aagm09', 'anat',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-anat_follower_ROI',       ['aegm09', 'epi',
                                   'SUMA/aparc.a2009s+aseg_REN_gmrois.nii']],
      ['-ROI_import',              ['BrodPijn', 'Brodmann_pijn_afni.nii.gz']],
      ['-ROI_import',              ['SchYeo7N', 'Schaefer_7N_400.nii.gz']],
      ['-blocks',                  ['ricor', 'tshift', 'align', 'tlrc',
                                    'volreg', 'mask', 'scale', 'regress']],
      ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
      ['-tcat_remove_first_trs',   ['4']],
      ['-ricor_regs',              ['physio/sub-005_rest_physio.slibase.1D']],
      ['-ricor_regs_nfirst',       ['4']],
      ['-ricor_regress_method',    ['per-run']],
      ['-align_unifize_epi',       ['local']],
      ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
      ['-tlrc_base',               ['MNI152_2009_template_SSW.nii.gz']],
      ['-tlrc_NL_warp',            []],
      ['-tlrc_NL_warped_dsets',    ['ssw/anatQQ.sub-005.nii',
                                    'ssw/anatQQ.sub-005.aff12.1D',
                                    'ssw/anatQQ.sub-005_WARP.nii']],
      ['-volreg_align_to',         ['MIN_OUTLIER']],
      ['-volreg_align_e2a',        []],
      ['-volreg_tlrc_warp',        []],
      ['-volreg_warp_dxyz',        ['3']],
      ['-volreg_compute_tsnr',     ['yes']],
      ['-mask_epi_anat',           ['yes']],
      ['-regress_apply_mot_types', ['demean', 'deriv']],
      ['-regress_motion_per_run',  []],
      ['-regress_censor_motion',   ['0.2']],
      ['-regress_censor_outliers', ['0.05']],
      ['-regress_make_corr_vols',  ['aegm09']],
      ['-regress_est_blur_epits',  []],
      ['-regress_est_blur_errts',  []],
      ['-regress_compute_tsnr_stats', ['BrodPijn', '7', '10', '12', '39',
                                       '107', '110', '112', '139']],
      ['-regress_compute_tsnr_stats', ['SchYeo7N', '161', '149', '7', '364',
                                       '367', '207']],
      ['-html_review_style',       ['pythonic']],
       ],
     ))

   examples.append( APExample('AP publish 3h',
     source='AP_paper/scripts_rest/do_38_ap_ex8_mesurf_oc.tcsh',
     descrip='do_38_ap_ex8_mesurf_oc.tcsh - multi-echo surface analysis.',
     moddate='2024.09.04',
     keywords=['complete', 'blip', 'ME', 'noshow', 'publish',
               'rest', 'surface'],
     header="""
              (recommended?  yes)

         This example is based on the APMULTI_Demo1_rest tree, to perform a
         resting state analysis on the surface with multi-echo data.

         This is the same as AP publish 3d, but it uses AFNI's OC echo combine
         method rather than tedana (-combine_method).  Because a mask is not
         required, the mask block and option have also been removed (though it
         is certainly fine to keep them).

         This is a surface-based resting state processing command, including:
            - slice timing correction (using wsinc9 interpolation)
            - distortion correction using reverse blip phase encoding
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - all registration transformations are concatenated, and
              based on echo 2 (as we did not specify), but applied to all
              echoes, and resampled using a wsinc9 interpolant
            - echos are combined and then "cleaned" by tedana
            - the EPI time series are then  projected onto the surface
              (a previously computed set of surfaces, registered to the
              current anat, making a new SurfVol_Alnd_Exp anat dset)
                - might have surf data gaps, due to coverage
            - (light) blurring _to_ of FWHM of 4 mm is applied on the surface
            - nodewise scaling to percent signal change
            - (light, since tedana) regression (projection) of:
                - per run motion and first differences
                - censor motion exceeding 0.2 ~mm from enorm time series,
                  or outliers exceeding 5% of brain 

            - QC options:
                -anat_follower (with skull), -radial_correlate_blocks,
                (-align_opts_aea) -check_flip, -volreg_compute_tsnr,
                -html_review_style

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',                  ['sub-005.ex8']],
        ['-dsets_me_run',             ['func/sub-005_rest_echo-1_bold.nii.gz',
                                       'func/sub-005_rest_echo-2_bold.nii.gz',
                                       'func/sub-005_rest_echo-3_bold.nii.gz']],
        ['-echo_times',               ['12.5', '27.6', '42.7']],
        ['-blip_forward_dset',        ['func/sub-005_blip-match.nii.gz[0]']],
        ['-blip_reverse_dset',        ['func/sub-005_blip-opp.nii.gz[0]']],
        ['-copy_anat',                ['ssw/anatSS.sub-005.nii']],
        ['-anat_has_skull',           ['no']],
        ['-anat_follower',            ['anat_w_skull', 'anat',
                                       'ssw/anatU.sub-005.nii']],
        ['-blocks',                   ['tshift', 'align', 'volreg',
                                       'combine', 'surf', 'blur', 'scale',
                                       'regress']],
        ['-radial_correlate_blocks',  ['tcat', 'volreg']],
        ['-tcat_remove_first_trs',    ['4']],
        ['-tshift_interp',            ['-wsinc9']],
        ['-align_unifize_epi',        ['local']],
        ['-align_opts_aea',           ['-cost', 'lpc+ZZ', '-giant_move',
                                       '-check_flip']],
        ['-volreg_align_to',          ['MIN_OUTLIER']],
        ['-volreg_align_e2a',         []],
        ['-volreg_warp_final_interp', ['wsinc5']],
        ['-volreg_compute_tsnr',      ['yes']],
        ['-combine_method',           ['OC']],
        ['-surf_anat',                ['SUMA/sub-005_SurfVol.nii']],
        ['-surf_spec',                ['SUMA/std.141.sub-005_lh.spec',
                                       'SUMA/std.141.sub-005_rh.spec']],
        ['-blur_size',                ['4']],
        ['-regress_apply_mot_types',  ['demean', 'deriv']],
        ['-regress_motion_per_run',   []],
        ['-regress_censor_motion',    ['0.2']],
        ['-regress_censor_outliers',  ['0.05']],
        ['-html_review_style',        ['pythonic']],

       ],
     ))

   examples.append( APExample('AP publish 3i',
     source='AP_paper/scripts_rest/do_39_ap_ex9_mevol_oc.tcsh',
     descrip='do_39_ap_ex9_mevol_oc.tcsh - ME volume rest analysis.',
     moddate='2024.08.27',
     keywords=['blip', 'complete', 'ME', 'publish', 'rest'],
     header="""
              (recommended?  yes, an example of resting state analysis)

         This example is based on the APMULTI_Demo1_rest tree, to perform a
         resting state analysis with a multi-echo time series.

         This is a multi-echo resting state processing command, including:
            - 1 run with 3 echoes of EPI time series data
            - reverse phase encoding distortion correction
            - slice timing correction
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - NL warp to MNI152_2009 template, as computed by sswarper2
            - apply 4 mm FWHM Gaussian blur, approx 1.5*voxel size,
              but lower because of multi-echo noise cancellation
            - all registration transformations are concatenated
            - combine echoes using the base OC (optimally combined) method
            - voxelwise scaling to percent signal change
            - regression (projection) of:
                - per run motion and first differences
                - censor motion exceeding 0.2 ~mm from enorm time series,
                  or outliers exceeding 5% of brain 
            - estimate data blur from the regression residuals and the
              regression input (separately) using the mixed-model ACF function

            - QC options:
                -anat_follower (with skull), -anat_follower_ROI (Brodmann
                 and Schaefer ROIs) for TSNR statistics
                -radial_correlate_blocks, (-align_opts_aea) -check_flip,
                -volreg_compute_tsnr, -html_review_style

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
      ['-subj_id',                 ['sub-005.ex9']],
      ['-dsets_me_run',            ['func/sub-005_rest_r1_e1_bold.nii.gz',
                                    'func/sub-005_rest_r1_e2_bold.nii.gz',
                                    'func/sub-005_rest_r1_e3_bold.nii.gz']],
      ['-echo_times',              ['12.5', '27.6', '42.7']],
      ['-blip_forward_dset',       ['func/sub-005_blip-match.nii.gz[0]']],
      ['-blip_reverse_dset',       ['func/sub-005_blip-opp.nii.gz[0]']],
      ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
      ['-anat_has_skull',          ['no']],
      ['-anat_follower',           ['anat_w_skull', 'anat',
                                   'ssw/anatU.sub-005.nii']],
      ['-ROI_import',              ['BrodPijn', 'Brodmann_pijn_afni.nii.gz']],
      ['-ROI_import',              ['SchYeo7N', 'Schaefer_7N_400.nii.gz']],
      ['-blocks',                  ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'combine', 'blur', 'scale',
                                    'regress']],
      ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
      ['-tcat_remove_first_trs',   ['4']],
      ['-align_unifize_epi',       ['local']],
      ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
      ['-tlrc_base',               ['MNI152_2009_template_SSW.nii.gz']],
      ['-tlrc_NL_warp',            []],
      ['-tlrc_NL_warped_dsets',    ['ssw/anatQQ.sub-005.nii',
                                    'ssw/anatQQ.sub-005.aff12.1D',
                                    'ssw/anatQQ.sub-005_WARP.nii']],
      ['-volreg_align_to',         ['MIN_OUTLIER']],
      ['-volreg_align_e2a',        []],
      ['-volreg_tlrc_warp',        []],
      ['-volreg_warp_dxyz',        ['3']],
      ['-volreg_compute_tsnr',     ['yes']],
      ['-mask_epi_anat',           ['yes']],
      ['-combine_method',          ['OC']],
      ['-blur_size',               ['4']],
      ['-regress_censor_motion',   ['0.2']],
      ['-regress_censor_outliers', ['0.05']],
      ['-regress_apply_mot_types', ['demean', 'deriv']],
      ['-regress_motion_per_run',  []],
      ['-regress_est_blur_epits',  []],
      ['-regress_est_blur_errts',  []],
      ['-regress_compute_tsnr_stats', ['BrodPijn', '7', '10', '12', '39',
                                       '107', '110', '112', '139']],
      ['-regress_compute_tsnr_stats', ['SchYeo7N', '161', '149', '7', '364',
                                       '367', '207']],
      ['-html_review_style',       ['pythonic']],
     ],
     ))

   examples.append( APExample('AP publish 3j',
     source='AP_paper/scripts_task/do_40_ap_ex10_task_bder.tcsh',
     descrip='do_40_ap_ex10_task_bder.tcsh - pamenc task analysis.',
     moddate='2024.02.20',
     keywords=['complete', 'noshow', 'publish', 'task'],
     header="""
              (recommended?  yes, for a volumetric task analysis)

         This example is based on the AFNI_demos/AFNI_pamenc data.

         This is the same as AP publish 3b, but it sets taskname as a uvar
         and outputs an extra bids derivative tree.

         This is a full analysis, including:
            - slice time correction (alt+z2 timing pattern)
            - EPI registration to MIN_OUTLIER vr_base volume
            - EPI/anat alignment, with -align_unifize_epi local
            - NL warp to MNI152_2009 template, as computed by @SSwarper
            - all registration transformations are concatenated
            - computing an EPI mask intersected with the anatomical mask
              for blurring and QC (-mask_epi_anat)
            - applying a 6 mm FWHM Gaussian blur, restricted to the EPI mask
            - voxelwise scaling to percent signal change
            - linear regression of task events using duration modulation with
              the BLOCK basis function (dmUBLOCK(-1)), where the ideal response
              height is unit for a 1 s event; stim_type AM1 is required here
            - censoring time points where motion exceeds 0.3 mm or the outlier
              fraction exceeds 5%
            - regression is performed by 3dREMLfit, accounting for voxelwise
              temporal autocorrelation in the noise
            - estimate data blur from the regression residuals using
              the mixed-model ACF function

            - QC options:
                -anat_follower (with skull), (-align_opts_aea) -check_flip,
                -radial_correlate_blocks, -volreg_compute_tsnr,
                -regress_make_ideal_sum, -html_review_style

         * input dataset names have been shortened

            """,
     trailer=""" """,
     olist = [
      ['-subj_id',                      ['sub-10506.ex10']],
      ['-uvar',                         ['taskname', 'pamenc']],
      ['-dsets',                        ['func/sub-10506_pamenc_bold.nii.gz']],
      ['-copy_anat',                    ['ssw/anatSS.sub-10506.nii']],
      ['-anat_has_skull',               ['no']],
      ['-anat_follower',                ['anat_w_skull', 'anat',
                                         'ssw/anatU.sub-10506.nii']],
      ['-blocks',                       ['tshift', 'align', 'tlrc', 'volreg',
                                         'mask', 'blur', 'scale', 'regress']],
      ['-radial_correlate_blocks',      ['tcat', 'volreg', 'regress']],
      ['-tcat_remove_first_trs',        ['0']],
      ['-tshift_opts_ts',               ['-tpattern', 'alt+z2']],
      ['-align_unifize_epi',            ['local']],
      ['-align_opts_aea',               ['-giant_move', '-cost', 'lpc+ZZ',
                                         '-check_flip']],
      ['-tlrc_base',                    ['MNI152_2009_template_SSW.nii.gz']],
      ['-tlrc_NL_warp',                 []],
      ['-tlrc_NL_warped_dsets',         ['ssw/anatQQ.sub-10506.nii',
                                         'ssw/anatQQ.sub-10506.aff12.1D',
                                         'ssw/anatQQ.sub-10506_WARP.nii']],
      ['-volreg_align_to',              ['MIN_OUTLIER']],
      ['-volreg_align_e2a',             []],
      ['-volreg_tlrc_warp',             []],
      ['-volreg_warp_dxyz',             ['3.0']],
      ['-volreg_compute_tsnr',          ['yes']],
      ['-mask_epi_anat',                ['yes']],
      ['-blur_size',                    ['6']],
      ['-blur_in_mask',                 ['yes']],
      ['-regress_stim_times',           ['timing/times.CONTROL.txt',
                                         'timing/times.TASK.txt']],
      ['-regress_stim_labels',          ['CONTROL', 'TASK']],
      ['-regress_stim_types',           ['AM1']],
      ['-regress_basis_multi',          ['dmUBLOCK(-1)']],
      ['-regress_opts_3dD',             ['-jobs', '8',
                                         '-gltsym', 'SYM: TASK -CONTROL',
                                         '-glt_label', '1', 'T-C',
                                         '-gltsym',
                                         'SYM: 0.5*TASK +0.5*CONTROL',
                                         '-glt_label', '2', 'meanTC']],
      ['-regress_motion_per_run',       []],
      ['-regress_censor_motion',        ['0.3']],
      ['-regress_censor_outliers',      ['0.05']],
      ['-regress_compute_fitts',        []],
      ['-regress_fout',                 ['no']],
      ['-regress_3dD_stop',             []],
      ['-regress_reml_exec',            []],
      ['-regress_make_ideal_sum',       ['sum_ideal.1D']],
      ['-regress_est_blur_errts',       []],
      ['-regress_run_clustsim',         ['no']],
      ['-html_review_style',            ['pythonic']],
      ['-bids_deriv',                   ['yes']],
       ],
     ))

   return examples

def egs_demo():
   """AP demo examples
   """

   examples =  []

   examples.append( APExample('AP demo 1a',
     source='ap_run_simple_rest.tcsh',
     descrip='for QC, ap_run_simple_rest.tcsh with EPI and anat',
     moddate='2024.02.20',
     keywords=['rest'],
     header="""
              (recommended?  yes, for quick quality control)

         This example was generated by running ap_run_simple_rest.tcsh,
         providing a single subject anat and (3 runs of) EPI.  It could
         be generated (and run) using the following:

            cd AFNI_data6/FT_analysis/FT
            ap_run_simple_rest.tcsh -subjid FT -run_proc \\
              -anat FT_anat+orig -epi FT_epi_r*.HEAD

         This is highly recommended as a tool for quick quality control to be
         run on all EPI data right out of the scanner.  It is fine to run on
         task data, but without worrying about the actual task regression.

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT']],
        ['-dsets',                 ['FT/FT_epi_r1+orig.HEAD',
                                    'FT/FT_epi_r2+orig.HEAD',
                                    'FT/FT_epi_r3+orig.HEAD']],
        ['-copy_anat',             ['FT_anat+orig']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
        ['-tcat_remove_first_trs', ['2']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template_SSW.nii.gz']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-volreg_compute_tsnr',   ['yes']],
        ['-mask_epi_anat',         ['yes']],
        ['-blur_size',             ['6']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.25']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-regress_make_ideal_sum',  ['sum_ideal.1D']],
        ['-html_review_style',     ['pythonic']],
       ],
     ))

   examples.append( APExample('AP demo 1b',
     source='ap_run_simple_rest.tcsh',
     descrip='for QC, ap_run_simple_rest.tcsh with no anat',
     moddate='2022.11.23',
     keywords=['rest'],
     header="""
              (recommended?  yes, for quick quality control of EPI)

         This example was generated by running ap_run_simple_rest.tcsh,
         providing only 3 runs of EPI data.  It could be generated (and run)
         using the following:

            cd AFNI_data6/FT_analysis/FT
            ap_run_simple_rest.tcsh -subjid FT -run_proc -epi FT_epi_r*.HEAD

         No anatomical volume is included, excluding many options from
         example simple_rest_QC.

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['FT']],
        ['-script',                ['proc.FT']],
        ['-out_dir',               ['FT.results']],
        ['-dsets',                 ['FT/FT_epi_r1+orig.HEAD',
                                    'FT/FT_epi_r2+orig.HEAD',
                                    'FT/FT_epi_r3+orig.HEAD']],
        ['-blocks',                ['tshift', 'volreg', 'mask',
                                    'blur', 'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg']],
        ['-tcat_remove_first_trs', ['2']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_compute_tsnr',   ['yes']],
        ['-blur_size',             ['6']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.25']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-regress_make_ideal_sum',  ['sum_ideal.1D']],
        ['-html_review_style',     ['pythonic']],
       ],
     ))

   examples.append( APExample('AP demo 1c',
     source='ap_run_simple_rest_me.tcsh',
     descrip='for QC, ap_run_simple_rest_me.tcsh with ME EPI and anat',
     moddate='2024.08.09',
     keywords=['rest', 'ME'],
     header="""
              (recommended?  yes, for quick quality control)

         This example was generated by running ap_run_simple_rest_me.tcsh,
         providing a single subject anat, EPI (1 run of 3 echoes), and
         the 3 echo times.

         It could be generated using the following, where the dataset names
         have been slightly truncated to save screen space:

            cd data_00_basic/sub-005/ses-01
            ap_run_simple_rest_me.tcsh                                  \\
                -run_ap                                                 \\
                -subjid      sub-005                                    \\
                -nt_rm       4                                          \\
                -anat        anat/sub-005*mprage_run-1_T1w.nii.gz       \\
                -epi_me_run  func/sub-005*rest*bold.nii.gz              \\
                -echo_times  12.5 27.6 42.7                             \\
                -template    MNI152_2009_template_SSW.nii.gz

         This is highly recommended as a tool for quick quality control to be
         run on all EPI data right out of the scanner.

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',               ['sub-005']],
        ['-dsets_me_run',          ['func/sub-005_rest_r1_e1_bold.nii.gz',
                                    'func/sub-005_rest_r1_e2_bold.nii.gz',
                                    'func/sub-005_rest_r1_e3_bold.nii.gz']],
        ['-echo_times',            ['12.5', '27.6', '42.7']],
        ['-reg_echo',              ['2']],
        ['-copy_anat',             ['anat/sub-005_mprage_run-1_T1w.nii.gz']],
        ['-blocks',                ['tshift', 'align', 'tlrc', 'volreg',
                                    'mask', 'combine', 'blur', 'scale',
                                    'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg', 'regress']],
        ['-tcat_remove_first_trs', ['4']],
        ['-tshift_interp',         ['-wsinc9']],
        ['-align_unifize_epi',     ['local']],
        ['-align_opts_aea',        ['-cost', 'lpc+ZZ', '-giant_move',
                                    '-check_flip']],
        ['-tlrc_base',             ['MNI152_2009_template_SSW.nii.gz']],
        ['-volreg_align_to',       ['MIN_OUTLIER']],
        ['-volreg_align_e2a',      []],
        ['-volreg_tlrc_warp',      []],
        ['-volreg_warp_final_interp', ['wsinc5']],
        ['-volreg_compute_tsnr',   ['yes']],
        ['-mask_epi_anat',         ['yes']],
        ['-combine_method',        ['OC']],
        ['-blur_size',             ['4']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run', []],
        ['-regress_censor_motion', ['0.25']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_est_blur_epits', []],
        ['-regress_est_blur_errts', []],
        ['-regress_make_ideal_sum',  ['sum_ideal.1D']],
        ['-html_review_style',     ['pythonic']],
       ],
     ))

   examples.append( APExample('AP demo 2a',
     source='APMULTI_Demo1_rest/scripts_desktop/do_20_ap_se.tcsh',
     descrip='do_20_ap_se.tcsh - one way to process rest data',
     moddate='2023.04.19',
     keywords=['complete', 'rest'],
     header="""
              (recommended?  somewhat, includes tissue-based regression)

         This example is part of the APMULTI_Demo1_rest tree, installable by
         running :

            @Install_APMULTI_Demo1_rest

         This is a sample rest processing command, including:

            - despike block for high motion subjects

            - QC options:
                -radial_correlate_blocks, (-align_opts_aea) -check_flip
                -volreg_compute_tsnr, -regress_make_corr_vols,
                -html_review_style, -anat_follower_ROI (some are for QC)

            - non-linear template alignment (precomputed warp is provided)

            - noise removal of:
                - motion and derivatives, per run
                - ventricle principal components (top 3 per run)
                - fast ANATICOR
                - censoring for both motion and outliers

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',                 ['sub-005']],
        ['-dsets',                   ['func/sub-005_rest_echo-2_bold.nii.gz']],
        ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
        ['-anat_has_skull',          ['no']],
        ['-anat_follower',           ['anat_w_skull', 'anat',
                                      'ssw/anatU.sub-005.nii']],
        ['-anat_follower_ROI',       ['aaseg', 'anat',
                                      'SUMA/aparc.a2009s+aseg_REN_all.nii.gz']],
        ['-anat_follower_ROI',       ['aeseg', 'epi',
                                      'SUMA/aparc.a2009s+aseg_REN_all.nii.gz']],
        ['-anat_follower_ROI',       ['FSvent', 'epi',
                                      'SUMA/fs_ap_latvent.nii.gz']],
        ['-anat_follower_ROI',       ['FSWe', 'epi', 'SUMA/fs_ap_wm.nii.gz']],
        ['-anat_follower_erode',     ['FSvent', 'FSWe']],
        ['-blocks',                  ['despike', 'tshift', 'align', 'tlrc',
                                      'volreg', 'mask', 'blur', 'scale',
                                      'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg']],
        ['-tcat_remove_first_trs',   ['4']],
        ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                      '-check_flip']],
        ['-tlrc_base',               ['MNI152_2009_template_SSW.nii.gz']],
        ['-tlrc_NL_warp',            []],
        ['-tlrc_NL_warped_dsets',    ['ssw/anatQQ.sub-005.nii',
                                      'ssw/anatQQ.sub-005.aff12.1D',
                                      'ssw/anatQQ.sub-005_WARP.nii']],
        ['-volreg_align_to',         ['MIN_OUTLIER']],
        ['-volreg_align_e2a',        []],
        ['-volreg_tlrc_warp',        []],
        ['-volreg_warp_dxyz',        ['3']],
        ['-volreg_compute_tsnr',     ['yes']],
        ['-mask_epi_anat',           ['yes']],
        ['-blur_size',               ['5']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run',  []],
        ['-regress_anaticor_fast',   []],
        ['-regress_anaticor_label',  ['FSWe']],
        ['-regress_ROI_PC',          ['FSvent', '3']],
        ['-regress_ROI_PC_per_run',  ['FSvent']],
        ['-regress_censor_motion',   ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_make_corr_vols',  ['aeseg', 'FSvent']],
        ['-regress_est_blur_epits',  []],
        ['-regress_est_blur_errts',  []],
        ['-html_review_style',       ['pythonic']],
       ],
     ))

   examples.append( APExample('AP demo 2b',
     source='APMULTI_Demo1_rest/scripts_desktop/do_44_ap_me_bTs.tcsh',
     descrip='do_44_ap_me_bTs.tcsh - ME surface rest with tedana',
     moddate='2024.01.04',
     keywords=['complete', 'blip', 'ME', 'rest', 'surface', 'tedana'],
     header="""
              (recommended?  yes)

         This example is based on the APMULTI_Demo1_rest tree, installable by
         running :

            @Install_APMULTI_Demo1_rest

         This is a sample rest processing command, including:

            - reverse phase encoding (blip) distortion correction
              (-blip_forward_dset, -blip_reverse_dset)
            - multi-echo EPI (-dsets_me_run, -echo_times)
            - MEICA-group tedana usage
              (-combine_method m_tedana, -volreg_warp_final_interp wsinc5)

            - surface-based analysis (-surf_anat, -surf_spec)

            - despike block for high motion subjects

            - QC options:
                -radial_correlate_blocks, -align_opts_aea -check_flip,
                -volreg_compute_tsnr, -regress_make_corr_vols,
                -anat_follower anat_w_skull, -anat_follower_ROI (some for QC),
                -html_review_style

            - noise removal of:
                - tedana
                - motion and derivatives, per run
                - censoring for both motion and outliers

         * input dataset names have been shortened to protect the margins

            """,
     trailer=""" """,
     olist = [
        ['-subj_id',                 ['sub-005']],
        ['-dsets_me_run',            ['func/sub-005_rest_echo-1_bold.nii.gz',
                                      'func/sub-005_rest_echo-2_bold.nii.gz',
                                      'func/sub-005_rest_echo-3_bold.nii.gz']],
        ['-echo_times',              ['12.5', '27.6', '42.7']],
        ['-blip_forward_dset',       ['func/sub-005_blip-match.nii.gz[0]']],
        ['-blip_reverse_dset',       ['func/sub-005_blip-opp.nii.gz[0]']],
        ['-copy_anat',               ['ssw/anatSS.sub-005.nii']],
        ['-anat_has_skull',          ['no']],
        ['-anat_follower',           ['anat_w_skull', 'anat',
                                      'ssw/anatU.sub-005.nii']],
        ['-anat_follower_ROI',       ['aaseg', 'anat',
                                      'SUMA/aparc.a2009s+aseg_REN_all.nii.gz']],
        ['-anat_follower_ROI',       ['aeseg', 'epi',
                                      'SUMA/aparc.a2009s+aseg_REN_all.nii.gz']],
        ['-anat_follower_ROI',       ['FSvent', 'epi',
                                      'SUMA/fs_ap_latvent.nii.gz']],
        ['-anat_follower_ROI',       ['FSWe', 'epi', 'SUMA/fs_ap_wm.nii.gz']],
        ['-anat_follower_erode',     ['FSvent', 'FSWe']],
        ['-blocks',                  ['despike', 'tshift', 'align', 'volreg',
                                      'mask', 'combine', 'surf', 'blur',
                                      'scale', 'regress']],
        ['-radial_correlate_blocks', ['tcat', 'volreg']],
        ['-tcat_remove_first_trs',   ['4']],
        ['-tshift_interp',           ['-wsinc9']],
        ['-align_unifize_epi',       ['local']],
        ['-align_opts_aea',          ['-cost', 'lpc+ZZ', '-giant_move',
                                      '-check_flip']],
        ['-volreg_align_to',         ['MIN_OUTLIER']],
        ['-volreg_align_e2a',        []],
        ['-volreg_warp_final_interp',['wsinc5']],
        ['-volreg_compute_tsnr',     ['yes']],
        ['-mask_epi_anat',           ['yes']],
        ['-combine_method',          ['m_tedana']],
        ['-surf_anat',               ['SUMA/sub-005_SurfVol.nii']],
        ['-surf_spec',               ['SUMA/std.141.sub-005_lh.spec',
                                      'SUMA/std.141.sub-005_rh.spec']],
        ['-blur_size',               ['4']],
        ['-regress_apply_mot_types', ['demean', 'deriv']],
        ['-regress_motion_per_run',  []],
        ['-regress_censor_motion',   ['0.2']],
        ['-regress_censor_outliers', ['0.05']],
        ['-regress_make_corr_vols',  ['aeseg', 'FSvent']],
        ['-html_review_style',       ['pythonic']],
       ],
     ))

   return examples

def egs_short():
   """AP short examples (examples of only partial processing)
   """

   examples =  []

   return examples

def find_eg(name):
   """try to find a matching ap_examples instance

      consider things like:
        if failure, and if 'name' is substring of only 1 example, return it

      return None on failure
   """
   global ap_examples
   populate_examples()

   if len(name) < 1: return None

   # use lower case for searching names
   nlist = [eg.name.lower() for eg in ap_examples]
   lname = name.lower()

   # if lname is exactly (no case) in nlist, return the respective example
   if lname in nlist:
      return ap_examples[nlist.index(lname)]

   # otherwise, try harder

   # If number (possibly with trailing a,b,c,...) search after 'example'.
   if lname[0].isdigit():
      tname = 'example ' + lname
      if tname in nlist:
         return ap_examples[nlist.index(tname)]

   # otherwise, just see if there is a unique substring match
   ind = unique_substr_name_index(lname, nlist)
   if ind >= 0:
      return ap_examples[ind]

   return None

def unique_substr_name_index(nsub, nlist, endswith=0):
   """search for nsub in nlist, where nsub can be a substring
      if endswith, use name.endswith(), rather than name.find()
      return : index >= 0 on success
             : -1, if no match is found
             : -2, if the match is not uniq
   """
   findex = -1
   for ind, name in enumerate(nlist):
      # first, check to see if nsub matches name
      found = 0
      if endswith:
         if name.endswith(nsub):
            found = 1
      elif name.find(nsub) >= 0:
         found = 1

      # if not, just move along
      if not found:
         continue

      # if so, fail on non-unique
      if findex >= 0:
         # not unique
         return -2

      # we have a match, keep looking for non-uniqueness
      findex = ind

   # we have either failed to find a match, or have a unique one
   # - either way, return findex
   return findex

def show_enames(verb=1):
   """list all ap_example names
   """
   global ap_examples
   populate_examples()

   if len(ap_examples) == 0:
      return

   nlist = [eg.name for eg in ap_examples]

   # basic: show list
   if verb == 0:
      print("%s" % nlist)
      return

   # nicer: show pretty list
   if verb == 1:
      istr = ' '*2
      jstr = '\n%s' % istr
      print("%s%s\n" % (istr, jstr.join(nlist)))
      return

   # so verb > 1
   maxn = max([len(name) for name in nlist])
   indent = ' '*2
   for eg in ap_examples:
      if verb <= 2: dstr = ''
      else:         dstr = '%s : ' % eg.moddate
      print("%s%-*s : %s%s" % (indent, maxn, eg.name, dstr, eg.descrip))

def compare_eg_pair(eg1, eg2, eskip=[], verb=1):
   """similar to compare(), above, but compare 2 known examples"""

   # set names and get an instance for the first entry, eg1
   if isinstance(eg1, APExample):
      n1 = eg1.name
      eg = eg1
   else:
      n1 = eg1
      eg = find_eg(eg1)
      if not isinstance(eg, APExample):
         print("** compare_eg_pair: failed to find example named '%s'" % n1)
         return

   # so eg is the instance of eg1, just use it
   return eg.compare(eg2, eskip=eskip, verb=verb)

def display_eg_all(aphelp=1, source='', verb=0):
   """display the examples array if someone wants it
         aphelp :  1   show only AP help examples
                   0   show only others
                  -1   show all (leaving higher numbers for other cases)
         source : if set, restrict to matching
         verb   : verbosity level, pass on to eg.display
   """
   global ap_examples

   for eg in ap_examples:
      # skip what we do not want to show
      if aphelp >= 0 and aphelp != eg.aphelp:
         continue
      if source != '' and source != eg.source:
         continue
      eg.display(verb=verb)

if __name__ == '__main__':
   print('** this is not a main module')
   sys.exit(1)

