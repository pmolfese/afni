#!/usr/bin/env python

# python3 status: started

# afni_util.py : general utilities for python programs

# ------------------------------------------------------
# no longer usable as a main: see afni_python_wrapper.py
# ------------------------------------------------------

import sys, os, math, copy
from afnipy import afni_base as BASE
from afnipy import lib_textdata as TD
from afnipy import lib_format_cmd_str as lfcs
import glob
import re
from   datetime   import datetime

# global lists for basis functions
basis_known_resp_l = ['GAM', 'BLOCK', 'dmBLOCK', 'dmUBLOCK', 'SPMG1',
                      'WAV', 'MION']
basis_one_regr_l   = basis_known_resp_l[:]
basis_one_regr_l.append('MION')
stim_types_one_reg = ['file', 'AM1', 'times']
g_valid_shells = ['csh','tcsh','sh','bash','zsh']
g_text_file_suffix_list = ['1D', 'json', 'niml', 'tsv', 'txt']
g_valid_slice_patterns = [ # synonymous pairs      # z2-types
                           'zero',  'simult',
                           'seq+z', 'seqplus',
                           'seq-z', 'seqminus',
                           'alt+z', 'altplus',     'alt+z2',    
                           'alt-z', 'altminus',    'alt-z2',    
                         ]
g_tpattern_irreg = 'irregular'


# this file contains various afni utilities   17 Nov 2006 [rickr]

def change_path_basename(orig, prefix='', suffix='', append=0):
    """given a path (leading directory or not) swap the trailing
       filename with the passed prefix and suffix
          e.g. C_P_B('my/dir/pickles.yummy','toast','.1D')
                 --> 'my/dir/toast.1D' 
       or with append...
          e.g. C_P_B('my/dir/pickles.yummy','toast','.1D', append=1)
                 --> 'my/dir/toastpickles.yummy.1D'
       or maybe we want another dot then
          e.g. C_P_B('my/dir/pickles.yummy','toast.','.1D', append=1)
                 --> 'my/dir/toast.pickles.yummy.1D'
    """
    if not orig: return ''
    if not prefix and not suffix: return orig

    (head, tail) = os.path.split(orig)
    if append: tail = '%s%s%s' % (prefix, tail, suffix)
    else:      tail = '%s%s' % (prefix, suffix)

    if head == '': return tail
    return "%s/%s" % (head, tail)

# write text to a file
def write_text_to_file(fname, tdata, mode='w', wrap=0, wrapstr='\\\n', exe=0,
                       method='rr'):
    """write the given text (tdata) to the given file
          fname   : file name to write (or append) to
          dtata   : text to write
          mode    : optional write mode 'w' or 'a' [default='w']
          wrap    : optional wrap flag
          wrapstr : optional wrap string: if wrap, apply this string
          exe     : whether to make file executable
          method  : either 'rr' or 'pt'

       return 0 on success, 1 on error
    """

    if not tdata or not fname:
        print("** WTTF: missing text or filename")
        return 1

    if wrap:
       tdata = add_line_wrappers(tdata, wrapstr, method=method, verb=1)
    
    if fname == 'stdout' or fname == '-':
       fp = sys.stdout
    elif fname == 'stderr':
       fp = sys.stderr
    else:
       try:
           fp = open(fname, mode)
       except:
           print("** failed to open text file '%s' for writing" % fname)
           return 1

    fp.write(tdata)

    if fp != sys.stdout and fp != sys.stderr:
       fp.close()
       if exe:
           try: code = eval('0o755')
           except: code = eval('0755')
           try:
               os.chmod(fname, code)
           except:
               print("** failed chmod 755 on %s" % fname)

    return 0

def wrap_file_text(infile='stdin', outfile='stdout', method='pt'):
   """make a new file with line wrappers                14 Mar 2014

      The default parameters makes it easy to process as a stream:

          cat INPUT | afni_python_wrapper.py -eval 'wrap_file_text()'
                or
          afni_python_wrapper.py -eval 'wrap_file_text()' < INPUT > OUTPUT
                or
          afni_python_wrapper.py -eval 'wrap_file_text("INPUT", "OUTPUT")'
                or
          afni_python_wrapper.py -eval "wrap_file_text('$f1', '$f2')"
   """

   tdata = read_text_file(fname=infile, lines=0, strip=0)
   if tdata != '': write_text_to_file(outfile, tdata, wrap=1, method=method)
   

def read_text_file(fname='stdin', lines=1, strip=1, nonl=0, noblank=0, verb=1):
   """return the text text from the given file as either one string
      or as an array of lines

        strip:   remove all surrounding whitespace from each line
        nonl:    remove only trailing '\n' characters (useless with strip)
        noblank: remove all completely blank lines
   """

   if fname == 'stdin' or fname == '-': fp = sys.stdin
   else:
      try: fp = open(fname, 'r')
      except:
        if verb: print("** read_text_file: failed to open '%s'" % fname)
        if lines: return []
        else:     return ''

   if lines:
      tdata = fp.readlines()
      if strip: tdata = [td.strip() for td in tdata]
      if nonl:  tdata = [td.replace('\n','') for td in tdata]
      if noblank: tdata = [td for td in tdata if td != '']
   else:
      tdata = fp.read()
      if strip: tdata.strip()

   fp.close()

   return tdata

def read_tsv_file(fname='stdin', strip=0, verb=1):
   """read a TSV file, and return a 2-D matrix of strings in row-major order
      (it must be rectangular)

      allow tab separators, else comma, else space
      if only 1 column, hard to guess
   """

   # get lines of text, omitting blank ones, and not including newlines
   tdata = read_text_file(fname, strip=strip, nonl=1, noblank=1, verb=verb)
   nlines = len(tdata)
   if verb > 1: print("-- TSV '%s' has %d lines" % (fname, nlines))
   if nlines == 0:
      return []

   # test for separators, require rectangular input
   sep = ''
   for tsep in ['\t', ',', ' ']:
     nsep = tdata[0].count(tsep)
     # is there anything?
     if nsep == 0:
        continue
     # if so, every line should have the same count
     matches = 1
     for tline in tdata:
       nn = tline.count(tsep)
       if nn != nsep:
          if verb > 2:
             print("-- sep '%s' mismatch: %d vs. %d, skipping" % (tsep,nn,nsep))
          matches = 0
          break
     if matches:
        # declare a winner
        sep = tsep

   if verb > 1:
      print("-- read_tsv_file: have sep '%s'" % sep)

   # if nothing found, assume one column, but each column must be a list
   if sep == '':
      table = [[tline] for tline in tdata]
   # otherwise, partition the table based on sep
   else:
      table = [tline.split(sep) for tline in tdata]

   # and make sure it is rectangular
   if len(table) == 0:
      return table

   ncols = len(table[0])
   for rind, row in enumerate(table):
      if len(row) != ncols:
         print("** table %s is not rectangular at line %d" % (fname, rind))
         return []

   if verb > 2:
      print("-- have %d x %d TSV data from %s" % (len(table), ncols, fname))

   return table

def tsv_get_vals_where_condition(fname, lab_val, where_val, verb=1):
   """from a TSV file, return value list where a condition is met

        fname     : (str) TSV file name to process, must have header row
        lab_val   : (str) column label to return values from
        where_val : (str) condition for when to include a value in return list
        verb      : (int) how much to chit-chat

      This function was written for tedana processing, where the input might
      be desc-tedana_metrics.tsv.  The metric file holds many details about
      each ICA component, where each component is on a row.  The columns to
      consider include
        'Component'      : the list of which we want to return
        'classification' : the decision for whether to return

      In this example,
        fname       = 'desc-tedana_metrics.tsv'
        lab_val     = 'Component'
        where_val   = 'classification=accepted' (or rejected)

      And the return status and value might be something like
        0, ['ICA_08', 'ICA_11', 'ICA_49']

      That is to say we might return a list of values from the 'Component'
      column where the 'classification' column value is 'accepted'.

      return status, value list
   """

   if verb > 1:
      print("-- using %s to report %s when %s" % (fname, lab_val, where_val))

   # parse inputs: where_val must currently be of the form A=B
   if '=' not in where_val:
      print("** TSV_GVWC: mal-formed where string '%s'" % where_val)
      return 1, []

   # might allow multiple where entries later, start with [0]
   where = where_val.split()[0].split('=')
   if len(where) != 2:
      print("** TSV_GVWC: bad where string '%s'" % where_val)
      return 1, []

   # read the file with the lab_val and 'where' column headers
   imat = read_tsv_file(fname, verb=verb)
   if len(imat) == 0: return 1, []  # error

   # the first row must be a header (with the 2 entries)
   ihead = imat.pop(0)
   if len(imat) == 0: return 0, []  # empty matrix, no worries

   # we must have columns lab_val and where[0] now
   if (lab_val not in ihead) or (where[0] not in ihead):
      print("** TSV_GVWC: missing header entries '%s', '%s' in %s" \
            % (lab_val, where[0], fname))
      return 1, []

   # okay, we should be ready to roll
   lab_ind = ihead.index(lab_val)
   wh_ind  = ihead.index(where[0])
   matlen  = len(imat)

   outvals = [irow[lab_ind] for irow in imat if irow[wh_ind] == where[1]]

   if verb > 1:
      print("++ TSV %s : '%s' when '%s'" % (fname, lab_val, where_val))
      print("     %d entries : %s\n" % (len(outvals), ','.join(outvals)))

   return 0, outvals


def read_top_lines(fname='stdin', nlines=1, strip=0, verb=1):
   """use read_text_file, but return only the first 'nlines' lines"""
   tdata = read_text_file(fname, strip=strip, verb=verb)
   if nlines != 0: tdata = tdata[0:nlines]
   return tdata

def read_text_dictionary(fname, verb=1, mjdiv=None, mndiv=None, compact=0,
                         qstrip=0):
   """this is the same as read_text_dict_list(), but it returns a dictionary

      if compact, collapse single entry lists
      if qstrip, strip any containing quotes
   """
   rv, ttable = read_text_dict_list(fname, verb=verb, mjdiv=mjdiv, mndiv=mndiv,
                                    compact=compact, qstrip=qstrip)
   if rv: return rv, {}
   
   rdict = {}
   for row in ttable:
      if len(row) != 2:
         print("** RT_dict: table '%s' has bad row length %d" \
               % (fname, len(row)))
         return 1, {}
      rdict[row[0]] = row[1]

   return 0, rdict

def read_text_dict_list(fname, verb=1, mjdiv=None, mndiv=None, compact=0,
                        qstrip=0):
   """read file as if in a plain dictionary format (e.g. LABEL : VAL VAL ...)

         mjdiv : major divider can be a single string (':') or a list of them
                 to try, in order
                 (if None, partitioning will be attempted over a few cases)
         mkdiv : minor divider should be a single string ('SPACE', ',,')

            if either divider is 'SPACE', the natural .split() will be used

         compact: collapse any single entry lists to return
                  DLIST enties of form [LABEL, VAL] or [LABEL, [V0, V1...]]

         qstrip:  strip any surrounding quotes

      return status, DLIST
             where status == 0 on success
                   DLIST is a list of [LABEL, [VAL,VAL...]]

      The list is to keep the original order.
   """
   tlines = read_text_file(fname, lines=1, strip=1, noblank=1, verb=verb)
   nrows = len(tlines)
   if nrows == 0: return 0, []

   # note the major divider
   if mjdiv == None:
      major_list = [':', '=', 'SPACE', '::']
   elif type(mjdiv) == list:
      major_list = mjdiv
   else:
      major_list = [mjdiv]

   # see if we have a major partitioning
   for mdiv in major_list:
      ttable = [] # keep the sorting
      good = 1
      for line in tlines:
         if mdiv == 'SPACE': entries = line.split()
         else:               entries = line.split(mdiv)
         if len(entries) == 2:
            entries = [e.strip() for e in entries]
            ttable.append(entries)
         else:
            good = 0
            break
      if good:
         if verb > 1: print("text_dict_list: using major delim '%s'" % mdiv)
         break

   # did we find a major partitioning?
   if not good:
      if verb:
         print("** failed read_text_dictionary('%s')" % fname)
         return 1, []

   # try to guess the minor divider
   if mndiv == None:
      has_cc = 0
      has_c  = 0
      has_s  = 0
      for tline in ttable:
         if ',,' in tline:
            has_cc += 1
            mndiv = ',,'
            # give double comma priority, since it should be uncommon
            break
         if ',' in tline:
            has_c += 1
   if mndiv == None:
      if has_cc: mndiv = ','
      else:      mndiv = 'SPACE'

   # now actually make a dictionary
   outtable = []
   for tline in ttable:
      if qstrip:
         tline[1] = tline[1].strip("'")
         tline[1] = tline[1].strip('"')
      if mndiv == 'SPACE': entries = tline[1].split()
      else:                entries = tline[1].split(mndiv)
      if compact and len(entries) == 1:
         entries = entries[0]
      outtable.append([tline[0], entries])

   return 0, outtable

def convert_table2dict(dlist):
   """try to convert a table to a dictionary
      each list entry must have exactly 2 elements

      return status, dictionary
   """
   if type(dlist) == dict:
      return 0, dlist
   if not type(dlist) == list:
      print("** convert_table2dict: unknown input of type %s" % type(dlist))
      return 1, None

   rdict = {}
   for dentry in dlist:
      if not type(dentry) == list:
         print("** convert_table2dict: unknown entry type %s" % type(dentry))
         return 1, None
      if len(dentry) != 2:
         print("** convert_table2dict: invalid entry length %d" % len(dentry))
         return 1, None

      # finally, the useful instruction
      rdict[dentry[0]] = dentry[1]

   return 0, rdict

def data_file_to_json(fin='stdin', fout='stdout', qstrip=0, sort=1, verb=1):
   """convert a data file to json format - this should be in a main program
      Ah, abids_json_tool.py does this.

            fin     : input  data file (can be - or stdin  for sys.stdin)
            fout    : output json file (can be - or stdout for sys.stdout)
            qstrip  : strip any containing quotes
            sort    : sort output dict

      Input should be of the form (such as):
        Label_0 : val
        Label_1 : v0 v1 v2
        ...

      Output is a dictionary form of this.
   """
   rv, tdict = read_text_dictionary(fname=fin, verb=verb, qstrip=qstrip)
   if rv: return
   rv = write_data_as_json(tdict, fname=fout, sort=sort)

def write_data_as_json(data, fname='stdout', indent=3, sort=1, newline=1,
                       table2dict=0):
   """dump to json file; check for stdout or stderr
      table2dict: if set, convert table Nx2 table into a dictionary
      return 0 on success
   """
   # possibly convert data to a dictionary
   if table2dict:
      rv, dtable = convert_table2dict(data)
      if rv:
         print("** write_data_as_json: failed to create file '%s'" % fname)
         return 1
      data = dtable

   # import locally, unless it is needed in at least a few functions
   # (will make default in afni_proc.py, so be safe)
   try:
      import json
   except:
      print("** afni_util.py: 'json' python library is missing")
      return 1

   if fname == 'stdout' or fname == '-':
      fp = sys.stdout
   elif fname == 'stderr':
      fp = sys.stderr
   else:
      try: fp = open(fname, 'w')
      except:
         print("** write_as_json: could not open '%s' for writing" % fname)
         return 1

   # actual write
   json.dump(data, fp, sort_keys=sort, indent=indent)

   if newline: fp.write('\n')
   if fp != sys.stdout and fp != sys.stderr:
      fp.close()

   return 0

def read_json_file(fname='stdin'):
   try:    import json
   except: return {}
   textdata = read_text_file(fname, lines=0)
   return json.loads(textdata)

def print_dict(pdict, fields=[], nindent=0, jstr=' ', verb=1):
   """print the contents of a dictionary, and possibly just a list of fields"""

   istr = ' '*nindent

   # allow for passing a simple string
   if type(fields) == str:
      fields = [fields]

   nfields = len(fields)
   for kk in pdict.keys():
      if nfields > 0 and kk not in fields:
         continue
      val = pdict[kk]
      if type(val) == dict:
         if verb: print('%s%s :' % (istr, kk))
         # recur
         print_dict(val, nindent=(nindent+3), verb=verb)
      elif type(val) == list:
         lstr = jstr.join('%s' % v for v in val)
         if verb: print('%s%s : %s' % (istr, kk, lstr))
         else:    print(lstr)
      else:
         if verb: print('%s%s : %s' % (istr, kk, val))
         else:    print(val)

def print_json_dict(fname='stdin', fields=[], verb=1):
   jdata = read_json_file(fname)
   print_dict(jdata, fields=fields, verb=verb)

def read_AFNI_version_file(vdir='', vfile='AFNI_version.txt', delim=', '):
   """read AFNI_version.txt from vdir (else executable_dir)
      return comma-delimited form
   """

   if vdir == '': vdir = executable_dir()
   if vdir == '': return ''

   vpath = '%s/%s' % (vdir, vfile)

   if not os.path.isfile(vpath): return ''

   vdata = read_text_file(vpath, verb=0)
   if vdata == '': return ''

   return delim.join(vdata)

def write_to_timing_file(data, fname='', nplaces=-1, verb=1):
   """write the data in stim_times format, over rows
      (this is not for use with married timing, but for simple times)"""

   if fname == '': return

   fp = open(fname, 'w')
   if not fp:
      print("** failed to open '%s' for writing timing" % fname)
      return 1

   if verb > 0:
      print("++ writing %d timing rows to %s" % (len(data), fname))

   fp.write(make_timing_data_string(data, nplaces=nplaces, flag_empty=1,
                                    verb=verb))
   fp.close()

   return 0

def make_timing_data_string(data, row=-1, nplaces=3, flag_empty=0,
                            mesg='', verb=1):
   """return a string of row data, to the given number of decimal places
      if row is non-negative, return a string for the given row, else
      return a string of all rows"""

   if verb > 2:
      print('++ make_data_string: row = %d, nplaces = %d, flag_empty = %d' \
            % (row, nplaces, flag_empty))

   if row >= 0:
      return make_single_row_string(data[row], row, nplaces, flag_empty)

   # make it for all rows
   if len(mesg) > 0: rstr = "%s :\n" % mesg
   else:             rstr = ''
   for ind in range(len(data)):
      rstr += make_single_row_string(data[ind], ind, nplaces, flag_empty)

   return rstr

def make_single_row_string(data, row, nplaces=3, flag_empty=0):
   """return a string of row data, to the given number of decimal places
      if row is non-negative, return a string for the given row"""

   rstr = ''

   # if flagging an empty run, use '*' characters
   if len(data) == 0 and flag_empty:
      if row == 0: rstr += '* *'
      else:        rstr += '*'

   for val in data:
      if nplaces >= 0: rstr += '%.*f ' % (nplaces, val)
      else:            rstr += '%g ' % (val)

   return rstr + '\n'

def quotize_list(inlist, opt_prefix='', skip_first=0, quote_wild=0,
                 quote_chars='', ok_chars=''):
    """given a list of text elements, return a new list where any existing
       quotes are escaped, and then if there are special characters, put the
       whole string in single quotes

       if the first character is '-', opt_prefix will be applied
       if skip_first, do not add initial prefix
       if quote_wild, quotize any string with '*' or '?', too

       add quote_chars to quote list, remove ok_chars
    """
    if not inlist or len(inlist) < 1: return inlist

    # okay, we haven't yet escaped any existing quotes...

    # default to ignoring wildcards, can always double-nest if needed
    if quote_wild: qlist = "[({*? "
    else:          qlist = "[({ "

    for c in quote_chars:
        if not c in qlist: qlist += c
    for c in ok_chars:
        posn = qlist.find(c)
        if posn >= 0: qlist = qlist[0:posn]+qlist[posn+1:]

    newlist = []
    first = 1   # ugly, but easier for option processing
    for qstr in inlist:
        prefix = ''
        if skip_first and first: first = 0       # use current (empty) prefix
        elif len(qstr) == 0: pass
        elif qstr[0] == '-': prefix = opt_prefix

        quotize = 0
        for q in qlist:
            if q in qstr:
                quotize = 1
                break
        if quotize: newlist.append("'%s%s'" % (prefix,qstr))
        else:       newlist.append("%s%s" % (prefix,qstr))

    return newlist

def args_as_command(args, prefix='', suffix=''):
    """given an argument list (such as argv), create a command string,
       including any prefix and/or suffix strings"""

    if len(args) < 1: return

    cstr = "%s %s" % (os.path.basename(args[0]),
                            ' '.join(quotize_list(args[1:],'')))
    fstr = add_line_wrappers('%s%s%s' % (prefix,cstr,suffix))

    return fstr

def show_args_as_command(args, note='command:'):
     """print the given argument list as a command
        (this allows users to see wildcard expansions, for example)"""

     print(args_as_command(args,
     "----------------------------------------------------------------------\n"
     "%s\n\n  " % note,
     "\n----------------------------------------------------------------------"
     ))

def exec_tcsh_command(cmd, lines=0, noblank=0, showproc=0):
    """execute cmd via: tcsh -cf "cmd"
       return status, output
          if status == 0, output is stdout
          else            output is stderr+stdout

       if showproc, show command and text output, and do not return any text

       if lines: return a list of lines
       if noblank (and lines): omit blank lines
    """

    # if showproc, show all output immediately
    if showproc: capture = 0
    else:        capture = 1

    # do not re-process .cshrc, as some actually output text
    cstr = 'tcsh -cf "%s"' % cmd
    if showproc:
       print("%s" % cmd)
       sys.stdout.flush()
    status, so, se = BASE.simple_shell_exec(cstr, capture=capture)

    if not status: otext = so
    else:          otext = se+so
    if lines:
       slines = otext.splitlines()
       if noblank: slines = [line for line in slines if line != '']
       return status, slines

    return status, otext

def limited_shell_exec(command, nlines=0):
   """run a simple shell command, returning the top nlines"""
   st, so, se = BASE.shell_exec2(command, capture=1)
   if nlines > 0:
      so = so[0:nlines]
      se = se[0:nlines]
   return st, so, se

def write_afni_com_history(fname, length=0, wrap=1):
   """write the afni_com history to the given file

      if length > 0: limit to that number of entries
   """
   com = BASE.shell_com('hi there')
   hist = com.shell_history(nhist=length)
   script = '\n'.join(hist)+'\n'
   write_text_to_file(fname, script, wrap=wrap)

def write_afni_com_log(fname=None, length=0, wrap=1):
   """write the afni_com log to the given file

      if length > 0: limit to that number of entries

      if no fname is given, simply print the log
   """
   com = BASE.shell_com('hi there')
   log = com.shell_log(nlog=length)
   
   # wrapping will occur *here*, if used
   log2   = proc_log(log, wrap=wrap)
   script = '\n'.join(log2)+'\n'

   if fname is None :
      print(script)
   else:
      # wrapping will have already occurred, above
      write_text_to_file(fname, script, wrap=0)

def proc_log(log, wid=78, wrap=1):
    """Process the log, which is a list of dictionaries (each of cmd,
status, so and se objects), and prepare it for string output.  The
output is a list of strings to be concatenated.

    """

    N = len(log)
    if not(N) :    return ''

    log2 = []
    for ii in range(N):
        D = log[ii]
        topline = not(ii)
        log2.extend( format_log_dict(D, wid=wid, wrap=wrap,
                                     topline=topline) )
    return log2

def format_log_dict(D, wid=78, wrap=1, topline=True):
    """Each dictionary contains the following keys: cmd, status, so, se.
Turn these into a list of strings, to be joined when displaying the log.

    """

    L = []
    if topline :
        L.append("="*wid)

    # cmd
    cmd = D['cmd'].strip()
    if cmd[0] == '#' :
       # comment dumping, essentially
       if wrap :
          cmd = add_line_wrappers(cmd, wrapstr='\\\n#')
    elif not(len(cmd.split()) > 3 and len(cmd) > 40) :
       # short commands
       if wrap :
          cmd = add_line_wrappers(cmd)
    else:
       # long/general commands
       ok, cmd = lfcs.afni_niceify_cmd_str(cmd)
    nline = cmd.count('\n') + 1
    L.append('cmd: ' + str(nline))
    L.append(cmd)
    L.append("-"*wid)

    # status
    L.append('stat: ' + str(D['status']))
    L.append("-"*wid)

    # so
    ooo = some_types_to_str(D['so'])
    if ooo :
       if wrap :
          ooo = add_line_wrappers(ooo, wrapstr='\\\n')
       nline = ooo.count('\n') + 1
       L.append('so: ' + str(nline))
       L.append(ooo)
    else:
       L.append('so: 0')
    L.append("-"*wid)

    # se
    eee = some_types_to_str(D['se'])
    if eee :
       if wrap :
          eee = add_line_wrappers(eee, wrapstr='\\\n')
       nline = eee.count('\n') + 1
       L.append('se: ' + str(nline))
       L.append(eee)
    else:
       L.append('se: 0')
    L.append("="*wid)

    return L

def some_types_to_str(x):
    """return a string form of a list, str or 'other' type,
    ... now implementing Reynoldsian Recursion!"""
    if type(x) == str :     return x
    elif type(x) == list :  return '\n'.join([some_types_to_str(v) for v in x])
    else:                   return str(x)

def get_process_depth(pid=-1, prog=None, fast=1):
   """print stack of processes up to init"""

   pstack = get_process_stack(pid=pid, fast=fast)

   if prog == None: return len(pstack)

   pids = [pp[0] for pp in pstack if pp[3] == prog]
   return len(pids)

# get/show_process_stack(), get/show_login_shell()   28 Jun 2013 [rickr]
def get_process_stack(pid=-1, fast=1, verb=1):
   """the stack of processes up to init

      return an array of [pid, ppid, user, command] elements
      --> element 0 should always be init
  """
   def get_pid_index(pids, plist, pid):
      try:
         pind = pids.index(pid)
      except:
         if verb > 1:
            print('** GPS pid %s not in pid list:\n' % pid)
            print('%s' % plist)
         return -1
      return pind

   def get_ancestry_indlist(pids, ppids, plist, pid=-1):
      """return bad status if index() fails"""
      if pid >= 0: mypid = pid
      else:        mypid = os.getpid()
      pind = get_pid_index(pids, plist, mypid)
      if pind < 0: return 1, []
      indtree = [pind]
      while mypid > 1:
         mypid = ppids[pind]
         pind = get_pid_index(pids, plist, mypid)
         if pind < 0: return 1, []
         indtree.append(pind)
      return 0, indtree

   if not fast:
      stack = get_process_stack_slow(pid=pid)
      stack.reverse()
      return stack

   if verb < 2: cmd = 'ps -eo pid,ppid,user,comm'
   else:        cmd = 'ps -eo pid,ppid,user,args'
   ac = BASE.shell_com(cmd, capture=1)
   ac.run()
   if ac.status:
      print('** GPS command failure for: %s\n' % cmd)
      print('error output:\n%s' % '\n'.join(ac.se))
      return []

   plist = [ll.split() for ll in ac.so]
   names = plist[0]
   plist = plist[1:]

   try:
      pids = [int(psinfo[0]) for psinfo in plist]
      ppids = [int(psinfo[1]) for psinfo in plist]
   except:
      print('** GPS type failure in plist\n%s' % plist)
      return []

   # maybe the ps list is too big of a buffer, so have a backup plan
   rv, indlist = get_ancestry_indlist(pids, ppids, plist, pid=pid)
   # if success, set stack, else get it from backup function
   if rv == 0: stack = [plist[i] for i in indlist]
   else:       stack = get_process_stack_slow(pid=pid)
   stack.reverse()

   return stack

def get_process_stack_slow(pid=-1, verb=1):
   """use repeated calls to get stack:
        ps h -o pid,ppid,user,comm -p PID

        if pid >= 0, use as seed, else use os.getpid()
   """
   if verb > 1: base_cmd = 'ps h -o pid,ppid,user,args -p'
   else:        base_cmd = 'ps h -o pid,ppid,user,comm -p'

   def get_cmd_entries(cmd):
      ac = BASE.shell_com(cmd, capture=1)
      ac.run()
      if ac.status:
         print('** GPSS command failure for: %s\n' % cmd)
         print('error output:\n%s' % '\n'.join(ac.se))
         return 1, []
      ss = ac.so[0]
      entries = ss.split()
      if len(entries) == 0: return 1, []
      
      return 0, entries

   def get_ppids(cmd, entries):
      try:
         pid = int(entries[0])
         ppid = int(entries[1])
      except:
         print('** bad GPSS entries for cmd: %s\n  %s' % (cmd, entries))
         return 1, -1, -1
      return 0, pid, ppid

   # get pid and ppid
   if pid >= 0: mypid = pid
   else:        mypid = os.getpid()
   cmd = '%s %s' % (base_cmd, mypid)
   rv, entries = get_cmd_entries(cmd)
   if rv:
      if mypid == pid: print('** process ID %d not found' % pid)
      else:            print('** getpid() process ID %d not found' % pid)
      return []
   rv, mypid, ppid = get_ppids(cmd, entries)
   if rv: return []

   stack = [entries] # entries is valid, so init stack
   # now some parents to straight to 0 without 1  [28 Feb 2024]
   while mypid > 1 and ppid > 0:
      cmd = '%s %s' % (base_cmd, ppid)
      rv, entries = get_cmd_entries(cmd)
      if rv: return []
      rv, mypid, ppid = get_ppids(cmd, entries)
      if rv: return []
      stack.append(entries)

   return stack

def show_process_stack(pid=-1,fast=1,verb=1):
   """print stack of processes up to init"""
   pstack = get_process_stack(pid=pid,fast=fast,verb=verb)
   if len(pstack) == 0:
      print('** empty process stack')
      return
   ulist = [pp[2] for pp in pstack]
   ml = max_len_in_list(ulist)
   header = '   PID   PPID  [USER]'
   dashes = '  ----   ----  ------'
   form = '%6s %6s  [%s]'
   ilen = len(form)+4+ml

   print('%-*s : %s' % (ilen, header, 'COMMAND'))
   print('%-*s   %s' % (ilen, dashes, '-------'))
   for row in pstack:
      ss = form % (row[0], row[1], row[2])
      if len(row) > 3: rv = ' '.join(row[3:])
      else:            rv = row[3]
      print('%-*s : %s' % (ilen, ss, rv))

def get_login_shell():
   """return the apparent login shell
      from get_process_stack(), search down s[3] until a shell is found
   """

   pstack = get_process_stack()
   if len(pstack) == 0:
      print('** cannot detect shell: empty process stack')
      return

   # start from init and work down to find first valid shell
   for pline in pstack:
      shell = shell_name_strip(pline[3])
      if shell: return shell

   return 'SHELL_NOT_DETECTED'

def get_current_shell():
   """return the apparent login shell
      from get_process_stack(), search down s[3] until a shell is found
   """

   pstack = get_process_stack()
   if len(pstack) == 0:
      print('** cannot detect shell: empty process stack')
      return
   pstack.reverse()

   # start from init and work down to find first valid shell
   for pline in pstack:
      shell = shell_name_strip(pline[3])
      if shell: return shell

   return 'SHELL_NOT_DETECTED'

def shell_name_strip(name):
   """remove any leading dash or path, then return if valid shell, else ''"""
   global g_valid_shells
   if len(name) < 2: return ''

   # strip any leading '-' or path
   if name[0] == '-': name = name[1:]
   if name[0] == '/': name = name.split('/')[-1]

   if name in g_valid_shells: return name

   return ''    # not a shell

def show_login_shell(verb=0):
   """print the apparent login shell

      from get_process_stack(), search down s[3] until a shell is found
   """
   shells = ['csh','tcsh','sh','bash','zsh']

   pstack = get_process_stack()
   if len(pstack) == 0:
      print('** cannot detect shell: empty process stack')
      return

   # start from init and work down to find first valid shell
   shell = ''
   for pline in pstack:
      shell = shell_name_strip(pline[3])
      if not shell: continue
      if verb: print('apparent login shell: %s' % shell)
      else: print('%s' % shell)
      break

   if shell == '':
      if verb:
         print('** failed to determine login shell, see process stack...\n')
         show_process_stack()
         return

   # in verbose mode, see if parent shell is different from login
   if verb:
      pstack.reverse()
      for pline in pstack:
         sh = shell_name_strip(pline[3])
         if not sh: continue

         if sh != shell: print('differs from current shell: %s' % sh)
         break

def get_unique_sublist(inlist, keep_order=1):
    """return a copy of inlist, but where elements are unique

       if keep_order, the order is not altered (first one found is kept)
          (easy to code, but maybe slow)
       else, sort (if needed), and do a linear pass

       tested with:
          llist = [3, 4, 7, 4, 5,5,5, 4, 7, 9]
          get_unique_sublist()
          get_unique_sublist(keep_order=0)
          llist.sort()
          get_unique_sublist(keep_order=0)
    """

    if len(inlist) == 0: return []

    # if keep_order, be slow
    if keep_order:
       newlist = []
       for val in inlist:
           if not val in newlist: newlist.append(val)
       return newlist

    # else, sort only if needed
    if vals_are_sorted(inlist):
       slist = inlist
    else:
       slist = inlist[:]
       slist.sort()

    newlist = slist[0:1]
    for ind in range(len(slist)-1):
       # look for new values
       if slist[ind+1] != slist[ind]: newlist.append(slist[ind+1])

    return newlist

def uniq_list_as_dsets(dsets, byprefix=0, whine=0):
    """given a list of text elements, create a list of afni_name elements,
       and check for unique prefixes"""

    if not dsets or len(dsets) < 2: return 1

    if type(dsets[0]) == str:
       anlist = [BASE.afni_name(dset) for dset in dsets]
    elif isinstance(dsets[0], BASE.afni_name):
       anlist = dsets
    else:
       print('** ULAD: invalid type for dset list, have value %s' % dsets[0])
       return 0

    if byprefix:
       plist = [an.prefix for an in anlist]
    else:
       plist = dsets[:]
    plist.sort()

    # iterate over dsets, searching for matches
    uniq = 1
    for ind in range(len(plist)-1):
       if anlist[ind].prefix == anlist[ind+1].prefix:
          uniq = 0
          break

    if not uniq and whine:
      print("-----------------------------------------------------------\n" \
          "** dataset names are not unique\n\n"                             \
          "   (#%d == #%d, '%s' == '%s')\n\n"                               \
          "   note: if using a wildcard, please specify a suffix,\n"        \
          "         otherwise datasets may be listed twice\n"               \
          "            e.g.  bad use:    ED_r*+orig*\n"                     \
          "            e.g.  good use:   ED_r*+orig.HEAD\n"                 \
          "-----------------------------------------------------------\n"   \
          % (ind, ind+1, anlist[ind].pve(), anlist[ind+1].pve()))

    return uniq

def uniq_list_as_dset_names(dsets, whine=0):
    """like as_dsets, but the input is a simlpe array of names"""
    asets = list_to_datasets(dsets, whine=whine)
    if not asets: return 0
    return uniq_list_as_dsets(asets, whine=whine)

def list_to_datasets(words, whine=0):
    """given a list, return the list of afni_name elements
         - the list can include wildcarding
         - they must be valid names of existing datasets
         - return None on error"""

    if not words or len(words) < 1: return []
    dsets = []
    wlist = []
    errs  = 0
    for word in words:
        glist = glob.glob(word)  # first, check for globbing
        if glist:
            glist.sort()
            wlist += glist
        else: wlist.append(word)
    # now process all words
    for word in wlist:
        dset = BASE.afni_name(word)
        if dset.exist():
            dsets.append(dset)
        else:
            if whine: print("** no dataset match for '%s'" % word)
            errs = 1

    if errs:
        if whine: print('') # for separation
        return None
    return dsets

def dset_prefix_endswith(dname, suffix):
    """return 1 if an afni_name dset based on dname has a prefix
       that ends with suffix"""
    aname = BASE.afni_name(dname)
    rv = aname.prefix.endswith(suffix)
    del(aname)
    return rv

def basis_has_known_response(basis, warn=0):
    """given a string, if the prefix is either GAM or BLOCK, then the basis
       function has a known response curve

       if warn, warn users about any basis function peculiarities"""
    if not basis: return 0

    if starts_with_any_str(basis, basis_known_resp_l): return 1
    return 0

def basis_is_married(basis):
    """if the given basis function is known to require married times, return 1
    """
    if not basis: return 0

    if starts_with(basis, 'dmBLOCK'): return 1
    else:                             return 0

def basis_has_one_reg(basis, st='times'):
    """if the given basis function is known to have 1 regressor, return 1
    """
    if not basis: return 0

    # only 'times', 'file' and 'AM1' are acceptable
    if not st in stim_types_one_reg: return 0

    if starts_with_any_str(basis, basis_one_regr_l): return 1

    return 0

def starts_with(word, sstr):
    """return 1 if word starts with sstr"""
    slen = len(sstr)
    if word[0:slen] == sstr: return 1
    return 0

def starts_with_any_str(word, slist):
    """return 1 if word starts with anything in slist"""
    for sstr in slist:
       slen = len(sstr)
       if word[0:slen] == sstr: return 1
    return 0

def get_default_polort(tr, reps):
    """compute a default run polort, as done in 3dDeconvolve
       (leave run_time as TR*reps, rather than TR*(reps-1), to match 3dD)
    """

    if tr <= 0 or reps <= 0:
        print("** cannot guess polort from tr = %f, reps = %d" % (tr,reps))
        return 2        # return some default

    return run_time_to_polort(tr*reps)

def run_time_to_polort(run_time):
    """direct computation: 1+floor(run_time/150)
       return as int
    """
    return int(1+math.floor(run_time/150.0))

def index_to_run_tr(index, rlens, rstyle=1, whine=1):
    """given index and a list of run lengths,
       return the corresponding run and TR index

       rstyle: 0/1 whether the run index is 0-based or 1-based

       any negative return indicates an error
    """
    if index < 0:
       if whine: print('** ind2run_tr: illegal negative index: %d' % index)
       return 1, index
    if len(rlens) == 0:
       if whine: print('** ind2run_tr: missing run lengths')
       return 1, index

    # if there is only 1 run and it is short, compute modulo
    rlengths = rlens
    if len(rlens) == 1 and index >= rlens[0]:
       rind = index // rlens[0]
       cind = index % rlens[0]
       if rstyle: return rind+1, cind
       else:      return rind, cind

    cind = index
    for rind, nt in enumerate(rlengths):
       if nt < 0:
          if whine:
             print('** ind2run_tr: have negative run length in %s' % rlengths)
          return -1, -1
       if cind < nt:
          if rstyle: return rind+1, cind
          else:      return rind, cind
       cind -= nt

    if whine: print('** ind2run_tr, index %d outside run list %s' \
                    % (index, rlengths))
    return rind, cind+nt


def get_num_warp_pieces(dset, verb=1):
    """return the number of pieces in the the WARP_DATA transformation
       for this dataset

       Note that 12 (30 float) pieces would imply manual tlrc, while 1 piece
       would suggest auto, +acpc, or perhaps an external transform.

       dset  : (string) name of afni dataset

       if len is a multiple of 30, return len(WARP_DATA)//30
       else return 0"""

    err, result = get_typed_dset_attr_list(dset, 'WARP_DATA', float, verb=verb)

    if err: return 0            # errors printed in called function

    nvals = len(result)
    npieces = nvals//30
    if npieces * 30 != nvals:
        print('** GNWP: invalid WARP_DATA length %d' % nvals)
        return 0

    if verb > 1: print('-- dset %s has a %d-piece warp' % (dset, npieces))

    del(result)

    return npieces

def get_typed_dset_attr_list(dset, attr, atype, verb=1):
    """given an AFNI dataset, return err (0=success), [typed attr list]

       dset  : (string) name of afni dataset
       attr  : (string) attribute name
       atype : (string) return type of attribute list
       verb  : verbose level

       This ends up running 3dAttribute for the given dataset name.
    """

    alist = BASE.read_attribute(dset, attr, verb=verb)
    if alist == None and verb > 0:
        print("** GTDAL: failed to read dset attr %s, dset = %s" % (attr,dset))
        return 1, []

    err = 0
    try: result = [atype(val) for val in alist]
    except:
        if verb > 0:
           print("** GTDAL: failed to convert attr %s to type %s for dset %s"%\
                 (attr, atype, dset))
        err = 1
        result = []

    return err, result

def get_dset_history(dname, lines=1):
   """run 3dinfo -history on 'dname' and return the list of lines"""
   command = '3dinfo -history %s' % dname
   status, olines = exec_tcsh_command(command, lines=lines, noblank=1)
   if status: return []
   return olines

def get_last_history_command(dname, substr):
   """run 3dinfo -history on 'dname' and return the last line
          containing 'substr'
      return one text line, or '' on failure
   """
   olines = get_dset_history(dname)
   olen = len(olines)
   if olen == 0: return ''

   # work backwards
   for ind in range(olen-1, -1, -1):
      if olines[ind].find(substr) >= 0: return olines[ind]

   return ''

def get_last_history_ver_pack(dname):
   """run 3dinfo -history on 'dname' and return the last line
          containing a valid version
      return status and the version and package strings
             status = 0 on success, 1 on error
   """
   olines = get_dset_history(dname)
   olen = len(olines)
   if olen == 0: return ''

   # work backwards and extract the first valid version found
   for ind in range(olen-1, -1, -1):
      st, vstrs = find_afni_history_version(olines[ind])
      if st == 0:
         return 0, vstrs[0], vstrs[1]

   return 1, '', ''

def get_last_history_version(dname):
   """get_last_history_version(DSET_NAME):
      return the last AFNI version from the HISTORY
      return an empty string on error
   """
   st, aver, package = get_last_history_ver_pack(dname)
   return aver

def get_last_history_package(dname):
   """get_last_history_package(DSET_NAME):
      return the last AFNI package from the HISTORY
      return an empty string on error
   """
   st, aver, package = get_last_history_ver_pack(dname)
   return package

def find_afni_history_version(av_str):
   """given {AFNI_A.B.C:PACKAGE},
      return the status and [AFNI_A.B.C, PACKAGE] pair as a list
      return 1, [] on error
   """
   re_format = r'{(AFNI_\d+\.\d+\.\d+):(.*)}'

   try:    match = re.search(re_format, av_str)
   except: return 1, []

   if match is None:
      return 1, []

   # we have a match, just verify that we found all 4 pieces
   rv = list(match.groups())

   if len(rv) != 2: return 1, []

   return 0, rv

def parse_afni_version(av_str):
   """given 'AFNI_18.2.10', return status 0 and the int list [18,2,10]
      return 1, [] on error
   """
   re_format = r'AFNI_(\d+)\.(\d+)\.(\d+)'

   try:    match = re.search(re_format, av_str)
   except: return 1, []

   if match is None:
      return 1, []

   # we have a match, convert to ints
   rv = list(match.groups())

   if len(rv) != 3: return 1, []

   # return what we have
   try: return 0, [int(val) for val in rv]
   except: return 1, []

def get_3dinfo(dname, lines=0, verb=0):
   """run 3dinfo, possibly -verb or -VERB
      if lines, splitlines
      return text, or None on failure (applies to string or lines)
   """
   vstr = ' '
   if verb == 1: vstr = ' -verb'
   elif verb > 1: vstr = ' -VERB'
   command = '3dinfo%s "%s"' % (vstr, dname)
   status, output = exec_tcsh_command(command, lines=lines, noblank=1)
   if status: return None

   return output

def get_3dinfo_nt(dname, verb=1):
   """run 3dinfo -nt

      return 0 on failure (>= 0 on success)
   """
   command = '3dinfo -nt "%s"' % dname
   status, output, se = limited_shell_exec(command, nlines=1)
   if status or len(output) == 0:
      if verb: print('** 3dinfo -nt failure: message is:\n%s%s\n' % (se,output))
      return 0

   output = output[0].strip()
   if output == 'NO-DSET' :
      if verb: print('** 3dinfo -nt: no dataset %s' % dname)
      return 0

   nt = 0
   try: nt = int(output)
   except:
      if verb: print('** 3dinfo -nt: cannot get NT from %s, for dset %s' \
                     % (output, dname))
      return 0

   return nt

def get_3dinfo_val(dname, val, vtype, verb=1):
   """run 3dinfo -val, and convert to vtype (also serves as a test)

      return vtype(0) on failure
   """
   command = '3dinfo -%s "%s"' % (val, dname)
   status, output, se = limited_shell_exec(command, nlines=1)
   if status or len(output) == 0:
      if verb:
         print('** 3dinfo -%s failure: message is:\n%s%s\n' % (val, se, output))
      return 0

   output = output[0].strip()
   if output == 'NO-DSET' :
      if verb: print('** 3dinfo -%s: no dataset %s' % (val, dname))
      return 0

   dval = 0
   try: dval = vtype(output)
   except:
      # allow conversion from float to int as a backup
      fail = 0
      if vtype == int:
         try:
            dval = float(output)
            dval = vtype(dval)
         except:
            fail = 1
      if verb and fail:
         print("** 3dinfo -%s: cannot get val from %s, for dset %s" \
               % (val, output, dname))
      if fail: return vtype(0)

   return dval

def get_3dinfo_val_list(dname, val, vtype, verb=1):
   """run 3dinfo -val, and convert to vtype (also serves as a test)

      return None on failure, else a list
   """
   command = '3dinfo -%s %s' % (val, dname)
   status, output, se = limited_shell_exec(command, nlines=1)
   if status or len(output) == 0:
      if verb:
         print('** 3dinfo -%s failure: message is:\n%s%s\n' % (val, se, output))
      return None

   output = output[0].strip()
   if output == 'NO-DSET' :
      if verb: print('** 3dinfo -%s: no dataset %s' % (val, dname))
      return None

   dlist = string_to_type_list(output, vtype)
   if dlist == None and verb:
      print("** 3dinfo -%s: cannot get val list from %s, for dset %s" \
            % (val, output, dname))

   return dlist

def dset_view(dname):
   """return the AFNI view for the given dset"""
   command = '3dinfo -av_space %s' % dname
   status, output = exec_tcsh_command(command)
   if status: return ''
   return output.replace('\n', '')

def get_3d_statpar(dname, vindex, statcode='', verb=0):
   """return a single stat param at the given sub-brick index
      if statcode, verify
      return -1 on error
   """
   ilines = get_3dinfo(dname, lines=1, verb=1)
   if ilines == None:
      print('** failed get_3dinfo(%s)' % dname)
      return -1

   N = len(ilines)

   # find 'At sub-brick #v' line
   sstr = 'At sub-brick #%d' % vindex
   posn = -1
   for ind, line in enumerate(ilines):
      posn = line.find(sstr)
      if posn >= 0: break

   if posn < 0:
      print('** 3d_statpar: no %s[%d]' % (dname, vindex))
      return -1       # failure

   # check statcode?
   lind = ind + 1
   if lind >= N:
      if verb > 1: print('** 3d_statpar: no space for statpar line')
      return -1

   sline = ilines[lind]
   plist = sline.split()
   if statcode: 
      olist = find_opt_and_params(sline, 'statcode', 2)
      if len(olist) < 3:
         print('** 3d_statpar: missing expected statcode')
         return -1
      code = olist[2]
      if code[-1] == ';': code = code[:-1]
      if code != statcode:
         print('** 3d_statpar: statcode %s does not match expected %s'\
               % (code, statcode))
         return -1
      if verb > 2: print('-- found %s' % olist)

   # now get something like "statpar = 32 x x"
   olist = find_opt_and_params(sline, 'statpar', 4)
   if len(olist) < 3:
      if verb: print('** 3d_statpar: missing expected statpar')
      if verb > 2: print('   found %s in %s' % (olist, sline))
      return -1 
   if verb > 2: print('-- found %s' % olist)

   par = -1
   try: par = int(olist[2])
   except:
      if verb: print('** 3d_statpar: bad stat par[2] in %s' % olist)
      return -1 

   return par

def converts_to_type(val, vtype):
   """  converts_to_type(val, vtype) :
        return whether vtype(val) succeeds"""
   rv = 1
   try: vtype(val)
   except: rv = 0

   return rv

def find_opt_and_params(text, opt, nopt=0):
   """given some text, return the option with that text, as well as
      the following 'nopt' parameters (truncated list if not found)"""
   tlist = text.split()

   if not opt in tlist: return []

   tind = tlist.index(opt)

   return tlist[tind:tind+1+nopt]

def get_truncated_grid_dim(dset, scale=1, verb=1):
    """return a new (isotropic) grid dimension based on the current grid
       - given md = min(DELTAS), return md truncated to 3 significant bits

            scale : scale dims by value so we do not miss something like 2.999
                  - this should likely be just greater than 1.0
                    (scale to make it independent of values)

       - return <= 0 on failure
    """
    err, dims = get_typed_dset_attr_list(dset, 'DELTA', float)
    if err: return -1
    if len(dims) < 1: return -1
    for ind in range(len(dims)):
        dims[ind] = abs(dims[ind])
    md = min(dims)
    # changed 2 -> 4  19 Mar 2010 
    if md >= 4.0: return math.floor(md)
    if md <= 0:
        print('** failed to get truncated grid dim from %s' % dims)
        return 0

    return truncate_to_N_bits(md, 3, method='r_then_t', scale=scale, verb=verb)

def truncate_to_N_bits(val, bits, method='trunc', scale=1, verb=1):
    """truncate the real value to most significant N bits
       allow for any real val and positive integer bits

       method   trunc           - truncate to 'bits' significant bits
                round           - round to 'bits' significant bits
                r_then_t        - round to 2*bits sig bits, then trunc to bits
    """

    # allow any real val
    if val == 0.0: return 0.0
    if val < 0.0: sign, fval = -1, -float(val)
    else:         sign, fval =  1,  float(val)

    if verb > 2:
       print('T2NB: applying sign=%d, fval=%g, scale=%g' % (sign,fval,scale))

    # if r_then_t, start by rounding to 2*bits, then continue to truncate
    meth = method
    if method == 'r_then_t':
        fval = truncate_to_N_bits(val, 2*bits, method='round', scale=scale,
                                  verb=verb)
        meth = 'trunc'

    if bits <= 0 or type(bits) != type(1):
        print("** truncate to N bits: bad bits = ", bits)
        return 0.0

    # possibly scale to just greater than 1
    if scale > 0:
        fval *= scale

    # find integer m s.t.  2^(bits-1) <= 2^m * fval < 2^bits
    log2 = math.log(2.0)
    m    = int(math.ceil(bits-1 - math.log(fval)/log2))
    pm   = 2**m

    # then (round or) truncate to an actual integer in that range
    # and divide by 2^m (cannot be r_then_t here)
    if meth == 'round': ival = round(pm * fval)
    else:               ival = math.floor(pm * fval)
    retval = sign*float(ival)/pm
    
    if verb > 2:
        print('-- T2NB: 2^%d <= 2^%d * %g < 2^%d' % (bits-1,m,fval,bits))
        print('         ival = %g, returning %g' % (ival,retval))

    return retval

def test_truncation(top=10.0, bot=0.1, bits=3, e=0.0001):
    """starting at top, repeatedly truncate to bits bits, and subtract e,
       while result is greater than bot"""

    print('-- truncating from %g down to %g with %d bits' % (top,bot,bits))
    val = top
    while val > bot:
        trunc = truncate_to_N_bits(val, bits, scale=1.0001)
        print(val, ' -> ', trunc)
        val = min(trunc-e, val-e)

def round_int_or_nsig(x, nsig=None, stringify=False):
    """if int(x) has at least nsig significant digits, then return
    round(x); otherwise, round and keep nsig significant digits at a
    minimum. If stringify=True, return as a string (since this
    function is often for preparing string output to report).

    if nsig=None, then just return the round(x)

    """

    # simple case: round to int
    if nsig is None :
       y = round(x)
       if stringify :  return "{:d}".format(y)
       else:           return y

    # some work: in sci notation, what is the exponent?
    nleft  = int(("{:e}".format(x)).split('e')[-1]) + 1
    nright = nsig - nleft
 
    # case: number of ints is >= num of sig figs, so round to int
    if nright <= 0 :
       y = round(x)
       if stringify :  return "{:d}".format(y)
       else:           return y

    # more work: round to appropriate number of decimals
    y = round(x, nright)

    if stringify : return "{:.{:d}f}".format(y, nright)
    else:          return y

def get_dset_reps_tr(dset, notr=0, verb=1):
    """given an AFNI dataset, return err, reps, tr

       if notr: use 3dinfo -nt

       err  = error code (0 = success, else failure)
       reps = number of TRs in the dataset
       tr   = length of TR, in seconds
    """

    # use 3dinfo directly, instead of TAXIS attributes  30 Jun 2016

    reps = get_3dinfo_val(dset, 'nt', int, verb=verb)
    tr = get_3dinfo_val(dset, 'tr', float, verb=verb)

    # store timing info in a list (to get reps and timing units)
    if reps == 0:
        print("** failed to find the number of TRs from dset '%s'" % dset)
        return 1, None, None

    if verb > 1: print('-- dset %s : reps = %d, tr = %ss' % (dset, reps, tr))

    return 0, reps, tr

def gaussian_width_to_fwhm(width, mode):
    """convert the given 'width' of gaussian 'mode' to FWHM
       (full width at half max)

       mode can be one of: fwhm, rms, sigma

            conversion based on valid mode is:
                rms   = 0.42466090 * fwhm
                sigma = 0.57735027 * rms

            implying:
                fwhm = 2.354820 * rms
                fwhm = 4.078668 * sigma

        return 0 on failure or error"""

    if width <= 0.0: return 0.0
    if mode == 'fwhm':  return width
    if mode == 'rms':   return width * 2.354820
    if mode == 'sigma': return width * 4.078668

    print("** GW2F: illegal mode '%s'" % mode)

    return 0.0

def attr_equals_val(object, attr, val):
    """return 1 of the object has attribute attr and it equals val"""

    rv = 0
    try:
       oval = getattr(object, attr)
       if oval == val: rv = 1
    except: pass

    return rv

def median(vals):
   """return the median of a list of values
      (dupe to sort without modifying input)
   """
   nvals = len(vals)

   # trivial cases
   if nvals == 0:
      return 0
   if nvals == 1:
      return vals[0]

   svals = sorted(vals)

   # truncate nvals/2, both as a test and as an index
   nby2 = int(nvals/2)

   # set based on parity of nvals
   # even: if nvals == 20, want ave(index 9 + index 10)
   if nby2 == nvals/2:
      med = (svals[nby2-1]+svals[nby2]) / 2.0
   # odd:  if nvals == 21, want index 10
   else:
      med = float(svals[nby2])

   del(svals)
   return med

def __mean_slice_diff(vals, verb=1):
   """return what seems to be the slice diff
      - get unique, sorted sublist, then diffs, then mean
   """
   unique = get_unique_sublist(vals)
   nunique = len(unique)
   # quick return - when there are no diffs
   if nunique < 2:
      return 0.0

   # sort unique sublist
   unique.sort()

   # get first diffs
   diffs = [unique[i+1]-unique[i] for i in range(nunique-1)]
   diffs.sort()

   if verb > 1:
      print("-- MSD: slice diffs: %s" % gen_float_list_string(diffs))
      print("      : min diff %g, max diff %g, diffdiff %g" \
            % (diffs[0], diffs[-1], diffs[-1]-diffs[0]))

   # return mean
   avediff = mean(diffs)
   del(diffs)

   return avediff

def numerical_resolution(vals):
   """return the apparent resolution of values expected to be on a grid
      (zero is good)

      input:  a list of real numbers
      output: the data resolution, largest minus smallest first diff

      The input vals are supposed to be multiples of some constant C, such
      that the sorted list of unique values should be:
             {0*C, 1*C, 2*C, ..., (N-1)*C}.
      In such a case, the first diffs would all be C, and the second diffs
      would be zero.  The returned resolution would be zero.

      If the first diffs are not all exactly some constant C, the largest
      difference between those diffs should implicate the numerical resolution,
      like a truncation error.  So return the largest first diff minus the
      smallest first diff.

      ** importantly: returns closer to zero are "better"

      - get unique, sorted sublist
      - take diffs
      - get unique, sorted sublist (of the diffs)
      - return last-first (so largest diff minus smallest)
   """
   unique = get_unique_sublist(vals)
   nunique = len(unique)
   # quick return - when there are no diffs
   if nunique < 2:
      return 0.0

   # sort unique sublist
   unique.sort()

   # get first diffs
   diffs = [unique[i+1]-unique[i] for i in range(nunique-1)]
   diffs.sort()

   return diffs[-1]-diffs[0]

def timing_to_slice_pattern(timing, rdigits=1, verb=1):
   """given an array of slice times, try to return multiband level and
      a pattern in g_valid_slice_patterns

      inputs:
         timing     : <float list> : slice times
         rdigits    : [1] <int>    : num digits to round to for required test
         verb       : [1] <int>    : verbosity level

      method:
         - detect timing grid (TR = tgrid*nslice_times (unique))
         - multiband = number of repeated time sets (ntimes/nunique)
         - round timing/tgrid
            - test as ints in {0..nunique-1}
            - detect timing pattern in this int list
            - check for multiband # of repeats

      return status/mb level (int), tpattern (string):
        status     -1   invalid timing
                    0   invalid multiband (not even 1)
                 >= 1   multiband level of timing (usually 1)
        tpattern        val in g_valid_slice_patterns or 'irregular'
   """

   # ----------------------------------------------------------------------
   # prep: get basic data
   #    slice diff (tgrid), TR, mblevel
   #    unique sublist (tunique), timing and unique lengths

   # default pattern and bad ones (to avoid random typos)
   defpat = 'simult'

   # first, estimate the slice diff (prefer mean over median for weak timing)
   tunique = get_unique_sublist(timing)
   tgrid = __mean_slice_diff(timing, verb=verb)

   ntimes = len(timing)
   nunique = len(tunique)

   # be sure there is something to work with
   if nunique <= 1:
      return 1, defpat

   # TR is slice time * num_slices_times
   TR = tgrid*nunique

   # note multiband level (number of repeated times)
   mblevel = int(round(ntimes/nunique))

   if verb > 2:
      print("-- TR +~ %g, MB %g, rdig %d, nunique %g, med slice diff: %g" \
            % (TR, mblevel, rdigits, nunique, tgrid))

   # if TR is not valid, we are out of here
   if TR < 0.0:
      return -1, defpat
   if TR == 0.0:
      return 1, defpat

   # ----------------------------------------------------------------------
   # scale timing: divide by tgrid to put in {0, 1, ..., nunique-1}
   scalar = 1.0/tgrid
   tscale = [t*scalar for t in timing]

   # get rounded unique sublist and sort, to compare against ints
   tround = get_unique_sublist([int(round(t,ndigits=rdigits)) for t in tscale])
   tround.sort()

   # chat
   if verb > 1:
      # print("++ t2sp: TR %s, min %s, max %s" % (TR, nzmin, nzmax))
      if verb > 2: print("-- times : %s" % timing)
      if verb > 3:
         print("== (sorted) tscale should be close to ints, tround must be ints")
         print("-- tscale: %s" % gen_float_list_string(sorted(tscale)))
         print("-- tround: %s" % gen_float_list_string(tround))

   # ----------------------------------------------------------------------
   # tests:

   # tround MUST be consecutive ints
   # (test that they round well, and are consecutive and complete)
   for ind in range(len(tround)):
      if ind != tround[ind]:
         if verb > 1:
            print("** timing is not multiples of expected %s" % tgrid)
         return 1, defpat

   # and tround must be the same length as nunique
   if len(tround) != nunique:
      if verb > 1:
         print("** have %d unique times, but %d unique scaled and rounded times" \
               % (nunique, len(tround)))
      return 1, defpat

   del(tround) # finished with this

   # does mblevel partition ntimes?
   if ntimes != mblevel * round(ntimes/mblevel,ndigits=(rdigits+1)):
      if verb > 1:
         print("** mblevel %d does not divide ntimes %d" % (mblevel, ntimes))
      return 1, defpat

   # check tscale timing, warn when not close enough to an int
   warnvec = []
   for ind in range(ntimes):
      tsval = tscale[ind]
      if round(tsval) != round(tsval, ndigits=2):
         warnvec.append(ind)

   # actually warn if we found anything
   badmults = 0
   if verb > 0 and len(warnvec) > 0:
      badmults = 1
      print("** warning: %d/%d slice times are only approx multiples of %g" \
            % (len(warnvec), ntimes, tgrid))
      if verb > 1:
         # make this pretty?
         ratios = [t/tgrid for t in timing]
         maxlen, strtimes = floats_to_strings(timing)
         maxlen, strratio = floats_to_strings(ratios)
         for ind in range(ntimes):
            print(" bad time[%2d] : %s  /  %g  =  %s" \
                  % (ind, strtimes[ind], tgrid, strratio[ind]))

   if verb > 1 or badmults:
      print("-- timing: max ind,diff = %s" % list(__max_int_diff(tscale)))

   # ----------------------------------------------------------------------
   # at this point, the sorted list has a regular (multiband?) pattern
   # so now we :
   #   - choose a pattern based on the first nunique entries
   #   - verify that the pattern repeats mblevel times

   # variables of importance: timing, tgrid, nunique, mblevel, tscale
   #   - convert timing to ints in range(nunique) (first nunique of tscale)
   #   - and then we can ignore tpattern and tgrid
   # then new vars of importance: tings, nunique, mblevel

   # round scaled times to be ints in {1..nunique-1} (but full length)
   tints = [int(round(t)) for t in tscale]
   ti0   = tints[0:nunique]

   # finally, the real step: try to detect a pattern in the first nunique
   tpat = _uniq_ints_to_tpattern(ti0)
   if verb > 1:
      if tpat == g_tpattern_irreg:
         print("** failed to detect known slice pattern from order:\n" \
               "   %s" % ti0)
      else:
         print("++ detected pattern %s in first slice set" % tpat)

   # pattern must match for each other mblevel's
   for bandind in range(1, mblevel):
      offset = bandind * nunique
      if not lists_are_same(ti0, tints[offset:offset+nunique]):
         # failure, not a repeated list
         return 0, g_tpattern_irreg

   return mblevel, tpat

def floats_to_strings(fvals):
    """return a corresponding list of floats converted to strings, such that:
          - they should align (all have same length)
          - the decimals are all in the same position (might have space pad)
          - using %g to convert

       return the new lengths (all the same) and the strign list
    """
    if len(fvals) == 0:
       return 0, []

    slist = ["%g" % v for v in fvals]
    for ind in range(len(slist)):
       if slist[ind].find('.') < 0:
          slist[ind] += '.'

    # now get max digits before and after decimal
    maxb = 0
    maxa = 0
    for ind in range(len(slist)):
       fs = slist[ind]

       # num before and after 
       numb = fs.find('.')
       numa = len(fs) - numb - 1
       if numb > maxb:
          maxb = numb
       if numa > maxa:
          maxa = numa

    # now modify slist by padding with needed spaces
    for ind in range(len(slist)):
       fs = slist[ind]

       # same num before and after, but subtracted from maxes
       numb = fs.find('.')
       newb = maxb - numb
       newa = maxa - (len(fs) - numb - 1)

       slist[ind] = ' '*newb + slist[ind] + ' '*newa

    return len(slist[0]), slist

def __max_int_diff(vals):
    """return the (index and) maximum difference from an int (max is 0.5)
       if all the same, return 0, diff[0]
    """
    nvals = len(vals)
    if nvals == 0: return 0, 0
  
    diff = abs(vals[0]-round(vals[0]))
    mdiff = diff
    mind = 0
    for ind in range(nvals):
       diff = abs(vals[ind]-round(vals[ind]))
       if diff > mdiff:
          mind = ind
          mdiff = diff

    return mind, mdiff

def _uniq_ints_to_tpattern(tints):
   """given a list of (unique) ints 0..N-1, try to match a timing pattern
        since uniq, do not test 'simult'
        test for : 'seq+z', 'seq-z', 'alt+z', 'alt-z', 'alt+z2', 'alt-z2'
        - for each test pattern:
            - compare with slice_pattern_to_timing()
        if no match, return 'irregular'

      return something in g_valid_slice_patterns or 'irregular'
   """
   nslices = len(tints)

   for tpat in ['seq+z', 'seq-z', 'alt+z', 'alt-z', 'alt+z2', 'alt-z2']:
      # get the expected list, compare and clean up
      ttimes = slice_pattern_to_timing(tpat, len(tints))
      rv = lists_are_same(tints, ttimes)
      del(ttimes)

      # did it match?
      if rv:
         return tpat

   # failure
   return g_tpattern_irreg

def slice_pattern_to_order(pattern, nslices):
   """return a list of slice order indices
      (a permuted list of 0..nslices-1)

        pattern : slice pattern (based on to3d - e.g. alt+z, simult)
        nslices : number of slices to apply the pattern to

      Given one of the g_valid_slice_patterns and a number of slices, return
      an array of slice indices, in the order they would be acquired.

      This assumes equal timing across the slices over the course of a TR.

      Note that zero/simult are not considered valid patterns here, since there
      is no sorted order in such a case (they are all equal).

      return : a list of slice indices, in order
             : None on error
   """

   if pattern not in g_valid_slice_patterns:
      print("** pattern_to_order, invalid pattern", pattern)
      return None
   if pattern in ['zero', 'simult']:
      print("** pattern_to_order, cannot make ordering from pattern", pattern)
      return None

   # sequential
   if pattern == 'seq+z' or pattern == 'seqplus':
      order = [ind for ind in range(nslices)]
   # reverse sequential
   elif pattern == 'seq-z' or pattern == 'seqminus':
      order = [(nslices-1-ind) for ind in range(nslices)]

   # alternating positive, get evens then odds
   elif pattern == 'alt+z' or pattern == 'altplus':
      order = list(range(0, nslices, 2))
      order.extend(list(range(1, nslices, 2)))
   # alternating negative, similar but top-down
   elif pattern == 'alt-z' or pattern == 'altminus':
      # start from final position (nslices-1) and work downward
      order =      list(range(nslices-1, -1, -2))
      order.extend(list(range(nslices-2, -1, -2)))

   # the z2 patterns are similar to alt, but odds come first
   elif pattern == 'alt+z2' :
      order = list(range(1, nslices, 2))
      order.extend(list(range(0, nslices, 2)))
   # alternating negative, similar but top-down
   elif pattern == 'alt-z2' :
      order =      list(range(nslices-2, -1, -2))
      order.extend(list(range(nslices-1, -1, -2)))

   else:
      print("** pattern_to_order, unhandled pattern", pattern)
      return None

   return order

def slice_pattern_to_timing(pattern, nslices, TR=0, mblevel=1, verb=1):
   """given tpattern, nslices, TR, and multiband level,
         return a list of slice times

      parameters:
         pattern    : (string) one of g_valid_slice_patterns :
                                  'zero',  'simult',
                                  'seq+z', 'seqplus',
                                  'seq-z', 'seqminus',
                                  'alt+z', 'altplus',     'alt+z2',    
                                  'alt-z', 'altminus',    'alt-z2',    
         nslices    : (int)    total number of output slice times
         TR         : (float)  total time to acquire all slices
         mblevel    : (int)    multiband level (number of repeated time sets)
         verb       : (int)    verbosity level

      special case: if TR == 0 (or unspecified)
         - do not scale (so output is int list, as if TR==nslices/mblevel)

      method:
         - verify that nslices is a multiple of mblevel
         - get result for ntimes = nslices/mblevel
            - get slice_pattern_to_order()
              - this is a list of slice indexes in the order acquired
            - attach the consecutive index list, range(nslices)
              - i.e, make list of [ [slice_index, acquisition_index] ]
            - sort() - i.e. by slice_index
              - so element [0] values will be the sorted list of slices
            - grab element [1] from each
              - this is the order the given slice was acquired in
            - scale all by TR/nslices
         - duplicate the result across other levels

      return a list of slice times, or an empty list on error
   """
   # ---------- sanity checks  ----------

   if nslices <= 0 or TR < 0.0 or mblevel <= 0:
      return []
   if nslices == 1:
      return [0]

   if pattern not in g_valid_slice_patterns:
      print("** slice_pattern_to_timing, invalid pattern", pattern)
      return []

   # if there is no time to partition or slices are simulaneous, return zeros
   if pattern in ['zero', 'simult']:
      return [0] * nslices

   # ---------- check for multiband  ----------

   ntimes = int(nslices/mblevel)
   if mblevel > 1:
      if nslices != ntimes*mblevel:
         print("** error: nslices (%d) not multiple of mblevel (%d)" \
               % (nslices, mblevel))
         return []

   # ---------- get result for ntimes ----------

   # first get the slice order
   order = slice_pattern_to_order(pattern, ntimes)
   if order is None:
      return []

   # attach index and sort
   slice_ordering = [ [order[ind], ind] for ind in range(ntimes)]
   slice_ordering.sort()

   # grab each element [1] and scale by TR/ntimes
   # (if TR == 0, do not scale)
   if TR == 0:
      stimes = [so[1]           for so in slice_ordering]
   else:
      stimes = [so[1]*TR/ntimes for so in slice_ordering]

   # ---------- duplicate results to mblevel ----------

   return stimes*mblevel

# ----------------------------------------------------------------------
# begin matrix functions

def num_cols_1D(filename):
    """return the number of columns in a 1D file"""
    mat = TD.read_1D_file(filename)
    if not mat or len(mat) == 0: return 0
    return len(mat[0])

def num_rows_1D(filename):
    """return the number of columns in a 1D file"""
    mat = TD.read_1D_file(filename)
    if not mat: return 0
    return len(mat)

def max_dim_1D(filename):
    """return the larger of the number of rows or columns"""
    mat = TD.read_1D_file(filename)
    if not mat: return 0
    rows = len(mat)
    cols = len(mat[0])
    if rows >= cols: return rows
    else:            return cols

def is_matrix_square( mat, full_check=False ):
    """return {0|1} about matrix being square"""

    if not(mat) or len(mat)==0 : return 0

    rows = len(mat)

    # can do complete check, or just [0]th ele check
    N = 1
    if full_check :
        N = rows
    
    for i in range(N):
        if len(mat[i]) != rows :
            return 0

    # have we survived to here? then -> square.
    return 1

def mat_row_mincol_maxcol_ragged_square(M):
    """check 5 properties of a matrix (list of lists) and return 5 ints:
      nrow
      min ncol
      max ncol  
      is_ragged 
      is_square

    """

    if not(M):  return 0,0,0,0,0

    is_square = 0             # just default; can change below

    nrow      = len(M)        
    all_clen  = [len(r) for r in M]
    ncolmin    = min(all_clen)
    ncolmax    = max(all_clen)

    if ncolmin == ncolmax :
        is_ragged = 0
        if ncolmin == nrow :
            is_square = 1
    else:
        is_ragged = 1

    return nrow, ncolmin, ncolmax, is_ragged, is_square
    
def transpose(matrix):
    """transpose a 2D matrix, returning the new one"""
    if matrix is None: return []

    rows = len(matrix)

    # handle trivial case
    if rows == 0: return []

    cols = len(matrix[0])
    newmat = []
    for c in range(cols):
        newrow = []
        for r in range(rows):
            newrow.append(matrix[r][c])
        newmat.append(newrow)
    return newmat

def derivative(vector, in_place=0, direct=0):
    """take the derivative of the vector, setting d(t=0) = 0

        in_place: if set, the passed vector will be modified
        direct:   if 0, use backward difference, if 1 use forward

       return v[t]-v[t-1] for 1 in {1,..,len(v)-1}, d[0]=0"""

    if type(vector) != type([]):
        print("** cannot take derivative of non-vector '%s'" % vector)
        return None

    if in_place: vec = vector    # reference original
    else:        vec = vector[:] # start with copy
    
    # count from the end to allow memory overwrite
    if direct:  # forward difference
       vlen = len(vec)
       for t in range(vlen-1):
           vec[t] = vec[t+1] - vec[t]
       vec[vlen-1] = 0
    else:       # backward difference
       for t in range(len(vec)-1, 0, -1):
           vec[t] -= vec[t-1]
       vec[0] = 0

    return vec

def matrix_multiply_2D(A, B, zero_dtype=''):
    '''Perform matrix multiplication of 2D lists, A and B, which can be of
arbitrary dimension (subject to Arow/Bcol matching).

    Each input must have 2 indices.  So, the following are valid inputs:
        x = [[1, 2, 3]]
    while these are NOT valid:
        x = [1, 2, 3] 

    Output a new array of appropriate dims, which will also always
    have 2 inds (unless one input had zero row or col; in this case,
    return empty arr).  So, output might look like:
       [[16], [10], [12]]
       [[2]] 
       []
    etc.

    The dtype of zeros used to initialize the matrix will be either:
    + whatever the [0][0]th element of the A matrix is, or
    + whatever is specified with 'zero_dtype' (which can be: bool,
      int, float, complex, ...)
    NB: any valid zero_dtype value will take precedence.

    '''

    # ------------ check dims/shapes ------------
    NrowA = len(A)
    if not(NrowA) :
        print("** ERROR: No rows in matrix A\n")
        sys.exit(4)
    NcolA = len(A[0])
    NrowB = len(B)
    if not(NrowB) :
        print("** ERROR: No rows in matrix B\n")
        sys.exit(4)
    NcolB = len(B[0])

    # Zero row/col case: return empty 1D arr
    if not(NcolA) or not(NcolB) or not(NrowA) or not(NrowB):
        return []

    if NcolA != NrowB :
        print("** ERROR: Ncol in matrix A ({}) != Nrow in matrix B ({})."
        "".format(NcolA, NrowB))
        sys.exit(5)

    # Allow for different types of matrices
    ZZ = calc_zero_dtype(A[0][0], zero_dtype)

    # Initialize output list of correct dimensions
    C = [[ZZ] * NcolB for row in range(NrowA)] 

    # Ye olde matrix multiplication, boolean style
    if type(ZZ) == bool :
        for i in range(NrowA):
            for j in range(NcolB):
                for k in range(NcolA):
                    C[i][j] = C[i][j] or (A[i][k] and B[k][j])
    else: 
        # Ye olde matrix multiplication
        for i in range(NrowA):
            for j in range(NcolB):
                for k in range(NcolA):
                    C[i][j]+= A[i][k] * B[k][j]

    return C

def matrix_sum_abs_val_ele_row(A, zero_dtype='' ):
    '''Input is a 2D list (matrix).

    Output is a list of sums of the absolute values of elements in
    rows (SAVERs!).

    '''

    N = len(A)
    if not(N) :
        return []
    Ncol = len(A[0]) # if this doesn't exist, someone will be unhappy

    # initialize for summation
    ZZ  = calc_zero_dtype(A[0][0], zero_dtype)
    out = [ZZ]*N

    # sum abs vals across rows
    if type(ZZ) == bool :
        for i in range(N):
            for j in range(len(A[i])) :
                out[i] = out[i] or abs(A[i][j])
    else:
        for i in range(N):
            for j in range(len(A[i])) :
                out[i]+= abs(A[i][j])

    return out

def calc_zero_dtype(x, zero_dtype=''):
    '''Calc an appropriate zero for initializing a list-matrix.

    '''

    # priority to an entered zero_dtype
    if    zero_dtype == int     : return 0
    elif  zero_dtype == float   : return 0.0
    elif  zero_dtype == bool    : return False
    elif  zero_dtype == complex : return 0j
    elif  type(x)    == int     : return 0
    elif  type(x)    == float   : return 0.0
    elif  type(x)    == bool    : return False
    elif  type(x)    == complex : return 0j
    else:
        print("+* Warning: unrecognized numeric type: {}"
              "".format(type(x)))
        print("   Initial matrix zeros will be type 'int'")
        return 0

def make_timing_string(data, nruns, tr, invert=0):
   """evaluating the data array as boolean (zero or non-zero), return
      non-zero entries in stim_times format

      data      : single vector (length must be multiple of nruns)
      nruns     : number of runs (must divide len(data))
      tr        : TR in seconds, with data viewed as per TR

      invert    : (optional) invert the boolean logic (write 0 times)

      return err code (0 on success), stim_times string"""

   if not data:
      print("** make_timing_string: missing data")
      return 1, ''
   if not type(data) == type([]):
      print("** make_timing_string: data is not a list")
      return 1, ''

   nvals = len(data)
   rlen  = nvals // nruns

   if nruns * rlen != nvals:
      print("** make_timing_str: nruns %d does not divide nvals %d"%(rlen,nvals))
      return 1, ''
   if tr <= 0.0:
      print("** make_timing_string: bad tr = %g" % tr)
      return 1, ''

   rstr = ''

   for run in range(nruns):
      bot = run*rlen
      if invert: rvals = [1*(data[i] == 0) for i in range(bot,bot+rlen)]
      else:      rvals = [1*(data[i] != 0) for i in range(bot,bot+rlen)]
      # if the run is empty, print 1 or 2 '*'
      nzero = rvals.count(0)
      if nzero == rlen:
         if run == 0: rstr += '* *'
         else:        rstr += '*'
      else:
         rstr += ' '.join(['%g'%(i*tr) for i in range(rlen) if rvals[i]])

      # if run0 and exactly 1 non-zero value, print a trailing '*'
      if run == 0 and nzero == rlen-1: rstr += ' *'
      rstr += '\n'

   return 0, rstr

def make_CENSORTR_string(data, nruns=0, rlens=[], invert=0, asopt=0, verb=1):
   """evaluating the data array as boolean (zero or non-zero), return
      non-zero entries in CENSORTR format

      data      : single vector (length must be multiple of nruns)
      nruns     : number of runs (must divide len(data))
                  (ignored if rlens is provided)
      rlens     : run lengths (required if run lengths vary)
      asopt     : if set, return as complete -CENSORTR option

      invert    : (optional) invert the boolean logic (write 0 TRs)

      return err code (0 on success), CENSORTR string"""

   if not data:
      print("** CENSORTR_str: missing data")
      return 1, ''
   if not type(data) == type([]):
      print("** CENSORTR_str: data is not a list")
      return 1, ''

   nvals = len(data)

   # we must have either nruns or a valid rlens list
   if nruns <= 0 and len(rlens) < 1:
      print('** make_CENSORTR_string: neither nruns nor rlens')
      return 1, ''

   if rlens:
      rlist = rlens
      runs  = len(rlist)
   else:
      rlist = [(nvals//nruns) for run in range(nruns)]
      runs  = nruns

   if verb > 1:
      print('-- CENSORTR: applying run lengths (%d) : %s' % (runs, rlist))

   if loc_sum(rlist) != nvals:
      print("** CENSORTR_str: sum of run lengths %d != nvals %d" \
            % (loc_sum(rlist),nvals))
      return 1, ''

   rstr = ''

   bot = 0
   for run in range(runs):
      rlen = rlist[run]
      if invert: rvals = [1*(data[i] == 0) for i in range(bot,bot+rlen)]
      else:      rvals = [1*(data[i] != 0) for i in range(bot,bot+rlen)]
      bot += rlen  # adjust bottom index for next run

      # if the run is empty, print 1 or 2 '*'
      nzero = rvals.count(0)
      if nzero == rlen: continue

      # make a ',' and '..' string listing TR indices
      estr = encode_1D_ints([i for i in range(rlen) if rvals[i]])

      # every ',' separated piece needs to be preceded by RUN:
      rstr += "%d:%s " % (run+1, estr.replace(',', ',%d:'%(run+1)))

   if asopt and rstr != '': rstr = "-CENSORTR %s" % rstr

   return 0, rstr

# ------------------
# funcs for specific kind of matrix: list of 2D list matrices

def check_list_2dmat_and_mask(L, mask=None):
    """A check to see if L is a list of 2D matrices.

    Also, check if any input mask matches the dimensions of the
    submatrices of L.

    This does not check if all submatrices have the same dimensions.
    Maybe someday it will.

    Returns 4 integers:
    + 'success' value   :1 if L is 2D mat and any mask matches
    + len(L)
    + nrow in L submatrix
    + ncol in L submatrix

    """

    # calc recursive list of embedded list lengths
    Ldims = get_list_mat_dims(L)

    # need a [N, nrow, ncol] here
    if len(Ldims) != 3 :   
        BASE.EP("Matrix fails test for being a list of 2D matrices;\n"
                "instead of having 3 dims, it has {}".format(len(Ldims)))

    if mask != None :
        mdims = get_list_mat_dims(mask)
        if mdims[0] != Ldims[1] or mdims[1] != Ldims[2] :
            BASE.EP("matrix dimensions don't match:\n"
                    "mask dims: {}, {}"
                    "each matrix in list has dims: {}, {}"
                    "".format(mdims[0], mdims[1], Ldims[1], Ldims[2]))

    return 1, Ldims[0], Ldims[1], Ldims[2]

def calc_list_2dmat_count_nonzero(L, mask=None, mode='count'):
    """L is a list of 2D list matrices.  Each submatrix must be the same
    shape.

    Calculate the number of nonzero entries across len(L) for each
    submatrix element.

    Return a matrix of a different style, depending on the specified
    'mode':
    'count'   :integer number of nonzero elements, i.e., max value is
               len(L); this is the default
    'frac'    :output is floating point: number of nonzero values at a 
               given matrix element, divided by len(L); therefore,  
               values here are in a range [0, 1]
    'all_nz'  :output is a matrix of binary vals, where 1 means that 
               all values across the length of L are nonzero for that 
               element, and 0 means that one or more value is zero
    'any_nz'  :output is a matrix of binary vals, where 1 means that 
               at least one value across the length of L is nonzero for
               that element, and 0 means that all values are zero

    A mask can be entered, so that all the above calcs are only done
    in the mask; outside of it, all values are 0 or False.

    if len(L) == 0:
        1 empty matrix is returned
    else:
        normal calcs

    """

    # check if L is a list of 2D mats, and if the mask matches
    # submatrices (if used)
    IS_OK, N, nrow, ncol = check_list_2dmat_and_mask(L, mask)

    mmat = [[0]*ncol for i in range(nrow)]  # [PT: June 9, 2020] fixed
    Nfl  = float(N)

    if mask == None :
        mask = [[1]*ncol for i in range(nrow)]  # [PT: June 9, 2020] fixed

    for ii in range(nrow) :
        for jj in range(ncol) :
            if mask[ii][jj] :
                for ll in range(N) :
                    if L[ll][ii][jj] :
                        mmat[ii][jj]+= 1

                if mode == 'frac':
                    mmat[ii][jj]/= Nfl
                elif mode == 'any_nz' :
                    mmat[ii][jj] = int(bool(mmat[ii][jj]))
                elif mode == 'all_nz' :
                    if mmat[ii][jj] == N :
                        mmat[ii][jj] = 1
                    else:
                        mmat[ii][jj] = 0

    return mmat

def calc_list_2dmat_mean_stdev_max_min(L, mask=None, ddof=1 ):
    """L is a list of 2D list matrices.  Each submatrix must be the same
    shape.

    Calculate four elementwise properties across the length of L, and
    return four separate matrices.  

    NB: by "list of 2D matrices," we mean that having L[0] = [[1,2,3]]
    would lead to actual calcs, but having L[0] = [1,2,3] would not.
    That is, each element of L must have both a nonzero number of rows
    and cols.

    For any shape of list matrix that isn't a "list of 2D matrices",
    we return all empty matrices.

    mask    :An optional mask can be be input; it will be treated as
             boolean.  Values outside the mask will simply be 0.
 
    ddof    :control the 'delta degrees of freedom' in the divisor of
             stdev calc; numpy's np.std() has ddof=1 by default, so the 
             denominator is "N"; here, default is ddof=1, so that 
             the denominator is N-1

    if len(L) == 0:
        4 empty matrices are returned
    elif len(L) == 1: 
        mean, min and max matrices are copies of input, and stdev is all 0s
    else:
        normal calcs

    """
    
    # check if L is a list of 2D mats, and if the mask matches
    # submatrices (if used)
    IS_OK, N, nrow, ncol = check_list_2dmat_and_mask(L, mask)

    if N == 1 :
        mmean  = copy.deepcopy(L[0])
        mstdev = [[0.]*ncol for i in range(nrow)]  # [PT: June 9, 2020] fixed
        mmax   = copy.deepcopy(L[0])
        mmin   = copy.deepcopy(L[0])
        if mask != None :
            for ii in range(nrow) :
                for jj in range(ncol) :
                    if not(mask[ii][jj]) :
                        mmean[ii][jj] = 0
                        mmin[ii][jj]  = 0
                        mmax[ii][jj]  = 0
        return mmean, mstdev, mmax, mmin

    Nfl   = float(N)
    Nstdev = Nfl - ddof  # ... and we know Nfl>1.0, if here

    # Initialize mats
    mmean  = [[0]*ncol for i in range(nrow)]  # [PT: June 9, 2020] fixed
    mstdev = [[0]*ncol for i in range(nrow)]  # [PT: June 9, 2020] fixed
    mmax   = copy.deepcopy(L[0]) 
    mmin   = copy.deepcopy(L[0])

    if mask == None :
        mask = [[1]*ncol for i in range(nrow)]  # [PT: June 9, 2020] fixed

    for ii in range(nrow) :
        for jj in range(ncol) :
            if mask[ii][jj] :
                for ll in range(N) :
                    x = L[ll][ii][jj]
                    mmean[ii][jj]+= x
                    mstdev[ii][jj]+= x*x
                    if x < mmin[ii][jj] :
                        mmin[ii][jj] = x
                    if x > mmax[ii][jj] :
                        mmax[ii][jj] = x
                print(mmean[ii][jj],Nfl)
                mmean[ii][jj]/= Nfl
                mstdev[ii][jj]-= N*mmean[ii][jj]*mmean[ii][jj]
                mstdev[ii][jj] = (mstdev[ii][jj]/Nstdev)**0.5

    return mmean, mstdev, mmax, mmin
    
def get_list_mat_dims(L, lengths=[]) :
    """Calc dims of list-style matrix.  Recursive.

    Returns a list of lengths of embedded lists; does not check that
    each list is the same length at a given level.

    """

    if type(L) == list:
        a = copy.deepcopy(lengths)
        a.append(len(L))
        return get_list_mat_dims(L[0], lengths=a)
    else:
        return lengths


# end matrix functions
# ----------------------------------------------------------------------
# index encode/decode functions

def encode_1D_ints(ilist):
   """convert a list of integers to a ',' and '..' separated string"""
   if not ilist: return ''
   if len(ilist) < 1: return ''

   text = '%d' % ilist[0]
   prev = ilist[0]
   ind  = 1
   while ind < len(ilist):
      ncontinue = consec_len(ilist, ind-1) - 1
      if ncontinue <= 1:     # then no '..' continuation, use ','
         text = text + ',%d' % ilist[ind]
         ind += 1
      else:
         text = text + '..%d' % ilist[ind+ncontinue-1]
         ind += ncontinue

   return text

def consec_len(ilist, start):
   """return the length of consecutive integers - always at least 1"""
   prev = ilist[start]
   length = len(ilist)
   for ind in range(start+1, length+1):
      if ind == length: break
      if ilist[ind] != prev + 1:
         break
      prev = ilist[ind]
   if ind == start:  length = 1
   else:             length = ind-start

   return length

def restrict_by_index_lists(dlist, ilist, base=0, nonempty=1, verb=1):
    """restrict elements of dlist by indices in ilist

        ilist    : can be string or list of strings
                  (require unique composite list)
        base     : can be 0 or 1 (0-based or 1-based)
        nonempty : if set, sub-lists are not allowed to be empty
        verb     : verbose level, default is to only report errors

       return status, sub-list
              status = 0 on success, 1 on error
    """

    # if either object is empty, there is nothing to do
    if not ilist or not dlist: return 0, []

    if type(ilist) == str: ilist = [ilist]

    if base not in [0,1]:
        if verb: print('** restrict_by_index_list: bad base = %d' % base)
        return 1, []

    # set imax to correctly imply '$' index
    if base: imax = len(dlist)          # 1-based
    else:    imax = len(dlist)-1        # 0-based

    composite = []
    for ind, istr in enumerate(ilist):
        if type(istr) != str:
            print('** RBIL: bad index selector %s' % istr)
            return 1, []
        curlist = decode_1D_ints(istr, imax=imax, verb=verb)
        if not curlist and nonempty:
            if verb: print("** empty index list for istr[%d]='%s'" % (ind,istr))
            return 1, []
        composite.extend(curlist)
        if verb > 3: print('-- index %d, ilist %s' % (ind, curlist))

    if not vals_are_unique(composite):
        if verb: print('** RBIL: composite index list elements are not unique')
        return 1, []

    cmin = min(composite)
    cmax = max(composite)
    if cmin < 0:
        if verb: print('** RBIL: cannot choose negative indices')
        return 1, []
    elif base and cmin == 0:
        if verb: print('** RBIL: 1-based index list seems 0-based')
        return 1, []
    elif cmax > imax:
        if verb: print('** RBIL: index value %d exceeds %d-based limit %d' \
                       % (cmax, base, imax))
        return 1, []

    # now convert 1-based to 0-based, if needed
    if base: clist = [v-1 for v in composite]
    else:    clist = composite

    # the big finish
    return 0, [dlist[ind] for ind in clist]

def decode_1D_ints(istr, imax=-1, labels=[], verb=1):
    """Decode a comma-delimited string of ints, ranges and A@B syntax,
       and AFNI-style sub-brick selectors (including A..B(C)).
       If the A..B format is used, and B=='$', then B gets 'imax'.
       If the list is enclosed in [], <> or ##, strip those characters.

          istr      : the int string to search through
          imax      : the max int (like final sub-brick index)
          labels    : array of labels, to possibly convert istr values to ints
          verb      : how chatty to be

       First split istr by ','.  Each returned element, can have the form of:
          - A           : a single integer (or label)
          - A@B         : (int) A entries of (int or label) B
          - A..B        : the inclusive values from (int/label) A to (I/L) B
          - A..B(C)     : similar, but with step of (int) C
          - wildcards   : like 'A', but it may contain '*' or '?'
                        - '*' matches any number of characters, '?' matches one
                        - like filename wildcard matching

       Examples:

          using ints for values:

               2,6..10      : return [2, 6, 7, 8, 9, 10]
               2,6..10(2)   : return [2, 6, 8, 10]
               2,4@8        : return [2, 8, 8, 8, 8]
               2,5..$(3)    : return [2, 5, 8, 11]   (given imax=11)

           with labels: replace labels by their index in the labels array,
           given labels = ['b0', 'b1', 'b2',
                           'mot01_roll', 'mot02_pitch', 'mot03_yaw',
                           'mot04_RL', 'mot05_AP', 'mot06_IS',
                           'ma', 'mb', 'mc']

               'mot04_RL,mot05_AP' : return [6, 7]
               'm*'           : return [3, 4, 5, 6, 7, 8, 9, 10, 11]
               'mot*_??'      : return [6, 7, 8]
               '??'           : return [0, 1, 2, 9, 10, 11]
               '5,??'         : return [5, 0, 1, 2, 9, 10, 11]
               'mot03*,??'    : return [5, 0, 1, 2, 9, 10, 11]

       - return a list of ints
    """

    newstr = strip_list_brackets(istr, verb)
    slist = newstr.split(',')
    if len(slist) == 0:
        if verb > 1: print("-- empty 1D_ints from string '%s'" % istr)
        return []
    elif verb > 3:
       print("-- decoding stripped list '%s' (nlabs = %d)" \
             % (newstr,len(labels)))
    ilist = []                  # init return list
    for s in slist:
        try:
            if s.find('@') >= 0:        # then expect "A@B"
                [N, val] = [n for n in s.split('@')]
                N = int(N)
                val = to_int_special(val, '$', imax, labels)
                ilist.extend([val for i in range(N)])
                if verb > 2: print("-- decode_1D_ints: @ special %s" % s)
            elif s.find('..') >= 0:     # then expect "A..B"
                pos = s.find('..')
                if s.find('(', pos) > 0:    # look for "A..B(C)"
                   [v1, v2] = [n for n in s.split('..')]
                   v1 = to_int_special(v1, '$', imax, labels)
                   [v2, step] = v2.split('(')
                   v2 = to_int_special(v2, '$', imax, labels)
                   # have start and end values, get step
                   step, junk = step.split(')')
                   step = int(step)
                   if   step > 0: inc = 1
                   elif step < 0: inc = -1
                   else:
                        print("** decode: illegal step of 0 in '%s'" % istr)
                        return []
                   ilist.extend([i for i in range(v1, v2+inc, step)])
                else:
                   [v1, v2] = [n for n in s.split('..')]
                   v1 = to_int_special(v1, '$', imax, labels)
                   v2 = to_int_special(v2, '$', imax, labels)
                   if v1 < v2 : step = 1
                   else:        step = -1
                   ilist.extend([i for i in range(v1, v2+step, step)])
                if verb > 2: print("-- decode_1D_ints: .. special %s" % s)
            elif '*' in s or '?' in s:
                lll = to_intlist_wild(s, labels)
                ilist.extend(lll)
                if verb > 2:
                   print("-- decode_1D_ints: wild special %s, list %s" %(s,lll))
            else:
                ilist.extend([to_int_special(s, '$', imax, labels)])
        except:
            print("** cannot decode_1D '%s' in '%s'" % (s, istr))
            return []
    if verb > 3: print('++ ilist: %s' % ilist)
    del(newstr)
    return ilist

def to_intlist_wild(cval, labels=[]):
   """return the index list of any labels that match cval, including wildcards

      Use '*' and '?' for wildcard matching.
      In the regular expression, replace '*' with '.*', and '?' with '.'.

        cval:   int as character string, or a label
        labels: labels to consider
   """

   cval = cval.replace('*', '.*')
   cval = cval.replace('?', '.')

   # look for any matching label
   ilist = []
   for lind, label in enumerate(labels):
      if re.fullmatch(cval, label):
         ilist.append(lind)
   return ilist

def to_int_special(cval, spec, sint, labels=[]):
   """basically return int(cval), but if cval==spec, return sint

      if labels were passed, allow cval to be one

        cval: int as character string, or a label
        spec: special value as string
        sint: special value as int"""
   if cval == spec:
      return sint
   elif cval[0].isalpha() and cval in labels:
      return labels.index(cval)
   else:
      return int(cval)

def extract_subbrick_selection(sstring):
   """search sstring for something like: [DIGITS*($|DIGITS)]
      - do not consider all DIGITS, '..', ',', '(DIGITS)' pieces,
        just let '*' refer to anything but another '['
   """
   import re
   res = re.search(r'\[\d+[^\[]*]', sstring)
   if res == None: return ''
   return res.group(0)

def strip_list_brackets(istr, verb=1):
   """strip of any [], {}, <> or ## surrounding this string
        - can only be surrounded by whitespace
        - assume only one pair
        - allow the trailing character to be missing
      return the remaining string"""

   if len(istr) < 1: return istr

   icopy = istr.strip()

   # strip any of these pairs
   for pairs in [ ['[',']'],  ['{','}'],  ['<','>'],  ['#','#'] ]:
      if icopy[0] == pairs[0]:
         ind1 = icopy.find(pairs[1], 1)
         if verb > 1: print('-- stripping %s%s at %d,%d in %s' % \
                            (pairs[0],pairs[1],0,ind1,icopy))
         if ind1 > 0: return icopy[1:ind1]
         else:        return icopy[1:]

   if verb > 2: print("-- nothing to strip from '%s'" % icopy)

   return icopy

def replace_n_squeeze(instr, oldstr, newstr):
   """like string.replace(), but remove all spaces around oldstr
      (so probably want some space in newstr)"""
   # while oldstr is found
   #   find last preceding keep posn (before oldstr and spaces)
   #   find next following keep posn (after oldstr and spaces)
   #   set result = result[0:first] + newstr + result[last:]
   newlen = len(newstr)
   result = instr
   posn = result.find(oldstr)
   while posn >= 0:
      rlen = len(result)
      start = posn-1
      while start >= 0 and result[start] == ' ': start -= 1
      if start >= 0: newres = result[0:start+1] + newstr
      else:          newres = newstr
      end = posn + newlen
      while end < rlen and result[end] == ' ': end += 1
      if end < rlen: newres += result[end:]

      result = newres
      posn = result.find(oldstr)

   return result

# ----------------------------------------------------------------------
# line wrapper functions

def list_to_wrapped_command(cname, llist, nindent=10, nextra=3, maxlen=76):
    """return a string that is a 'cname' command, indented by
         nindent, with nextra indentation for option continuation

       This function takes a command and a list of options with parameters,
       and furnishes a wrapped command, where each option entry is on its
       own line, and any option entry line wraps includes nextra indentation.

           - wrap each option line without indentation
           - split across \\\n; ==> have good indentation
           - keep all in array
           
           - join array with indentation and \\\n
           - finally, align the \\\n wrappers
    """
    ntot = nindent + nextra
    isind = ' '*nindent
    istot = ' '*ntot
    jstr  = ' \\\n' + istot

    clist = [cname]
    for line in llist:
       wline = add_line_wrappers(line, maxlen=(maxlen-ntot))
       wlist = wline.split('\\\n')
       clist.extend(wlist)

    return align_wrappers(isind + jstr.join(clist))


# MAIN wrapper: add line wrappers ('\'), and align them all
def add_line_wrappers(commands, wrapstr='\\\n', maxlen=78, verb=1,
                      method='rr'):
    """wrap long lines with 'wrapstr' (probably '\\\n' or just '\n')
       if '\\\n', align all wrapstr strings

       method can be rr or pt
    """
    new_cmd = ''
    posn = 0

    while needs_wrapper(commands, maxlen, posn):
            
        end = find_command_end(commands, posn)

        if not needs_wrapper(commands, maxlen, posn, end): # command is okay
            if end < 0: new_cmd = new_cmd + commands[posn:]
            else      : new_cmd = new_cmd + commands[posn:end+1]
            posn = end+1
            continue

        new_line = insert_wrappers(commands, posn, end, wstring=wrapstr,
                                   maxlen=maxlen, method=method, verb=verb)

        new_cmd += new_line
            
        posn = end + 1     # else, update posn and continue

    result = new_cmd + commands[posn:]

    # wrappers are in, now align them (unless method == 'pt')
    if wrapstr == '\\\n' and method != 'pt':
       return align_wrappers(result)
    else:
       return result

def align_wrappers(command):
    """align all '\\\n' strings to be the largest offset
       from the previous '\n'"""

    # first, find the maximum offset
    posn = 0
    max  = -1
    while 1:
        next = command.find('\n',posn)
        if next < 0: break
        if next > posn and command[next-1] == '\\':  # check against max
            width = next - 1 - posn
            if width > max: max = width
        posn = next + 1 # look past it

    if max < 0: return command  # none found

    # repeat the previous loop, but adding appropriate spaces
    new_cmd = ''
    posn = 0
    while 1:
        next = command.find('\n',posn)
        if next < 0: break
        if next > posn and command[next-1] == '\\':  # check against max
            width = next - 1 - posn
            if width < max:     # then insert appropriate spaces
                new_cmd += command[posn:next-1] + ' '*(max-width) + '\\\n'
                posn = next + 1
                continue

        # just duplicate from the previous posn
        new_cmd += command[posn:next+1]
        posn = next + 1 # look past it

    if posn < len(command): new_cmd += command[posn:]

    return new_cmd

def insert_wrappers(command, start=0, end=-1, wstring='\\\n',
                    maxlen=78, method='rr', verb=1):
    """insert any '\\' chars for the given command
         - insert between start and end positions
         - apply specified wrap string wstring

       return a new string, in any case
    """

    global wrap_verb

    if end < 0: end = len(command) - start - 1

    nfirst = num_leading_line_spaces(command,start,1) # note initial indent
    prefix = get_next_indentation(command,start,end)
    sskip  = nfirst             # number of init spaces expected
    plen   = len(prefix)
    newcmd = ''
    cur    = start

    if verb > 1: print("+d insert wrappers: nfirst=%d, prefix='%s', plen=%d" \
                       % (nfirst, prefix, plen))

    # if P Taylor special, short circuit the rest
    if method == 'pt':
        cline = command[start:end+1]
        clist = cline.replace('\\\n', ' ').split()
        cline = ' '.join(clist)
        short_pre = prefix[:-4]
        rv, cline = lfcs.afni_niceify_cmd_str(cline, comment_start=short_pre)
        return cline + '\n'

    # rewrite: create new command strings after each wrap     29 May 2009
    # (only care where new wrappers are needed)
    while needs_new_wrapper(command,maxlen,cur,end):
        endposn = command.find('\n',cur)

        if needs_new_wrapper(command,maxlen,cur,endposn):
            lposn = find_last_space(command, cur+sskip, endposn, maxlen-sskip)

            # if the last space is farther in than next indent, wrap
            # (adjust initial skip for any indent)
            if sskip+cur < lposn:   # woohoo, wrap away (at lposn)
                newcmd = newcmd + command[cur:lposn+1] + wstring
                # modify command to add prefix, reset end and cur
                command = prefix + command[lposn+1:]
                end = end + plen - (lposn+1)
                sskip = nfirst + plen   # now there is a prefix to skip
                cur = 0
                continue

        # no change:
        # either the line does not need wrapping, or there is no space to do it
        if endposn < 0: endposn = end     # there may not be a '\n'
        newcmd += command[cur:endposn+1]
        cur = endposn + 1

    if cur <= end: newcmd += command[cur:end+1]   # add remaining string

    return newcmd

def get_next_indentation(command,start=0,end=-1):
    """get any '#' plus leading spaces, from beginning or after first '\\\n'"""
    if end < 0: end = len(command) - start - 1

    spaces = num_leading_line_spaces(command,start,1)
    prefix = command[start:start+spaces]+'    ' # grab those spaces, plus 4
    # now check for an indentation prefix
    posn = command.find('\\\n', start)
    pn = command.find('\n', start)      # but don't continue past current line
    if posn >= 0 and posn < pn:
        spaces = num_leading_line_spaces(command,posn+2,1)
        if posn > start and spaces >= 2:
            prefix = command[posn+2:posn+2+spaces] # grab those spaces

    return prefix

def needs_wrapper(command, maxlen=78, start=0, end=-1):
    """does the current string need line wrappers

       a string needs wrapping if there are more than 78 characters between
       any previous newline, and the next newline, wrap, or end"""

    if end < 0: end_posn = len(command) - 1
    else:       end_posn = end

    # cur_posn should always point to the beginning of a line
    cur_posn = start
    remain = end_posn - cur_posn

    # find end of current line (\\\n does not count)
    # (newend will point to '\n', if it exists)
    while cur_posn < end_posn:
        newend = find_command_end(command, cur_posn, check_cmnt=0)
        # print("\n== new length %d, command:\n%s\n===\n\n" \
        #       % (newend-cur_posn,command[cur_posn:newend+1]))
        # if it is long, we want to wrap
        if newend - cur_posn >= maxlen:
           return 1
        # otherwise, see if we are done
        if newend >= end_posn:
           return 0

        # we are not done, adjust cur_posn and continue
        cur_posn = newend+1

    return 0

def needs_new_wrapper(command, maxlen=78, start=0, end=-1):
    """does the current string need NEW line wrappers
       (different to needing ANY)

       a string needs wrapping if there are more than 78 characters between
       any previous newline, and the next newline, wrap, or end"""

    if end < 0: end_posn = len(command) - 1
    else:       end_posn = end

    # cur_posn should always point to the beginning of a line
    cur_posn = start
    remain = end_posn - cur_posn

    while remain > maxlen:
        # find next '\\\n'
        posn = command.find('\\\n', cur_posn)
        if 0 <= posn-cur_posn <= maxlen: # adjust and continue
            cur_posn = posn + 2
            remain = end_posn - cur_posn
            continue

        # else find next '\n'
        posn = command.find('\n', cur_posn)
        if 0 <= posn-cur_posn <= maxlen: # adjust and continue
            cur_posn = posn + 1
            remain = end_posn - cur_posn
            continue

        # otherwise, space means wrap, else not
        if find_next_space(command, cur_posn, 1) > cur_posn: return 1
        return 0

    return 0        # if we get here, line wrapping is not needed

def find_command_end(command, start=0, check_cmnt=1):
    """find the next '\n' that is not preceded by '\\', or return the
       last valid position (length-1)"""

    length = len(command)
    end = start-1   # just to re-init start
    while 1:
        start = end + 1
        end = command.find('\n',start)

        if end < 0: return length-1   # not in command
        elif end > start and command[end-1] == '\\':
            if check_cmnt:
                if length > end+1 and command[start] == '#'   \
                                  and command[end+1] != '#':
                    return end      # since comments cannot wrap
            continue 
        return end              # found

def num_leading_line_spaces(istr,start,pound=0):
    """count the number of leading non-whitespace chars
       (newline chars are not be counted, as they end a line)
       if pound, skip any leading '#'"""

    length = len(istr)
    if start < 0: start = 0
    if length < 1 or length <= start: return 0
    posn = start
    if pound and istr[posn] == '#': posn += 1

    while posn < length and istr[posn].isspace() and istr[posn] != '\n':
        posn += 1

    if posn == length: return 0   # none found
    return posn-start             # index equals num spaces from start

def find_next_space(istr,start,skip_prefix=0):
    """find (index of) first space after start that isn't a newline
       (skip any leading indentation if skip_prefix is set)
       return -1 if none are found"""

    length = len(istr)
    index  = start
    if skip_prefix: index += num_leading_line_spaces(istr,start,1)
    
    while 1:
        if index >= length: break
        if istr[index] != '\n' and istr[index].isspace(): break
        index += 1

    if index >= length : return -1
    return index

def find_last_space(istr,start,end,max_len=-1,stretch=1):
    """find (index of) last space in current line range that isn't a newline
       if stretch and not found, search towards end
       return start-1 if none are found"""

    if end < 0: end = len(istr) - 1
    if max_len >= 0 and end-start >= max_len: index = start+max_len-1
    else:                                     index = end

    posn = index        # store current position in case of stretch
    
    while posn >= start and (istr[posn] == '\n' or not istr[posn].isspace()):
        posn -= 1

    if posn < start and stretch:       # then search towards end
        posn = index
        while posn <= end and (istr[posn] == '\n' or not istr[posn].isspace()):
            posn += 1
        if posn > end: posn = start-1 # still failed

    return posn   # for either success or failure

def nuke_final_whitespace(script, skipchars=[' ', '\t', '\n', '\\'],
                                append='\n\n'):
    """replace final white characters with newlines"""
    slen = len(script)
    ind = slen-1
    while ind > 0 and script[ind] in skipchars: ind -= 1

    return script[0:ind+1]+append

# end line_wrapper functions
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# other functions

# 17 May, 2008 [rickr]
def vals_are_multiples(num, vals, digits=4):
    """decide whether every value in 'vals' is a multiple of 'num'
       (vals can be a single float or a list of them)

       Note, 'digits' can be used to specify the number of digits of accuracy
       in the test to see if a ratio is integral.  For example:
           vals_are_multiples(1.1, 3.3001, 3) == 1
           vals_are_multiples(1.1, 3.3001, 4) == 0

       return 1 if true, 0 otherwise (including error)"""

    if num == 0.0: return 0

    try:
        l = len(vals)
        vlist = vals
    except:
        vlist = [vals]

    for val in vlist:
        rat = val/num
        rem = rat - int(rat)

        if round(rem,digits) != 0.0: return 0

    return 1

def vals_are_constant(vlist, cval=None):
   """vals_are_constant(vlist, cval=None)

      determine whether every value in vlist is equal to cval
      (if cval == None, use vlist[0])"""

   if vlist == None: return 1
   if len(vlist) < 2: return 1

   if cval == None: cval = vlist[0]

   for val in vlist:
      if val != cval: return 0
   return 1

def vals_are_positive(vlist):
   """determine whether every value in vlist is positive"""
   for val in vlist:
      if val <= 0: return 0
   return 1

def vals_are_0_1(vlist):
   """determine whether every value in vlist is either 0 or 1"""
   for val in vlist:
      if val != 0 and val != 1: return 0
   return 1

def vals_are_sorted(vlist, reverse=0):
   """determine whether values non-decreasing (or non-inc if reverse)"""
   if vlist == None: return 1
   if len(vlist) < 2: return 1

   rval = 1
   try:
      for ind in range(len(vlist)-1):
         if reverse:
            if vlist[ind] < vlist[ind+1]:
               rval = 0
               break
         else:
            if vlist[ind] > vlist[ind+1]:
               rval = 0
               break
   except:
      print("** failed to detect sorting in list: %s" % vlist)
      rval = 0
      
   return rval

def vals_are_increasing(vlist, reverse=0):
   """determine whether values strictly increasing (or dec if reverse)"""
   if vlist == None: return 1
   if len(vlist) < 2: return 1

   rval = 1
   try:
      for ind in range(len(vlist)-1):
         if reverse:
            if vlist[ind] <= vlist[ind+1]:
               rval = 0
               break
         else: # verify increasing
            if vlist[ind] >= vlist[ind+1]:
               rval = 0
               break
   except:
      print("** failed to detect sorting in list: %s" % vlist)
      rval = 0
      
   return rval

def vals_are_unique(vlist, dosort=1):
   """determine whether (possibly unsorted) values are unique
      - use memory to go for N*log(N) speed"""

   if vlist == None: return 1
   if len(vlist) < 2: return 1

   # copy and sort
   dupe = vlist[:]
   if dosort: dupe.sort()

   rval = 1
   try:
      for ind in range(len(dupe)-1):
         if dupe[ind] == dupe[ind+1]:
            rval = 0
            break
   except:
      print("** uniq: failed to compare list elements in %s" % vlist)
      rval = 0

   del(dupe)
      
   return rval

def lists_are_same(list1, list2, epsilon=0, doabs=0):
   """return 1 if the lists have similar values, else 0

      similar means difference <= epsilon
   """
   if not list1 and not list2: return 1
   if not list1: return 0
   if not list2: return 0
   if len(list1) != len(list2): return 0

   for ind in range(len(list1)):
      if doabs:
         v1 = abs(list1[ind])
         v2 = abs(list2[ind])
      else:
         v1 = list1[ind]
         v2 = list2[ind]
      if v1 != v2: return 0
      if epsilon:
         if abs(v1-v2) > epsilon: return 0

   return 1

def list_intersect(listA, listB, sort=1):
   """return a list of the elements that are in both lists

      if sort, sort the result
   """

   # if either is empty, the match list is as well
   if not listA or not listB: return []

   rlist = [v for v in listA if v in listB]

   if sort: rlist.sort()

   return rlist

def list_diff(listA, listB, dtype='A-B', sort=1):
   """return a list of the elements differ between the lists

      dtype:    A-B     : return elements in listA that are not list B
                B-A     : return elements in listB that are not list A
                all     : return all diffs (same as A-B and B-A)

      return a list of newly created elements
   """

   # keep logic simple

   if dtype == 'A-B':
      rlist = [v for v in listA if v not in listB]
   elif dtype == 'B-A':
      rlist = [v for v in listB if v not in listA]
   else:
      rlist = [v for v in listA if v not in listB]
      rlist.extend([v for v in listB if v not in listA])

   if sort:
      rlist.sort()

   return rlist

def string_to_float_list(fstring):
   """return a list of floats, converted from the string
      return None on error
   """

   if type(fstring) != str: return None
   slist = fstring.split()

   if len(slist) == 0: return []

   try: flist = [float(sval) for sval in slist]
   except: return None

   return flist

def string_to_type_list(sdata, dtype=float):
   """return a list of dtype, converted from the string
      return None on error
   """

   if type(sdata) != str: return None
   slist = sdata.split()

   if len(slist) == 0: return []

   # if going to int, use float as an intermediate step
   if dtype == int:
      try: slist = [float(sval) for sval in slist]
      except: return None

   try: dlist = [dtype(sval) for sval in slist]
   except: return None

   return dlist

def float_list_string(vals, mesg='', nchar=7, ndec=3, nspaces=2, left=0, sep=' '):
   """return a string to display the floats:
        vals    : the list of float values
        mesg    : []  message to precede values
        nchar   : [7] number of characters to display per float
        ndec    : [3] number of decimal places to print to
        nspaces : [2] number of spaces between each float

        left    : [0] whether to left justify values
        sep     : [ ] separator for converted strings
   """

   if left: form = '%-*.*f%*s'
   else:    form = '%*.*f%*s'

   svals = [form % (nchar, ndec, val, nspaces, '') for val in vals]
   istr = mesg + sep.join(svals)

   return istr

def gen_float_list_string(vals, mesg='', nchar=0, left=0, sep=' '):
   """mesg is printed first, if nchar>0, it is min char width"""

   istr = mesg

   if left: form = '%-'
   else:    form = '%'

   if nchar > 0:
      form += '*g'
      svals = [form % (nchar, val) for val in vals]
   else:
      form += 'g'
      svals = [form % val for val in vals]

   istr = mesg + sep.join(svals)

   return istr

def int_list_string(ilist, mesg='', nchar=0, sepstr=' '):
   """like float list string, but use general printing
      (mesg is printed first, if nchar>0, it is min char width)"""

   istr = mesg

   if nchar > 0: slist = ['%*d' % (nchar, val) for val in ilist]
   else:         slist = ['%d' % val for val in ilist]
   istr += sepstr.join(slist)

   return istr

def invert_int_list(ilist, top=-1, bot=0):
   """invert the integer list with respect to bot and top
      i.e. return a list of integers from bot to top that are not in
           the passed list
   """
   if top < bot:
      print('** invert_int_list requires bot<=top (have %d, %d)' % (bot, top))
      return []

   return [ind for ind in range(bot, top+1) if not ind in ilist]

def is_valid_int_list(ldata, imin=0, imax=-1, whine=0):
   """check whether:
        o  ldata is a of type []
        o  values are of type int
        o  values are in within imin..imax (only if imin <= imax)
        o  if whine: complain on error
      return 1 on true, 0 on false"""

   if not ldata or type(ldata) != type([]):
      if whine: print("** not valid as a list: '%s'" % ldata)

   for ind in range(len(ldata)):
      val = ldata[ind]
      if type(val) != type(0):
         if whine: print("** non-int value %d in int list (@ %d)" % (val,ind))
         return 0
      if imin <= imax: # then also test bounds
         if val < imin:
            if whine: print("** list value %d not in [%d,%d]" %(val,imin,imax))
            return 0
         elif val > imax:
            if whine: print("** list value %d not in [%d,%d]" %(val,imin,imax))
            return 0
   return 1

def data_to_hex_str(data):
   """convert raw data to hex string in groups of 4 bytes"""

   if not data: return ''

   dlen = len(data)             # total length in bytes
   groups = (dlen+3) // 4       # number of 4-byte blocks to create
   remain = dlen
   retstr = ''  # return string

   for group in range(groups):
      if group > 0: retstr += ' '
      retstr += '0x'
      if remain >= 4: llen = 4
      else:           llen = remain

      for ind in range(llen):
         retstr += '%02x' % data[dlen-remain+ind]

      remain -= llen

   return retstr

def section_divider(hname='', maxlen=74, hchar='=', endchar=''):
    """return a title string of 'hchar's with the middle chars set to 'hname'
       if endchar is set, put at both ends of header
       e.g. section_divider('volreg', endchar='##') """
    if len(hname) > 0: name = ' %s ' % hname
    else:              name = ''

    if endchar != '': maxlen -= 2*len(endchar)
    rmlen = len(name)
    if rmlen >= maxlen:
        print("** S_DIV_STR, rmlen=%d exceeds maxlen=%d" % (rmlen, maxlen))
        return name
    prelen  = (maxlen - rmlen) // 2     # basically half the chars
    postlen = maxlen - rmlen - prelen   # other 'half'

    return endchar + prelen*hchar + name + postlen*hchar + endchar

def get_command_str(args=[], preamble=1, comment=1, quotize=1, wrap=1):
    """return a script generation command

        args:           arguments to apply
        preamble:       start with "script generated by..."
        comment:        have text '#' commented out
        quotize:        try to quotize any arguments that need it
        wrap:           add line wrappers
    """

    if args == []: args = sys.argv

    if preamble: hstr = '\n# %s\n# script generated by the command:\n#\n' \
                        % section_divider()
    else:        hstr = ''

    if comment: cpre = '# '
    else:       cpre = ''

    # note command and args
    cmd = os.path.basename(args[0])
    if quotize: args = ' '.join(quotize_list(args[1:],''))
    else:       args = ' '.join(args[1:],'')

    cstr = '%s%s%s %s\n' % (hstr, cpre, cmd, args)

    if wrap: return add_line_wrappers(cstr)
    else:    return cstr

def max_len_in_list(vlist):
    mval = 0
    try:
       for v in vlist:
          if len(v) > mval: mval = len(v)
    except:
       print('** max_len_in_list: cannot compute lengths')
    return mval

def get_rank(data, style='dense', reverse=0, uniq=0):
   """return the rank order of indices given values,
      i.e. for each value, show its ordered index
      e.g. 3.4 -0.3 4.9 2.0   ==>   2 0 3 1

      Sort the data, along with indices, then repeat on the newly
      sorted indices.  If not uniq, set vals to minimum for repeated set.
      Note that sorting lists sorts on both elements, so first is smallest.

      styles:  (e.g. given input -3 5 5 6, result is ...)

         competition:  competition ranking (result 1 2 2 4 )
         dense:        dense ranking (result 1 2 2 3 )
                       {also default from 3dmerge and 3dRank}

      return status (0=success) and the index order
   """

   dlen = len(data)

   # maybe reverse the sort order
   if reverse:
      maxval = max(data)
      dd = [maxval-val for val in data]
   else: dd = data

   # sort data and original position
   dd = [[dd[ind], ind] for ind in range(dlen)]
   dd.sort()

   # invert position list by repeating above, but with index list as data
   # (bring original data along for non-uniq case)
   dd = [[dd[ind][1], ind, dd[ind][0]] for ind in range(dlen)]

   # deal with repeats: maybe modify d[1] from ind, depending on style
   if not uniq:
      if style == 'dense':
         cind = dd[0][1] # must be 0
         for ind in range(dlen-1):      # compare next to current
            if dd[ind+1][2] == dd[ind][2]:
               dd[ind+1][1] = cind
            else:
               cind += 1
               dd[ind+1][1] = cind
      elif style == 'competition':
         for ind in range(dlen-1):      # compare next to current
            if dd[ind+1][2] == dd[ind][2]:
               dd[ind+1][1] = dd[ind][1]
      else:
         print("** UTIL.GR: invalid style '%s'" % style)
         return 1, []

   dd.sort()

   return 0, [dd[ind][1] for ind in range(dlen)]

# ----------------------------------------------------------------------
# wildcard construction functions
# ----------------------------------------------------------------------

def first_last_match_strs(slist):
   """given a list of strings, return the first and last consistent strings
      (i.e. all strings have the form first*last)

        e.g. given ['subj_A1.txt', 'subj_B4.txt', 'subj_A2.txt' ]
             return 'subj_' and '.txt'
   """

   if type(slist) != list:
      print('** FL match strings requires a list')
      return '', ''

   if not slist: return '', ''

   maxlen = len(slist[0])
   hmatch = maxlen              # let them shrink
   tmatch = maxlen
   for sind in range(1, len(slist)):
      if slist[0] == slist[sind]: continue

      hmatch = min(hmatch, len(slist[sind]))
      tmatch = min(tmatch, len(slist[sind]))

      # find first left diff
      i = 0
      while i < hmatch:
         if slist[sind][i] != slist[0][i]: break
         i += 1
      hmatch = min(hmatch, i)

      # find first right diff (index from 1)
      i = 1
      while i <= tmatch:
         if slist[sind][-i] != slist[0][-i]: break
         i += 1
      tmatch = min(tmatch, i-1)

   if hmatch+tmatch > maxlen:           # weird, but constructable
      tmatch = maxlen - hmatch          # so shrink to fit

   # list[-0:] is not empty but is the whole list
   if tmatch > 0: tstr = slist[0][-tmatch:]
   else:          tstr = ''

   return slist[0][0:hmatch], tstr

def list_files_by_glob(L, sort=False, exit_on_miss=False) :
    """Input a list L of one or more strings to glob (fiiine, the input L
    can also be a string, which will just be made into a list
    immediately, L = [L]), and then glob each one, sort that piece's
    results, and append it to the list of results.  Return the full
    list of files.

    The program can be made to exist if any piece fails to find a
    result, by setting 'exit_on_miss=True'.

    Sorting of the full/final list is NOT done by default (can be
    turned on with 'sort=True').

    """

    if not(L): return []

    if   type(L) == str :      L = [L]
    elif type(L) != list :     BASE.EP("Input a list (or str)")

    N        = len(L)
    cbad     = 0
    out      = []

    for ii in range(N):
        sublist = L[ii].split()
        for item in sublist:
            gout = glob.glob(item)
            if not(gout) :
                cbad+= 1 
                BASE.WP("No files for this part of list: {}".format(item))
            else:
                gout.sort()
                out.extend(gout)

    if cbad and exit_on_miss :
        BASE.EP("No findings for {} parts of this list.  Bye.".format(cbad))
    
    if sort :
        out.sort()

    return out

def glob2stdout(globlist):
   """given a list of glob forms, print all matches to stdout

      This is meant to be a stream workaround to shell errors
      like, "Argument list too long".

      echo 'd1/*.dcm' 'd2/*.dcm'            \\
        | afni_python_wrapper.py -listfunc glob2stdout -
   """
   for gform in globlist:
      for fname in glob.glob(gform):
         print(fname)

def glob_form_from_list(slist):
   """given a list of strings, return a glob form

        e.g. given ['subjA1.txt', 'subjB4.txt', 'subjA2.txt' ]
             return 'subj*.txt'

      Somewhat opposite list_minus_glob_form().
   """

   if len(slist) == 0: return ''
   if vals_are_constant(slist): return slist[0]

   first, last = first_last_match_strs(slist)
   if not first and not last: return '' # failure
   globstr = '%s*%s' % (first,last)

   return globstr

def glob_form_matches_list(slist, ordered=1):
   """given a list of strings, make a glob form, and then test that against
      the actual files on disk

      if ordered: files must match exactly (i.e. slist must be sorted)
      else:       slist does not need to be sorted
   """

   slen = len(slist)

   # check trivial cases of lengths 0 and 1
   if slen == 0: return 1
   if slen == 1:
      if os.path.isfile(slist[0]): return 1
      else:                        return 0

   globstr = glob_form_from_list(slist)
   glist = glob.glob(globstr)
   glist.sort()

   # quick check: lengths must match
   if len(glist) != slen: return 0

   if ordered:
      inlist = slist
   else: 
      inlist = slist[:]
      inlist.sort()

   # now files must match exactly (between inlist and glist)
   for ind in range(slen):
      if glist[ind] != inlist[ind]: return 0

   # they must match
   return 1
   

def list_minus_glob_form(inlist, hpad=0, tpad=0, keep_dent_pre=0, strip=''):
   """given a list of strings, return the inner part of the list that varies
      (i.e. remove the consistent head and tail elements)

        e.g. given ['subjA1.txt', 'subjB4.txt', 'subjA2.txt' ]
             return [ 'A1', 'B4', 'A2' ]

      hpad NPAD         : number of characters to pad at prefix
      tpad NPAD         : number of characters to pad at suffix
      keep_dent_pre     : possibly keep directory entry prefix
                          0 : never
                          1 : keep entire prefix from directory entry
                          2 : do it if dir ent starts with sub
      strip             : one of ['', 'dir', 'file', 'ext', 'fext']

      If hpad > 0, then pad with that many characters back into the head
      element.  Similarly, tpad pads forward into the tail.

        e.g. given ['subjA1.txt', 'subjB4.txt', 'subjA2.txt' ]
             if hpad = 926 (or 4 :) and tpad = 1,
             return [ 'subjA1.', 'subjB4.', 'subjA2.' ]

      If keep_dent_pre is set, then (if '/' is found) decrement hlen until 
      that '/'.  If '/' is not found, start from the beginning.

        e.g. given ['dir/subjA1.txt', 'dir/subjB4.txt', 'dir/subjA2.txt' ]
                -> return = [ 'A1.', 'B4.', 'A2.' ]
             with keep_dent_pre == 1:
                -> return = [ 'subjA1.', 'subjB4.', 'subjA2.' ]

      Somewhat opposite glob_form_from_list().
   """

   if len(inlist) <= 1: return inlist

   # init with original
   slist = inlist

   # maybe make a new list of stripped elements
   stripnames = ['dir', 'file', 'ext', 'fext']
   if strip != '' and strip not in stripnames:
      print('** LMGF: bad strip %s' % strip)
      strip = ''

   if strip in stripnames:
      ss = []
      for inname in inlist:
         if strip == 'dir':
            dname, fname = os.path.split(inname)
            ss.append(fname)
         elif strip == 'file':
            dname, fname = os.path.split(inname)
            ss.append(dname)
         elif strip == 'ext':
            fff, ext = os.path.splittext(inname)
            ss.append(fff)
         elif strip == 'fext':
            fff, ext = os.path.splittext(inname)
            ss.append(fff)
         else:
            print('** LMGF: doubly bad strip %s' % strip)
            break
      # check for success
      if len(ss) == len(slist): slist = ss

   if hpad < 0 or tpad < 0:
      print('** list_minus_glob_form: hpad/tpad must be non-negative')
      hpad = 0 ; tpad = 0

   # get head, tail and note lengths
   head, tail = first_last_match_strs(slist)
   hlen = len(head)
   tlen = len(tail)

   # adjust by padding, but do not go negative
   if hpad >= hlen: hlen = 0
   else:            hlen -= hpad
   if tpad >= tlen: tlen = 0
   else:            tlen -= tpad

   # apply directory entry prefix, if requested
   if keep_dent_pre:
      s = slist[0]
      posn = s.rfind('/', 0, hlen)
      # if found, start at position to right of it
      # otherwise, use entire prefix
      if posn >= 0: htmp = posn + 1
      else:         htmp = 0

      # apply unless KDP == 2 and not 'subj'
      if keep_dent_pre != 2 or s[htmp:htmp+3] == 'sub':
         hlen = htmp

   # and return the list of center strings
   if tlen == 0: return [ s[hlen:]      for s in slist ]
   else:         return [ s[hlen:-tlen] for s in slist ]

def glob_list_minus_pref_suf(pref, suf):
   """just strip the prefix and suffix from string list elements
      (for now, assume they are there)
   """
   glist = glob.glob('%s*%s' % (pref, suf))

   plen = len(pref)
   slen = len(suf)

   return [d[plen:-slen] for d in glist]

def list_minus_pref_suf(slist, pref, suf, stripdir=1):
   """just strip the prefix and suffix from string list elements

      if stripdir, remove leading directories

      return status, stripped list

      status =  0 : all strings have prefix and suffix
                1 : not all do
               -1 : on error
   """

   plen = len(pref)
   slen = len(suf)

   # possibly strip of directory names
   if stripdir:
      flist = []
      for sname in slist:
         dd, ff = os.path.split(sname)
         flist.append(ff)
   else: flist = slist

   
   rv = 0
   rlist = []
   for fname in flist:
      if fname.startswith(pref): poff = plen
      else:                      poff = 0

      if fname.endswith(suf): soff = slen
      else:                   soff = 0

      if soff: rlist.append(fname[poff:-soff])
      else:    rlist.append(fname[poff:])

      if not poff or not soff: rv = 1

   return rv, rlist

def okay_as_lr_spec_names(fnames, verb=0):
   """check that names are okay as surface spec files, e.g. for afni_proc.py
        - must be 1 or 2 file names
        - must contain lh and rh, respectively
        - must otherwise be identical

        if verb, output why failure occurs
        return 1 if okay, 0 otherwise
   """
   nfiles = len(fnames)
   if nfiles == 0: return 1     # no problems, anyway
   if nfiles > 2:
      if verb: print('** only 1 or 2 spec files allowed (have %d)' % nfiles)
      return 0

   if nfiles == 1:
      if fnames[0].find('lh')>=0 or fnames[0].find('rh')>=0: return 1 # success
      # failure
      if verb: print("** spec file '%s' missing 'lh' or 'rh'" % fnames[0])
      return 0

   # so we have 2 files, get the varying part

   # - tpad=1 to include following 'h'
   # - do not return dir entry prefix (e.g. avoid sub-*)
   hlist = list_minus_glob_form(fnames, tpad=1, keep_dent_pre=0)

   for h in hlist:
      if h != 'rh' and h != 'lh':
         if verb: print('** multiple spec files must only differ in lh vs. rh')
         return 0

   return 1

def make_spec_var(fnames, vname='hemi'):
   """return a spec file variable and a list of replaced hemispheres

        e.g. make_spec_var(['surface.lh.spec', 'surface.rh.spec']) returns
                surface.${hemi}.spec, ['lh', 'rh']
      given 1 or 2 spec file names, return a single variable that
      represents them with $vname replacing the lh or rh

      return '' on failure
   """
   if not okay_as_lr_spec_names(fnames): return '', []

   nfiles = len(fnames)
   if nfiles == 0 or nfiles > 2: return '', []

   sfile = fnames[0]

   if nfiles == 1:
      # just find lh or rh and replace it
      hh = 'lh'
      posn = sfile.find(hh)
      if posn < 0:
         hh = 'rh'
         posn = sfile.find(hh)
      if posn < 0: return '', [] # should not happen

      return sfile[0:posn] + '${%s}'%vname + sfile[posn+2:], [hh]

   # so nfiles == 2, use glob

   head, tail = first_last_match_strs(fnames)
   hlen = len(head)

   hemi = sfile[hlen:hlen+2]
   if hemi != 'lh' and hemi != 'rh':
      print('** MSV: bad lh/rh search from spec files: %s' % fnames)
      return '', []

   return sfile[0:hlen] + '${%s}'%vname + sfile[hlen+2:], ['lh', 'rh']

def parse_as_stim_list(flist):
   """parse filename list as PREFIX.INDEX.LABEL.SUFFIX, where the separators
        can be '.', '_' or nothing (though ignore PREFIX and SUFFIX, as well
        as those separators)

        - strip PREFIX and SUFFIX (those are garbage)
          (if SUFFIX has a '.' after position 0, adjust the SUFFIX)
        - strip any leading digits as INDEXes

      return Nx2 table where if one column entry is filled, they all are
             (None on failure for form a complete table)
             (note: blank INDEX or LABEL is okay, but they must all be)
   """

   if len(flist) < 1: return []

   # first get PREFIX and SUFFIX
   prefix, suffix = first_last_match_strs(flist)

   # if suffix contains an extension, make the suffix into the extension
   dot = suffix.find('.')
   if dot < 0: dot = 0

   # strip prefix, suffix: might include part of 'suffix' in label
   inner_list = list_minus_glob_form(flist, tpad=dot)

   # then make table of the form <NUMBER><SEP><LABEL>
   s_table = [list(_parse_leading_int(name)) for name in inner_list]

   # if any number does not exist, just return inner_list as LABELs
   for entry in s_table:
      if entry[0] < 0: return [[-1, label] for label in inner_list]

   # return INDEX and LABEL (no SEP)
   return [[entry[0], entry[2]] for entry in s_table]

def _parse_leading_int(name, seplist=['.','_','-']):
   """assuming name is a string starting with digits, parse name into
      val, sep, suffix

        val    = -1 if name does not start with a digit
        sep    = one of {'.', '_', '-', ''}
        suffix = whatever remains after 'sep'
   """

   nlen = len(name)

   if nlen < 1: return -1, '', ''

   # first strip of any leading (non-negative) integer
   posn = 0     # count leading digits
   for char in name:
      if char.isdigit(): posn += 1
      else:              break

   if posn == 0: val = -1
   else:
      try: val = int(name[0:posn])
      except:
         print("** _parse_leading_int: can't parse int from %s" % name)
         return

   # if only a number, we're outta here
   if posn == nlen: return val, '', ''

   # note any separator
   if name[posn] in seplist:
      sep = name[posn]
      posn += 1
   else:
      sep = ''

   # aaaaaand, we're done
   return val, sep, name[posn:]

def glob_form_has_match(form):
   """see if anything at all exists according to this glob form"""
   glist = glob.glob(form)
   glen = len(glist)
   del(glist)
   if glen > 0: return 1
   return 0

def executable_dir(ename=''):
   """return the directory whre the ename program is located
      (by default, use argv[0])"""
   if ename == '': ee = sys.argv[0]
   else:           ee = ename

   dname = os.path.dirname(ee)
   return os.path.abspath(dname)

def common_dir(flist):
   """return the directory name that is common to all files (unless trivial)"""
   dname, junk = first_last_match_strs(flist)
   if len(dname) > 0 and dname[-1] == '/': dname = dname[0:-1]
   # enter the JR scenario:
   #   if dname IS a directory, make sure the flist items are under it
   #   - consider: data/s1001/... data/s1002/... data/s1003/...
   #   - so dname = data/s100
   #   - BUT, what if dname exists as a separate subject directory?!?
   if not os.path.isdir(dname):
      dname = os.path.dirname(dname)
   else:
      # if we are in JR scenario, we still want to use dirname()
      dlen = len(dname)
      # if dname ends in / or flist[0] continues with /, we are okay
      if dname[-1] != '/' and len(flist[0]) > dlen:
         if flist[0][dlen] != '/':
            # aha!
            dname = os.path.dirname(dname)

   if is_trivial_dir(dname): return ''
   return dname

def common_parent_dirs(flists):
   """return parent directories

      flists = lists of file names (each element should be a list)

      return:
         top_dir    (common to all parents (files), '' if not used)
         parent_dir (for each flist, common parent)
         short_dir  (for each flist, common parent under top_dir)
         short_name (for each flist, file names under parent dirs)

      if top_dir has at least 2 levels, use it
   """
   if type(flists) != list:
      print('** common_parent_dirs: bad flists type')
      return None, None, None, None
   for ind in range(len(flists)):
      flist = flists[ind]
      if type(flist) != list:
         print('** common_parent_dirs: bad flist[%d] type' % ind)
         return None, None, None, None

   # get top_dir and parents
   all_pars    = []
   par_dirs    = []
   short_names = []
   for flist in flists:
      # track parent dirs
      parent = common_dir(flist)
      if parent == '/' or is_trivial_dir(parent): parent = ''
      par_dirs.append(parent)

      # and make short names
      plen = len(parent)
      if plen > 0: start = plen+1
      else:        start = 0
      short_names.append([fname[start:] for fname in flist])

   # top is common to all parents
   top_dir = common_dir(par_dirs)
   if top_dir.count('/') <= 1: top_dir = ''

   # now get all short dir names, under top dir
   if top_dir == '': short_dirs = par_dirs
   else: short_dirs = [child_dir_name(top_dir, pdir) for pdir in par_dirs]

   return top_dir, par_dirs, short_dirs, short_names

def child_dir_name(parent, child):
   """return the child directory name truncated under the parent"""
   if parent == '' or child == '': return child
   plen = len(parent)
   clen = len(child)

   if child[0:plen] != parent: return child     # not a proper child

   # return everything after separator
   if clen < plen + 2: return '.'               # trivial as child
   else:               return child[plen+1:]    # remove parent portion

def is_trivial_dir(dname):
   """input a string
      return 1 if dname is empty or '.'
      else return 0
   """
   if dname == None: return 1
   if dname == '' or dname == '.' or dname == './' : return 1

   return 0

def flist_to_table_pieces(flist):
   """dissect a file list
      input: list of file names
      output:
        - common directory name
        - short name list (names after common directory)
        - glob string from short names
      note: short names will be new data, never just references to input
   """
   if len(flist) == 0: return '', [], ''

   ddir = common_dir(flist)
   dirlen = len(ddir)
   if dirlen > 0: snames = [dset[dirlen+1:] for dset in flist]
   else:          snames = [dset[:]         for dset in flist]

   globstr = glob_form_from_list(snames)

   return ddir, snames, globstr

def get_ids_from_dsets(dsets, prefix='', suffix='', hpad=0, tpad=0, verb=1):
   """return a list of subject IDs corresponding to the datasets

      Make sure we have afni_name objects to take the prefix from.
      If simple strings, turn into afni_names.

      Try list_minus_glob_form on the datasets.  If that fails, try
      on the directories.

      prefix, suffix: attach these to the resulting IDs
      hpad, tpad:     padding to pass to list_minus_glob_form

      return None on failure
   """
   if hpad < 0 or tpad < 0:
      print('** get_ids_from_dsets: will not apply negative padding')
      hpad, tpad = 0, 0

   if len(dsets) == 0: return None

   # be more aggressive, use dataset prefix names
   # dlist = [dset.split('/')[-1] for dset in dsets]
   if type(dsets[0]) == str: 
      nlist = [BASE.afni_name(dset) for dset in dsets]
   elif isinstance(dsets[0], BASE.afni_name):
      nlist = dsets
   else:
      print('** GIFD: invalid type for dset list, have value %s' % dsets[0])
      return None

   dlist = [dname.prefix for dname in nlist]

   # if nothing to come from file prefixes, try the complete path names
   if vals_are_constant(dlist): dlist = dsets

   slist = list_minus_glob_form(dlist, hpad, tpad, keep_dent_pre=2)

   # do some error checking
   for val in slist:
      if '/' in val:            # no directories
         if verb > 0: print('** GIFD: IDs would have directories')
         return None

   if len(slist) != len(dsets): # appropriate number of entries
      if verb > 0: print('** GIFD: length mismatch getting IDs')
      return None

   if not vals_are_unique(slist):
      if verb > 0: print('** GIFD: final IDs are not unique')
      return None

   return slist

def insensitive_word_pattern(word):
   """replace each alphabetical char with [Ul], an upper/lower search pair
      use of 'either' suggested by G Irving at stackoverflow.com
   """
   def either(c):
      if c.isalpha: return '[%s%s]'%(c.lower(),c.upper())
      else:         return c
   return ''.join(map(either,word))
   
def insensitive_glob(pattern):
   """return glob.glob, but where every alphabetic character is
      replaced by lower/upper pair tests
   """
   return glob.glob(insensitive_word_pattern(pattern))


def search_path_dirs(word, mtype=0, casematch=1):
   """return status and list of matching files

      Could be more efficient, esp. with mtype=exact and casematch set, but
      I will strive for simplicity and consistency and see how it goes.

        mtype     : 0 = match any sub-word (i.e. look for DIR/*word*)
                    1 = exact match (i.e. no wildcard, look for DIR/word)
                    2 = prefix match (i.e. look for DIR/word*)
        casematch : flag: if set, case must match
                          else, 'word' letters can be either case
   """
   try:
      plist = os.environ['PATH'].split(os.pathsep)
   except:
      print('** search_path_dirs: no PATH var')
      return 1, []

   # if no casematch, look for upper/lower pairs
   if casematch: wpat = word
   else:         wpat = insensitive_word_pattern(word)

   # if not exact, surround with wildcard pattern
   if   mtype == 0: form = '%s/*%s*'    # any sub-word
   elif mtype == 1: form = '%s/%s'      # exact match
   elif mtype == 2: form = '%s/%s*'     # prefix match
   else:
      print('** search_path_dirs: illegal mtype = %s' % mtype)
      return 1, []

   # now just search for matches
   rlist = []
   for pdir in plist:
      glist = glob.glob(form % (pdir, wpat))
      glist.sort()
      if len(glist) > 0: rlist.extend(glist)

   # make a new list based on os.path.realpath, to avoid links
   rlist = [os.path.realpath(pfile) for pfile in rlist]

   return 0, get_unique_sublist(rlist)

def which(pname):
   """like unix which command: return the first 'pname' in PATH
      (this is a simplified version of search_path_dirs())

      return a simple string, empty on error
   """
   try:
      plist = os.environ['PATH'].split(os.pathsep)
   except:
      print('** which_prog: no PATH var')
      return ''

   for pdir in plist:
      # accept pname having a path
      search = os.path.join(pdir, pname)
      if os.path.isfile(search) and os.access(search, os.X_OK):
         return search

   return ''
   

def num_found_in_path(word, mtype=0, casematch=1):
   """a simple wrapper to print search_path_dirs results

      Search for given 'word' in path, and print out list elements
      with element prefix of 'indent'.

        mtype     : 0 = match any sub-word (i.e. look for DIR/*word*)
                    1 = exact match (i.e. no wildcard, look for DIR/word)
                    2 = prefix match (i.e. look for DIR/word*)
        casematch : flag: if set, case must match
                          else, 'word' letters can be either case
        indent    : prefix/separator string for list elements
   """
   rv, rlist = search_path_dirs(word, mtype=mtype, casematch=casematch)
   if rv: return 0
   return len(rlist)

def show_found_in_path(word, mtype=0, casematch=1, indent='\n   '):
   """a simple wrapper to print search_path_dirs results

      Search for given 'word' in path, and print out list elements
      with element prefix of 'indent'.

        mtype     : 0 = match any sub-word (i.e. look for DIR/*word*)
                    1 = exact match (i.e. no wildcard, look for DIR/word)
                    2 = prefix match (i.e. look for DIR/word*)
        casematch : flag: if set, case must match
                          else, 'word' letters can be either case
        indent    : prefix/separator string for list elements
   """
   rv, rlist = search_path_dirs(word, mtype=mtype, casematch=casematch)
   if not rv: print(indent+indent.join(rlist))

# ----------------------------------------------------------------------
# mathematical functions:
#    vector routines: sum, sum squares, mean, demean
# ----------------------------------------------------------------------

def loc_sum(vals):
   """in case 'sum' does not exist, such as on old machines"""

   try: tot = sum(vals)
   except:
      tot = 0
      for val in vals: tot += val
   return tot

def sumsq(vals):
   """return the sum of the squared values"""
   ssq = 0
   for val in vals: ssq += (val*val)
   return ssq

def euclidean_norm(vals):
   """name is toooooo long"""

   if len(vals) < 1: return 0.0
   return math.sqrt(sumsq(vals))

def L2_norm(vals):
   return euclidean_norm(vals)

def weighted_enorm(vals, weights):

   if len(vals) < 1: return 0.0
   if len(vals) != len(weights): return 0.0

   sum = 0.0
   for ind in range(len(vals)):
      vv = vals[ind]*weights[ind]
      sum += vv*vv
   return math.sqrt(sum)

def dotprod(v1,v2):
   """compute the dot product of 2 vectors"""
   try: dsum = loc_sum([v1[i]*v2[i] for i in range(len(v1))])
   except:
      print('** cannot take dotprod() of these elements')
      dsum = 0
   return dsum

def affine_to_params_6(avec, verb=1):
   """convert rotation/shift affine "matrix" to 6 parameters
      (e.g. as in 3dvolreg 1Dmatrix format to 1Dfile format)

      matvec: length 12+ vector (row major order)
      return: length 6 param vector:
        roll, pitch, yaw, dx, dy, dz
   """

   rvec = [0.0]*6

   if len(avec) < 12:
      print('** affine_to_params_6: requires length 12+ vector, have %d' \
            % len(avec))
      return rvec

   # rotations
   rvec[0] = 180.0/math.pi * math.atan2(avec[9], avec[10])
   rvec[1] = 180.0/math.pi *-math.asin (avec[8])
   rvec[2] = 180.0/math.pi * math.atan2(avec[4], avec[0])

   # deltas
   rvec[3] = avec[3]
   rvec[4] = avec[7]
   rvec[5] = avec[11]

   return rvec

def maxabs(vals):
   """convenience function for the maximum of the absolute values"""
   if len(vals) == 0: return 0
   return max([abs(v) for v in vals])

def ndigits_lod(num, base=10):
   """return the number of digits to the left of the decimal"""
   anum = abs(num)
   if base == 10: return 1+int(math.log10(anum))
   else:          return 1+int(math.log10(anum)/math.log10(base))

# almost identical to demean, but just return the mean
def mean(vec, ibot=-1, itop=-1):
    """return the vector mean, from index ibot to index itop

        if ibot == -1, ibot will be 0
        if itop == -1, itop will be len-1"""

    if not vec: return 0.0
    if ibot > itop:
        print('** afni_util.mean: ibot (%d) > itop (%d)' % (ibot, itop))
        return 0.0

    vlen = len(vec)

    if ibot < 0: ibot = 0
    if ibot > vlen-1: ibot = vlen-1
    if itop < 0: itop = vlen-1
    if itop > vlen-1: itop = vlen-1

    tot = 0.0
    tot = loc_sum(vec[ibot:itop+1])

    return tot/(itop-ibot+1)


# convert from degrees to chord length
def deg2chordlen(theta, radius=1.0):
   """deg2chord_length(theta, radius=1.0):

      Given theta in degrees (0<=theta<=180) and a radius (>=0), return the
      length of the chord with an arc that subtends central angle theta.
      (So the chord corresponds with an arc of a circle, and the arc subtends
      central angle theta.)
      This might help estimate motion distance due to a rotation.

      For a circle of radius R and a central angle T in degrees, compute the
      resulting chord length (distance between points on the circle for the
      corresponding arc).  If the chord has endpoints A and B, we are looking
      for the length of the segment (AB).

      Note that if a perpendicular is dropped from the circle's center to AB,
      cutting it in half (call the length h), then we have:

        sin(T/2) = h/R, so      h = R * sin(T/2)

      return 2*h (to get the entire chord length)
   """

   # put theta in [0,180]
   if theta < 0.0: theta = abs(theta)
   if theta > 360: theta = theta % 360
   if theta > 180: theta = 180 - theta

   # ignore a negative radius
   if radius <= 0.0: return 0.0

   # math.tan takes input in radians
   theta = theta * math.pi / 180.0

   return 2.0 * radius * math.sin(theta/2.0)

# ----------------------------------------------------------------------
# vector manipulation functions
# ----------------------------------------------------------------------

# almost identical to mean, but subtract the mean instead of returning it
def demean(vec, ibot=-1, itop=-1):
    """demean the vector (in place), from index ibot to index itop

        if ibot == -1, ibot will be 0
        if itop == -1, itop will be len-1
    
       return 0 on success, 1 on error"""

    if not vec: return 0
    if ibot > itop:
        print('** afni_util.demean: ibot (%d) > itop (%d)' % (ibot, itop))
        return 1

    vlen = len(vec)

    if ibot < 0: ibot = 0
    if ibot > vlen-1: ibot = vlen-1
    if itop < 0: itop = vlen-1
    if itop > vlen-1: itop = vlen-1

    # first compute the mean
    tot = 0.0
    for ind in range(ibot,itop+1):
       tot += vec[ind]
    mm = tot/(itop-ibot+1)

    # now subtract it
    for ind in range(ibot,itop+1):
       vec[ind] -= mm

    return vec

def lin_vec_sum(s1, vec1, s2, vec2):
   """return s1*[vec1] + s2*[vec2]
      note: vec2 can be None"""

   if vec2 == None:
      return [s1*vec1[i] for i in range(len(vec1))]

   l1 = len(vec1)
   l2 = len(vec2)
   if l1 != l2:
      print('** LVC: vectors have different lengths (%d, %d)' % (l1, l2))
      return []

   return [s1*vec1[i]+s2*vec2[i] for i in range(l1)]

def proj_onto_vec(v1, v2, unit_v2=0):
   """return vector v1 projected onto v2

      unit_v2: flag denoting whether v2 is a unit vector

      return <v1,v2>/<v2,v2> * v2"""

   if unit_v2: scalar = dotprod(v1,v2)
   else:
      len2 = L2_norm(v2,v2)
      if len2 == 0:
         print('** cannot project onto <0> vector')
         return []
      scalar = dotprod(v1,v2) * 1.0 / len2

   return lin_vec_sum(scalar, v2, 0, None)

def proj_out_vec(v1, v2, unit_v2=0):
   """return vector v1 with v2 projected out
      (return the component of v1 that is orthogonal to v2)

      unit_v2: flag denoting whether v2 is a unit vector

      return v1 - (v1 projected onto v2)

      Note: y - proj(y,p)
          = y - <y,p>/<p,p> * pT        = y - yTp/pTp * pT
          = y - <y,p>/|p|^2 * pT
          = y - <y,p>*(pTp)^-1 * pT     (where (pTp)^-1*pT = pseudo-inverse)
          = (I - p (pTp)^-1 * pT) * y
   """

   return lin_vec_sum(1, v1, -1, proj_onto_vec(v1, v2, unit_v2))

# ----------------------------------------------------------------------
# statistical routines - stdev, variance, ttest
# ----------------------------------------------------------------------

def stat_mean_abs_dev(data):
    """return the mean absolute deviation"""

    if not data: return 0
    length = len(data)
    # length == 1 has MAD 0
    if length <=  1: return 0

    if type(data[0]) == str:
       try: dd = [float(val) for val in data]
       except:
          print('** bad data for min_mean_max_stdev')
          return 0, 0, 0, 0
    else: dd = data

    meanval = loc_sum(dd)/float(length)
    dsum = 0.0
    for val in dd:
       dsum += abs(val-meanval)
    return 1.0*dsum/length

def min_mean_max_stdev(data):
    """return 4 values for data: min, mean, max, stdev (unbiased)"""

    if not data: return 0,0,0,0
    length = len(data)
    if length <  1: return 0,0,0,0

    if type(data[0]) == str:
       try: dd = [float(val) for val in data]
       except:
          print('** bad data for min_mean_max_stdev')
          return 0, 0, 0, 0
    else: dd = data
    if length == 1: return dd[0], dd[0], dd[0], 0.0

    minval  = min(dd)
    maxval  = max(dd)
    meanval = loc_sum(dd)/float(length)

    return minval, meanval, maxval, stdev_ub(dd)

def interval_offsets(times, dur):
    """given a list of times and an interval duration (e.g. TR), return
       the offsets into respective intervals"""

    if not times or dur <= 0:
        print("** interval offsets: bad dur (%s) or times: %s" % (dur, times))
        return []

    length = len(times)
    if length <  1: return []

    fdur = float(dur)   # to make sure (e.g. avoid int division)

    try: offlist = [val % fdur for val in times]
    except:
        print("** interval offsets 2: bad dur (%s) or times: %s" % (dur, times))
        return []
   
    return offlist

def fractional_offsets(times, dur):
    """given a list of times and an interval duration (e.g. TR), return
       the fractional offsets into respective intervals

       i.e. similar to interval offsets, but times are divided by dur"""

    # rely on i_o for error checking
    olist = interval_offsets(times, dur)
    if len(olist) < 1 or dur <= 0: return []

    dur = float(dur)
    for ind, val in enumerate(olist):
        olist[ind] = val/dur

    return olist

def stdev_ub(data):
    """unbiased standard deviation (divide by len-1, not just len)
              stdev_ub = sqrt( (sumsq - N*mean^2)/(N-1) )
    """

    return math.sqrt(variance_ub(data))

def stdev(data):
    """(biased) standard deviation (divide by len, not len-1)
       standard deviation = sqrt(variance)
    """

    return math.sqrt(variance(data))

def variance_ub(data):
    """unbiased variance (divide by len-1, not just len)

       variance = mean squared difference from the mean
                = sum(x-mean)^2 / N

     * unbiased variance
                = sum(x-mean)^2 / (N-1)
                = (sumsq - N*mean^2)/(N-1)
    """

    length = len(data)
    if length <  2: return 0.0

    meanval = loc_sum(data)/float(length)
    # compute variance
    ssq = 0.0
    for val in data: ssq += val*val
    val = (ssq - length*meanval*meanval)/(length-1.0)

    # watch for truncation artifact
    if val < 0.0 : return 0.0
    return val

def variance(data):
    """(biased) variance (divide by len, not len-1)

       variance = mean squared difference from the mean
                = sum(x-mean)^2 / N
                = (sumsq - N*mean^2)/N
    """

    length = len(data)
    if length <  2: return 0.0

    meanval = loc_sum(data)/float(length)
    # compute variance
    ssq = 0.0
    for val in data: ssq += val*val
    val = (ssq - length*meanval*meanval)/length

    # watch for truncation artifact
    if val < 0.0 : return 0.0
    return val

def covary(x, y):
    """return the raw covariance:
       sum[(xi - mean_x)*(yi - mean_y)] / (N-1)
    """

    ll = len(x)
    mx = mean(x)
    my = mean(y)

    vv = loc_sum([(x[i]-mx)*(y[i]-my) for i in range(ll)])

    if ll > 1: return 1.0 * vv / (ll-1.0)
    else:      return 0.0

def r(vA, vB, unbiased=0):
    """return Pearson's correlation coefficient

       for demeaned and unit vectors, r = dot product
       for unbiased correlation, return r*(1 + (1-r^2)/2N)

       note: correlation_p should be faster
    """
    la = len(vA)

    if len(vB) != la:
        print('** r (correlation): vectors have different lengths')
        return 0.0
    ma = mean(vA)
    mb = mean(vB)
    dA = [v-ma for v in vA]
    dB = [v-mb for v in vB]
    sA = L2_norm(dA) # is float
    sB = L2_norm(dB)
    dA = [v/sA for v in dA]
    dB = [v/sB for v in dB]

    r = dotprod(dA,dB)

    if unbiased: return r * (1.0 + (1-r*r)/(2.0*la))
    return r

def linear_fit(vy, vx=None):
   """return slope and intercept for linear fit to data

      if vx is not provided (i.e. no scatterplot fit), then return
      fit to straight line (i.e. apply as vx = 1..N, demeaned)

      slope = N*sum(xy) - (sum(x)*sum(y)]
              ---------------------------
              N*sum(x*x) - (sum(x))^2

      inter = 1/N * (sum(y) - slope*sum(x))

      note: we could compute slope = covary(x,y)/covary(x,x)
   """

   N = len(vy)
   mn = (N-1.0)/2

   # maybe use demeaned, slope 1 line
   if vx == None: vx = [i-mn for i in range(N)]
   else:
      if len(vx) != N:
         print('** cannot fit vectors of different lengths')
         return 0.0, 0.0

   sx   = loc_sum(vx)
   sy   = loc_sum(vy)
   ssx  = dotprod(vx, vx)
   ssxy = dotprod(vy, vx)

   slope = (1.0 * N * ssxy - sx * sy) / (N * ssx - sx*sx )
   inter = 1.0 * (sy - slope * sx) / N

   return slope, inter


def eta2(vA, vB):
    """return eta^2 (eta squared - Cohen, NeuroImage 2008

                        SUM[ (a_i - m_i)^2 + (b_i - m_i)^2 ]
         eta^2 =  1  -  ------------------------------------
                        SUM[ (a_i - M  )^2 + (b_i - M  )^2 ]

         where  a_i and b_i are the vector elements
                m_i = (a_i + b_i)/2
                M = mean across both vectors

    """

    length = len(vA)
    if len(vB) != length:
        print('** eta2: vectors have different lengths')
        return 0.0
    if length < 1: return 0.0

    ma = mean(vA)
    mb = mean(vB)
    gm = 0.5*(ma+mb)

    vmean = [(vA[i]+vB[i])*0.5 for i in range(length)]

    da = [vA[i] - vmean[i] for i in range(length)]
    db = [vB[i] - vmean[i] for i in range(length)]
    num = sumsq(da) + sumsq(db)

    da = [vA[i] - gm       for i in range(length)]
    db = [vB[i] - gm       for i in range(length)]
    denom = sumsq(da) + sumsq(db)

    if num < 0.0 or denom <= 0.0 or num >= denom:
        print('** bad eta2: num = %s, denom = %s' % (num, denom))
        return 0.0
    return 1.0 - float(num)/denom

def correlation_p(vA, vB, demean=1, unbiased=0):
    """return the Pearson correlation between the 2 vectors
       (allow no demean for speed)
    """

    la = len(vA)
    if len(vB) != la:
        print('** correlation_pearson: vectors have different lengths')
        return 0.0

    if la < 2: return 0.0

    if demean:
       ma = mean(vA)
       mb = mean(vB)
       dA = [v-ma for v in vA]
       dB = [v-mb for v in vB]
    else:
       dA = vA
       dB = vB

    sAB = dotprod(dA, dB)
    ssA = sumsq(dA)
    ssB = sumsq(dB)

    if demean: del(dA); del(dB)

    if ssA <= 0.0 or ssB <= 0.0: return 0.0
    else:
       r = sAB/math.sqrt(ssA*ssB)
       if unbiased: return r * (1 + (1-r*r)/(2*la))
       return r

def ttest(data0, data1=None):
    """just a top-level function"""

    if data1: return ttest_2sam(data0, data1)
    return ttest_1sam(data0)

def ttest_1sam(data, m0=0.0):
    """return (mean-m0) / (stdev_ub(data)/sqrt(N)),

              where stdev_ub = sqrt( (sumsq - N*mean^2)/(N-1) )

       or faster, return: (sum-N*m0)/(sqrt(N)*stdev_ub)

       note: move 1/N factor to denominator
    """

    # check for short length
    N = len(data)
    if N < 2: return 0.0

    # check for division by 0
    sd = stdev_ub(data)
    if sd <= 0.0: return 0.0

    # and return, based on any passed expected mean
    if m0: t = (loc_sum(data) - N*m0)/(math.sqrt(N)*sd)
    else:  t =  loc_sum(data)        /(math.sqrt(N)*sd)

    return t

def ttest_paired(data0, data1):
    """easy: return 1 sample t-test of the difference"""

    N0 = len(data0)
    N1 = len(data1)
    if N0 < 2 or N1 < 2: return 0.0
    if N0 != N1:
        print('** ttest_paired: unequal vector lengths')
        return 0.0

    return ttest_1sam([data1[i] - data0[i] for i in range(N0)])

def ttest_2sam(data0, data1, pooled=1):
    """if not pooled, return ttest_2sam_unpooled(), otherwise

       return (mean1-mean0)/sqrt(PV * (1/N0 + 1/N1))

              where PV (pooled_variance) = ((N0-1)*V0 + (N1-1)*V1)/(N0+N1-2)

       note: these lists do not need to be of the same length
       note: the sign is as with 3dttest (second value(s) minus first)
    """

    if not pooled: return ttest_2sam_unpooled(data0, data1)

    N0 = len(data0)
    N1 = len(data1)
    if N0 < 2 or N1 < 2: return 0.0

    m0 = loc_sum(data0)/float(N0)
    v0 = variance_ub(data0)

    m1 = loc_sum(data1)/float(N1)
    v1 = variance_ub(data1)

    pv = ((N0-1)*v0 + (N1-1)*v1) / (N0+N1-2.0)
    if pv <= 0.0: return 0.0

    return (m1-m0)/math.sqrt(pv * (1.0/N0 + 1.0/N1))

def ttest_2sam_unpooled(data0, data1):
    """return (mean1-mean0)/sqrt(var0/N0 + var1/N1)

       note: these lists do not need to be of the same length
       note: the sign is as with 3dttest (second value(s) minus first)
    """

    N0 = len(data0)
    N1 = len(data1)
    if N0 < 2 or N1 < 2: return 0.0

    m0 = loc_sum(data0)/float(N0)
    v0 = variance_ub(data0)

    m1 = loc_sum(data1)/float(N1)
    v1 = variance_ub(data1)

    if v0 <= 0.0 or v1 <= 0.0: return 0.0

    return (m1-m0)/math.sqrt(v0/N0 + v1/N1)


def p2q(plist, do_min=1, verb=1):
    """convert list of p-value to a list of q-value, where
         q_i = minimum (for m >= i) of N * p_m / m
       if do min is not set, simply compute q-i = N*p_i/i

       return q-values in increasing significance
              (i.e. as p goes large to small, or gets more significant)
    """

    q = plist[:]
    q.sort()
    N = len(q)

    # work from index N down to 0 (so index using i-1)
    min = 1
    for i in range(N,0,-1):
       ind = i-1
       q[ind] = N * q[ind] / i
       if do_min:
          if q[ind] < min: min = q[ind]
          if min < q[ind]: q[ind] = min

    # and flip results
    q.reverse()

    return q

def argmax(vlist, absval=0):
   if len(vlist) < 2: return 0
   if absval: vcopy = [abs(val) for val in vlist]
   else:      vcopy = vlist

   mval = vcopy[0]
   mind = 0
   for ind, val in enumerate(vlist):
      if val > mval:
         mval = val
         mind = ind

   return mind

def argmin(vlist, absval=0):
   if len(vlist) < 2: return 0
   if absval: vcopy = [abs(val) for val in vlist]
   else:      vcopy = vlist

   mval = vcopy[0]
   mind = 0
   for ind, val in enumerate(vlist):
      if val < mval:
         mval = val
         mind = ind

   return mind

def gaussian_at_hwhm_frac(frac):
   """gaussian_at_hwhm_frac(frac):

      return the magnitude of a (unit, 0-centered) Gaussian curve at a
      fractional offset of the HWHM (half width at half maximum) = FWHM/2

         return h(f) = 2^[-f^2]

      HWHM is a logical input, as it is a radius, while FWHM is a diameter.

      So gaussian_at_hwhm_frac(1) = 0.5, by definition.

        - the return value should be in (0,1], and == 1 @ 0
        - the return value should be 0.5 @ 0.5 (half max @ FWHM radius)
          (i.e. since FWHM is a diameter, FWHM/2 is the radius)
        - if frac < 0, whine and return 0 (it is undefined)

      Gaussian curves have the form: G(x) = a*e^-[ (x-b)^2 / (2*c^2) ],
      where     a = scalar, maybe 1/(c*sqrt(2*PI)), to be unit integral
                b = expected value, the central point of symmetry
                c = standard deviation

      A unit height, zero-centered curve has a=1, b=0: g(x)=e^-[x^2 / (2*c^2)]

      To find (the HWHM) where g(x) = 1/2, solve: g(w) = 1/2 for w.

         w = sqrt(c^2 * 2*ln(2))    {just use positive}

      We want an equation for g(x), but where x is a fraction of the HWHM.
      Rewrite g(x) in terms of w, by solving the above for c:

        c = w / sqrt(2 * ln2)

      and substitute back into g(x).  g(x) = e^-[x^2 * ln2 / w^2]

      and finally write in terms of f, where f is a fraction of w.

      Define h(f) = g(fw) = e^-[f^2*w^2 * ln2 / w^2]  =  e^-[f^2 * ln2]

      Now we can cancel the 2 and ln:

        h(f) = e^-ln[2^(f^2)] = e^ln[2^(-f^2)] = 2^(-f^2)

      return h(f) = 2^(-f^2)
   """
   if frac < 0:
      print("** gaussian_at_hwhm_frac: illegal frac < 0 of %s", frac)
      return 0

   return 2.0 ** -(frac*frac)

def gaussian_at_fwhm(x, fwhm):
   """gaussian_at_fwhm(x, fwhm):

      return the magnitude of unit, zero-centered Gaussian curve at the given x

      The Gaussian curve is defined by the given FWHM value, e.g.

         g(x) = e^-[x^2 * ln2 / FWHM^2]

      This actually returns gaussian_at_hwhm_frac(x/(fwhm/2)), or of (2x/fwhm).
   """
   if fwhm <= 0:
      print("** gaussian_at_fwhm: illegal fwhm <= 0 of %s", fwhm)
      return 0

   return gaussian_at_hwhm_frac(2.0*x/fwhm)

def convolve(data, kernel, length=0):
   """simple convolution of data with a kernel

      data      : array of values to convolve with kernel
      kernel    : convolution kernel (usually shorter)
      length    : if > 0: defines output length (else len(data))

      return convolved array
   """
   klen = len(kernel)
   if length == 0: rlen = len(data)
   else:           rlen = length

   if len(data) == 0 or klen == 0:
      return []

   res = [0]*rlen
   for dind, dval in enumerate(data):
      if dind >= rlen:
         break

      for kind, kval in enumerate(kernel):
         if dind+kind >= rlen:
            break
         res[dind+kind] += dval * kval

   return res

# ----------------------------------------------------------------------
# random list routines: shuffle, merge, swap, extreme checking
# ----------------------------------------------------------------------

def swap_2(vlist, i1, i2):
    if i1 != i2:
       val = vlist[i2]
       vlist[i2] = vlist[i1]
       vlist[i1] = val

def shuffle(vlist, start=0, end=-1):
    """randomize the order of list elements, where each permutation is
       equally likely

       - akin to RSFgen, but do it with equal probabilities
         (search for swap in [index,N-1], not in [0,N-1])
       - random.shuffle() cannot produce all possibilities, don't use it
       - start and end are indices to work with
    """

    # if we need random elsewhere, maybe do it globally
    import random

    vlen = len(vlist)

    # check bounds
    if start < 0 or start >= vlen: return
    if end >= 0 and end <= start:  return

    # so start is valid and end is either < 0 or big enough
    if end < 0 or end >= vlen: end = vlen-1

    nvals = end-start+1

    # for each index, swap with random other towards end
    for index in range(nvals-1):
        rind = int((nvals-index)*random.random())
        swap_2(vlist, start+index, start+index+rind)

    # return list reference, though usually ignored
    return vlist

def shuffle_blocks(vlist, bsize=-1):
    """like shuffle, but in blocks of given size"""

    vlen = len(vlist)

    if bsize < 0 or bsize >= vlen:
       shuffle(vlist)
       return

    if bsize < 2: return

    nblocks = vlen // bsize
    nrem    = vlen  % bsize

    boff = 0
    for ind in range(nblocks):
       shuffle(vlist, boff, boff+bsize-1)
       boff += bsize
    shuffle(vlist, boff, boff+nrem-1)
        
    return vlist

def random_merge(list1, list2):
    """randomly merge 2 lists (so each list stays in order)

       shuffle a list of 0s and 1s and then fill from lists
    """

    # if we need random elsewhere, maybe do it globally
    import random

    mlist = [0 for i in range(len(list1))] + [1 for i in range(len(list2))]
    shuffle(mlist)

    i1, i2 = 0, 0
    for ind in range(len(mlist)):
        if mlist[ind] == 0:
            mlist[ind] = list1[i1]
            i1 += 1
        else:
            mlist[ind] = list2[i2]
            i2 += 1

    return mlist

def show_sum_pswr(nT, nR):
    cp = 0.0
    prev = 0
    for r in range(nR+1):
       # already float, but be clear
       p = float(prob_start_with_R(nT,nR,r))
       cp += p
       # print 'prob at %3d = %g (cum %g)' % (r, p, cp)
       if prev == 0: prev = p
       print(p, p/prev)
       prev = p
    print('cum result is %g' % cp)


def prob_start_with_R(nA, nB, nS):
    """return the probability of starting nS (out of nB) class B elements
       should equal: choose(nB, nS)*nS! * nA *(nB+nA-nS-1)! / (nA+nB)!
       or: factorial(nB, init=nB-nS+1) * nA / fact(nA+nB, init=nA+nB-nS)

       or: choose(nB,nS)/choose(nA+nB,nS) * nA/(nA+nB-nS)
       
    """
    return 1.0 * nA * factorial(nB,    init=nB-nS+1) \
                    / factorial(nA+nB, init=nA+nB-nS)

def choose(n,m):
    """return n choose m = n! / (m! * (n-m)!)"""
    # integral division (or use floats, to better handle overflow)
    return factorial(n,init=n-m+1) / factorial(m)

def factorial(n, init=1):
    prod = 1
    val = init
    while val <= n:
       prod *= val
       val += 1
    return prod

def swap2(data):
    """swap data elements in pairs"""
    
    size  = 2
    nsets = len(data)//size
    if nsets <= 0: return

    for ind in range(nsets):
        off = ind*size
        v           = data[off]     # swap d[0] and d[1]
        data[off]   = data[off+1]
        data[off+1] = v

def swap4(data):
    """swap data elements in groups of 4"""
    
    size  = 4
    nsets = len(data)//size
    if nsets <= 0: return

    for ind in range(nsets):
        off = ind*size
        v           = data[off]     # swap d[0] and d[3]
        data[off]   = data[off+3]
        data[off+3] = v
        v           = data[off+1]   # swap d[1] and d[2]
        data[off+1] = data[off+2]
        data[off+2] = v

def vec_extremes(vec, minv, maxv, inclusive=0):
   """return a integer array where values outside bounds are 1, else 0

      if inclusive, values will also be set if they equal a bound

      return error code, new list
             success: 0, list
             error  : 1, None"""

   if not vec: return 1, None

   if minv > maxv:
      print('** extremes: minv > maxv (', minv, maxv, ')') 
      return 1, None

   if inclusive:
      elist = [1*(vec[t]>=maxv or vec[t]<=minv) for t in range(len(vec))]
   else:
      elist = [1*(vec[t]> maxv or vec[t]< minv) for t in range(len(vec))]

   return 0, elist

def vec_moderates(vec, minv, maxv, inclusive=1):
   """return a integer array where values inside bounds are 1, else 0

      if inclusive, values will also be set if they equal a bound

      return error code, new list
             success: 0, list
             error  : 1, None"""

   if not vec: return 1, None

   if minv > maxv:
      print('** moderates: minv > maxv (', minv, maxv, ')') 
      return 1, None

   if inclusive:
      elist = [1*(vec[t]>=minv and vec[t]<=maxv) for t in range(len(vec))]
   else:
      elist = [1*(vec[t]> minv and vec[t]< maxv) for t in range(len(vec))]

   return 0, elist

def vec_range_limit(vec, minv, maxv):
   """restrict the values to [minv, maxv]

      This function modifies the past vector.

      return 0 on success, 1 on error"""

   if not vec: return 0

   if minv > maxv:
      print('** range_limit: minv > maxv (', minv, maxv, ')')
      return 1

   for ind in range(len(vec)):
      if   vec[ind] < minv: vec[ind] = minv
      elif vec[ind] > maxv: vec[ind] = maxv

   return 0

# for now, make 2 vectors and return their correlation
def test_polort_const(ntrs, nruns, verb=1):
    """compute the correlation between baseline regressors of length ntrs*nruns
       - make vectors of 11...10000...0 and 00...011...100..0 that are as the
         constant polort terms of the first 2 runs
       - return their correlation

       - note that the correlation is easily provable as -1/(N-1) for N runs
    """

    if ntrs <= 0 or nruns <= 2: return -1  # flag

    # lazy way to make vectors
    v0 = [1] * ntrs + [0] * ntrs + [0] * (ntrs * (nruns-2))
    v1 = [0] * ntrs + [1] * ntrs + [0] * (ntrs * (nruns-2))

    if verb > 1:
        print('++ test_polort_const, vectors are:\n' \
              '   v0 : %s \n'                        \
              '   v1 : %s' % (v0, v1))

    return correlation_p(v0, v1)

# for now, make 2 vectors and return their correlation
def test_tent_vecs(val, freq, length):
    a = []
    b = []
    for i in range(length):
        if (i%freq) == 0:
            a.append(val)
            b.append(1-val)
        elif ((i-1)%freq) == 0:
            a.append(0.0)
            b.append(val)
        else:
            a.append(0.0)
            b.append(0.0)

    return correlation_p(a,b)


# -----------------------------------------------------------------------
# [PT: June 8, 2020] for matching str entries in list (for the FATCAT
# -> MVM and other group analysis programs)
## [PT: June 8, 2020] updated to allow cases where only a subset of
## either A or B is matched; but will still give errors if something
## matches more than 1 item in the other list.

def match_listA_str_in_listB_str(A, B):
    """Input: two lists (A and B), each of whose elements are strings.  A
    and B don't have to have the same length.

    See if each string in A is contained is contained within one (and
    only one) string in list B.  If yes, return:
      1 
      the dictionary of matches, A->B
      the dictionary of matches, B->A
    elif only a subset of A is contained in B or vice versa, return:
      0
      the dictionary of matches, A->B
      the dictionary of matches, B->A
    otherwise, error+exit.

    The primary/first usage of this program is for the following case:
    matching subject IDs from a CSV file (the list of which form A)
    with path names for *.grid|*.netcc files (the list of which form
    B), for combining all that information.  We want to find 1 subject
    ID from the CSV file with 1 (and only 1) matrix file, for each
    subject ID.

    NB: in the above case, the users are alternatively able to enter
    the list for matching subj ID with matrix file, if the above
    name-matching wouldn't work.

    """

    if type(A) != list or type(B) != list :
        BASE.EP("Both inputs A and B must be lists, not {} and {}, "
              "respectively".format(type(A), type(B)))

    na = len(A)
    nb = len(B)

    if not(na and nb) :
        BASE.EP("One of the sets is empty: len(A) = {}; len(B) = {}"
              "".format(na, nb))
 
    # check that each ele is a str
    ta = [type(x)!=str for x in A]
    tb = [type(x)!=str for x in B]
    if max(ta) or max(tb) :
        BASE.EP("All elements of A and B must be str, but that is not true\n"
              "for a least one list:  for A, it's {};  for B it's {}"
              "".format(not(max(ta)), not(max(tb))))

    checklist_a = [0]*na
    checklist_b = [0]*nb
    matchlist_a = [-1]*na
    matchlist_b = [-1]*nb

    for ii in range(nb):
        for jj in range(na):
            if B[ii].__contains__(A[jj]) :
                checklist_a[jj]+= 1
                checklist_b[ii]+= 1
                matchlist_a[jj] = ii
                matchlist_b[ii] = jj
                # *should* be able to break here, but will check
                # *all*, to verify there are no problems/ambiguities

    CHECK_GOOD = True   # determines if we 'return' anything
    FULL_MATCH = 1      # flags type of matching when returning

    # --------- now check the outcomes
    if min(checklist_a) == 1 and max(checklist_a) == 1 :
        # all matches found, all singletons
        BASE.IP("Found single matches for each element in A")
    else:
        FULL_MATCH = 0
        if min(checklist_a) == 0 :
            # all found matches are singletons, but there are gaps;
            # not a fatal error
            BASE.WP("Some elements of A are unmatched:")
            for ii in range(na):
                if not(checklist_a[ii]) :
                    print("\t unmatched: {}".format(A[ii]))
        if max(checklist_a) > 1 :
            # ambiguities in matching --> badness
            CHECK_GOOD = False
            BASE.WP("Some elements of A are overmatched:")
            for ii in range(na):
                if checklist_a[ii] > 1 :
                    print("\t overmatched: {}".format(A[ii]))
        if not(max(checklist_a)) :
            CHECK_GOOD = False
            BASE.WP("No elements in A matched??:")

    if min(checklist_b) == 1 and max(checklist_b) == 1 :
        # all matches found, all singletons
        BASE.IP("Found single matches for each element in B")
    else:
        FULL_MATCH = 0
        if min(checklist_b) == 0 :
            # all found matches are singletons, but there are gaps;
            # not a fatal error
            BASE.WP("Some elements of B are unmatched:")
            for ii in range(nb):
                if not(checklist_b[ii]) :
                    print("\t unmatched: {}".format(B[ii]))
        if max(checklist_b) > 1 :
            # ambiguities in matching --> badness
            CHECK_GOOD = False
            BASE.WP("Some elements of B are overmatched:")
            for ii in range(nb):
                if checklist_b[ii] > 1 :
                    print("\t overmatched: {}".format(B[ii]))
        if not(max(checklist_b)) :
            CHECK_GOOD = False
            BASE.WP("No elements in B matched??:")

    if not(CHECK_GOOD) :
        BASE.EP('Exiting')

    # if we made it here, things are good
    da = {}
    for ii in range(na):
        if checklist_a[ii] :
            da[ii] = matchlist_a[ii]
    db = {}
    for ii in range(nb):
        if checklist_b[ii] :
            db[ii] = matchlist_b[ii]

    na_keys = len(da.keys())
    nb_keys = len(db.keys())

    # final check on consistency
    if na_keys != nb_keys :
        BASE.EP("Mismatch in number of keys in output lists?\n"
                "There are {} for A and {} for B"
                "".format(na_keys, nb_keys))

    # e.g., 
    #for key in da: 
    #    print(A[key] + '  <---> ' +  B[da[key]]) 
    #for key in db: 
    #    print(A[db[key]] + '  <---> ' +  B[key]) 

    return FULL_MATCH, da, db

def invert_dict(D):
    """Take a dictionary D of key-value pairs and return a dictionary
    Dinv, which the keys of D as values and the (matched) values of D
    as keys.

    """

    if not(D) :      return {}

    Dinv = {}
    for k, v in D.items():     
        Dinv[v] = k

    return Dinv


# ----------------------------------------------------------------------
# [PT: Jan 13, 2020] Pieces for dealing with files of AFNI seed points

class afni_seeds:

    def __init__(self, file_line=[]):
        '''Make an object for seed-based correlation.

        Use the optional 'file_line' input to read a (non-split) line
        from a file easily.

        '''

        self.xyz        = []
        self.space      = ''
        self.roi_label  = ''
        self.netw       = None

        # these for future/other use;  cannot be changed at the moment
        self.rad    = None
        self.shape  = 'sphere'

        # use this to read a (non-split) line from a file easily
        if file_line :
            self.set_all_info_from_file(file_line)

    def set_space(self, x):
        self.space = x

    def set_xyz(self, xlist):
        '''Input a list of strings or floats to be a list of float coords.
        '''
        if len(xlist) != 3:
            print("** ERROR: need 3 items input, not {}:\n"
                  "      '{}'".format(len(xlist), xlist))
            sys.exit(3)
        self.xyz = [float(x) for x in xlist]

    def set_roi_label(self, x):
        self.roi_label = x

    def set_netw(self, x):
        self.netw = x

    def set_rad(self, x):
        self.rad = float(x)

    def set_shape(self, x):
        self.shape = x

    def set_all_info_from_file(self, line):
        '''Read in an unsplit line from a file, and save all the information
        appropriately.  Column order matters!

        '''
        y    = line.split()
        Ny   = len(y)
        if Ny < 5 :
            print("** ERROR: too few columns ({}) in this line:\n"
                  "   {}\n"
                  "   e.g., did you forget something?".format(Ny, line))
            sys.exit(2)

        elif Ny > 6 :
            print("** ERROR: too many columns ({}) in this line:\n"
                  "   {}\n"
                  "   e.g., did you put spaces in names?".format(Ny, line))
            sys.exit(2)

        xyz  = [ float(y[j]) for j in range(3) ]

        self.set_xyz(xyz)
        self.set_space(y[3])
        self.set_roi_label(y[4])

        if len(y) == 6 :
            self.set_netw(y[5])

# [PT: Jan 13, 2020] for APQC (and maybe other purposes later), read
# in a simple text file that contains both numbers and strings.
# See text file:  afni_seeds_per_space.txt
def read_afni_seed_file(fname, only_from_space=None):
    '''Read in a file of columns of info, and return list of seed objs.

    Input
    -----
    fname           : (req) text file of seed point information.
    only_from_space : (opt) can choose to only attach seeds from a given space
                      (which is entered here); otherwise, all seeds are output.
    
    Output: list of seed objects.

    '''

    # expected/allowed cols ('netw' is opt)
    cols  = ['x', 'y', 'z', 'space', 'roi_label', 'netw']
    Ncols = len(cols)

    # input data, ignore blank lines
    x   = read_text_file( fname, noblank=1 )
    Nx  = len(x)

    dat = []

    for i in range(Nx):
        row = x[i]
        if row[0] != '#': 
            seed_obj = afni_seeds(file_line=row)
            if not(only_from_space) :
                dat.append( seed_obj )
            else:
                if seed_obj.space == only_from_space :
                    dat.append( seed_obj )

    return dat

# [PT: June 5, 2023] For APQC HTML generation, and likely other
# things.
def rename_label_safely(x, only_hash=False):
    """Make safe string labels that can be used in filenames (so no '#')
and NiiVue object names (so no '-', '+', etc.; sigh).  

For example, 'vstat_V-A_GLT#0_Coef' -> 'vstat_V__A_GLT_0_Coef'.

The mapping rules are in the text of this function.  This function
might (likely) update over time.

Parameters
----------
x : str
    a name
only_hash : bool
    simpler renaming, only replacing the hash symbol (that appears
    in stat outputs)

Returns
-------
y : str
    a name that has (hopefully) been made safe by various letter 
    substitutions.

    """

    y = x.replace('#', '_')

    if only_hash :
        return y

    y = y.replace('-', '__')
    y = y.replace('+', '___')
    y = y.replace('.', '____')

    return y

# ----------------------------------------------------------------------
# [PT: 2024-05-08] when output is a new dir, this is useful to control
# behavior (happens in APQC generation, and now more)

# control overwriting/backing up any existing dirs
dict_ow_modes = {
    'simple_ok'   : 'make new dir, ok if pre-exist (mkdir -p ..)',
    'shy'         : 'make new dir only if one does not exist',
    'overwrite'   : 'purge old dir and make new dir in its vacant place',
    'backup'      : 'move existing dir to dir_<time>; then make new dir',
}

list_ow_modes = list(dict_ow_modes.keys())
list_ow_modes.sort()
str_ow_modes  = ', '.join([x for x in list_ow_modes])
hstr_ow_modes = '\n'.join(['{:12s} -> {}'.format(x, dict_ow_modes[x]) \
                           for x in list_ow_modes])

def is_valid_ow_mode(ow_mode):
    """Simple check about whether input ow_mode is a valid one. Return
True if valid, and False otherwise."""

    is_valid = ow_mode in list_ow_modes

    return is_valid

def make_new_odir(new_dir, ow_mode='backup', bup_dir=''):
    """When outputting to a new directory new_dir, just create it if it
doesn't exist already; but if a dir *does* pre-exist there, then do
one of the following behaviors, based on keyword values of ow_mode
(and with consideration of bup_dir value):
  'simple_ok'   : make new dir, ok if pre-exist (mkdir -p ..)
  'overwrite'   : remove pre-existing new_dir and create empty one in
                  its place
  'backup' and bup_dir != '' : move that pre-existing dir to bup_dir
  'backup' and bup_dir == '' : move that pre-existing dir to new_dir_<TIMESTAMP>
  'shy'         : make new_dir only if one does not pre-exist.

Parameters
----------
new_dir : str
    name of new dir to make
ow_mode : str
    label for the overwrite mode behavior of replacing or backing up
    an existing new_dir (or a file with the same name as the dir)
bup_dir : str
    possible name for backup directory    

Returns
----------
num : int
    return 0 up on success, or a different int if failure

    """

    do_cap = True

    # valid ow_mode?
    if not( is_valid_ow_mode(ow_mode) ) :
        print("** ERROR: illegal ow_mode '{}', not in the list:\n"
              "   {}".format(ow_mode, str_ow_modes))
        sys.exit(11)

    # check if the main QC dir exists already
    DIR_EXISTS = os.path.exists(new_dir)

    if DIR_EXISTS :
        if ow_mode=='shy' or ow_mode==None :
            print("** ERROR: output dir exists already: {}\n"
                  "   Exiting.".format(new_dir))
            print("   Check out '-ow_mode ..' for overwrite/backup opts.")
            sys.exit(10)
        
        elif ow_mode=='backup' :
            if not(bup_dir) :
                # make our own backup dir with timestamp
                now     = datetime.now()          # current date and time
                bup_dir = now.strftime( new_dir + "_%Y-%m-%d-%H-%M-%S")

            print("+* WARN: output dir exists already: {}\n"
                  "   -> backing it up to: {}".format(new_dir, bup_dir))
            cmd  = '''\\mv {} {}'''.format(new_dir, bup_dir)
            com  = BASE.shell_com(cmd, capture=do_cap)
            stat = com.run()

        elif ow_mode=='overwrite' :
            print("+* WARN: output dir exists already: {}\n"
                  "   -> overwriting it".format(new_dir))
            cmd    = '''\\rm -rf {}'''.format(new_dir)
            com    = BASE.shell_com(cmd, capture=do_cap)
            stat   = com.run()

        elif ow_mode=='simple_ok' :
            # just leads to essentially doing 'mkdir -p ..' with the new_dir
            print("++ OK, output dir exists already: {}".format(new_dir))

    # Now make the new output dir
    cmd    = '''\\mkdir -p {}'''.format(new_dir)
    com    = BASE.shell_com(cmd, capture=do_cap)
    stat   = com.run()

    return 0

# ----------------------------------------------------------------------

if __name__ == '__main__':
   print('afni_util.py: not intended as a main program')
   print('              (consider afni_python_wrapper.py)')
   sys.exit(1)

