#!/usr/bin/env python

# python3 status: started

# Timing library, accessed via timing_tool.py.
#
# - offers class AfniTiming
#
# R Reynolds    December, 2008

import sys
from afnipy import module_test_lib
g_testlibs = ['sys', 'math', 'copy']
if module_test_lib.num_import_failures(g_testlibs): sys.exit(1)
   
# import libraries
import math
import copy
from afnipy import afni_util as UTIL
from afnipy import lib_afni1D as LD

g_marry_AM_methods = ['lin_run_fraq', 'lin_event_index']
g_tsv_def_labels   = ['onset', 'duration', 'trial_type']

# ----------------------------------------------------------------------
# AfniTiming class - for stim times that are married with
#                    either time intervals or magnitudes

class AfniTiming(LD.AfniData):
   def __init__(self, filename="", dur=-1, mdata=None, fsl_flist=[], verb=1):
      """AFNI married stimulus timing class"""

      super(AfniTiming, self).__init__(filename, mdata=mdata, 
                                       fsl_flist=fsl_flist, verb=verb)

      # initialize duration data
      if self.ready: self.init_durations(dur)

      if self.verb > 3: self.show()

   # old accessor functions from AfniTiming
   def add_rows(self, newdata):
      """add newdata.nrows of data (additional runs), return 0 on success"""
      if not self.ready or not newdata.ready:
         print('** AMTiming elements not ready for adding rows (%d,%d)' % \
               (self.ready, newdata.ready))
         return 1

      if self.verb > 1: print('++ Timing: adding %d rows' % newdata.nrows)

      if self.mtype != newdata.mtype:
         print('** add rows: mismatch mtypes (%d vs. %d)' \
               % (self.mtype, newdata.mtype))

      # get data and mdata
      edata = copy.deepcopy(newdata.data)
      self.data.extend(edata)

      edata = copy.deepcopy(newdata.mdata)
      self.mdata.extend(edata)

      self.nrows += newdata.nrows

      if self.row_lens and newdata.row_lens:
         self.row_lens.extend(newdata.row_lens)

      return 0

   def init_durations(self, dur=-1):
      """initialize ddata as format [start_time, end_time]
      """

      # start with general format
      if not self.ready: return

      # if dur is passed, first apply it
      if dur >= 0:
         for row in self.mdata:
            for entry in row: entry[2] = dur

      self.ddata = []
      for row in self.mdata:
         self.ddata.append([[val[0],val[0]+val[2]] for val in row])

      if   dur  > 0: self.mtype |= LD.MTYPE_DUR # set bit
      elif dur == 0: self.mtype &= LD.MTYPE_AMP # only allow AMP to go through
      # else leave as is

      self.dur_len = dur

   def set_dur_len(self):
      self.dur_len = self.get_duration()

   def extend_rows(self, newdata):
      """extend each row by the corresponding row of brows"""

      if self.extend_data_rows(newdata): return 1

      # if interval lengths are not equal, must clear it
      if self.dur_len != newdata.dur_len: self.dur_len = -1

      return 0

   def get_start_end_timing(self, sort=0):
      """return a 2D array of [start, finish] times (so 3D object)"""
      sftimes = []
      for row in self.mdata:
         # if married, include modulators
         if self.married:
            times = [[e[0],e[0]+e[2], e[1]] for e in row]
         else:
            times = [[e[0],e[0]+e[2]] for e in row]
         if sort: times.sort()
         sftimes.append(times)
      return sftimes

   def partition(self, part_file, prefix):
      """partition timing based on part_file labels, write each part to
         prefix_label.1D"""

      if not self.ready:
         print('** Timing element not ready for partitioning')
         return 1

      labels = read_value_file(part_file)
      if labels == None:
         print("** failed to read partition label file '%s'" % part_file)
         return 1

      nlabr = len(labels)

      if self.verb > 3: print('-- partition: %d labels: %s' % (nlabr, labels))

      # first test nrows, then test lengths per row
      if self.nrows != nlabr:
         print('** Timing nrows differs from partition nrows (%d, %d)' % \
               (self.nrows,nlabr))
         return 1
      for ind in range(self.nrows):
         if len(self.data[ind]) != len(labels[ind]):
            print("** timing and label row lengths differ at line %d"%(ind+1))
            print("   (%d != %d)" % (len(self.data[ind]), len(labels[ind])))
            return 1

      # make unique label list
      ulabs = []
      for line in labels:
         ulabs.extend(line)
      ulabs = UTIL.get_unique_sublist(ulabs)
      if self.verb > 2: print('-- partition: unique label list: %s' % ulabs)
      if ulabs.count('0'): ulabs.remove('0')

      if self.verb > 1: print('++ Timing: partitioning with %s' % part_file)

      # ------------------------------------------------------------
      # do the work, starting with copy:
      # for each label, extract those times and write out as timing file
      dupe = self.copy()        # keep results in new class instance
      for lab in ulabs:
         # extract timing for this label 'lab'
         mdata = []
         for r in range(nlabr):
            drow = []           # make one row of times for this label
            for c in range(len(labels[r])):
               if labels[r][c] == lab: drow.append(self.mdata[r][c])
            mdata.append(drow)  # and append the new row
         del(dupe.mdata)        # out with the old,
         dupe.mdata = mdata     # and in with the new
         dupe.write_times('%s_%s.1D' % (prefix,lab))    # and write, yay

      del(dupe)                 # nuke the temporary instance

      return 0

   def global_to_local(self, run_len):
      """convert global times to local, based in run_len array
         return 0 on success, 1 on any error"""

      if not self.ready:
         print('** global timing not ready')
         return 1

      if len(run_len) == 0:
         print('** global_to_local requires -run_len')
         return 1

      if not self.is_rect():
         print('** global timing is not rectangular')
         return 1

      rlen = len(self.mdata[0])         # note row length

      if rlen > 1:
         print('** global timing is not a single column')
         return 1

      if rlen < 1: return 0             # nothing to do
      if self.nrows < 2: return 0       # nothing to do

      # make just one row and sort
      self.transpose()
      self.sort()

      if self.verb > 2:
         self.show('global_to_local')
         print('-- run lengths : %s' % run_len)

      if self.verb > 4:
         print('-- global start time matrix %s' % self.data)

      # now walk through runs and insert times as we go
      newdata = []
      newmdata = []
      endtime = 0.0
      sind = 0
      stimes = self.mdata[0]
      ntimes = len(stimes)
      for etime in run_len:
         starttime = endtime
         endtime += etime
         if sind >= ntimes:  # only empty runs left
            newdata.append([])
            newmdata.append([])
            continue
         # times are left, decide which go for this run
         last = sind
         while last < ntimes and stimes[last][0] < endtime: last += 1
         mtimes = stimes[sind:last]
         for tentry in mtimes: tentry[0] -= starttime
         newdata.append([t[0] for t in mtimes])
         newmdata.append(mtimes)
         sind = last

      # insert any remaining times at end of last run(and warn user)
      if sind < ntimes:
         if self.verb > 0:
            print('** global to local: %d times after last run' % (ntimes-sind))
         mtimes = stimes[sind:]
         for tentry in mtimes: tentry[0] -= starttime
         newdata[-1].extend([t[0] for t in mtimes])
         newmdata[-1].extend(mtimes)

      del(self.data)
      del(self.mdata)
      self.data = newdata
      self.mdata = newmdata
      self.nrows = len(self.mdata)

      if self.verb > 4:
         print('-- local time matrix %s' % newdata)

      return 0

   def local_to_global(self, run_len):
      """convert local times to global times, based in run_len array
         return 0 on success, 1 on any error"""

      if not self.ready:
         print('** local timing not ready')
         return 1

      if len(run_len) != self.nrows:
         print('** local_to_global: have %d run times but %d data rows' \
               % (len(run_len), self.nrows))
         return 1

      # make sure it is sorted for this
      self.sort()

      if self.verb > 2:
         self.show('local_to_global')
         print('-- run lengths : %s' % run_len)

      if self.verb > 4:
         print('-- local time matrix %s' % self.data)

      # now walk through runs and insert times as we go
      newdata = []
      newmdata = []
      runstart = 0.0
      for rind, rtime in enumerate(run_len):
         # each new time is a new row
         mrow = copy.deepcopy(self.mdata[rind])
         for ind, mtime in enumerate(mrow):
            mtime[0] += runstart
            newdata.append([mtime[0]])
            newmdata.append([mtime])
         runstart += rtime      # last one is useless

      del(self.data)
      self.data = newdata
      del(self.mdata)
      self.mdata = newmdata
      self.nrows = len(self.data)

      if self.verb > 4:
         print('-- global time matrix %s' % newdata)

      return 0

   def add_val(self, val):
      """add the given value to each element"""
      if not self.ready: return 1

      if type(val) == str:
         try: val = float(val)
         except:
            print("** invalid value to add to timing: '%s'" % val)
            return 1

      if self.verb > 1: print('-- Timing: adding %f to data...' % val)

      for row in self.data:
         for ind in range(len(row)):
            row[ind] += val

      for row in self.mdata:
         for ind in range(len(row)):
            row[ind][0] += val

      return 0

   def shift_to_offset(self, offset=0):
      """shift all run times to start at the given offset
         (offset should be <= first values)"""

      if not self.ready: return 1

      if type(offset) == str:
         try: offset = float(offset)
         except:
            print("** invalid offset to add to timing: '%s'" % offset)
            return 1

      if self.verb > 1: print('-- timing: setting offset to %f ...' % offset)

      # make sure it is sorted for this
      self.sort()

      del(self.data)
      self.data = []

      for rind, row in enumerate(self.mdata):
         if len(row) < 1:
            self.data.append([])
            continue
         diff = row[0][0] - offset
         if diff < 0:
            print('** offset shift to %f too big for run %d' % (offset, rind))
            return 1
         for ind in range(len(row)):
            row[ind][0] -= diff
         self.data.append([e[0] for e in row])

      return 0

   def scale_val(self, val):
      """multiply the given value into each element"""
      if not self.ready: return 1

      if type(val) == type('hi'):
         try: val = float(val)
         except:
            print("** invalid value to scale into timing: '%s'" % val)
            return 1

      if self.verb > 1: print('-- Timing: scaling data by %f ...' % val)

      for row in self.data:
         for ind in range(len(row)):
            row[ind] *= val

      for row in self.mdata:
         for ind in range(len(row)):
            row[ind][0] *= val

      return 0

   def round_times(self, tr, round_frac=1.0):
      """round/truncate times to multiples of the TR

         round_frac : fraction of TR required to "round up"
                      e.g. 0.5  : normal rounding
                           1.0  : never round up, floor()
                           0.0  : always round up, ceil()
      """
      if not self.ready: return 1
      if tr <= 0:
         print("** truncate_times: invalid tr %s" % tr)
         return 1

      # convert to fraction to add before truncation (and be sure of type)
      try: rf = 1.0 - round_frac
      except:
         print("** truncate_times: invalid round_frac '%s'" % round_frac)
         return 1

      try: tr = float(tr)
      except:
         print("** truncate_times: invalid tr '%s'" % tr)
         return 1

      if rf < 0.0 or rf > 1.0:
         print("** truncate_times: bad round_frac outside [0,1]: %s" % rf)
         return 1

      if self.verb > 1:
         print('-- Timing: round times to multiples of %s (frac = %s)'%(tr,rf))

      # fr to use for floor and ceil (to assist with fractional TRs)
      if tr == math.floor(tr): tiny = 0.0
      else:                    tiny = 0.0000000001

      # scale mdata and re-create data
      del(self.data)
      self.data = []
      for row in self.mdata:
         for ind in range(len(row)):
            # start with TR index
            tind = row[ind][0]/tr
            # note that rf = 0 now means floor and 1 means ceil

            # add/subtract a tiny fraction even for truncation
            if rf == 1.0   :
               if tind == 0: val = 0.0  # to avoid tiny negatives
               else:         val = math.ceil(tind-tiny) * tr
            elif rf == 0.0 : val = math.floor(tind+tiny) * tr
            else           : val = math.floor(tind+rf) * tr
            row[ind][0] = val
         self.data.append([e[0] for e in row])

      return 0

   def marry_AM(self, mtype, rlens=[], nplaces=-1):
      """add modulator of given type

         lin_run_fraq     : [0,1] modulator = event time/run dur
         lin_run_fraq_off : [0,1] modulator = event time from first/run dur
         lin_event_index  : 0, 1, 2, ... = event index

         if rlens is needed should match len(mdata)
      """
      if not self.ready: return 1
      # g_marry_AM_methods = ['lin_run_fraq', 'lin_event_index']

      if mtype not in g_marry_AM_methods:
         print("** marry_AM: invalid mod type '%s'" % mtype)
         return 1

      if len(rlens) == 1 and len(self.mdata) > 1:
         rlens = rlens*len(self.mdata)

      # most types need the run lengths
      if mtype != 'lin_event_index':
         nruns = len(self.mdata)
         nlens = len(rlens)
         if nlens == 0:
            print('** marry_AM needs run lengths')
            return 1
         if nlens != nruns:
            print('** marry_AM: have %d runs but %d run lengths'%(nruns, nlens))
            return 1

      # append the modulator to each event in mdata
      for rind, mrun in enumerate(self.mdata):
         if mtype == 'lin_event_index': rlen = 0
         else:                          rlen = rlens[rind]

         for ind, event in enumerate(mrun):
            if mtype == 'lin_event_index':
               mval = ind+1
            elif mtype == 'lin_run_fraq':
               mval = event[0]/float(rlen)
               if nplaces >= 0:
                  power = 10.0**nplaces
                  mval = int(power*mval)/power
            event[1].append(mval)

      # and be sure it is married
      self.married = 1
      self.mtype |= LD.MTYPE_AMP

      return 0

   def write_times(self, fname='', nplaces=-1, mplaces=-1, force_married=0):
      """write the current M timing out, with nplaces right of the decimal

         if force_married, force to married timing, if appropriate
      """

      # inherited from lib_afni1D
      simple = 1
      if force_married:
         self.write_dm = 1
         simple = 0
         if self.verb > 2:
            print("-- forcing write timing as married")
      self.write_as_timing(fname, nplaces=nplaces, mplaces=mplaces,
                                  check_simple=simple)

      return 0

   def timing_to_1D(self, run_len, tr, min_frac, per_run=0, allow_warns=0,
                    write_mods=0):
      """return a 0/1 array of trs surviving min_frac cutoff

                run_len   : list of run durations, in seconds
                            (each must encompass the last stimulus, of course)
                tr        : time resolution of output 1D file
                min_frac  : minimum TR fraction occupied by stimulus required
                            for a set value in output array
                            (must be in [0.0, 1.0])
                per_run   : if set, result is list of runs, else catenated
                write_mods: if set, do not write 1.0, but the modulator

         return an error string and the 0/1 array
                (on success, the error string is '')
      """

      if self.verb > 2:
         print("-- t21D: per_run %s, allow_wars %s, write_mods %s" \
               % (per_run, allow_warns, write_mods) )

      # maybe the user provided only one run length
      if self.nrows > 0 and len(run_len) == 1:
         run_len = [run_len[0]] * self.nrows

      errstr, result, modres = self.timing_to_tr_frac(run_len, tr, per_run,
                               allow_warns=allow_warns, write_mods=write_mods)

      if errstr != '' or len(result) < 1: return errstr, result

      if per_run:
         new_res = []
         for rind, row in enumerate(result):
            thr_res = self.get_threshold_list(row, min_frac)
            # if write_mods, convert binary to modulator
            if write_mods:
               mrow = modres[rind]
               if len(mrow) == len(thr_res):
                  for ind in range(len(thr_res)): thr_res[ind] *= mrow[ind]
            new_res.append(thr_res)
      else:
         new_res = self.get_threshold_list(result, min_frac)
         # if write_mods, convert binary to modulator
         if write_mods and len(modres) == len(new_res):
            for ind in range(len(new_res)): new_res[ind] *= modres[ind]

      del(result)
      del(modres)

      return '', new_res

   def get_threshold_list(self, data, min_val):
      """return a 0/1 copy of 'data', with each value set iff >= min_val"""
      result = [0 for ind in range(len(data))]
      for ind, val in enumerate(data):
         if val >= min_val: result[ind] = 1
      return result

   def timing_to_tr_frac(self, run_len, tr, per_run=0, allow_warns=0,
                         write_mods=0):
      """return an array of stim fractions, where is value is the fraction
         of the current TR occupied by stimulus

         --> if per_run is set, output will be one row per run

         The output will be of length sum(TRs per run) and will be in [0,1.0].
         Note that a single stimulus can span multiple TRs.

         ** save this comment for elsewhere
         ** save min_tr_frac for elsewhere
                min_tr_frac     : minimum fraction of a TR required to set it
                                  (must be in (0.0, 1.0])
         Note that if they are not TR-locked and min_tr_frac <= 0.5, then
         one stimulus lasting one TR but occurring on the half TR can result
         in a pair of consecutive 1s.
         ** end save

                run_len         : list of run durations, in seconds
                                  (each must encompass the last TR, of course)
                tr              : time resolution of output 1D file
                per_run         : make a list of results (else all one)
                allow_warns     : make some issues non-fatal
                write_mods      : return mod values along with tr fractions

         On success, the error string should be empty and stim_list should not.

         return error string, stim list, mod list (if write_mods)

         ** testing **

         import math, copy
         from afnipy import afni_util as UTIL, lib_timing as LT
         reload LT
         t = LT.AfniTiming('ch_fltr.txt')
         t.timing_to_tr_frac(run_len, 2.5)
         t.timing_to_1D(run_len, 2.5, 0.3)
      """

      if not self.ready:
         return '** M Timing: nothing to compute ISI stats from', [], []

      #if not self.mtype & LD.MTYPE_DUR:
      #   return '** M Timing: cannot compute stats without stim duration', []

      if self.nrows != len(self.data):
         return '** bad MTiming, nrows=%d, datalen=%d, failing...' % \
               (self.nrows, len(self.data)), [], []

      if self.nrows != len(run_len):
         return '** run_len list is %d of %d runs in timing_to_1D: %s'   \
               % (len(run_len), self.nrows, run_len), [], []

      if tr <= 0.0:
         return '** timing_to_tr, illegal TR <= 0: %g' % tr, [], []

      if write_mods and not self.married:
         print("** to1D: asking for modulatros, but times are not married")
         write_mods = 0

      # make a sorted copy of format run x stim x [start,end], i.e. is 3-D
      tdata = self.get_start_end_timing(sort=1)

      if self.verb > 1:
         if write_mods: mstr = ', will write mods'
         else:          mstr = ''
         print('-- timing_to_tr_fr, tr = %g, nruns = %d%s' \
               % (tr,len(run_len),mstr))

      # need to check each run for being empty
      for ind, data in enumerate(tdata):
          if len(data) < 1: continue
          if data[-1][1] > run_len[ind] or run_len[ind] < 0:
             entry = data[-1]
             if entry[0] > run_len[ind] or run_len[ind] < 0:
                return '** run %d, stim starts after end of run' % (ind+1),[],[]
             elif not allow_warns:
                return '** run %d, stim ends after end of run' % (ind+1), [],[]
             else:
                # entry[1] > run_len[ind]
                print('** WARNING: run %d stim ends after end of run'%(ind+1))
                print('            onset %g, offset %g, end of run %g' \
                      % (entry[0], entry[1], run_len[ind]))
                print('        --> truncating end time')
                entry[1] = entry[0] + 0.99*(run_len[ind]-entry[0])
          
      result = []
      modres = []   # for modulated results
      # process one run at a time, first converting to TR indices
      for rind, data in enumerate(tdata):
         if self.verb > 4:
            print('\n++ stimulus on/off times, run %d :' % (rind+1))
            print(data)

         for tind in range(len(data)):  # convert seconds to TRs
            data[tind][0] = round(data[tind][0]/float(tr),3)
            data[tind][1] = round(data[tind][1]/float(tr),3)

            if tind > 0 and data[tind][0] < data[tind-1][1]:
               estr = '(event times %g and %g)' \
                      % (tr*data[tind-1][0], tr*data[tind][0])
               emesg = '** run %d, index %d %s, stimulus overlap' \
                         % (rind+1, tind, estr)
               if allow_warns:
                  print(emesg)
                       
               else:
                  return emesg, [], []

         if self.verb > 4:
            print('++ stimulus on/off TR times, run %d :' % (rind+1))
            print(data)
         if self.verb > 3:
            print('++ tr fractions, run %d :' % (rind+1))
            print([data[tind][1]-data[tind][0] for tind in range(len(data))])

         # do the real work, for each stimulus, fill appropriate tr fractions
         # init frac list with TR timing (and enough TRs for run)
         num_trs = int(math.ceil(run_len[rind]/float(tr)))
         if self.verb>2: print('-- TR frac: have %d TRs and %d events over run'\
                               % (num_trs, len(data)))
         rdata = [0] * num_trs
         mdata = [0] * num_trs
         for sind in range(len(data)):
            # note the first and last time indices
            start  = data[sind][0]      # initial TR (fractional) time
            end    = data[sind][1]      # end TR time
            startp = int(start)         # initial TR index
            endp   = int(end)           # end TR index

            # and decide on any modval, in case of write_mods
            modval = 0.0
            if write_mods:
               if len(data[sind]) < 3:
                  print("** married but missing mods to write in %d: %s" \
                        % (sind, data[sind]))
                  write_mods = 0
               elif len(data[sind][2]) == 0:
                  print("** married but no mods to write in %d: %s" \
                        % (sind, data[sind]))
                  write_mods = 0
               else:
                  modval = data[sind][2][0]

            # deal with easy case : quick process of single TR result
            if endp == startp:
               rdata[startp] = end-start
               mdata[startp] = modval
               continue

            # otherwise, fill startp, intermediate, and endp fractions
            # (start and end are probably < 1, intermediates are 1)
            rdata[startp] = round(1-(start-startp),3)
            for tind in range(startp+1,endp):
               rdata[tind] = 1.0
            rdata[endp] = round(end-endp, 3)

            # if mods, just write all startp through endp as modval
            if write_mods:
               for tind in range(startp,endp+1):
                  mdata[tind] = modval

         if self.verb > 3:
            print('\n++ timing_to_tr_fr, result for run %d:' % (rind+1))
            print(' '.join(["%g" % r for r in rdata]))
            if write_mods:
               print(' '.join(["%g" % r for r in mdata]))

         if per_run:
            result.append(rdata)
            modres.append(mdata)
         else:
            result.extend(rdata)
            modres.extend(mdata)

      del(tdata)
      del(rdata)
      del(mdata)

      return '', result, modres

   def show_isi_stats(self, mesg='', run_len=[], tr=0, rest_file=''):
      """display ISI timing statistics

            mesg        : display the user message first
            run_len     : can be empty, length 1 or length nrows
            tr          : if > 0: show mean/stdev for stimuli within TRs
                          (so 0 <= mean < tr)
            rest_file   : if set, save all rest durations

         display these statistics:
            - total time, total stim time, total rest time

            - total time per run
            - total stim time per run
            - total rest time per run

            - pre-stim rest per run
            - post-stim response rest per run (if run_len is given)
            - total ISI rest per run

            - min/mean/max (stdev) of stim duration
            - min/mean/max (stdev) of ISI rest
      """

      if not self.ready:
         print('** M Timing: nothing to compute ISI stats from')
         return 1

      if not (self.mtype & LD.MTYPE_DUR):
         print('** warning: computing stats without duration')

      if self.nrows != len(self.data):
         print('** bad MTiming, nrows=%d, datalen=%d, failing...' % \
               (self.nrows, len(self.data)))
         return 1

      # make a sorted copy
      scopy = self.copy()
      scopy.sort()

      # make a copy of format run x stim x [start,end], i.e. is 3-D
      tdata = scopy.get_start_end_timing()

      # make an updated run lengths list
      if len(run_len) == 0:
         rlens = [0 for rind in range(self.nrows)]
      elif len(run_len) == 1:
         rlens = [run_len[0] for rind in range(self.nrows)]
      elif len(run_len) == self.nrows:
         rlens = run_len
      else:     # failure
         print('** invalid len(run_len)=%d, must be one of 0,1,%d' % \
               (len(run_len), self.nrows))
         return 1

      if self.verb > 1:
         print('-- show_isi_stats, run_len = %s, tr = %s, rest_file = %s' \
               % (run_len, tr, rest_file))

      if self.verb > 3:
         print(scopy.make_data_string(nplaces=1, flag_empty=0, check_simple=0,
                        mesg='scopy data'))
      
      all_stim  = []    # all stimulus durations (one row per run)
      all_isi   = []    # all isi times (one row per run)
      pre_time  = []    # pre-stim, per run
      post_time = []    # pose-stim, per run
      run_time  = []    # total run time, per run
      errs      = 0     # allow a few errors before failing
      max_errs  = 10
      for rind in range(self.nrows):
         run  = tdata[rind]
         rlen = rlens[rind]

         if len(run) == 0:      # empty run
            all_stim.append([])
            all_isi.append([])
            pre_time.append(rlen)
            post_time.append(0.0)
            run_time.append(rlen)
            continue

         # if no run len specified, use end of last stimulus
         if rlen == 0.0: rlen = run[-1][1]
         elif rlen < run[-1][1]:
            print('** run %d: given length = %s, last stim ends at %s' % \
                  (rind+1, rlen, run[-1][1]))
            errs += 1
            if errs > max_errs:
               print('** bailing...')
               return 1

         # pre- and post-stim times are set
         pre = run[0][0]
         post = rlen - run[-1][1]

         if pre < 0:
            print('** ISI error: first stimulus of run %d at negative time %s'%\
                  (rind+1, run[0][0]))
            errs += 1
            if errs > max_errs:
               print('** bailing...')
               return 1

         # init accumulation vars
         stimes = [run[0][1] - run[0][0]]
         itimes = []

         # for each following index, update stim and isi times
         # (check for bad overlap)
         tot_olap = 0
         n_olap = 0
         for sind in range(1, len(run)):
            stime = run[sind][1]-run[sind][0]
            itime = run[sind][0]-run[sind-1][1]
            olap  = -itime # for clarity
            # there may be float/ascii reason for a tiny overlap...
            if olap > 0.01:
               if self.verb > 1:
                  print('** ISI error: stimuli overlap at run %d, time %s,' \
                        ' overlap %s' % (rind+1, run[sind][0], olap))
            if olap > 0.0001:
               # get rid of the overlap (note that itime < 0)
               n_olap += 1
               # fix itime and stime
               itime = 0
               stime -= olap
               tot_olap += olap
            stimes.append(stime)
            itimes.append(itime)

         if n_olap > 0:
            print('** run %d, adjusted for %d overlaps in stim, total = %g s' \
                  % (rind, n_olap, tot_olap))

         # store results
         all_stim.append(stimes)
         all_isi.append(itimes)
         pre_time.append(pre)
         post_time.append(post)
         run_time.append(rlen)

      if errs > 0: return 1

      # tally the results
      rtot_stim = [] ; rtot_isi = [] ; rtot_rest = []
      stim_list = [] ; isi_list = [] ; nstim_list = []
      for rind in range(self.nrows):
         rtot_stim.append(UTIL.loc_sum(all_stim[rind]))
         rtot_rest.append(pre_time[rind] + UTIL.loc_sum(all_isi[rind]) +
                          post_time[rind])
         rtot_isi.append(UTIL.loc_sum(all_isi[rind]))
         stim_list.extend(all_stim[rind])
         isi_list.extend(all_isi[rind])
         nstim_list.append(len(all_stim[rind]))

      if mesg: mstr = '(%s) ' % mesg
      else:    mstr = ''
      print('\nISI statistics %s:\n' % mstr)

      print('                        total      per run')
      print('                       ------      ------------------------------')
      print('    total time         %6.1f     %s'   % \
            (UTIL.loc_sum(run_time), float_list_string(run_time, ndec=1)))
      print('    total time: stim   %6.1f     %s'   % \
            (UTIL.loc_sum(rtot_stim),float_list_string(rtot_stim,7,ndec=1)))
      print('    total time: rest   %6.1f     %s'   % \
            (UTIL.loc_sum(rtot_rest),float_list_string(rtot_rest,7,ndec=1)))
      print('')
      print('    rest: total isi    %6.1f     %s'   % \
            (UTIL.loc_sum(rtot_isi), float_list_string(rtot_isi,7,ndec=1)))
      print('    rest: pre stim     %6.1f     %s'   % \
            (UTIL.loc_sum(pre_time), float_list_string(pre_time,7,ndec=1)))
      print('    rest: post stim    %6.1f     %s'   % \
            (UTIL.loc_sum(post_time),float_list_string(post_time,7,ndec=1)))
      print('')
      print('    num stimuli      %6d     %s'   % \
            (UTIL.loc_sum(nstim_list), float_list_string(nstim_list,7,ndec=0)))
      print('\n')

      print('                         min      mean     max     stdev')
      print('                       -------  -------  -------  -------')

      print('    rest: pre-stim     %7.3f  %7.3f  %7.3f  %7.3f' % \
            (UTIL.min_mean_max_stdev(pre_time)))
      print('    rest: post-stim    %7.3f  %7.3f  %7.3f  %7.3f' % \
            (UTIL.min_mean_max_stdev(post_time)))
      print('')

      for ind in range(self.nrows):
         m0, m1, m2, s = UTIL.min_mean_max_stdev(all_isi[ind])
         print('    rest: run #%d ISI   %7.3f  %7.3f  %7.3f  %7.3f' % \
                    (ind, m0, m1, m2, s))

      print('')
      print('    all runs: ISI      %7.3f  %7.3f  %7.3f  %7.3f' % \
            (UTIL.min_mean_max_stdev(isi_list)))
      print('    all runs: stimuli  %7.3f  %7.3f  %7.3f  %7.3f' % \
            (UTIL.min_mean_max_stdev(stim_list)))
      print('')

      # and possibly print out offset info
      if tr > 0: self.show_TR_offset_stats(tr, '')

      # maybe write all rest durations
      if rest_file:
         all_rest = copy.deepcopy(all_isi)
         for run, rest in enumerate(all_rest):
             rest[:0] = [pre_time[run]]
             rest.append(post_time[run])
         UTIL.write_to_timing_file(all_rest, rest_file)

      # clean up, just to be kind
      del(all_stim); del(all_isi); del(pre_time); del(post_time); del(run_time)
      
      del(rtot_stim); del(rtot_isi); del(rtot_rest); del(stim_list)
      del(isi_list); del(nstim_list)

      del(scopy)
      del(tdata)

   def get_TR_offsets(self, tr):
      """create a list of TR offsets (per-run and overall)

            tr : must be positive

         return: all offsets
      """

      if not self.ready:
         print('** M Timing: nothing to compute ISI stats from')
         return []

      if self.nrows != len(self.data):
         print('** bad MTiming, nrows=%d, datalen=%d, failing...' % \
               (self.nrows, len(self.data)))
         return []

      if tr < 0.0:
         print('** show_TR_offset_stats: invalid TR %s' % tr)
         return []

      tr = float(tr) # to be sure

      # make a copy of format run x stim x [start,end], i.e. is 3-D
      tdata = self.get_start_end_timing()

      offsets   = []    # stim offsets within given TRs
      for rind in range(self.nrows):
         run  = tdata[rind]
         if len(run) == 0: continue

         roffsets = UTIL.interval_offsets([val[0] for val in run], tr)
         offsets.extend(roffsets)

      return offsets

   def get_TR_offset_stats(self, tr):
      """create a list of TR offsets (per-run and overall)

            tr : must be positive

         return: 8 values in a list:
                    min, mean, maxabs, stdev of absolute and fractional offsets
                 empty list on error
      """

      if not self.ready:
         print('** M Timing: nothing to compute ISI stats from')
         return []

      if self.nrows != len(self.data):
         print('** bad MTiming, nrows=%d, datalen=%d, failing...' % \
               (self.nrows, len(self.data)))
         return []

      if tr < 0.0:
         print('** show_TR_offset_stats: invalid TR %s' % tr)
         return []

      tr = float(tr) # to be sure

      # make a copy of format run x stim x [start,end], i.e. is 3-D
      tdata = self.get_start_end_timing()

      offsets   = []    # stim offsets within given TRs
      for rind in range(self.nrows):
         run  = tdata[rind]
         if len(run) == 0: continue

         roffsets = UTIL.interval_offsets([val[0] for val in run], tr)
         offsets.extend(roffsets)

      if len(offsets) < 1: return []

      # get overall stats (absolute and fractional)

      # absolute
      m0, m1, m2, s = UTIL.min_mean_max_stdev(offsets)
      offmn = m0; offm = m1; offs = s
      mn = abs(min(offsets))
      offmax = abs(max(offsets))
      if mn > offmax: offmax = mn       

      # fractional
      for ind, val in enumerate(offsets):
         offsets[ind] = val/tr
      m0, m1, m2, s = UTIL.min_mean_max_stdev(offsets)

      del(offsets)

      return [offmn, offm, offmax, offs, m0, m1, offmax/tr, s]
      
   def show_TR_offset_stats(self, tr, mesg='', wlimit=0.4):
      """display statistics regarding within-TR offsets of stimuli

            tr          : show mean/stdev for stimuli within TRs
                          (so 0 <= mean < tr)
            mesg        : display the user message in the output
      """

      rv, rstr = self.get_TR_offset_stats_str(tr, mesg=mesg, wlimit=wlimit)
      print(rstr)

   def get_TR_offset_stats_str(self, tr, mesg='', wlimit=0.4):
      """return a string to display statistics regarding within-TR
                offsets of stimuli

            tr          : show mean/stdev for stimuli within TRs
                          (so 0 <= mean < tr)
            mesg        : display the user message in the output

         return status, stats string

                status > 0 : success, warnings were issued
                       = 0 : success, no warnings
                       < 0 : errors
      """

      if not self.ready:
         return 1, '** M Timing: nothing to compute ISI stats from'

      if self.nrows != len(self.data):
         return 1, '** bad MTiming, nrows=%d, datalen=%d, failing...' % \
                   (self.nrows, len(self.data))

      if tr < 0.0:
         return 1, '** show_TR_offset_stats: invalid TR %s' % tr

      off_means = []    # ... means per run
      off_stdev = []    # ... stdevs per run
      for rind in range(self.nrows):
         run  = self.data[rind]
         if len(run) == 0: continue

         # start with list of time remainders (offsets) within each TR
         roffsets = UTIL.interval_offsets([val for val in run], tr)

         m0, m1, m2, s = UTIL.min_mean_max_stdev(roffsets)
         off_means.append(m1)
         off_stdev.append(s)

      # if no events ere found, we're outta here
      if len(off_means) == 0:
         print('file %s: no events?' % self.name)
         return 0, ''

      # and get overall stats (absolute and fractional)
      offs = self.get_TR_offset_stats(tr)

      # print out offset info
      if mesg: mstr = '(%s) ' % mesg
      else:    mstr = ''

      rstr = '\nwithin-TR stimulus offset statistics %s:\n' % mstr

      hdr1 = '    overall:     '
      hdr2 = '    fractional:  '
      shdr = '                 '
      rv   = 0
      if self.nrows > 1:
         rstr += '                       per run\n'                        \
                 '                       ------------------------------\n' \
                 '    offset means       %s\n'                             \
                 '    offset stdevs      %s\n'                             \
                 '\n'                                                      \
                 % (float_list_string(off_means,ndec=3),
                    float_list_string(off_stdev,ndec=3))
      else: hdr1 = '    one run:     '

      rstr += '%smean = %.3f  maxoff = %.3f  stdev = %.4f\n'   \
              % (hdr1, offs[1], offs[2], offs[3])
      rstr += '%smean = %.3f  maxoff = %.3f  stdev = %.4f\n' \
              % (hdr2, offs[5], offs[6], offs[7])

      # a warning may be issued if the min is positive and the max is small
      if offs[6] == 0: rstr += '\n%s(stimuli are TR-locked)\n' % shdr
      elif wlimit > 0.0 and offs[4] > 0 and offs[6] < wlimit:
         rstr += '\n'                                                         \
            '%s** WARNING: small maxoff suggests (almost) TR-locked stimuli\n'\
            '%s   consider: timing_tool.py -round_times (if basis = TENT)\n'  \
            %(shdr,shdr)
         rv = 1

      # clean up, just to be kind
      del(off_means); del(off_stdev)

      return rv, rstr

   def modulator_stats_str(self, perrun=1):
      """return a string to display detailed statistics regarding
                amplitude modulators within a stimulus file

            perrun      : if more than 1 run, also include per-run results

         return status, stats string

                status > 0 : success, warnings were issued
                       = 0 : success, no warnings
                       < 0 : errors
      """
      if not self.ready:
         return 1, '** Timing: nothing to compute ISI stats from'
      
      num_am = self.num_amplitudes()
      if num_am == 0:
         return 1, '** no amplitude modulators for file %s' % self.name

      rv = 0

      # print out offset info
      if self.nrows > 1: rstr = '%d runs' % self.nrows
      else:              rstr = '1 run'
      rstr = '\namplitude modulator statistics (%s - %s):\n\n' \
             % (self.fname, rstr)

      # ------------------------------------------------------------
      # if more than one run, display per run stats
      if self.nrows > 1 and perrun:

         # for each modulator (index)
         for amindex in range(num_am):
            # mmms values across runs : [ [m,m,m,s], [m,m,m,s], ... ]
            rlist = []
            for rindex in range(self.nrows):
               amlist = self.get_am_list(amindex, rindex)
               rlist.append(list(UTIL.min_mean_max_stdev(amlist)))

            pstr = 'mod_%d' % amindex
            rstr += self._mms_per_run_str(rlist, pmesg=pstr, mtype='')

      # ------------------------------------------------------------
      # either way (1 run or more), display similar overall stats

      # for each modulator, per run lists of min, mean, max, stdev
      # [ MOD_0, MOD_1, ...], MOD_0 = [ [m,m,m,s], [m,m,m,s] ... runs ]
      mmms_lists = []
      lablist = []
      for amindex in range(num_am):
         amlist = self.get_am_list(amindex, 0)
         mmms = list(UTIL.min_mean_max_stdev(amlist))
         # each as a 1-run list of mmms values
         mmms_lists.append([mmms])
         lablist.append("mod_%d" % amindex)

      rstr += self._mms_overall_str(mmms_lists, lablist)

      return rv, rstr

   def detailed_TR_offset_stats_str(self, tr, perrun=1):
      """return a string to display detailed statistics regarding within-TR
                offsets of stimuli

            tr          : show mean/stdev for stimuli within TRs
                          (so 0 <= mean < tr)
            perrun      : if more than 1 run, also include per-run results

         This is like detailed_TR_offset_stats_str, but is more detailed.

         return status, stats string

                status > 0 : success, warnings were issued
                       = 0 : success, no warnings
                       < 0 : errors
      """
      if not self.ready:
         return 1, '** Timing: nothing to compute ISI stats from'
      
      if self.nrows != len(self.data):
         return 1, '** bad Timing, nrows=%d, datalen=%d, failing...' % \
                   (self.nrows, len(self.data))

      if tr < 0.0:
         return 1, '** show_TR_offset_stats: invalid TR %s' % tr

      rv = 0

      # print out offset info
      rstr = '\nwithin-TR stimulus offset statistics (%s):\n\n' % self.fname

      # ------------------------------------------------------------
      # if more than one run, display per run stats
      if self.nrows > 1 and perrun:

         # per run lists of min, mean, max, stdev
         # off_mmms : mmms for offsets, in seconds
         # fr_mmms  : mmms for TR fractional offsets, in [0,1.0)
         # frd_mmms : mmms for diffs of fr offsets, in [0,1]
         # fru_mmms : mmms for diffs of unique fr offsets, in [0,1]
         off_mmms, fr_mmms, frd_mmms, fru_mmms \
            = self.get_offset_mmms_lists(self.data, tr)

         # if no events are found, we're outta here
         if len(off_mmms) == 0:
            print('file %s: no events?' % self.name)
            return 0, ''

         rstr += self._mms_per_run_str(off_mmms, 'in seconds')
         rstr += self._mms_per_run_str(fr_mmms,  'in TR fractions', c=1)
         rstr += self._mms_per_run_str(frd_mmms, 'diffs of frac onsets')
         rstr += self._mms_per_run_str(fru_mmms, 'diffs of unique fracs')

      # ------------------------------------------------------------
      # either way (1 run or more), display similar overall stats

      allruns = []
      for drun in self.data:
         allruns.extend(drun)

      off_mmms, fr_mmms, frd_mmms, fru_mmms \
         = self.get_offset_mmms_lists([allruns], tr)

      # if no events ere found, we're outta here
      if len(off_mmms) == 0:
         print('file %s: no events?' % self.name)
         return 0, ''

      rstr += self._mms_overall_str([off_mmms, fr_mmms, frd_mmms, fru_mmms],
                                    ['offsets', 'frac', 'diffs', 'unique'])

      # comments are now in the overall str

      # clean up, just to be kind
      del(off_mmms); del(fr_mmms); del(frd_mmms); del(fru_mmms)

      return rv, rstr

   def _mms_overall_str(self, mmms_list, name_list, mtype='global',
                        ndec=3, do_c=1):
      """
            mmms_list : array of (single element array) mmms lists
                      : [ [[m,m,m,s]], [[m,m,m,s]], ... ]
                      : (where each mmms_list *could* have had multiple runs)
            name_list : names for each mmms list
            mtype     : type of mmms values to show, max 10 chars
            ndec      : number of decimal places to print
            do_c      : flag: do we show comment strings
      """

      lmin  = [m[0][0] for m in mmms_list]
      lmean = [m[0][1] for m in mmms_list]
      lmax  = [m[0][2] for m in mmms_list]
      lstd  = [m[0][3] for m in mmms_list]

      # make a single formatted string from the names
      nlist = ['%8s' % n for n in name_list]
      nstr  = ' '.join(nlist)

      # header:      each of the names in name_list
      # global vals: the mmms values for each of the mmms lists
      # print using 7.3 + 2 char spacing
      tstr = '%10s' % mtype
      mstr = '\n'                                                       \
             '                   %s\n'                                  \
             '                    ----------------------------------\n' \
             '%s min      %s\n'                                         \
             '%s mean     %s\n'                                         \
             '%s max      %s\n'                                         \
             '%s stdev    %s\n\n'                                       \
             % (nstr,
                tstr, float_list_string(lmin  ,ndec=ndec),
                tstr, float_list_string(lmean ,ndec=ndec),
                tstr, float_list_string(lmax  ,ndec=ndec),
                tstr, float_list_string(lstd  ,ndec=ndec))

      # generate any comment string, if there is something to say
      # this might grow : e.g. pass 'frac' and maybe 'unique' to get comment
      if do_c:
         cstr = ''
         if 'frac' in name_list:
            ind = name_list.index('frac')
            mmms = mmms_list[ind][0]

            c = self._mmms_comment(mmms)
            if c != '':
               mstr += '           comment: %s\n\n' % c

      return mstr

   def _mms_per_run_str(self, mmms, pmesg='', mtype='offset', ndec=3, c=0):
      """
            mmms    : array of mmms values, one set per run
                      [min, mean, max, stdev], [min, mean, max, stdev], ... ]
            pmesg   : parenthetical message to print
            mtype   : type of mms values to show, max 10 chars
            ndec    : number of decimal places to print
            c       : flag: do we show comment strings
      """

      lmin  = [m[0] for m in mmms]
      lmean = [m[1] for m in mmms]
      lmax  = [m[2] for m in mmms]
      lstd  = [m[3] for m in mmms]

      # generate and comment string, if there is something to say
      cstr = ''
      if c:
         cl = [self._mmms_comment(m) for m in mmms]
         ml = max([len(cc) for cc in cl])
         if ml > 0:
            # if we have something, pad it out
            cl = ['%7s' % cc for cc in cl]
            cstr = '    any comments...   %s\n' % '  '.join(cl)

      if pmesg != '':
         pstr = ' (%s)' % pmesg
      else:
         pstr = pmesg

      # print using 7.3 + 2 char spacing
      mstr = '                    per run%s\n'                      \
             '                    ------------------------------\n' \
             '%10s min      %s\n'                                   \
             '%10s mean     %s\n'                                   \
             '%10s max      %s\n'                                   \
             '%10s stdev    %s\n'                                   \
             '%s\n'                                                 \
             % (pstr,
                mtype, float_list_string(lmin  ,ndec=ndec),
                mtype, float_list_string(lmean ,ndec=ndec),
                mtype, float_list_string(lmax  ,ndec=ndec),
                mtype, float_list_string(lstd  ,ndec=ndec),
                cstr)

      return mstr

   def _mmms_comment(self, mmms):
      """return any apparent comment about the values
         - if stdev == 0, either tr-locked or constant
         - if min and max are close to but not 0 or 1,
           we may be 'almost' tr-locked
      """
      if mmms[0] == 0 and mmms[3] == 0:
         return 'tr-lock'
      # check close to lock before const
      if mmms[0] > 0 and mmms[2] < 0.1:
         return '~lock'
      if mmms[0] > 0.9:
         return '~lock'
      if mmms[0] > 0 and mmms[3] == 0:
         return 'const'
      return ''

   def get_offset_mmms_lists(self, etimes, tr):
      """given a list (stimulus, presumably) event onset times per run and
         a TR, return four lists of : min, mean, max, stdev of

         offset: for any event time, this is the time *within* a TR
                 - like (etime modulo TR)
                   e.g. etime 73.4, TR 1 -> offset 1.4
                 - an offset will real and in [0, tr)

         per-run mmms lists to return:

            raw   : raw within-TR offset times
            fr    : TR fractional offsets (so val/TR for each val in 'raw')
            frd   : from the diffs of sorted 'fr' values, including 0.0 and 1.0
                    e.g. fr [0.2, 0.3, 0.4] does not just give [0.1, 0.1],
                         but [0.2, 0.1, 0.1, 0.6] to include TR boundaries
            fru   : same as frd, but where 'fr' values are unique
      """

      # per run lists of [min, mean, max, stdev]
      raw_mmms = []     # mmms: for raw times
      fr_mmms  = []     # mmms: for TR fractional offsets, in [0,1.0)
      frd_mmms = []     # mmms: for diffs of fr offsets, in [0,1]
      fru_mmms = []     # mmms: for diffs of unique fr offsets

      # for each run of etimes
      for ron in etimes:
         # skip any empty run
         if len(ron) == 0: continue

         # raw offsets
         roffsets = UTIL.interval_offsets([val for val in ron], tr)
         # self._add_to_mmms_list(roffsets, raw_mmms)
         raw_mmms.append(list(UTIL.min_mean_max_stdev(roffsets)))

         # fractional offsets
         foffsets = [v/tr for v in roffsets]
         fr_mmms.append(list(UTIL.min_mean_max_stdev(foffsets)))
         # self._add_to_mmms_list(foffsets, fr_mmms)

         # diffs of fractional timing offsets
         # for diffs, require foffsets to include 0 and 1 (1, if not locked)
         foffsets.sort()
         if 0.0 not in foffsets: foffsets.insert(0,0.0)
         if (1.0 not in foffsets) and max(foffsets) > 0: foffsets.append(1.0)
         # and get diffs of these
         doffsets = [foffsets[i+1]-foffsets[i] for i in range(len(foffsets)-1)]
         frd_mmms.append(list(UTIL.min_mean_max_stdev(doffsets)))
         # self._add_to_mmms_list(doffsets, frd_mmms)

         # same, but restricted to a list of unique foffsets
         uoffsets = UTIL.get_unique_sublist(foffsets)
         doffsets = [uoffsets[i+1]-uoffsets[i] for i in range(len(uoffsets)-1)]
         fru_mmms.append(list(UTIL.min_mean_max_stdev(doffsets)))
         # self._add_to_mmms_list(doffsets, fru_mmms)

      return raw_mmms, fr_mmms, frd_mmms, fru_mmms

   def _add_to_mmms_list(self, vals, mmms):
      """for the given vals, append mmms results to the 4 mmms lists

         vals  : array of real values
         mmms  : list of 4 lists, for min, mean, max, stdev results
      """
      m0, m1, m2, ss = UTIL.min_mean_max_stdev(vals)
      mmms.append(list(UTIL.min_mean_max_stdev(vals)))

def float_list_string(vals, nchar=7, ndec=3, nspaces=2):
   str = ''
   for val in vals: str += '%*.*f%*s' % (nchar, ndec, val, nspaces, '')

   return str

def read_multi_ncol_tsv(flist, hlabels=None, def_dur_lab=None, 
                        show_only=0, tsv_int=None, verb=1):
   """Read a set of N column tsv (tab separated value) files
         - one file per run
         - each with a list of events for all classes
      and convert it to list of AfniTiming instances.

      A 3 column tsv file should have an optional header line,
      followed by rows of VAL VAL LABEL, separated by tabs.

      Allow 4+ column files that include trailing amplitudes.
      For now, there might be other event labels included,
         so go after trailing floats, until a more robust way
         can be defined.

      Use the labels to set name and possibly fname fields.

         show_only  : just show header details and return
   """

   tlist = []   # all AfniTiming instances to return

   h0 = []      # original header, once set
   cdict = {}   # dictionary of events per class type
                #   - array of event lists, per run
   elist = []   # temporary variable, events for 1 run at a time
   for rind, fname in enumerate(flist):
      if show_only: print("\nparsing TSV file : %s\n" % fname)
      nvals, header, elist = parse_Ncol_tsv(fname, hlabels=hlabels,
                     def_dur_lab=def_dur_lab, show_only=show_only,
                     tsv_int=tsv_int, verb=verb)
      if show_only: continue
      if nvals <= 0: return 1, tlist

      # store original header, else check for consistency
      if not h0:
         h0 = header
         if verb > 1: print('-- RMNCT: header = %s' % header)
      elif h0 != header and verb:
         print('** inconsistent column headers in N column tsv file %s' % fname)
         print('   orig:    %s' % ' '.join(h0))
         print('   current: %s' % ' '.join(header))

      # update list of class names, as they are found,
      # and add to cdict (including empty runs)
      for event in elist:
         cname = event[2]
         if cname not in list(cdict.keys()):
            cdict[cname] = [[]]*rind
            if verb > 4:
               print('++ RMNCT: init cdict[%s] with %s' % (cname, cdict[cname]))

      # partition elist per known class (should be complete, as all class
      # names were added to dict) - okay if empty
      skeys = list(cdict.keys())
      skeys.sort()
      for cname in skeys:
         # there might be an amplitude
         # if nvals > 3:
         cevents = [e for e in elist if e[2] == cname]
         # check for consisency:
         #   note whether AMs exist, and that use is constant
         numam = 0
         if len(cevents) > 0:
            le = len(cevents[0])
            numam = len(cevents[0][3])
            for e in cevents:
               if len(e) != le:
                  print("** inconsistent modulators for condition %s" % cname)
                  return 1, tlist
               if numam != len(e[3]):
                  print("** inconsistent num mods for condition %s" % cname)
                  return 1, tlist

         cevents = [[e[0], e[3], e[1]] for e in cevents]
         cdict[cname].append(cevents)

         if verb > 4:
            print('++ RM3CT: append cdict[%s] with %s' % (cname, cevents))

   # if just showing, we are done
   if show_only: return 0, []

   # now convert to AfniTiming instances
   # (sorted, to provide consistency across timing patterns)
   skeys = list(cdict.keys())
   skeys.sort()
   for cname in skeys:
      mdata = cdict[cname]
      timing = AfniTiming(mdata=cdict[cname], verb=verb)
      # init name and fname based on label, consider ability to change
      timing.name = cname
      timing.fname = 'times.%s.txt' % cname
      tlist.append(timing)
      if verb > 3: timing.show(mesg=('have timing for %s'%cname))

   return 0, tlist

def parse_Ncol_tsv(fname, hlabels=None, 
                   def_dur_lab=None, show_only=0, tsv_int=None, verb=1):
   """Read one N column tsv (tab separated value) file, and return:
        - ncol: -1 on error, else >= 0
        - header list (length ncol)
        - list of onset, duration, label values, amplitudes?

      If all hlabels exist in header labels, then extract those columns.
      If all hlabels are integers, extract those 0-based columns.
      The format for hlabels should indicate :
         onset, duration, type and amplitudes, if any:

         onset time, duration, trial type, amp1, amp2, ... ampA

      If hlabels is None, use defaults, plus any mod_*.

      def_dur_lab  None     : use default = 'missed'
                   'missed' : add event as "missed_LABEL" class
                   LABEL    : alternate column label for missed events

      show_only : show labels, then return

      tsv_int   : if set, write a smaller tsv
                  (columns of interest - union of col_inds and cols_alt)

      An N column tsv file should have an optional header line,
      followed by rows of VAL VAL LABEL, separated by tabs.
   """
   lines = UTIL.read_text_file(fname, lines=1)
   if len(lines) < 1:
      print("** failed parse_Ncol_tsv for '%s'" % fname)
      return -1, [], []

   # if show_only, be verbose
   if show_only and verb < 4: verb = 4

   # pare lines down to useful ones
   newlines = []
   norig = 0
   for lind, line in enumerate(lines):
      vv = line.split('\t')
      llen = len(line)
      vlen = len(vv)
      if vlen == 0 or llen == 0:
         if verb > 2: print('** skipping empty line %d of 3col tsv file %s' \
                            % (lind, fname))
         continue

      # possibly initialize norig
      if norig == 0: norig = vlen

      # require consistency
      if vlen != norig:
         if verb:
            print("** line %d ncol=%d, mismatch with orig ncol %d, v0 = '%s'" \
                  % (lind, vlen, norig, vv[0]))
            print('** skipping bad line of len %d: %s' % (len(line),line))
         continue

      newlines.append(vv)

   lines = newlines

   if len(lines) < 1:
      if verb: print("** parse_Ncol_tsv for '%s' is empty" % fname)
      return -1, [], []
   if norig < 3:
      if verb: print("** parse_Ncol_tsv: bad ncols = %d in %s" % (norig,fname))
      return -1, [], []

   # ----------------------------------------
   # decide on column extraction indices, based on hlabels and lines[0:2]

   # if nothing passed, set to defaults
   # - let tsv_hlabels_to_col_list append any mod_*
   # - if appending cols to col_inds, hlables will also be appended to
   if hlabels == None:
      defmods = 1
      hlabels = g_tsv_def_labels[:]
   else:
      defmods = 0

   col_inds = tsv_hlabels_to_col_list(hlabels, lines, defmods=defmods,
                                      verb=verb)
   if len(col_inds) < 3:
      if verb: print("** failed to make tsv column index list in %s" % fname)
      return -1, [], []

   # ----------------------------------------
   # decide how to handle n/a duration
    
   # default to 'missed'
   if def_dur_lab == None: def_dur_lab = 'missed'

   if def_dur_lab == 'missed':
      col_dur_alt = -1
   else:
      # replace duration label with alt, and see if it works
      h_alt = hlabels[:]
      h_alt[1] = def_dur_lab
      if verb > 1: print("\n-- processing for def_dur_label, %s" % def_dur_lab)
      cols_alt = tsv_hlabels_to_col_list(h_alt, lines, verb=verb)
      if verb > 2: print("++ have cols_alt: %s" \
                         % UTIL.int_list_string(cols_alt,sepstr=', '))
      if len(cols_alt) < 3:
         if verb:
           print("** could not find tsv column '%s' in %s"%(def_dur_lab,fname))
         return -1, [], []
      col_dur_alt = cols_alt[1]

   # perhaps we want to write out cols of interest
   if tsv_int is not None:
      # first be sure any col_dur_alt is in col_inds, then write
      if verb > 1:
         print("== writing to tsv %s" % tsv_int)
      csub = col_inds[:]
      if col_dur_alt >= 0 and col_dur_alt not in csub:
         csub.append(col_dur_alt)
      write_tsv_cols(lines, csub, ofile=tsv_int)

   # if show_only, we are done
   if show_only: return 0, [], []

   # ----------------------------------------
   # set header list, if a header exists
   header = []
   l0 = lines[0]
   try:
      onset = float(l0[col_inds[0]])
      dur   = float(l0[col_inds[1]])
      lab   = l0[col_inds[2]].replace(' ', '_') # convert spaces to underscores
   except:
      l0 = lines.pop(0)
      header = [l0[col_inds[i]] for i in range(len(col_inds))]

   # decide whether there are amplitudes
   ainds = col_inds[3:]

   # ----------------------------------------
   # now lines should be all: onset, duration, label and possibly amplitudes
   slist = []
   oind = col_inds[0]
   dind = col_inds[1]
   lind = col_inds[2]
   for line_no, line in enumerate(lines):

      # allow for alternate duration, another column or as a MISSING event
      missing_event = 0
      dur_txt = line[dind]
      if dur_txt in ['na', 'n/a', '']:
         if col_dur_alt >= 0:
            dur_txt = line[col_dur_alt]
         else: # then code as MISSING
            missing_event = 1
            dur_txt = '0'
         if verb > 2:
           print("-- line %d, swap dur txt [col %d] '%s' with [col %d] '%s'"\
                 % (line_no, dind, line[dind], col_dur_alt, dur_txt))

      try:
         onset = float(line[oind])
         dur = float(dur_txt)
         lab = line[lind].replace(' ', '_')   # convert spaces to underscores
      except:
         if verb:
            print('** bad line Ncol tsv file %s:\n   %s' \
                  % (fname, ' '.join(line)))
            print("   dur_txt = '%s'" % dur_txt)
         return -1, [], []

      # in modulators, check for na
      amps = []
      if len(ainds) > 0:
         # if any na exists, ignore mods
         avals = [line[aind] for aind in ainds]
         if not has_na(avals):
             amps = [float(v) for v in avals]
         elif verb > 3:
             print("-- ignoring mods due to n/a")

      # append new event, possibly with a 'MISSED' label
      if missing_event:
         use_lab = 'MISSED_%s' % lab
      else:
         use_lab = lab

      slist.append([onset, dur, use_lab, amps])

   nuse = len(col_inds)

   return nuse, header, slist

def has_na(vals):
   """return 1 if there are any na-type vals in the list"""
   if len(vals) < 1:
      return 0
   for val in ['na', 'NA', 'n/a', 'N/A']:
      if val in vals:
         return 1
   return 0

def tofloat(val,verb=1):
   """convert to float, but allow na, NA, n/a, N/A"""
   if val in ['na', 'NA', 'n/a', 'N/A']:
      if verb > 3: print("-- converting %s to 0.0" % val)
      return 0.0
   return float(val)

def write_tsv_cols(table, cols, ofile='stdout'):

   if ofile == '-' or ofile == 'stdout':
      fp = sys.stdout
   else:
      try:
         fp = open(ofile, 'w')
      except:
         print("** failed to open '%s' for writing tsv cols" % ofile)
         return 1

   for line in table:
      fp.write('\t'.join(['%s' % line[c] for c in cols]))
      fp.write('\n')

   if fp != sys.stdout:
      fp.close()

   return 0

def tsv_hlabels_to_col_list(hlabs, linelists,
                            defmods=0, mod_prefix='mod_', verb=1):
   """return a list of columns to extract, in order of:
        onset, duration, label, amplitude1, amplitude2 ...

      if defmods: append any columns with headers starting with mod_prefix
      mod_prefix: prefix for any defmods (default modulator columns)

      ** note: if defmods (and any are found), 'hlabs' will be appended to

      if hlabs are integers, extract those columns
      else, try to find in linelists[0]
         if linelists[0] is not text, return 0, 1, 2, ...
         else if labels do not match
            if ncols matches, return 0, 1, 2, ...
            else: fail
   """

   line0 = linelists[0]
   nlabs = len(hlabs)
   ncols = len(line0)

   if nlabs < 3:
      if verb: print("** tsv hlabs has only %d entries: %s" % (nlabs, hlabs))
      return []
   if ncols < 3:
      if verb: print("** tsv file has only %d columns" % ncols)
      return []

   # first, decide whether hlabs are text labels or integers
   try:
      lints = [int(entry) for entry in hlabs]
   except:
      lints = []
   have_label_ints = (len(lints) == len(hlabs))

   if verb > 2:
      print("-- TSV labels, wanted columns: %s" % ' '.join(hlabs))
      print("   enumerated header cols:")
      lfound = '(chosen) '
      lmissing = ' ' * len(lfound)
      for cind, label in enumerate(line0):
         if have_label_ints and cind in lints:
            lstr = lfound
         elif label in hlabs:
            lstr = lfound
         else:
            lstr = lmissing
         print("      col %02d %s: %s" % (cind, lstr, label))

   # if they are integers, we are done
   if len(lints) >= 3:
      if verb > 2: print("-- tsv labels is via index list: %s" % lints)
      return lints


   # decide whether line0 is text (count floats in line[0])
   nfloat = ntext = 0
   for entry in line0:
      try:
         fval = tofloat(entry, verb=verb)
         nfloat += 1
      except:
         ntext += 1

   # if at least 2 floats, assume no header and return
   if nfloat >= 2:
      if verb > 1: print("-- tsv file has no header, using first columns")
      # so hlabs must now be taken sequential, just check the length
      if nlabs <= ncols:
         return [i for i in range(nlabs)]
      else:
         return [i for i in range(ncols)]

   # so we have text in the header, check whether all hlabs are found
   lints = []
   for label in hlabs:
      if label in line0:
         lints.append(line0.index(label))
      elif verb > 0:
         print("** missing chosen label '%s'" % label)

   # if we did not find them all, require same ncols or fail
   if len(lints) < nlabs:
      if nlabs == ncols:
         if verb: print("-- tsv has unexpected labels, assuming sequential")
         return [i for i in range(nlabs)]
      else:
         if verb > 0:
            print("** tsv file has unexpected labels, failing")
            if verb < 3: print("   (consider -verb 3)")
         return []

   # automatically append any columns starting with mod_prefix
   if defmods and mod_prefix != '':
      for lind, label in enumerate(line0):
         if label.startswith(mod_prefix):
            lints.append(lind)
            hlabs.append(label)
       
   # all are found, yay!
   if verb > 1: print("++ found expected labels in tsv file")
   if verb > 2:
       print("-- labels     : %s" % ', '.join(hlabs))
       print("-- index list : %s" % UTIL.int_list_string(lints, sepstr=', '))

   return lints

def read_value_file(fname):
   """read value file, returning generic values in a matrix (no comments)"""
   try: fp = open(fname, 'r')
   except:
      print("** failed to open value file '%s'" % fname)
      return None

   data = []         # data lines

   lind = 0
   for line in fp.readlines():
      lind += 1
      lary = line.split()
      if len(lary) == 0: continue
      if lary[0] == '#': continue

      data.append([x for x in lary if x != '*'])

   fp.close()

   return data

