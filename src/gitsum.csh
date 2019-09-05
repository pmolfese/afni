#!/bin/tcsh

#################################################################################
## This script counts up who is to 'blame' for lines of code/text in the AFNI
## files. It runs for a long time (30+ minutes), running 'git blame' over
## 1000+ files, and grepping out counts from each one for over a dozen
## co-conspirators.
## Output files:
##  gitsum.out.txt    = author line counts (also cat-ed to stdout)
##  gitsum.unkown.txt = lines that had unknown authors (for further research)
##
## -- RWCox -- Sep 2019
#################################################################################

# set list of source files to query
# stuff from outside sources (e.g., eispack, sonnets.h) is not included here

if ( 1 ) then
  set flist = ( af*.[ch] mri*.[ch] thd*.[ch] cs*.[ch] 3d*.[ch] edt*.[ch] suma*.[ch]     \
                plug*.c niml/niml*.[ch] coxplot/*.[ch] SUMA/*.[ch] 1d*.[ch] model*.c    \
                display.[ch] imseq.[ch] bbox.[ch] xim.[ch] xutil.[ch] xutil_webber.[ch] \
                nifti_statlib.[ch] `find nifti -name '*.[ch]'`                          \
                rickr/*.[ch] ptaylor/*.[ch] gifti/*.[ch] svm/*.[ch]                     \
                `git ls-tree --name-only -r master | grep scripts_install`              \
                `find discoraj -type f` `find shiny -type f -name '*.R'`                \
                python_scripts/afni_python/*.py R_scripts/*.R                           \
                ../tests/scripts/*.py ../tests/scripts/utils/*.py                       \
                parser*.c powell_int.c dicom_hdr.c cat_matvec.c plugout*.c file_tool.c  \
                rtfeedme.c DTIStudioFibertoSegments.c im*.c                               )
else
# for quicker testing
  set flist = ( afni.c imseq.c suma_datasets.c )
endif

# run the Count Lines Of Code script, if present (this is fast)

which cloc-1.64.pl >& /dev/null
if ( $status == 0 ) then
  cloc-1.64.pl --quiet $flist
endif

# list of authors needing only one alias (not case sensitive)
# anyone whose alias has spaces in it is out of luck

set alist = ( Cox Craddock discoraj Froehlich Gang        \
              Gaudes Glen Hammett Kaczmarzyk LeeJ3        \
              Laconte Lisinski Clark Johnson Julia        \
              Molfese Oosterhof Rick Schwabacher          \
              Vincent Warren Markello                       )

# list of authors needing two aliases (i.e., troublemakers)
# anyone who has three aliases is out of luck

set blist1 = ( Nielson      Saad Taylor  afniHQ )
set blist2 = ( shotgunosine ziad mrneont Ubuntu )

# tsum = total sum of lines thus far
set tsum = 0

# nn = number of files processed thus far
set nn   = 0

# setup counts for the alist
set anum = $#alist
set aqq  = ( `count -dig 1 1 $anum` )
set asum = ( )
foreach uuu ( $alist )
  set asum = ( $asum 0 )
end

# setup counts for the blist
set bnum = $#blist1
set bqq  = ( `count -dig 1 1 $bnum` )
set bsum = ( )
foreach uuu ( $blist1 )
  set bsum = ( $bsum 0 )
end

# grep command option to remove known authors, to count unknowns

set gunk = ( -v -i )
foreach uuu ( $alist $blist1 $blist2 )
  set gunk = ( $gunk -e $uuu )
end

# will acccumulate all unknown lines in to one file, for later research
if( -f gitsum.unknown.txt ) \rm -f gitsum.unknown.txt
touch gitsum.unknown.txt

# loop over source files, plus README documents

printf "\nblaming "

set glist = ( $flist ../doc/README/README.* )
foreach fff ( $glist )

 # get the list of blamees for this file (grep out blank lines)
  git blame $fff | grep -v '[0-9]) $' > gitsum.junk.txt

 # count lines in this file
  set aa = `wc -l < gitsum.junk.txt` ; @ tsum += $aa

 # loop over the alist and grep out count for each author
  foreach qq ( $aqq )
    set aa = `grep -i -e "$alist[$qq]" gitsum.junk.txt | wc -l` ; @ asum[$qq] += $aa
  end

 # loop over the blist and get their counts
  foreach qq ( $bqq )
    set aa = `grep -i -e "$blist1[$qq]" -e "$blist2[$qq]" gitsum.junk.txt | wc -l` ; @ bsum[$qq] += $aa
  end

 # accumulate the lines with unknown authors into a separate file
  grep $gunk gitsum.junk.txt >> gitsum.unknown.txt

 # print a progress pacifier
  @ nn ++ ; if( $nn % 20 == 0 ) printf "%d/%d " $nn $#glist

end

# cleanup after loop over files

printf "... total line count = %d \n" $tsum
\rm -f gitsum.junk.txt
touch gitsum.junk.txt

# count total number of unknown lines now

set aa = `wc -l < gitsum.unknown.txt` ; @ unksum = $aa

# format lines for the final report into a temp file

foreach qq ( $aqq )
  if ( $asum[$qq] > 0 ) then
    set perc  = `ccalc "100*$asum[$qq]/$tsum"`
    printf " %12s  %6s  %5.2f%%\n" "$alist[$qq]" $asum[$qq] $perc >> gitsum.junk.txt
  endif
end

foreach qq ( $bqq )
  if ( $bsum[$qq] > 0 ) then
    set perc  = `ccalc "100*$bsum[$qq]/$tsum"`
    printf " %12s  %6s  %5.2f%%\n" "$blist1[$qq]" $bsum[$qq] $perc >> gitsum.junk.txt
  endif
end

if ( $unksum > 0 ) then
  set perc  = `ccalc "100*$unksum/$tsum"`
  printf " %12s  %6s  %5.2f%%\n" Unknown $unksum $perc >> gitsum.junk.txt
endif

# Put header lines into the final report

echo " Contributor   Lines   %-age"    > gitsum.out.txt
echo " ------------  ------  ------"  >> gitsum.out.txt

# sort output lines by second column, put in final report

sort -n -r --key=2 gitsum.junk.txt      >> gitsum.out.txt

# let the user see the results

echo
cat gitsum.out.txt

# toss out the junk

\rm -f gitsum.junk.txt

exit 0
