# echo "Setting KKpipi 00-00-01 in /afs/ihep.ac.cn/users/m/martintat/workdir/6.6.4.p02"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/contrib/CMT/v1r20p20081118
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=KKpipi -version=00-00-01 -path=/afs/ihep.ac.cn/users/m/martintat/workdir/6.6.4.p02  -no_cleanup $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

