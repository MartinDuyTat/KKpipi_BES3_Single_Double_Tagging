# echo "Setting KKpipi 00-00-01 in /afs/ihep.ac.cn/users/m/martintat/workdir/6.6.4.p02"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/contrib/CMT/v1r20p20081118; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=KKpipi -version=00-00-01 -path=/afs/ihep.ac.cn/users/m/martintat/workdir/6.6.4.p02  -no_cleanup $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

