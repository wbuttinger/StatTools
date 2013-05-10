# echo "cleanup StatTools StatTools-00-00-01 in /var/clus/usera/will/testareas/RootDev"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/atlas-nightlies.cern.ch/repo/sw/nightlies/x86_64-slc5-gcc43-opt/17.X.0/rel_0/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtStatToolstempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtStatToolstempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=StatTools -version=StatTools-00-00-01 -path=/var/clus/usera/will/testareas/RootDev  $* >${cmtStatToolstempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=StatTools -version=StatTools-00-00-01 -path=/var/clus/usera/will/testareas/RootDev  $* >${cmtStatToolstempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtStatToolstempfile}
  unset cmtStatToolstempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtStatToolstempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtStatToolstempfile}
unset cmtStatToolstempfile
return $cmtcleanupstatus

