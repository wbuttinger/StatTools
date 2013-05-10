# echo "cleanup StatTools StatTools-00-00-01 in /var/clus/usera/will/testareas/RootDev"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/atlas-nightlies.cern.ch/repo/sw/nightlies/x86_64-slc5-gcc43-opt/17.X.0/rel_0/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtStatToolstempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtStatToolstempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=StatTools -version=StatTools-00-00-01 -path=/var/clus/usera/will/testareas/RootDev  $* >${cmtStatToolstempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=StatTools -version=StatTools-00-00-01 -path=/var/clus/usera/will/testareas/RootDev  $* >${cmtStatToolstempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtStatToolstempfile}
  unset cmtStatToolstempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtStatToolstempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtStatToolstempfile}
unset cmtStatToolstempfile
exit $cmtcleanupstatus

