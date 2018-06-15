#!/usr/bin/perl

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------

use Env;

#PG lettura dei parametri da cfg file
#PG --------------------------------
print "reading ".$ARGV[0]."\n" ;

open (USERCONFIG,$ARGV[0]) ;

while (<USERCONFIG>)
{
    chomp; 
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
#    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences{$var} = $value;
}

$BASEDir          = $User_Preferences{"BASEDir"};
$INPUTdir         = $User_Preferences{"INPUTdir"};
$LISTOFSamples    = $User_Preferences{"LISTOFSamples"} ;
$JOBCfgTemplate   = $User_Preferences{"JOBCfgTemplate"} ;
$JOBCfiTemplate   = $User_Preferences{"JOBCfiTemplate"} ;
$OUTPUTSAVEPath   = $User_Preferences{"OUTPUTSAVEPath"} ;
$SCRAM_ARCH       = $User_Preferences{"SCRAM_ARCH"} ;
$JOBModulo        = $User_Preferences{"JOBModulo"} ;
$QUEUE            = $User_Preferences{"QUEUE"};


print "BASEDir = "          .$BASEDir."\n" ;
print "INPUTdir = "         .$INPUTdir."\n" ;
print "LISTOFSamples = "    .$LISTOFSamples."\n" ;
print "JOBCfgTemplate = "   .$JOBCfgTemplate."\n" ;
print "JOBCfiTemplate = "   .$JOBCfiTemplate."\n" ;
print "OUTPUTSAVEPath = "   .$OUTPUTSAVEPath."\n" ;
print "JOBModulo = "        .$JOBModulo."\n" ;
print "SCRAM_ARCH = "       .$SCRAM_ARCH."\n" ;
print "QUEUE  = "           .$QUEUE."\n" ;


#####################################################
# PG prepare the array containing the root files list
#####################################################

open (LISTOFSamples,$LISTOFSamples) ;
while (my $line = <LISTOFSamples>)
{
 
 if (index($line,"#") != -1)
 {
     next;
 }

 my @name = split / /, $line;	

 $sampleJobListFile = "./lancia.sh";
 open(SAMPLEJOBLISTFILE, ">", $sampleJobListFile);
 
 $JOBS_dir = $name[0];
 
 system("cd ".$BASEDir."\n");
    
 chomp($_);
    
 print("Sample: ".$JOBS_dir."\n") ;  

 system ("rm -r ".$JOBS_dir."\n") ;
 system ("mkdir ".$JOBS_dir."\n") ;
 system ("rm -rf ".$OUTPUTSAVEPath."/".$JOBS_dir."_ntuples");
    
 
 $LISTOFFiles = "./list_.txt" ;
 system ("find /eos/cms".$INPUTdir."/".$JOBS_dir." | grep root > ".$LISTOFFiles."\n") ;
   
 $totNumber = 0;
 $jobNumber = 0;
  
 open (LISTOFFiles,$LISTOFFiles) ;
 while (<LISTOFFiles>)
 {
	++$totNumber;
 }

 $jobNumber = int($totNumber/$JOBModulo);
 if( $totNumber%$JOBModulo != 0)
 {
	$jobNumber = $jobNumber+1;
 }
    
 print "NumberOfJobs = ".$jobNumber."\n";
    

  
    ################
    # loop over jobs 
    ################
    
 for($jobIt = 1; $jobIt <= $jobNumber; ++$jobIt)
 { 
	$currDir = `pwd` ;
	chomp ($currDir) ;
    
	$jobDir = $currDir."/".$JOBS_dir."/JOB_".$jobIt ;
	system ("mkdir ".$jobDir." \n") ;

        system ("cp mc_infos.dat ".$jobDir." \n") ;
    
	$tempBjob = $jobDir."/bjob_".$jobIt.".sh" ;
	$command = "touch ".$tempBjob ;
	system ($command) ;
	$command = "chmod 777 ".$tempBjob ;
	system ($command) ;
    

	$tempo1 = "./tempo1" ;
	system ("cat ".$JOBCfgTemplate."   | sed -e s%OUTPUTNAME%".$JOBS_dir."_".$jobIt.
		                       "%g > ".$tempo1) ;
    
	$it = 0;
	$JOBLISTOFFiles;

	open (LISTOFFiles2,$LISTOFFiles) ;
	while (<LISTOFFiles2>)
	{
	    chomp; 
	    s/#.*//;                # no comments
	    s/^\s+//;               # no leading white
	    s/\s+$//;               # no trailing white
	    $file = $_ ;
	    
	    if( ($it >= ($jobIt - 1)*$JOBModulo) && ($it < ($jobIt)*$JOBModulo) )
	    { 
                #print $JOBLISTOFFiles."APICE".$file."APICE,";
		#$JOBLISTOFFiles = "APICE".$INPUTSAVEPath."/".$sample."/".$file."APICE,";
                $file = substr($file,8);
		$JOBLISTOFFiles = $JOBLISTOFFiles."APICE".$file."APICE,";
	    }
	    ++$it;
	}
	
        
	$tempo2 = "./tempo2" ;    
	system ("cat ".$tempo1." | sed -e s%LISTOFFILES%".$JOBLISTOFFiles."%g > ".$tempo2) ;
	$JOBLISTOFFiles = "" ;

	$tempo3 = "tempo3" ;
	system ("cat ".$tempo2." | sed -e s%APICE%\\'%g > ".$tempo3) ;

        $JOBCfiFile = "FlashGXAnalysis_cfi" ;
        $tempo4 = "tempo4" ;
	system ("cat ".$tempo3." | sed -e s%NAMECfi%".$JOBCfiFile."%g > ".$tempo4) ;
        $JOBCfgFile = $jobDir."/FlashGXAnalysis_cfg.py" ;
        system ("mv ".$tempo4." ".$JOBCfgFile) ;
        
        $tempo5 = "tempo5" ;
	system ("cat ".$JOBCfiTemplate." | sed -e s%NAMESample%".$JOBS_dir."%g > ".$tempo5) ;
        system ("mv ".$tempo5." ".$jobDir."/".$JOBCfiFile.".py") ;

        system ("rm tempo*") ;    
    
  
    ######################
    # make job files
    ######################    
    
	open (SAMPLEJOBFILE, ">", $tempBjob) or die "Can't open file ".$tempBjob;

	$command = "#!/bin/tcsh" ;
	print SAMPLEJOBFILE $command."\n";

	$command = "cd ".$jobDir ;
	print SAMPLEJOBFILE $command."\n";

	#$command = "setenv SCRAM_ARCH ".$SCRAM_ARCH ;
	#print SAMPLEJOBFILE $command."\n";
    
	$command = "eval `scramv1 ru -csh`" ;
	print SAMPLEJOBFILE $command."\n";
  
	$command = "eos mkdir ".$OUTPUTSAVEPath."/".$JOBS_dir."_ntuples";
	print SAMPLEJOBFILE $command."\n";

	$command = "cmsRun FlashGXAnalysis_cfg.py";
	print SAMPLEJOBFILE $command."\n";

	$command = "eos cp ".$JOBS_dir."_".$jobIt.".root root://eoscms.cern.ch/".$OUTPUTSAVEPath."/".$JOBS_dir."_ntuples/";
	print SAMPLEJOBFILE $command."\n";
        
        $command = "rm *.root";
	print SAMPLEJOBFILE $command."\n";

	
	############
	# submit job
	############
	
        $command = "bsub -cwd ".$jobDir." -q ".$QUEUE." ".$tempBjob."\n" ; 
	print SAMPLEJOBLISTFILE $command."\n";
    
 }

 system ("rm ".$LISTOFFiles) ;
 system ("sh lancia.sh");
 system ("rm lancia.sh");
}

