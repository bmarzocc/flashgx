
import sys, math, os, collections, getopt, subprocess, time, shutil, datetime
from ROOT import *
from array import array

def hadd(ls_dump):
   with open(ls_dump) as f_dump:
      data_dump = f_dump.read()
      lines_dump = data_dump.splitlines()  
      for pos,x in enumerate(lines_dump):
         command = os.system("rm -rf "+str(x)+".root")
         command = os.system("hadd -f "+str(x)+".root "+str(x)+"/*.root")

def usage():
    print "Usage: python hadd_samples.py --dir=[dir]"

#---------------------------------------------------------- MAIN ----------------------------------------------------------
def main():
    
 try:
     opts, args = getopt.getopt(sys.argv[1:], "", ["dir=","help"])

 except getopt.GetoptError:
     #* print help information and exit:*
     usage()
     sys.exit(2)


 dir = ""
 help = False

 for opt, arg in opts:
    
     if opt in ("--dir"):
        dir = arg
     if opt in ("--help"):
        help = True     

 if(help == True):
   usage()
   sys.exit(2)

 if(dir == ""):
   usage()
   sys.exit(2)

 if(dir != ""):
   print "dir = ",dir

 command = os.system("ls -d /eos/cms"+str(dir)+"/*_ntuples > ls_dump")
 hadd("ls_dump")
 command = os.system("rm ls_dump")

if __name__ == "__main__":
    main()

