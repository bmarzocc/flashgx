
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
    print "Usage: python check_events.py --file=[file]"

#---------------------------------------------------------- MAIN ----------------------------------------------------------
def main():
    
 try:
     opts, args = getopt.getopt(sys.argv[1:], "", ["file=","help"])

 except getopt.GetoptError:
     #* print help information and exit:*
     usage()
     sys.exit(2)


 file = ""
 help = False

 for opt, arg in opts:
    
     if opt in ("--file"):
       file = arg
     if opt in ("--help"):
        help = True     

 if(help == True):
   usage()
   sys.exit(2)

 if(file == ""):
   usage()
   sys.exit(2)

 if(file != ""):
   print "file = ",file

 run = 0
 lumi = 0
 event = 0

 inFile = TFile.Open(file)
 tree = inFile.Get("flashgxanalysis/FlashGXTree")
 tree.SetBranchAddress("run",run)
 tree.SetBranchAddress("lumi",lumi)
 tree.SetBranchAddress("event",event)

 for i in range(0, tree.GetEntries()):
    tree.GetEntry(i)
    print run,lumi,event

if __name__ == "__main__":
    main()

