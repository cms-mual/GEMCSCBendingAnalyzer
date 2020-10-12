import os

output_path = "condor_MWGR4_1011/output"

filelist = os.listdir(output_path)


nCSCSegs = 0
nGE11 = 0

for filename in filelist:
  f = open(output_path+"/"+filename)
  for line in f:
    if "Muons with cscSegs =" in line:
      #print line[-2]
      nCSCSegs += int(line[-2])
    if "Muons with prop to GE11 =" in line:
      #print line[-2]
      nGE11 += int(line[-2])


print "total with CSCSegs = ", nCSCSegs

print "total with GE11 prop = ", nGE11

print "ratio of segs propagated = ", float(nGE11)/float(nCSCSegs)
