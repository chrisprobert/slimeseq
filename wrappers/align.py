import sys
import os
import subprocess

for filename in os.listdir('./aligned') :
  if "out" not in filename : continue
  root = filename.split('.')[0]
  cmd = ["./pal2nal.pl", "./aligned/" + filename, "./aligned/" + root + ".fasta", "-output", "paml", "-nogap", "-codontable", "4", ">", root + ".paml"]
  subprocess.Popen(cmd).wait()
