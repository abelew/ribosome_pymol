from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkColorChooser
import tkFileDialog
import os, sys, string, re
import urllib
import gzip
from pymol import stored, cmd, selector
## The function 'fetch_then_chains' was taken with very minor changes
## from remote_pdb_load.py.  The copyright notice is at the bottom of
## this file as well as the README

## Set up the directory where the data files live and some global variables
## which will be used to store the names of the molecules and helices from the
## PDB files.
helices_dir = os.path.dirname(__file__)
datadir=helices_dir + "/data/"
global molecule_list
molecule_list = []
global original_list
original_list = []
global helices_list
helices_list = []
global default_colors
default_colors = {
  'default' : 'gray10',
  'helix' : 'red',
  'RACK' : 'cyan',
  'SSU_protein' : 'cyan',
  'LSU_protein' : 'skyblue',
  'unknown' : 'green',
  'tRNA' : 'slate',
  'mRNA' : 'forest',
  'SSU_RNA' : 'gray20',
  'LSU_RNA' : 'gray30',
  '5S_RNA' : 'gray40',
  '5.8S_RNA' : 'gray50',
  'RNA' : 'gray60',
  'other' : 'red'
}
## The delete_ functions do just that, delete the various molecules associated
## with them.  The protein and helices functions do some very simplistic
## string.find calls to identify the helices and proteins.
def delete_all_helices():
  delete_ssu_helices
  delete_lsu_helices

def delete_all_rna():
  delete_ssu_rna
  delete_lsu_rna
  delete_mrna
  delete_trna

def delete_all_protein():
  delete_lsu_protein
  delete_ssu_protein

def delete_ssu_helices():
  for helix in helices_list:
    helix = helix.lstrip('/')
    if helix.find('SSU_h') > -1:
      cmd.delete(helix)

def delete_lsu_helices():
  for helix in helices_list:
    helix = helix.lstrip('/')
    helix = str(helix)
    if helix.find('LSU_H') > -1:
      cmd.delete(helix)

def delete_ssu_protein():
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('S_S') > -1:
      cmd.delete(mol)
    elif re.compile('^S').search(mol) is not None:
      cmd.delete(mol)

def delete_lsu_protein():
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('S_L') > -1:
      cmd.delete(mol)
    elif re.compile('^L').search(mol) is not None:
      cmd.delete(mol)

def delete_ssu_rna():
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('16S_RRNA') > -1:
      cmd.delete(mol)

def delete_lsu_rna():
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('23S_RRNA') > -1:
      cmd.delete(mol)
    elif mol.find('5S_RR') > -1:
      cmd.delete(mol)
    elif mol.find('5.8S') > -1:
      cmd.delete(mol)

def delete_trna():
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('TRNA') > -1:
      cmd.delete(mol)

def delete_mrna():
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('MRNA') > -1:
      cmd.delete(mol)

def delete_original():
  for original_molecule in original_list:
    cmd.delete(original_molecule)

def delete_lsuh():
  counter = 0
  while (counter <= 104):
      counter = counter + 1
      string = "LSU_H", counter
      cmd.delete(string)

def delete_ssuh():
    counter = 0
    while (counter <= 45):
        counter = counter + 1
        string = "SSU_h", counter
        cmd.delete(string)

## I changed like 2 lines from remote_load_pdb.py
## The main change is at the end of fetch()
class fetch_then_chains:
  def __init__(self, app):
    import tkSimpleDialog
    import tkMessageBox
    pdbCode = tkSimpleDialog.askstring('PDB Loader Service',
                                       'Please enter a 4-digit pdb code:',
                                       parent=app.root)
    if pdbCode: # None is returned for user cancel
      pdbCode = string.upper(pdbCode)
      fetch(pdbCode,"")


def fetch(pdb, splitp):
  import urllib
  import gzip
  import os
  import string
  import tkMessageBox
  try:
    filename = urllib.urlretrieve('http://www.rcsb.org/pdb/files/' + pdb + '.pdb.gz')[0]
  except:
    tkMessageBox.showerror('Connection Error',
                           'Can not access to the PDB database.\n'+
                           'Please check your Internet access.',)
  else:
    if (os.path.getsize(filename) > 0): # If 0, then pdb code was invalid
      fpin = gzip.open(filename)
      outputname = os.path.dirname(filename) + os.sep + pdb + '.pdb'
      fpout = open(outputname, 'w')
      fpout.write(fpin.read()) # Write pdb file
      fpin.close()
      fpout.close()
      cmd.load(outputname,quiet=0) # Load the fresh pdb
      ## This is the change from Trey, a callout to random_chains()
      random_chains(outputname,splitp)
    else:
      tkMessageBox.showerror('Invalid Code', 'You entered an invalid pdb code:' + pdb,)
      os.remove(filename) # Remove tmp file (leave the pdb)

## chain_color will color arbitrary bases/residues with the colors
## specified in data/color_definitions.txt
## The residues chosen should be in a text file as per
## data/modifications.txt
def chain_color(bases):
  ## The next 8 lines attempts to figure out what file
  ## to use to define the residues to color
  input_file = ""
  if bases == "modified":
    input_file = datadir + 'modifications.txt'
  else:
    tmp_filename = tkFileDialog.askopenfile(title="Open a session")
    if not tmp_filename: return
    input_file = tmp_filename.name
  comment = ''

  ## From here until 'if input_file:' the color definitions
  ## are specified
  colors_file = datadir + 'color_definitions.txt'
  colors = dict({None : 'gray',})
  if colors_file:
    color_lines = file(colors_file).readlines()
  for color_line in color_lines:
    color_datum = color_line.split()
    colors[color_datum[0]] = color_datum[1]

  ## Now read the input file and select the appropriate residues
  ## and color them according to the rules in the colors dictionary
  if input_file:
    lines = file(input_file).readlines()
  else:
    lines = sys.stdin.readlines()
  subunit = ''
  chain = ''
  for line in lines:
    chain = ""
    if re.compile('^#').search(line) is not None: # skip commented lines
      tmpre = re.compile('^#')
      tmpre = tmpre.sub('', line)
      try:
        (subunit, chain) = tmpre.split()
      except:
        subunit = tmpre.strip()
    else:
      datum = line.split()
      try:
        num = datum[0].strip()
        if chain == "":
          selection_string = '/' + subunit + '///' + num
          selection_name = subunit + '_' + num
        else:
          selection_string = '/' + subunit + '//' + chain + '/' + num
          selection_name = subunit + '_' + chain + '_' + num
        color_choice = datum[1].strip()
        color_name = ''
        try:
          color_name = colors[color_choice]
        except:
          color_name = color_choice
        try:
          cmd.color(color_name, selection_string)
        except:
          cmd.color(colors[None], selection_string)
      except:
        print "Cannot find your selection, perhaps you must split the chains first"


def make_pretty():
  ## These are some settings our professor prefers.
  cmd.bg_color("white")
  cmd.show("cartoon")
  cmd.set("cartoon_ring_mode", 3)

## This function is the toplevel function to make pretty helices
## Change the default_colors['helix'] to whatever color you prefer.
def helices():
  make_chains(organism, 'sticks', default_colors['helix'])

def make_chains(chains, showastype, showascolor):
  ## Start out figuring out the data file to specify the helices
  ## Currently I just have a stupid if/elif chain for the few species
  ## I have annotated.
  if chains == 'wtf':
    print "WTF"
    chains_filenames = [datadir + 'wtf.txt',]
  elif chains == 'ssh':
    chains_filenames = [datadir + 'helices_ssu.txt',]
  elif chains.find('escherichia_coli') > -1:
    chains_filenames = [datadir + 'escherichia_coli.txt',]
  elif chains == 'thermomyces_lanuginosus':
    chains_filenames = [datadir + 'saccharomyces_helices.txt',]
  elif chains == 'thermus_thermophilus':
    chains_filenames =  [datadir + 'thermus_thermophilus.txt',]
  elif chains == 'haloarcula_marismortui':
    chains_filenames = [datadir + 'haloarcula_marismortui.txt',]
  elif chains == 'saccharomyces':
    chains_filenames = [datadir + 'saccharomyces_helices.txt',]
  elif chains == 'saccharomyces_cerevisiae':
    chains_filenames = [datadir + 'saccharomyces_helices.txt',]
  else:
    print "Could not understand the argument:" + chains +", using the wtf file"
    chains_filenames = [datadir + 'wtf.txt',]
  ## Once the species has been decided, open the appropriate file and start
  for chains_filename in chains_filenames:
    if chains_filename:
      chains_lines = file(chains_filename).readlines()
      for ch in chains_lines:
        if re.compile('^#').search(ch) is not None:
          continue
        name = ''
        location = ''
        ## Each line of the file is a name, pymol_specification
        ## so just split by comma and run with it
        (name, location) = ch.split(',')
        try:
          cmd.create(name, location)
          cmd.disable(name)
          new_selection = "/" + name
          if showastype:
            helices_list.append(new_selection)
            cmd.show(showastype, new_selection)
            if showascolor:
              cmd.color(showascolor, new_selection)
        except:
          print "There was an error."
    ## Zoom to something sane
    cmd.zoom("all")


def load_session(filename):
    filename = tkFileDialog.askopenfile(title="Open a session")
    if not filename: return
    file_path = filename.name
    cmd.load(file_path)

## This function will split apart ribosomal PDB files into 
## the individual pieces by reading the header and attempting
## to choose sane names from the information there.
def random_chains(pdb_file, splitp):
  cmd.bg_color("white")
  import re
  if pdb_file is None:
    pdb_file = tkFileDialog.askopenfile(title="Open a session")
  if not pdb_file: return
  ## pdb_filename is the full filename
  ## pdb_shortname is the 2WGD or whathaveyou
  ## pdb_basename is the path it lives in
  ## pdb_file is the file object which has all the attributes etc
  pdb_filename = str(pdb_file)
  pdb_basename = os.path.basename(pdb_filename)
  pdb_shortname = os.path.splitext(pdb_basename)
  pdb_shortname = pdb_shortname[0] 
  original_list.append(pdb_shortname)
  cmd.load(pdb_filename)
  pdb_lines = file(pdb_filename).readlines()
  chain = ''
  source_count = 0
  global organism
  organism = str('saccharomyces_cerevisiae')
  ## The following lines are attempting to properly decide
  ## when to stop reading a PDB file
  for pdb_line in pdb_lines:
    if re.compile("^HEADER").search(pdb_line) is not None:
      continue
    if re.compile("^TITLE").search(pdb_line) is not None:
      continue
    ## If a PDB file has a SPLIT entry, then it is part of
    ## a group.  So make a list of all entries in the group
    ## and fetch/split them all.
    if re.compile("^SPLIT").search(pdb_line) is not None:
      if splitp == "":
        chains_list = pdb_line.split()
        chains_list.pop(0)
        for pdb_id in chains_list:
          if (pdb_id != pdb_shortname):
            fetch(pdb_id,"1")
      else:
        continue
    if re.compile("^CAVEAT").search(pdb_line) is not None:
      continue
    ## The SOURCE stanza contains the species and so will be useful
    ## for figuring out the helices later if need be.
    if re.compile("^SOURCE").search(pdb_line) is not None:
      source_count = source_count + 1
      if (source_count > 3):
        break
      else:
        line_array = pdb_line.split()
        line_type = line_array[0]
        num = line_array[1]
        chain_mol = line_array[2]
        if (chain_mol == 'ORGANISM_SCIENTIFIC:'):
          (pre, organism) = pdb_line.split(": ")
          org = str(organism)
          org = org.rstrip()
          org = org.rstrip(';')
          org = org.lower()
          org = org.replace(' ', '_')
          organism = "%s" %(org)

    ## The surviving lines should proved a means to find the
    ## name of every chain and molecule of the PDB file.
    ## A little minor sterilizing of the molecule's name
    ## will likely be required, and then pull them apart,
    ## select them, and color them, voila.
    line_array = pdb_line.split()
    line_type = line_array[0]
    num = line_array[1]
    chain_mol = line_array[2]
    if (chain_mol == 'MOLECULE:'):
      (pre, mol_name) = pdb_line.split(": ")
    if (chain_mol == 'CHAIN:'):
      (pre, chain_name) = pdb_line.split(": ")
      mol_name = str(mol_name)
      chain_name = str(chain_name)
      mol_name = mol_name.rstrip()
      mol_name = mol_name.rstrip(';')
      mol_name = mol_name.replace(')','')
      mol_name = mol_name.replace('(','')
      mol_name = mol_name.replace('*','')
      mol_name = mol_name.replace('\'','')
      mol_name = mol_name.replace('\`','')
      mol_name = mol_name.replace('\\','')
      mol_name = mol_name.replace(' ', '_')
      color = default_colors['default']
      if (mol_name.find('PROTEIN') > -1):
        if (mol_name.find('RACK') > -1):
          color = default_colors['RACK']
        elif (mol_name.find('30S') > -1):
          color = default_colors['SSU_protein']
        elif (mol_name.find('50S') > -1):
          color = default_colors['LSU_protein']
        elif (mol_name.find('40S') > -1):
          color = default_colors['SSU_protein']
        elif (mol_name.find('60S') > -1):
          color = default_colors['LSU_protein']
        else:
          color = default_colors['unknown']
      elif (mol_name.find('RNA') > -1):
        if (mol_name.find('TRNA') > -1):
          color = default_colors['tRNA']
        elif (mol_name.find('MRNA') > -1 or mol_name.find('MESSENGER') > -1):
          color = default_colors['mRNA']
        elif (mol_name.find('RRNA') > -1 or mol_name.find('S_RNA') > -1 or mol_name.find('RIBOSOMAL') > -1):
          if (mol_name.find('16S') > -1):
            color = default_colors['SSU_RNA']
          elif (mol_name.find('23S') > -1):
            color = default_colors['LSU_RNA']
          elif (mol_name.find('5S') > -1):
            color = default_colors['5S_RNA']
          elif (mol_name.find('5.8S') > -1):
            color = default_colors['5.8S_RNA']
          else:
            color = default_colors['RNA']
      else:
        color = default_colors['other']

      ## Now that the various colors have been chosen, find the chains and make them
      mol_name = mol_name.replace('RIBOSOMAL_','')
      mol_name = mol_name.replace('S_RN', 'S_RRN')
      mol_name = mol_name.replace('PROTEIN_','')
      chain_name = chain_name.rstrip()
      chain_name = chain_name.rstrip(';')
      if (chain_name.find(',') > -1):
        chain_array = chain_name.split(',')
        num = 0
        for single_chain_name in chain_array:
          single_chain_name = single_chain_name.replace(' ','')
          selection_string = '/' + pdb_shortname + '//' + single_chain_name
          mol_name = mol_name + str(num)
          define_chain(selection_string, color, mol_name)
          num = num + 1
      else:
        selection_string = '/' + pdb_shortname + '//' + chain_name
        define_chain(selection_string, color, mol_name)
  cmd.disable(pdb_shortname)

## check_names is intended to avoid having duplicate names
## for those cases when there are crystals containing multiple
## ribosomes and/or multiple tRNAs in the same ribosome
## If that happens, an '_' is just appended to the molecule name.
## If more are found, more '_'s are added
def check_names(current, used_names=[]):
  current = str(current)
  found = 0
  for used_mol_name in used_names:
    if current == used_mol_name:
      current = current + '_'
      found = found + 1
  if (found > 0):
    return check_names(current)
  else:
    used_names.append(current)
    return current

def define_chain(selection_string, color, mol_name):
  mol_name = check_names(mol_name)
  cmd.create(mol_name, selection_string)
  cmd.show("cartoon", mol_name)
  cmd.color(color, mol_name)
  cmd.disable(mol_name)
  cmd.zoom("all")
  molecule_list.append(mol_name)


cmd.extend("chain_color", chain_color)
cmd.extend("make_pretty", make_pretty)
cmd.extend("make_chains", make_chains)
cmd.extend("load_session", load_session)
cmd.extend("split_pdb", random_chains("", ""))
cmd.extend("helices", helices)
cmd.extend("delete_original",delete_original)
cmd.extend("delete_mrna", delete_mrna)
cmd.extend("delete_trna", delete_trna)
cmd.extend("delete_lsu_rna", delete_lsu_rna)
cmd.extend("delete_ssu_rna", delete_ssu_rna)
cmd.extend("delete_lsu_protein", delete_lsu_protein)
cmd.extend("delete_ssu_protein", delete_ssu_protein)
cmd.extend("delete_lsu_helices", delete_lsu_helices)
cmd.extend("delete_ssu_helices", delete_ssu_helices)
cmd.extend("delete_all_helices", delete_all_helices)

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2004 by Charles Moad <cmoad@indiana.edu>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
