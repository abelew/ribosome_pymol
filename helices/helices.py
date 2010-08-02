### Stealing horribly from remote_pdb_load.py
from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkColorChooser
import tkFileDialog
import os, sys, string, re
import urllib
import gzip
from pymol import stored, cmd, selector

helices_dir = os.path.dirname(__file__)
datadir=helices_dir + "/data/"
global molecule_list
molecule_list = []
global original_list
original_list = []
global helices_list
helices_list = []

#global organism

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
    print "TESTME: " + helix
    if helix.find('LSU_H') > -1:
      cmd.delete(helix)
    else:
      print "DID NOT FIND"

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
    # filename = urllib.urlretrieve('http://www.rcsb.org/pdb/cgi/export.cgi/' +
    # pdbCode + '.pdb.gz?format=PDB&pdbId=' +
    # pdbCode + '&compression=gz')[0]
    filename = urllib.urlretrieve('http://www.rcsb.org/pdb/files/' + pdb + '.pdb.gz')[0]
  except:
    tkMessageBox.showerror('Connection Error',
                           'Can not access to the PDB database.\n'+
                           'Please check your Internet access.',)
#                           parent=app.root)
  else:
    if (os.path.getsize(filename) > 0): # If 0, then pdb code was invalid
      # Uncompress the file while reading
      fpin = gzip.open(filename)
        # Form the pdb output name
      outputname = os.path.dirname(filename) + os.sep + pdb + '.pdb'
      fpout = open(outputname, 'w')
      fpout.write(fpin.read()) # Write pdb file
      fpin.close()
      fpout.close()
      cmd.load(outputname,quiet=0) # Load the fresh pdb
#      print "TESTME about to run: " + filename
      random_chains(outputname,splitp)
    else:
      tkMessageBox.showerror('Invalid Code', 'You entered an invalid pdb code:' + pdb,)
#                             parent=app.root)
          
      os.remove(filename) # Remove tmp file (leave the pdb)


def chain_color(bases):
    input_file = ""
    if bases == "modified":
        input_file = datadir + 'modifications.txt'
    else:
        tmp_filename = tkFileDialog.askopenfile(title="Open a session")
        if not tmp_filename: return
        input_file = tmp_filename.name
    comment = ''

    colors_file = datadir + 'color_definitions.txt'
    colors = dict({None : 'gray',})
    if colors_file:
        color_lines = file(colors_file).readlines()
    for color_line in color_lines:
        color_datum = color_line.split()
        colors[color_datum[0]] = color_datum[1]
#        print "Setting ", color_datum[0], " to ", color_datum[1]
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
#                print "Working with " + subunit + chain
            except:
                subunit = tmpre.strip()
#                print "Working with " + subunit
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
#                print 'Coloring: ', selection_string, ' ', color_name
                try:
                    cmd.color(color_name, selection_string)
#                    cmd.show("sticks", selection_string)
                except:
                    cmd.color(colors[None], selection_string)
            except:
                print "Cannot find your selection, perhaps you must split the chains first"


def make_pretty():
#    print "Setting the background to black, change this with 'bg_color <color>'"
    cmd.bg_color("white")
    print "Setting the protein of the large subunit to yelloworange"
    print "Change this with 'color yelloworange, 3JYW'"
    cmd.color("yelloworange", "/3JYW")
    print "Setting the large subunit RNA to gray60"
    print "Change this with 'color gray60, /3JYX'"
    cmd.color("gray60", "/3JYX")
    print "Setting the entire small subunit to palecyan"
    print "Including the proteins, tRNA, and RNA"
    cmd.color("palecyan", "/3JYV")
    print "Setting the tRNA to slate"
    cmd.color("slate", "/3JYV//7")
    print "Setting the small subunit RNA to gray70"
    cmd.color("gray70", "/3JYV//A")      ## //A is the small subunit RNA
    print "Setting the cartoon settings"
    cmd.show("cartoon")
    cmd.set("cartoon_ring_mode", 3)

def helices():
#  print "TESTME organism " + organism
  make_chains(organism, 'sticks', 'red')

def make_chains(chains, showastype, showascolor):
#  print "TESTME chain, type, color " + chains + showastype + showascolor
  if chains == 'wtf':
    print "WTF"
    chains_filenames = [datadir + 'wtf.txt',]
  elif chains == 'ssh':
    chains_filenames = [datadir + 'helices_ssu.txt',]
  elif chains == 'frank':
    load_frank_ribosome()
    chains_filenames = [datadir + 'frank_chains.txt',]
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
  elif chains == 'lsh':
    chains_filenames = [datadir + 'helices_lsu.txt',]
  elif chains == 'rna':
    chains_filenames = [datadir + 'chains_rna.txt',]
  elif chains == 'ssu':
    chains_filenames = [datadir + 'chains_ssu.txt',]
  elif chains == 'lsu':
    chains_filenames = [datadir + 'chains_lsu.txt',]
  elif chains == 'proteins':
    chains_filenames = [datadir + 'chains_lsu.txt', datadir + 'chains_ssu.txt']
  elif chains == 'all':
    chains_filenames = [datadir + 'chains_rna.txt', datadir + 'chains_ssu.txt', datadir + 'chains_lsu.txt',]
  else:
    print "Could not understand the argument:" + chains +", using the wtf file"
    chains_filenames = [datadir + 'wtf.txt',]
  for chains_filename in chains_filenames:
#    print "Working with " + chains_filename
    if chains_filename:
      chains_lines = file(chains_filename).readlines()
      for ch in chains_lines:
        if re.compile('^#').search(ch) is not None:
          continue
        name = ''
        location = ''
        (name, location) = ch.split(',')
#        print "Creating " + name
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
          if chains == 'frank':
            load_frank_ribosome()
          elif chains == 'saccharomyces':
            load_frank_ribosome()
            make_chains("frank", "", "")
    ## Zoom to something sane
        if chains == 'frank':
          delete_frank()
          cmd.zoom("all")


def load_session(filename):
    filename = tkFileDialog.askopenfile(title="Open a session")
    if not filename: return
    file_path = filename.name
    cmd.load(file_path)

def random_chains(pdb_file, splitp):
  cmd.bg_color("white")
  import re
  if pdb_file is None:
    pdb_file = tkFileDialog.askopenfile(title="Open a session")
  if not pdb_file: return
  pdb_filename = str(pdb_file)
#  print "TESTME " + pdb_filename
  pdb_basename = os.path.basename(pdb_filename)
  pdb_shortname = os.path.splitext(pdb_basename)
  pdb_shortname = pdb_shortname[0] 
  original_list.append(pdb_shortname)

#  print "TESTME " + pdb_filename
    ## pdb_filename is the full filename
    ## pdb_shortname is the 2WGD or whathaveyou
    ## pdb_basename is the path it lives in
    ## pdb_file is the file object which has all the attributes etc
  cmd.load(pdb_filename)
  pdb_lines = file(pdb_filename).readlines()
  chain = ''
  source_count = 0
  global organism
  organism = str('saccharomyces_cerevisiae')
  for pdb_line in pdb_lines:
    if re.compile("^HEADER").search(pdb_line) is not None:
      continue
    if re.compile("^TITLE").search(pdb_line) is not None:
      continue
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
      color = 'gray10'
      if (mol_name.find('PROTEIN') > -1):
        if (mol_name.find('RACK') > -1):
          color = 'cyan'
        elif (mol_name.find('30S') > -1):
          color = 'cyan'
        elif (mol_name.find('50S') > -1):
          color = 'skyblue'
        elif (mol_name.find('40S') > -1):
          color = 'cyan'
        elif (mol_name.find('60S') > -1):
          color = 'skyblue'
        else:
          color = 'green'
      elif (mol_name.find('RNA') > -1):
        if (mol_name.find('TRNA') > -1):
          color = 'slate'
        elif (mol_name.find('MRNA') > -1 or mol_name.find('MESSENGER') > -1):
          color = 'forest'
        elif (mol_name.find('RRNA') > -1 or mol_name.find('S_RNA') > -1 or mol_name.find('RIBOSOMAL') > -1):
          if (mol_name.find('16S') > -1):
            color = 'gray20'
          elif (mol_name.find('23S') > -1):
            color = 'gray30'
          elif (mol_name.find('5S') > -1):
            color = 'gray40'
          elif (mol_name.find('5.8S') > -1):
            color = 'gray50'
          else:
            color = 'gray60'
      else:
        color = 'red'

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
