from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkColorChooser
import tkFileDialog
import sys, string, re, os, csv
import urllib
import gzip
import subprocess
import platform
from pymol import stored, cmd, selector
pymol_data_path = os.getenv("PYMOL_DATA")
helices_path = pymol_data_path + "/pmg_tk"
sys.path.append(helices_path)
#helices_dir = os.path.dirname(__file__)
datadir=helices_path + "/helices_data/"
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
## Small subunit should have full subunit, bacterial then eukaryotic
## Then the RNA bacterial/euk
## Large subunit should be full bac/euk, then the big RNA, then little
global small_subunit_rnas
global large_subunit_rnas
global small_subunit_prot
global large_subunit_prot
small_subunit_rnas = ['18S', '16S']
large_subunit_rnas = ['25S', '23S', '5.8S', '5S']
small_subunit_prot = ['40S', '30S']
large_subunit_prot = ['60S', '50S', 'RACK']

## The function 'fetch_then_chains' was taken with very minor changes
## from remote_pdb_load.py.  The copyright notice is at the bottom of
## this file as well as the README

## Set up the directory where the data files live and some global variables
## which will be used to store the names of the molecules and helices from the
## PDB files.

## This function loads into the pymol menu system and creates 'Ribosome' menu.
def __init__(self):
    cmd.unset("ignore_case")
    self.menuBar.addcascademenu('Plugin','Ribosome')

    self.menuBar.addmenuitem('Ribosome','command','Load from PDB',
                             label = 'Load from PDB',
                             command = lambda s=self : fetch_then_chains(s))
    self.menuBar.addmenuitem('Ribosome','command','Helices',label='Helices',command = lambda: helices())
    self.menuBar.addmenuitem('Ribosome', 'command', 'Load another session',
                             label = "Load another session",
                             command = lambda:  load_session(""))
    self.menuBar.addcascademenu('Ribosome','Colors')
    self.menuBar.addmenuitem('Colors', 'command', 'Color by attribute', label = 'Color by attribute',
                             command = lambda: color_by_aa_residue_type())
    self.menuBar.addmenuitem('Colors','command','Color by amino acid', label = 'Color by amino acid',
                             command = lambda: color_by_amino_acid())
    self.menuBar.addmenuitem('Colors', 'command', 'Modified Bases',
                             label = 'Modified Bases',
                             command = lambda: chain_color("modified"))
    self.menuBar.addmenuitem('Colors', 'command', 'Custom file',
                             label = 'Custom file',
                             command = lambda: chain_color("custom"))
    self.menuBar.addcascademenu('Ribosome','Delete objects')
    self.menuBar.addmenuitem('Delete objects', 'command', 'Original', label='Original',
                             command = lambda: delete_original())
    self.menuBar.addmenuitem('Delete objects', 'command', 'ALL_RNA', label='ALL_RNA',
                             command = lambda: delete_all_rna())
    self.menuBar.addmenuitem('Delete objects', 'command', 'LSU_RNA', label='LSU_RNA',
                             command = lambda: delete_lsu_rna())
    self.menuBar.addmenuitem('Delete objects', 'command', 'SSU_RNA', label='SSU_RNA',
                             command = lambda: delete_ssu_rna())
    self.menuBar.addmenuitem('Delete objects', 'command', 'All_protein', label='All_protein',
                             command = lambda: delete_all_protein())
    self.menuBar.addmenuitem('Delete objects', 'command', 'LSU_protein', label='LSU_protein',
                             command = lambda: delete_lsu_protein())
    self.menuBar.addmenuitem('Delete objects', 'command', 'SSU_protein', label='SSU_protein',
                             command = lambda: delete_ssu_protein())
    self.menuBar.addmenuitem('Delete objects', 'command', 'All_helices', label='All_helices',
                             command = lambda: delete_all_helices())
    self.menuBar.addmenuitem('Delete objects', 'command', 'LSU_helices', label='LSU_helices',
                             command = lambda: delete_lsu_helices())
    self.menuBar.addmenuitem('Delete objects', 'command', 'SSU_helices', label='SSU_helices',
                             command = lambda: delete_ssu_helices())
    self.menuBar.addmenuitem('Ribosome','command','2dHelices',label='2dHelices',command = lambda: twod_helices())
    self.menuBar.addmenuitem('Ribosome','command','Get_Sequence',label='Get_Sequence',command = lambda: get_seq())
    self.menuBar.addmenuitem('Ribosome','command','Edit_Ribosomes',label='Edit_Ribosomes',command = lambda: edit_ribosomes())
    specific_ribosome_menu(self)

def specific_ribosome_menu(self):
    """
    specific_ribosome_menu:
    Called by init to load ribosomes by year/species/author using a CSV spreadsheet
    in the 'data' directory.
    """
    infile = datadir + "structures.csv"
    ribosome_species = dict({})
    ribosome_years = dict({})
    ribosome_authors = dict({})
    csvfile = open(infile)
    dialect = csv.Sniffer().sniff(csvfile.read())
    csvfile.seek(0)
    reader = csv.reader(csvfile, dialect)

    for datum in reader:
        species = datum[0]
        author = datum[1]
        year = datum[2]
        accession = datum[3]
        title = datum[4]
        if species not in ribosome_species:
            ribosome_species[species] = [(species,author,year,accession,title)]
        else:
            ribosome_species[species].append((species,author,year,accession,title))

        if year not in ribosome_years:
            ribosome_years[year] = [(species,author,year,accession,title)]
        else:
            ribosome_years[year].append((species,author,year,accession,title))

        if author not in ribosome_authors:
            ribosome_authors[author] = [(species,author,year,accession,title)]
        else:
            ribosome_authors[author].append((species,author,year,accession,title))

    self.menuBar.addcascademenu('Ribosome','Ribosomes by Species')
    for spec in sorted(ribosome_species.keys()):
        self.menuBar.addcascademenu('Ribosomes by Species', spec)
        entry_list = ribosome_species[spec]
        for entry in entry_list:
            entry_name = entry[1] + "-" + entry[2] + "-" + entry[3]
            self.menuBar.addmenuitem(spec, 'command', entry_name, label=entry_name, command = lambda s=entry : check_fetch(s))

    self.menuBar.addcascademenu('Ribosome','Ribosomes by Year')
    for year in sorted(ribosome_years.keys()):
        self.menuBar.addcascademenu('Ribosomes by Year', year)
        year_list = ribosome_years[year]
        for entry in year_list:
            entry_name = entry[1] + "-" + entry[0] + "-" + entry[3]
            self.menuBar.addmenuitem(year, 'command', entry_name, label=entry_name, command = lambda s=entry: check_fetch(s))

    self.menuBar.addcascademenu('Ribosome','Ribosomes by Author')
    for author in sorted(ribosome_authors.keys()):
        self.menuBar.addcascademenu('Ribosomes by Author', author)
        author_list = ribosome_authors[author]
        for entry in author_list:
            entry_name = entry[2] +  "-" + entry[0] + "-" + entry[3]
            self.menuBar.addmenuitem(author, 'command', entry_name, label=entry_name, command = lambda s=entry: check_fetch(s))


def edit_ribosomes():
    """
    edit_ribosomes: Opens the ribosome database in openoffice/excel
    """
    my_type = platform.system()
    ribosomes_path = datadir + "structures.csv"
    open_command = ""
    if my_type == "Linux":
        open_command = "openoffice.org"
    elif my_type == "MacOS":
        open_command = "open"
    elif my_type == "Windows":
        open_command = "explorer"
    os.system(open_command + " " + ribosomes_path + " &")
## End of edit_ribosomes


def delete_all_helices():
    """
    delete_all_helices
    Attempts to delete all helices to save memory
    """
    delete_ssu_helices
    delete_lsu_helices
## End of delete_all_helices


def delete_all_rna():
    """
    delete_all_rna
    Attempts to delete all the RNA molecules to save memory
    """
    delete_ssu_rna
    delete_lsu_rna
    delete_mrna
    delete_trna
## End of delete_all_rna


def delete_all_protein():
    """
    delete_all_protein
    Attempts to delete all the proteins to save memory
    """
    delete_lsu_protein
    delete_ssu_protein
## End of delete_all_protein


def delete_ssu_helices():
    """
    delete_ssu_helices
    Attempts to delete the small subunit helices to save memory
    """
    for helix in helices_list:
        helix = helix.lstrip('/')
        if helix.find('SSU_h') > -1:
            cmd.delete(helix)
## End of delete_ssu_helices


def delete_lsu_helices():
    """
    delete_lsu_helices
    Attempts to delete the large subunit helices to save memory
    """
    for helix in helices_list:
        helix = helix.lstrip('/')
        helix = str(helix)
        if helix.find('LSU_H') > -1:
            cmd.delete(helix)
## End of delete_lsu_helices


def delete_ssu_protein():
    """
    delete_ssu_protein
    Attempts to delete the small subunit proteins to save memory
    """
    for mol in molecule_list:
        mol = mol.lstrip('/')
        for ssu in small_subunit_prot:
            ssu_name = ssu + '_S'
            if mol.find(ssu_name) > -1:
                cmd.delete(mol)
## End of delete_ssu_protein


def delete_lsu_protein():
    """
    delete_lsu_protein
    Attempts to delete large subunit proteins to save memory
    """
  for mol in molecule_list:
    mol = mol.lstrip('/')
    for lsu in large_subunit_prot:
        lsu_mol = lsu + '_L'
        if mol.find(lsu_mol) > -1:
            cmd.delete(mol)
## End of delete_lsu_protein


def delete_ssu_rna():
    """
    delete_ssu_rna
    Attempts to delete small subunit rna to save memory
    """
  for mol in molecule_list:
    mol = mol.lstrip('/')
    for ssu in small_subunit_rnas:
        ssu_mol = ssu + '_RRNA'
        if mol.find(ssu_mol) > -1:
            cmd.delete(mol)
## End of delete_ssu_rna


def delete_lsu_rna():
    """
    delete_lsu_rna
    Attempts to delete the large subunit rna to save memory
    """
  for mol in molecule_list:
    mol = mol.lstrip('/')
    for lsu in large_subunit_rnas:
        lsu_mol = lsu + '_RRNA'
        if mol.find(lsu_mol) > -1:
            cmd.delete(mol)
## End of delete_lsu_rna


def delete_trna():
    """
    delete_trna
    Attempts to delete any tRNAs to save memory
    """
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('TRNA') > -1:
      cmd.delete(mol)
## End of delete_trna


def delete_mrna():
    """
    delete_mrna
    Attempts to delete any mRNA molecules to save memory
    """
  for mol in molecule_list:
    mol = mol.lstrip('/')
    if mol.find('MRNA') > -1:
      cmd.delete(mol)
## End of delete_mrna


def delete_original():
  for original_molecule in original_list:
    cmd.delete(original_molecule)
## End of delete_original


def delete_lsuh():
  counter = 0
  while (counter <= 104):
      counter = counter + 1
      string = "LSU_H", counter
      cmd.delete(string)
## End of delete_lsuh


def delete_ssuh():
    counter = 0
    while (counter <= 45):
        counter = counter + 1
        string = "SSU_h", counter
        cmd.delete(string)
## End of delete_ssuh

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
## End of fetch_then_chains


def check_fetch(information):
  mymessage = "The pdb " + information[3] + ", species: " + str(information[0]) + " came from the " + str(information[1]) + " lab in " + str(information[2]) + " described by:\n" + str(information[4]) + "\nClick 'yes' if you wish to view this pdb file."
  response = tkMessageBox.askyesno(title=information[3],message=mymessage)
  if response:
    fetch(information[3],"")
## End of check_fetch


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
## End of fetch



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
## End of chain_color


def make_pretty():
  ## These are some settings our professor prefers.
  cmd.bg_color("white")
  cmd.show("cartoon")
  cmd.set("cartoon_ring_mode", 3)
## End of make_pretty


## This function is the toplevel function to make pretty helices
## Change the default_colors['helix'] to whatever color you prefer.
def helices(new_organism = "saccharomyces_cerevisiae"):
    if organism is None:
        make_chains(new_organism, 'sticks', default_colors['helix'])
    else:
        make_chains(organism, 'sticks', default_colors['helix'])
## End of helices


## This should ask for the relevant data file and call the
## cheater perl scripts I wrote
def twod_helices():
  subprocess.Popen([r"2d/make_color_ps.pl"]).wait()
## End of twod_helices
  

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
  elif chains.find('thermomyces_lanuginosus') > -1:
    chains_filenames = [datadir + 'thermomyces_languinosus.txt',]
  elif chains.find('thermus_thermophilus') > -1:
    chains_filenames =  [datadir + 'thermus_thermophilus.txt',]
  elif chains.find('haloarcula_marismortui') > -1:
    chains_filenames = [datadir + 'haloarcula_marismortui.txt',]
  elif chains.find('saccharomyces') > -1:
    chains_filenames = [datadir + 'saccharomyces_helices.txt',]
  elif chains.find('saccharomyces_cerevisiae') > -1:
    chains_filenames = [datadir + 'saccharomyces_helices.txt',]
  elif chains.find('saccharomyces_cerevisiae_s288c') > -1:
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
## End of make_chains

def load_session(filename):
    filename = tkFileDialog.askopenfile(title="Open a session")
    if not filename: return
    file_path = filename.name
    cmd.load(file_path)
## End of load_session


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
    ## Go recursion!
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
    try:
      line_type = line_array[0]
    except:
      line_type = "UNKNOWN"
      print "problem with line_type on:" + pdb_line
    try:
      num = line_array[1]
    except:
      num = 0
      print "problem with num on:" + pdb_line
    try:
      chain_mol = line_array[2]
    except:
      chain_mol = "UNKNOWN"
      print "problem with chain_mol on:" + pdb_line
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
      ## This section needs to be generalized using
      ## large_subunit_rnas and the similar globals
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
## End of random_chains


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
## End of check_names


def get_seq(selection_string):
  from pymol import stored
  sequence=[]
  residues={}
  one_letter={
    "ALA":"A",
    "ARG":"R",
    "ASN":"N",
    "ASP":"D",
    "ASX":"B",
    "CYS":"C",
    "CYH":"C",#protonated cys
    "CYX":"C",#cystein
    "GLN":"Q",
    "GLU":"E",
    "GLY":"G",
    "GLX":"Z",
    "HIS":"H",
    "HIP":"H",#protonated
    "HID":"H",#h on delta
    "HIE":"H",#h on epsilon
    "ILE":"I",
    "LEU":"L",
    "LYS":"K",
    "MET":"M",
    "MSE":"M",#Seleno- Methionin
    "PHE":"F",
    "PRO":"P",
    "SER":"S",
    "THR":"T",
    "TRP":"W",
    "TYR":"Y",
    "VAL":"V",
    "UNK":"X",
    "A":"A",
    "C":"C",
    "G":"G",
    "T":"T",
    "U":"U",
    ## Below are residues which can't be added
    "MG":"",
    "HOH":"",
    "OHX":"",
    }
  stored.residue_names=[]
  cmd.iterate(selection_string , "stored.residue_names.append(resn)")
  stored.residue_numbers=[]
  cmd.iterate(selection_string , "stored.residue_numbers.append(resi)")

  for index in range(len(stored.residue_names)):
    residues[stored.residue_numbers[index]] = stored.residue_names[index]

#  for key in sorted(residues.keys(),lambda x,y: 1 if x[1] > y[1] else -1):
#    print "key: %s value: %s" % (key, residues[key])

  result = ""
  for k in sorted(residues.keys(), cmp=_compare_keys):
    result += one_letter[residues[k]]

  print result
## End of get_seq



def _compare_keys(x, y):
  try:
    x = int(x)
  except ValueError:
    xint = False
  else:
    xint = True
  try:
    y = int(y)
  except ValueError:
    if xint:
      return -1
    return cmp(x.lower(), y.lower())
    # or cmp(x, y) if you want case sensitivity.
  else:
    if xint:
      return cmp(x, y)
    return 1
## End of _compare_keys


def define_chain(selection_string, color, mol_name):
  if color is None:
    color = "black"
  mol_name = check_names(mol_name)
  try:
    cmd.create(mol_name, selection_string)
    cmd.show("cartoon", mol_name)
    cmd.color(color, mol_name)
    cmd.disable(mol_name)
    cmd.zoom("all")
  except:
    print "There was an error with:" + mol_name
  molecule_list.append(mol_name)
## End of define_chain


## The following functions were mostly taken from 
## Robert L. Campbell's 2004 color_by_restype.py
## which appears to no longer work
## I am changing it to (hopefully) work with pymol 1.3
def color_by_aa_residue_type(selection):
    if selection is None:
        selections = cmd.get_names(enabled_only = 1)
        for sel in selections:
            color_by_aa_residue_type(sel)
    residue_abbrev = {
        'A': 'ALA', 
        'C': 'CYS', 
        'D': 'ASP', 
        'E': 'GLU', 
        'F': 'PHE', 
        'G': 'GLY', 
        'H': 'HIS', 
        'I': 'ILE', 
        'K': 'LYS', 
        'L': 'LEU', 
        'M': 'MET', 
        'N': 'ASN', 
        'P': 'PRO', 
        'Q': 'GLN', 
        'R': 'ARG', 
        'S': 'SER', 
        'T': 'THR', 
        'V': 'VAL', 
        'W': 'TRP', 
        'Y': 'TYR',
        }
    abbrev_to_type = {
        'A': 'hydrophobic',
        'ALA' : 'hydrophobic',
        'C': 'cysteine',
        'CYS' : 'cysteine',
        'D': 'negative',
        'ASP' : 'negative',
        'E': 'negative',
        'GLU' : 'negative',
        'F': 'aromatic',
        'PHE' : 'aromatic',
        'G': 'hydrophobic',
        'GLY' : 'hydrophobic',
        'H': 'polar',
        'HIS' : 'polar',
        'I': 'hydrophobic',
        'ILE' : 'hydrophobic',
        'K': 'positive',
        'LYS' : 'positive',
        'L': 'hydrophobic',
        'LEU' : 'hydrophobic',
        'M': 'hydrophobic',
        'MET' : 'hydrophobic',
        'N': 'polar',
        'ASN' : 'polar',
        'P': 'proline',
        'PRO' : 'proline',
        'Q': 'polar',
        'GLN' : 'polar',
        'R': 'positive',
        'ARG' : 'positive',
        'S': 'polar',
        'SER' : 'polar',
        'T': 'polar',
        'THR' : 'polar',
        'V': 'hydrophobic',
        'VAL' : 'hydrophobic',
        'W': 'aromatic',
        'TRP' : 'aromatic',
        'Y': 'aromatic',
        'TYR' : 'aromatic',
        }
    type_colors = {
        'hydrophobic' : 'grey70',
        'aromatic' : 'lightmagenta',
        'polar' : 'teal',
        'positive' : 'blue',
        'negative' : 'red',
        'cysteine' : 'yellow',
        'proline' : 'green',
        }
    for aa in residue_abbrev:
        amino_acid = residue_abbrev[aa]
        sel = selection + " and resn %s" % amino_acid
        chosen_color = type_colors[abbrev_to_type[aa]]
        cmd.select("temp", sel)
        cmd.color(chosen_color, "temp")
    for type in type_colors:
        print "Coloring " + type + " residues " + type_colors[type]
## End of color_by_aa_residue_type

def color_by_amino_acid(selection):
    if selection is None:
        selections = cmd.get_names(enabled_only = 1)
        for sel in selections:
            color_by_amino_acid(sel)        
    residue_colors = {
        'ALA':'gray40',
        'CYS':'yellow',
        'ASP':'red', 
        'GLU':'tv_red',
        'PHE':'deepblue', 
        'GLY':'gray80',
        'HIS': 'lightblue',
        'ILE':'forest',
        'LYS':'blue',
        'LEU':'green',
        'MET':'tv_yellow',
        'ASN':'cyan',
        'PRO':'yelloworange',
        'GLN':'palecyan',
        'ARG':'tv_blue',
        'SER':'orange',
        'THR':'tv_orange',
        'VAL':'splitpea',
        'TRP':'purple',
        'TYR':'density',
        'A':'blue',
        'G':'gray50',
        'C':'red',
        'U':'green',
        'T':'green',
        }
    for aa in residue_abbrev:
        amino_acid = residue_abbrev[aa]
        sel = selection + " and resn %s" % amino_acid
        chosen_color = type_colors[abbrev_to_type[aa]]
        cmd.select("temp", sel)
        cmd.color(chosen_color, "temp")
    for type in type_colors:
        print "Coloring " + type + " residues " + type_colors[type]
## End of color_by_amino_acid


def find_neighbors(sel_one, sel_two, distance = 3, selection_name = "neighbors", pairs_mode = 0):
    cmd.select(selection_name, "name c")
    print "Looking for pairs using a distance of " + str(distance) + ".  " + sel_one + " " + sel_two
    pairs_list = cmd.find_pairs("(" + sel_one + ")", "(" + sel_two + ")", mode = pairs_mode, cutoff = float(distance))
    for pairs in pairs_list:
        print "\"\",\"" + sel_one + "\",\"",
        cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]), 'print "%s%s %s" % (resi, resn, name),')
        print "\",\"" + sel_two + "\",\"",
        cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]), 'print "%s%s %s" % (resi, resn, name),')
        print "\",\"",
        print "%.2f" % cmd.distance(selection_name, "%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1])),
        print "\""
        cmd.show("dashes", selection_name)
#  cmd.show("labels", selection_name)
        cmd.color("magenta", selection_name)
        cmd.zoom(selection_name)
## End of find_neighbors

def rename_chain(old_name, new_name, tmp_name = 'tmp_'):
    old_obj = "(" + old_name + ")"
    cmd.select(tmp_name, old_obj)
    cmd.create(new_name, tmp_name)
    cmd.delete(old_name)
    cmd.disable(new_name)

def switch_chains(old):
    new_name = old + "_"
    tmp = new_name + "tmp"
    rename_chain(new_name, tmp)
    rename_chain(old, new_name)
    rename_chain(tmp, old)

def rename_on_bridge_dist(lsu_nucleotide = "/25S_RRNA///1024" , ssu_nucleotide = "/18S_RRNA///1240" , switch_prot = "40S_S" , switch_rna = "18S_RRNA"):
    b1a_distance = cmd.distance("B1A_test", lsu_nucleotide, ssu_nucleotide)
    ## The default values for the B1A bridge are from Yusupov's 2010 paper about the yeast ribosome
    if (b1a_distance < 0):
        ## These values for B1A come from Noller 2001 and Yusupov 2010
        rename_on_bridge_dist("/23S_RRNA///886", "/30S_S13///93", "30S_S", "16S_RRNA")
    else:
        
        if (b1a_distance > 100):
            mol_list = cmd.get_names("all")
            for mol in mol_list:
                if re.compile('_$').search(mol) is None:
                    if mol.find(switch_prot) > -1:
                        switch_chains(mol)
                    elif mol.find(switch_rna) > -1:
                        switch_chains(mol)

            for mol in mol_list:
                if re.compile('UNASSIGNED').search(mol):
                    print "Skipping unassigned"
                elif re.compile('_$').search(mol):
                    name_two = mol + "2"
                    rename_chain(mol, name_two)

def search_interactions_helices(distance = 3 , species = "saccharomyces_cerevisiae"):
    lsu_protein_list = []
    ssu_protein_list = []
    rna_list = []
    molecule_list = cmd.get_names("all")
    for mol in molecule_list:
        for lsu_prot in large_subunit_prot:
            if mol.find(lsu_prot) > -1:
                print "Found " + mol
                lsu_protein_list.append(mol)
        for ssu_prot in small_subunit_prot:
            if mol.find(ssu_prot) > -1:
                print "Found " + mol
                ssu_protein_list.append(mol)

    helices_filename = datadir + species + ".txt"
    chains_lines = file(helices_filename).readlines()
    for ch in chains_lines:
        if re.compile('^#').search(ch) is not None:
            continue
        name = ''
        location = ''
        (name, location) = ch.split(',')
        rna_list.append(name)

    for prot in lsu_protein_list:
        for rna in rna_list:
            print "Searching " + prot + " against " + rna
            find_neighbors(prot, rna, distance)

    for prot in ssu_protein_list:
        for rna in rna_list:
            print "Searching " + prot + " against " + rna
            find_neighbors(prot, rna, distance)


# enabled_objs = cmd.get_names("all",enabled_only=1)
def search_interactions(distance = 3):
  rna_list = []
  lsu_protein_list = []
  ssu_protein_list = []
  molecule_list = cmd.get_names("all")
  for mol in molecule_list:
    if mol.find('RRNA') > -1:
      rna_list.append(mol)
  for lsu_mol in large_subunit_prot:
      if mol.find(lsu_mol) > -1:
          lsu_protein_list.append(mol)
  for ssu_mol in small_subunit_prot:
      if mol.find(ssu_mol) > -1:
          ssu_protein_list.append(mol)

  for prot in lsu_protein_list:
      for rna in rna_list:
          print "Searching " + prot + " against " + rna
          find_neighbors(prot, rna, distance)
  
  for prot in ssu_protein_list:
      for rna in rna_list:
          print "Searching " + prot + " against " + rna
          find_neighbors(prot, rna, distance)

  for rna_one in rna_list:
      for rna_two in rna_list:
          if rna_one != rna_two:
              find_neighbors(rna_one , rna_two)

  for lsu_prot in lsu_protein_list:
      for ssu_prot in ssu_protein_list:
          find_neighbors(lsu_prot , ssu_prot)

  for lsu_one in lsu_protein_list:
      for lsu_two in lsu_protein_list:
          if lsu_one != lsu_two:
              find_neighbors(lsu_one , lsu_two)

  for ssu_one in ssu_protein_list:
      for ssu_two in ssu_protein_list:
          if ssu_one != ssu_two:
              find_neighbors(ssu_one , ssu_two)
## End of search_interactions
cmd.extend("switch_chains", switch_chains)
cmd.extend("rename_chain", rename_chain)
cmd.extend("rename_on_bridge_dist", rename_on_bridge_dist)
cmd.extend("search_interactions_helices", search_interactions_helices)
cmd.extend("search_interactions", search_interactions)
cmd.extend("find_neighbors", find_neighbors)
cmd.extend("color_by_amino_acid", color_by_amino_acid)
cmd.extend("color_by_aa_residue_type", color_by_aa_residue_type)
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
cmd.extend("get_seq", get_seq)

## Below is Charles Moad's copyright notice for fetch.py
## Below that, I added the GPLv2, which if it does not conflict
## I would like to use for the code I wrote.


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






# (gpl-text
#         "  This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
 
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# ")
