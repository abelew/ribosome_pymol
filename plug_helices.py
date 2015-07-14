from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkColorChooser
import tkFileDialog
import sys
import string
import re
import os
import csv
import urllib
import gzip
import platform
from pymol import stored, cmd, selector, movie

stored.organism = "saccharomyces_cerevisiae"

helices_path = None
try:
    helices_path = os.environ["HELICES_HOME"]
except:
    helices_path = os.environ["PYMOL_DATA"]
sys.path.append(helices_path)
datadir = None
if platform.system() == "Windows":
    datadir = "C:\\Program Files\\PyMOL\\ribosome_pymol\\helices_data\\"
else:
    datadir = helices_path + "/helices_data/"

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
large_subunit_rnas = ['26S', '25S', '23S', '5.8S', '5S']
small_subunit_prot = ['40S', '30S']
large_subunit_prot = ['60S', '50S']

## The function 'fetch_then_chains' was taken with very minor changes
## from remote_pdb_load.py.  The copyright notice is at the bottom of
## this file as well as the README

## Set up the directory where the data files live and some global variables
## which will be used to store the names of the molecules and helices from the
## PDB files.

## This function loads into the pymol menu system and creates 'Ribosome' menu.
def __init__(self):
    cmd.unset("ignore_case")
    cmd.set("orthoscopic", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("depth_cue", 1)
    cmd.set("ray_trace_fog", 1)
    cmd.set("antialias", 1.0)
    cmd.set("cartoon_ring_mode", 3)
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
    ## Added for Suna to make some bases opaque, but the rest transparent
    self.menuBar.addmenuitem('Colors', 'command', 'Custom file trans.',
                             label = 'Custom file trans.',
                             command = lambda: chain_color("trans"))
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
    elif my_type == "Darwin":
        open_command = "open"
    elif my_type == "Windows":
        open_command = "explorer"
    os.system(open_command + " " + ribosomes_path + " &")
## End of edit_ribosomes

def del_enabled():
    mols = cmd.get_names(enabled_only = 1)
    for mol in mols:
        cmd.delete(mol)

def thick_lines_enabled(width):
    """
    thick_lines_enabled
    Attempts to set the width of the enabled molecules.
    """
    mols = cmd.get_names(enabled_only = 1)
    for mol in mols:
        cmd.set("line_width", width, mol)

def delete_all_helices():
    """
    delete_all_helices
    Attempts to delete all helices to save memory
    """
    delete_ssu_helices()
    delete_lsu_helices()
## End of delete_all_helices


def delete_all_rna():
    """
    delete_all_rna
    Attempts to delete all the RNA molecules to save memory
    """
    delete_ssu_rna()
    delete_lsu_rna()
    delete_mrna()
    delete_trna()
## End of delete_all_rna


def delete_all_protein():
    """
    delete_all_protein
    Attempts to delete all the proteins to save memory
    """
    delete_lsu_protein()
    delete_ssu_protein()
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
    """
    delete_original
    Delete the original pdb entries to save memory
    """
    for original_molecule in original_list:
        cmd.delete(original_molecule)
## End of delete_original


def delete_lsuh():
    """
    delete_lsuh
    Delete the helices of the large subunit
    """
    counter = 0
    while (counter <= 104):
        counter = counter + 1
        string = "LSU_H", counter
        cmd.delete(string)
## End of delete_lsuh


def delete_ssuh():
    """
    delete_ssuh
    Delete the helices of the small subunit
    """
    counter = 0
    while (counter <= 45):
        counter = counter + 1
        string = "SSU_h", counter
        cmd.delete(string)
## End of delete_ssuh

## I changed like 2 lines from remote_load_pdb.py
## The main change is at the end of fetch()
class fetch_then_chains:
    """
    fetch_then_chains
    Take in a PDB accession, download the file, parse its header
    and display the pieces.
    """
    def __init__(self, app):
        pdbCode = tkSimpleDialog.askstring('PDB Loader Service',
                                           'Please enter a 4-digit pdb code:',
                                           parent=app.root)
        if pdbCode: # None is returned for user cancel
            pdbCode = string.upper(pdbCode)
            fetch(pdbCode,"")
## End of fetch_then_chains


def check_fetch(information):
    """
    check_fetch
    Print some information about a ribosomal pdb before fetching it.
    This information lies in helices_data/structures.csv
    """
    mymessage = "The pdb " + information[3] + ", species: " + str(information[0]) + " came from the " + str(information[1]) + " lab in " + str(information[2]) + " described by:\n" + str(information[4]) + "\nClick 'yes' if you wish to view this pdb file."
    response = tkMessageBox.askyesno(title=information[3],message=mymessage)
    if response:
        fetch(information[3],"")
## End of check_fetch

def color_saccharomyces():
    """
    Color the yusupov 2011 ribosome according to Dr. Dinman's preferences.
    """
    protein_colors = {
        'P1_ALPHA' : 'blue',
        'P2_BETA' : 'hotpink',
        '60S_L2' : 'blue',
        '60S_L3' : 'green',
        '60S_L4' : 'hotpink',
        '60S_L5' : 'limon',
        '60S_L6' : 'forest',
        '60S_L7' : 'palegreen',
        '60S_L8' : 'tv_green',
        '60S_L9' : 'limegreen',
        '60S_L10' : 'red',
        '60S_L11' : 'cyan',
        '60S_L13' : 'orange',
        '60S_L14' : 'tv_blue',
        '60S_L15' : 'lightmagenta',
        '60S_L16' : 'magenta',
        '60S_L17' : 'smudge',
        '60S_L18' : 'slate',
        '60S_L19' : 'marine',
        '60S_L20' : 'brightorange',
        '60S_L21' : 'purpleblue',
        '60S_L22' : 'purple',
        '60S_L23' : 'slate',
        '60S_L24' : 'yellow',
        '60S_L25' : 'violet',
        '60S_L26' : 'teal',
        '60S_L27' : 'sand',
        '60S_L28' : 'chocolate',
        '60S_L29' : 'blue',
        '60S_L30' : 'red',
        '60S_L31' : 'warmpink',
        '60S_L32' : 'marine',
        '60S_L33' : 'tv_green',
        '60S_L34' : 'orange',
        '60S_L35' : 'limon',
        '60S_L36' : 'magenta',
        '60S_L37' : 'chocolate',
        '60S_L38' : 'limegreen',
        '60S_L39' : 'chocolate',
        'UBIQUITIN-60S_L40' : 'chocolate',
        '60S_ACIDIC_P0' : 'chocolate',
        '60S_L41' : 'red',
        '60S_L42' : 'wheat',
        '60S_L43' : 'smudge',
        '40S_S0' : 'cyan',
        '40S_S1' : 'purpleblue',
        '40S_S2' : 'teal',
        '40S_S3' : 'yellow',
        '40S_S4' : 'forest',
        '40S_S5' : 'chartreuse',
        '40S_S6' : 'orange',
        '40S_S7' : 'chartreuse',
        '40S_S8' : 'red',
        '40S_S9' : 'yellow',
        '40S_S10' : 'red',
        '40S_S11' : 'marine',
        '40S_S12' : 'tv_green',
        '40S_S13' : 'wheat',
        '40S_S14' : 'tv_red',
        '40S_S15' : 'orange',
        '40S_S16' : 'slate',
        '40S_S17' : 'red',
        '40S_S18' : 'blue',
        '40S_S19' : 'violetpurple',
        '40S_S20' : 'deepolive',
        '40S_S21' : 'red',
        '40S_S22' : 'orange',
        '40S_S23' : 'blue',
        '40S_S24' : 'density',
        '40S_S25' : 'red',
        '40S_S26' : 'yellow',
        '40S_S27' : 'deepblue',
        '40S_S28' : 'orange',
        '40S_S29' : 'marine',
        '40S_S30' : 'raspberry',
        'UBIQUITIN-40S_S31' : 'magenta',
        'GUANINE_NUCLEOTIDE-BINDING_SUBUNIT_BETA-LIKE' : 'chocolate',
        'SUPPRESSOR_STM1' : 'olive',    
        }
    for item in (protein_colors.keys()):
        string = '/' + item
        try:
            cmd.color(protein_colors[item], string)
        except:
            print "Failed " + string

def fetch(pdb, splitp):
    """
    fetch
    As per load_pdb, the only difference is that it calls to the following
    function, random_chains
    """
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
    """
    chain_color
    read over helices_data/color_definitions.txt to get an idea
    of likely colors for chains and residues.
    """
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

    if bases == "trans":
        objects = cmd.get_names()
        for o in objects:
            for v in ["cartoon_ring_transparency" , "cartoon_transparency" , "stick_transparency"]:
                cmd.set(v, 1.0, o)

  ## From here until 'if input_file:' the color definitions
  ## are specified
    colors_file = datadir + 'color_definitions.txt'
    colors = dict({None : 'gray',})
    transp = dict({None : '0.0',})
    if colors_file:
        color_lines = file(colors_file).readlines()
    for color_line in color_lines:
        color_datum = color_line.split()
        colors[color_datum[0]] = color_datum[1]
        try:
            transp[color_datum[0]] = color_datum[2]
        except:
            transp[color_datum[0]] = "0.0"

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
                color_name = colors[color_choice]
                print "Setting color to: " + color_name
                try:
                    print "In the try."
                    cmd.color(color_name, selection_string)
                except:
                    print "In the except."
                    cmd.color(colors[None], selection_string)

                print "TESMTE: " + transp[color_choice]
                if str(transp[color_choice]) != "0.0":
                    print "Setting trans to: " + my_trans
                    cmd.set("cartoon_ring_transparency", my_trans, selection_string)
                    cmd.set("cartoon_transparency", my_trans, selection_string)
                    cmd.set("stick_transparency", my_trans, selection_string)
            except:
                print "Cannot find your selection, perhaps you must split the chains first"
## End of chain_color


def make_pretty():
  ## These are some settings our professor prefers.
  cmd.bg_color("white")
  cmd.show("cartoon")
#  cmd.set("cartoon_ring_mode", 3)
## End of make_pretty


## This function is the toplevel function to make pretty helices
## Change the default_colors['helix'] to whatever color you prefer.
def helices(new_organism = stored.organism):
    """
    helices
    read a file in helices_data/ which corresponds to the species
    of this ribosome.  These contain definitions for every ribosomal
    helix.  Create individual pymol objects for every helix.
    """
    make_chains(stored.organism, 'sticks', default_colors['helix'])
## End of helices


## This should ask for the relevant data file and call the
## cheater perl scripts I wrote
def twod_helices():
    """
    twod_helices
    Ask the user for an input file named either:
    18S_rRNA.txt, 25S_rRNA_3p.txt, or 25S_rRNA_5p.txt
    These files should contain base numbers followed by
    an integer 'color.'  This will then color a 2d representation
    of the Saccharomyces cerevisiae ribosome with these colors.
    """
    print "Provide the text file you wish to use to color, the filename"
    print "Should be one of: 18S_rRNA.txt, 25S_rRNA_3p.txt, 25S_rRNA_5p.txt"
    print "Otherwise this is not smart enough to understand what you want."
    twod_textfile = tkFileDialog.askopenfile(parent=app.root,mode='rb',title='Choose a file')
    new_twod = twod_filename
    new_twod = re.sub('txt$', 'ps', str(twod_textfile))
    old_twod = os.path.basename(new_twod)
    old_twod = datadir + str(old_twod)
    file_new_twod = open(new_twod, 'w')
    file_old_twod = open(old_twod, 'r')
    file_twod_tex = open(twod_textfile, 'r')
    
    ## First get the numbers from the text file.
    if file(file_twod_text) is not None:
        twod_text_lines = file(file_twod_text).readlines()
        color_list = []
        for li in twod_text_lines:
            (num, col) = li.split()
            color_list.append(col)
        file.close(file_twod_text)

    if file(file_old_twod) is not None:
        test_string = ""
        if str(file_old_twod) == '18S_rRNA.ps':
            test_string = "290.00 -105.33 290.00 -98.67 lwline"
        elif str(fold_old_twod) == '25S_rRNA_3p.ps':
            test_string = "-148.33 -1010.00 -141.67 -1010.00 lwline"
        elif str(fold_old_twod) == '25S_rRNA_5p.ps':
            test_string = "360.00 0.00 1.00 1.00 1.00 431.01 154.00 lwfarc"
        count = None
        list_count = 0
        for li in file_old_twod:
            file_new_twod.write(li)
            if li == test_string:
                count = 0
            elif count == 0:
                count = count + 1
            elif count == 1:
                count = count - 1
                chosen_color = color_list[list_count]  ## A number from the input file
                if (chosen_color == 0):  ## black
                    file_new_twod.write("0 0 0 setrgbcolor\n")
                elif (chosen_color == 10):  ## gray
                    file_new_twod.write("0.6 0.6 0.6 setrgbcolor\n")
                elif (chosen_color == 11):   ## neon pink
                    file_new_twod.write("0.85 0.30 0.64 setrgbcolor\n")
                elif (chosen_color == -4):  ## purple
                    file_new_twod.write("0.36 0.18 0.64 setrgbcolor\n")
                elif (chosen_color == -3):  ## blue
                    file_new_twod.write("0.08 0.25 1.0 setrgbcolor\n")
                elif (chosen_color == -2):  ## greenblue
                    file_new_twod.write("0.25 0.90 0.92 setrgbcolor\n")
                elif (chosen_color == -1):  ## green
                    file_new_twod.write("0.1 0.90 0.1 setrgbcolor\n")
                elif (chosen_color == 1):  ## yellow
                    file_new_twod.write("0.9 0.9 0.15 setrgbcolor\n")
                elif (chosen_color == 2):  ## yelloworange
                    file_new_twod.write("0.90 0.60 0.10 setrgbcolor\n")
                elif (chosen_color == 3):  ## orangered
                    file_new_twod.write("0.92 0.34 0.08 setrgbcolor\n")
                elif (chosen_color == 4):  ## red
                    file_new_twod.write("0.92 0.10 0.10 setrgbcolor\n")
                else:
                    file_new_twod.write("0.57 0.08 0.32 setrgbcolor\n")
    file.close(file_old_twod)
    file.close(file_new_twod)
## End of twod_helices
  

def make_chains(chains, showastype, showascolor):
    """
    make_chains
    Running make_chains should have pymol read a file in helices_data
    which contains specifications of every ribosomal helix and some
    special features (the PTC for instance)
    Pymol will use this information to create objects corresponding to
    every object.
    """
  ## Start out figuring out the data file to specify the helices
  ## Currently I just have a stupid if/elif chain for the few species
  ## I have annotated.
    cmd.set("auto_zoom", "off")
    cmd.set("auto_show_selections", "off")
    cmd.set("cartoon_fancy_helices", 1)

    test_chains = datadir + '/' + str(chains) + '/helices.txt'
    if file(test_chains) is not None:
        chains_filenames = [ test_chains , ]
    else:
        chains_filenames = [ datadir + 'wtf.txt', ]
#        if chains == 'wtf':
#            print "WTF"
#            chains_filenames = [datadir + 'wtf.txt',]
#        elif chains.find('escherichia_coli') > -1:
#            chains_filenames = [datadir + 'escherichia_coli.txt',]
#        elif chains.find('thermomyces_lanuginosus') > -1:
#            chains_filenames = [datadir + 'thermomyces_languinosus.txt',]
#        elif chains.find('thermus_thermophilus') > -1:
#            chains_filenames =  [datadir + 'thermus_thermophilus.txt',]
#        elif chains.find('haloarcula_marismortui') > -1:
#            chains_filenames = [datadir + 'haloarcula_marismortui.txt',]
#        elif chains.find('saccharomyces') > -1:
#            chains_filenames = [datadir + 'saccharomyces_helices.txt',]
#        elif chains.find('saccharomyces_cerevisiae') > -1:
#            chains_filenames = [datadir + 'saccharomyces_helices.txt',]
#        elif chains.find('saccharomyces_cerevisiae_s288c') > -1:
#            chains_filenames = [datadir + 'saccharomyces_helices.txt',]
#        else:
#            print "Could not understand the argument:" + chains +", using the wtf file"
#            chains_filenames = [datadir + 'wtf.txt',]

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
                if re.compile('-[A-Z]$').search(name) is not None:
                    name = re.sub('-[A-Z]$', '', name)
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
    """
    random_chains
    Using random_chains will parse the pdb file's
    header in order to discover a few things:
    a)  Is there more than 1 pdb file in this complete image?
    b)  What species is this from?
    c)  What are the proteins, RNAs, and ligands in this image?
    With this information it will attempt to create individual
    objects for every chain in the pdb file which have helpful
    colors and names.
    """
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
                    (pre, org) = pdb_line.split(": ")
                    org = str(org)
                    org = org.rstrip()
                    org = org.rstrip(';')
                    org = org.lower()
                    org = org.replace(' ', '_')
                    stored.organism = "%s" %(org)
                    
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
        #      mol_name = re.sub(r')|(|*|\'|\`|\\)' , '' , mol_name)
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
                matched = 0
                if (mol_name.find('RACK') > -1):
                    color = default_colors['RACK']
                    matched = 1
                for ssu in small_subunit_prot:
                    if (mol_name.find(ssu) > -1):
                        color = default_colors['SSU_protein']
                        matched = 1
                for lsu in large_subunit_prot:
                    if (mol_name.find(lsu) > -1):
                        color = default_colors['LSU_protein']
                        matched = 1
                if matched == 0:
                    color = default_colors['unknown']

            elif (mol_name.find('RNA') > -1):
                if (mol_name.find('TRNA') > -1):
                    color = default_colors['tRNA']
                elif (mol_name.find('MRNA') > -1 or mol_name.find('MESSENGER') > -1):
                    color = default_colors['mRNA']
                elif (mol_name.find('RRNA') > -1 or mol_name.find('S_RNA') > -1 or mol_name.find('RIBOSOMAL') > -1):
                    matched = 0
                    for lsr in large_subunit_rnas:
                        if (mol_name.find(lsr) > -1):
                            color = default_colors['LSU_RNA']
                            matched = 1
                    for ssr in small_subunit_rnas:
                        if (mol_name.find(ssr) > -1):
                            color = default_colors['SSU_RNA']
                            matched = 1
                    if matched == 0:
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
    make_pretty()
## End of random_chains


## check_names is intended to avoid having duplicate names
## for those cases when there are crystals containing multiple
## ribosomes and/or multiple tRNAs in the same ribosome
## If that happens, an '_' is just appended to the molecule name.
## If more are found, more '_'s are added
## Changing this now to a number
def check_names(current, increment = 1):
    current = str(current)
    found = 0
    for used_mol_name in molecule_list:
        if current == used_mol_name:
            found = 1
            if increment == 1:
                current = current +  '_' + str(increment)
            else:
                current = current.rstrip(str(increment - 1))
                current = current + str(increment)
            increment = increment + 1

    if (found > 0):
        return check_names(current, increment)
    else:
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
        print "There was an error with:" + str(mol_name)
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
        sel = selecton + " and resn %s" % amino_acid
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
    

def switch_chains(old, delete_second = 0):
    new_name = old + "_1"
    tmp = new_name + "tmp"
    rename_chain(new_name, tmp)
    if str(delete_second) == "1":
        cmd.delete(old)
    else:
        rename_chain(old, new_name)
    rename_chain(tmp, old)

def check_ratcheted():
    ## Use B7a to figure out if any given ribosomal crystal structure is in 
## the ratcheted state or not.  
    print "Not yet implemented."

def rename_on_bridge_dist(delete_second = 0, lsu_nucleotide = "/25S_RRNA///1024" , ssu_nucleotide = "/18S_RRNA///1240" , switch_prot = "40S_S" , switch_rna = "18S_RRNA"):
    b1a_distance = 0
    try:
        b1a_distance = cmd.distance("B1A_test", lsu_nucleotide, ssu_nucleotide)
    except:
        lsu_nucleotide = "/23S_RRNA///886"
        ssu_nucleotide = "/30S_S13///93"
        switch_prot = "30S_S"
        switch_rna = "16S_RRNA"
        b1a_distance = cmd.distance("B1A_test", lsu_nucleotide, ssu_nucleotide)

    print "The distance across the current B1A is: " + str(b1a_distance)
    if (b1a_distance <= 0):
        lsu_nucleotide = "/23S_RRNA///886"
        ssu_nucleotide = "/30S_S13///93"
        switch_prot = "30S_S"
        switch_rna = "16S_RRNA"
        b1a_distance = cmd.distance("B1A_test", lsu_nucleotide, ssu_nucleotide)

    cmd.disable("B1A_test")
    print "The distance across the current B1A is: " + str(b1a_distance)
    ## The default values for the B1A bridge are from Yusupov's 2010 paper about the yeast ribosome
        
    if (b1a_distance > 100):
        mol_list = cmd.get_names("all")
        for mol in mol_list:
            if re.compile('_1$').search(mol) is None:
                if mol.find(switch_prot) > -1:
                    switch_chains(mol, delete_second)
                elif mol.find(switch_rna) > -1:
                    switch_chains(mol, delete_second)

        for mol in mol_list:
            if re.compile('UNASSIGNED').search(mol):
                print "Skipping unassigned"
            elif re.compile('_$').search(mol):
                name_two = mol + "2"
                rename_chain(mol, name_two)

def find_all_neighbors(distance = 3, species = "saccharomyces_cerevisiae"):
    all_mol = cmd.get_names("all")
    for mol in all_mol:
        for second_mol in all_mol:
            find_neighbors(mol, second_mol, distance)

def search_interactions_helices(distance = 3 , species = "saccharomyces_cerevisiae"):
    lsu_protein_list = []
    ssu_protein_list = []
    rna_list = []
    mol_list = cmd.get_names("all")
    for mol in mol_list:
        for lsu_prot in large_subunit_prot:
            if mol.find(lsu_prot) > -1:
                print "Found " + mol
                lsu_protein_list.append(mol)
    for ssu_prot in small_subunit_prot:
        if mol.find(ssu_prot) > -1:
            print "Found " + mol
            ssu_protein_list.append(mol)

    helices_filename = datadir + species + "/helices.txt"
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
    mol_list = cmd.get_names("all")
    for mol in mol_list:
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

def movie_stitch(png_dir=None,bitrate=2400000,opts="msmpeg4v2:dia=2:predia=2:qns=3"):
    """
    movie_stitch
    Choose a directory containing your pymol movie png files
    and mencoder will stitch them into a movie in the same
    directory named 'output.wmv'
    Alternate mencoder commands:
    mencoder mf://*.png -o mjpeg.avi -ovc lavc -lavcopts vcodec=mjpeg:vhq:psnr -noskip
    """
    if png_dir == None:
        png_dir = tkFileDialog.askdirectory(title="Pick the directory with your movie files.")
    if png_dir != None:
        mencoder_command = "cd " + png_dir + " && /usr/bin/mencoder mf://*.png -o output.wmv -of lavf -ovc lavc -lavcopts vcodec=" + opts + ":vbitrate=" + str(bitrate)
        print "Running: " + mencoder_command
        s = os.system(mencoder_command)
        print "Your movie lives in: " + png_dir + "/" + "output.wmv"
    else:
        print "Cannot run without a directory."


def delete_enabled():
    """
    delete_enabled
    Delete anything which is currently enabled.
    """
    en = cmd.get_names(enabled_only = 1)
    for e in en:
        cmd.delete(e)

def transparent_enabled(tr = 0.7):
    """
    transparent_enabled
    Make the currently enabled items in pymol transparent
    If you pass this a variable, it will change the transparency
    0 is completely opaque, 1 is completely transparent.
    The default is 0.7
    """
    en = cmd.get_names(enabled_only = 1)
    for e in en:
        for v in ["cartoon_ring_transparency" , "cartoon_transparency" , "stick_transparency"]:
            print "Setting " + v + " to 0.7 for " + e
            cmd.set(v, tr, e)
#            cmd.set(v, tr, e)
#            cmd.set(v, tr, e)
    print "Consider also setting the following variables:"
    print "cartoon_oval_width , cartoon_tube_radius , line_width"
    print "cartoon_loop_radius , cartoon_rect_width"

def crown_view(species = "saccharomyces_cerevisiae"):
    """
    Return the (currently yeast) ribosome to the crown view.
    """
    if (species == "saccharomyces_cerevisiae"):
        cmd.set_view("0.003838378, -0.220511124, -0.975378811, 0.732650876, -0.663220644, 0.152820781, -0.680590451, -0.715204060, 0.159013703, -0.000904173, -0.000284255, -898.858032227, 10.964979172, 21.375329971, 75.672447205, -387776.812500000, 389574.437500000, -20.000000000")
    else:
        print "I haven't yet set the coordinates for " + species
        print "Try this: "
        cmd.set_view("0.003838378, -0.220511124, -0.975378811, 0.732650876, -0.663220644, 0.152820781, -0.680590451, -0.715204060, 0.159013703, -0.000904173, -0.000284255, -898.858032227, 10.964979172, 21.375329971, 75.672447205, -387776.812500000, 389574.437500000, -20.000000000")
#
#"0.030628815 , 0.004798603 , -0.999527216 , 0.482913435 , -0.875605762 , 0.010592541 , -0.875139892 , -0.483004808 , -0.029141756 , 0.000000000 , 0.000000000 , -259.756744385 , 12.284235001 , -28.817157745 , 10.823875427 , 204.794189453 , 314.719299316 , -20.000000000")


## End of search_interactions
cmd.extend("color_saccharomyces",color_saccharomyces)
cmd.extend("movie_stitch", movie_stitch)
cmd.extend("transparent_enabled", transparent_enabled)
cmd.extend("crown_view", crown_view)
cmd.extend("del_enabled", del_enabled)
cmd.extend("switch_chains", switch_chains)
cmd.extend("rename_chain", rename_chain)
cmd.extend("rename_on_bridge_dist", rename_on_bridge_dist)
cmd.extend("search_interactions_helices", search_interactions_helices)
cmd.extend("search_interactions", search_interactions)
cmd.extend("find_neighbors", find_neighbors)
cmd.extend("find_all_neighbors", find_all_neighbors)
cmd.extend("color_by_amino_acid", color_by_amino_acid)
cmd.extend("color_by_aa_residue_type", color_by_aa_residue_type)
cmd.extend("chain_color", chain_color)
cmd.extend("make_pretty", make_pretty)
cmd.extend("make_chains", make_chains)
cmd.extend("load_session", load_session)
cmd.extend("split_pdb", random_chains)
cmd.extend("helices", helices)
cmd.extend("delete_enabled", delete_enabled)
cmd.extend("delete_original", delete_original)
cmd.extend("delete_mrna", delete_mrna)
cmd.extend("delete_trna", delete_trna)
cmd.extend("delete_lsu_rna", delete_lsu_rna)
cmd.extend("delete_ssu_rna", delete_ssu_rna)
cmd.extend("delete_lsu_protein", delete_lsu_protein)
cmd.extend("delete_ssu_protein", delete_ssu_protein)
cmd.extend("delete_lsu_helices", delete_lsu_helices)
cmd.extend("delete_ssu_helices", delete_ssu_helices)
cmd.extend("delete_all_helices", delete_all_helices)
cmd.extend("thick_lines_enabled", thick_lines_enabled)
cmd.extend("get_seq", get_seq)
cmd.extend("random_chains", random_chains)

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
