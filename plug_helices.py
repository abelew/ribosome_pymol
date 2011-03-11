from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkColorChooser
import tkFileDialog
import sys, string, re, os, csv
pymol_data_path = os.getenv("PYMOL_DATA")
helices_path = pymol_data_path + "/pmg_tk"
sys.path.append(helices_path)

from pymol import cmd
from helices.helices import *
## This file loads into the pymol menu system and creates 'Ribosome' menu.

def __init__(self):
    self.menuBar.addcascademenu('Plugin','Ribosome')

    self.menuBar.addmenuitem('Ribosome','command','Load from PDB',
                             label = 'Load from PDB',
                             command = lambda s=self : fetch_then_chains(s))
    self.menuBar.addmenuitem('Ribosome','command','Helices',label='Helices',command = lambda: helices())
    self.menuBar.addmenuitem('Ribosome', 'command', 'Load another session',
                             label = "Load another session",
                             command = lambda:  load_session(""))
    self.menuBar.addcascademenu('Ribosome','Color positions')
    self.menuBar.addmenuitem('Color positions', 'command', 'Modified Bases',
                             label = 'Modified Bases',
                             command = lambda: chain_color("modified"))
    self.menuBar.addmenuitem('Color positions', 'command', 'Custom file',
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
    specific_ribosome_menu(self)

def specific_ribosome_menu(self):
    infile = datadir + "/structures.csv"
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
