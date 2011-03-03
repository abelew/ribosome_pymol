from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkColorChooser
import tkFileDialog
import sys, string, re, os
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
    self.menuBar.addcascademenu('Ribosome','Specific Ribosomes')
    self.menuBar.addcascademenu('Specific Ribosomes','Saccharomyces cerevisiae')
    self.menuBar.addmenuitem('Saccharomyces cerevisiae', 'command', 'Beckmann 2010', label='Beckmann 2010', command = lambda: fetch("3IZS",""))
    self.menuBar.addmenuitem('Saccharomyces cerevisiae', 'command', 'Beckmann 2009', label='Beckmann 2009', command = lambda: fetch("2WW9",""))
    self.menuBar.addmenuitem('Saccharomyces cerevisiae', 'command', 'Yusupov 2010', label='Yusupov 2010', command = lambda: fetch("3O2Z",""))
    self.menuBar.addmenuitem('Saccharomyces cerevisiae', 'command', 'Doudna 2009', label='Doudna 2009', command = lambda: fetch("3FRX",""))
    self.menuBar.addmenuitem('Saccharomyces cerevisiae', 'command', 'Frank 2008_eEF2', label='Frank 2008_eEF2', command = lambda: fetch("3DNY",""))
    self.menuBar.addmenuitem('Saccharomyces cerevisiae', 'command', 'Spahn 2007_IRES', label='Spahn 2007_IRES', command = lambda: fetch("2NOQ",""))
    self.menuBar.addcascademenu('Specific Ribosomes','Thermus thermophilus')
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2010-1', label='Ramakrishnan 2010-1', command = lambda: fetch("2XQD",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2010-2', label='Ramakrishnan 2010-2', command = lambda: fetch("2XFZ",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2010-3', label='Ramakrishnan 2010-3', command = lambda: fetch("2X9R",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2009-1', label='Ramakrishnan 2009-1', command = lambda: fetch("3KIQ",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2009-2', label='Ramakrishnan 2009-2', command = lambda: fetch("2WRI",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2009-3', label='Ramakrishnan 2009-3', command = lambda: fetch("2WRN",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2009-4', label='Ramakrishnan 2009-4', command = lambda: fetch("2WH1",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2009-5', label='Ramakrishnan 2009-5', command = lambda: fetch("2WDG",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Spahn 2009', label='Spahn 2009', command = lambda: fetch("3FIC",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2007-1', label='Ramakrishnan 2007-1', command = lambda: fetch("2UXD",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2007-2', label='Ramakrishnan 2007-2', command = lambda: fetch("2V46",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2007-3', label='Ramakrishnan 2007-3', command = lambda: fetch("2UUC",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2006-1', label='Ramakrishnan 2006-1', command = lambda: fetch("2J00",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2006-2', label='Ramakrishnan 2006-2', command = lambda: fetch("2B64",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2006-3', label='Ramakrishnan 2006-3', command = lambda: fetch("2B9M",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2006-4', label='Ramakrishnan 2006-4', command = lambda: fetch("2B9P",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2004-1', label='Ramakrishnan 2004-1', command = lambda: fetch("1XM0",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Ramakrishnan 2004-2', label='Ramakrishnan 2004-2', command = lambda: fetch("1XNR",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Jogl 2010', label='Jogl 2010', command = lambda: fetch("3OTO",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Noller 2008', label='Noller 2008', command = lambda: fetch("3F1E",""))
    self.menuBar.addmenuitem('Thermus thermophilus', 'command', 'Noller 2007', label='Noller 2007', command = lambda: fetch("2OW8",""))
    self.menuBar.addcascademenu('Specific Ribosomes','Escherichia coli')
    self.menuBar.addmenuitem('Escherichia coli', 'command', 'Frank 2009', label='Frank 2009', command = lambda: fetch("3F1H",""))
    self.menuBar.addcascademenu('Specific Ribosomes','Thermomyces languinosus')
    self.menuBar.addmenuitem('Thermomyces languinosus', 'command', 'Frank 2009', label='Frank 2009', command = lambda: fetch("3JYV",""))

    self.menuBar.addmenuitem('Ribosome','command','2dHelices',label='2dHelices',command = lambda: twod_helices())

