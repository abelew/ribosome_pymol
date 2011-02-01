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
    self.menuBar.addmenuitem('Ribosome','command','2dHelices',label='2dHelices',command = lambda: twod_helices())

