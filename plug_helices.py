from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkColorChooser
import tkFileDialog
import sys, string, re, os
from pymol import cmd
from helices.helices import *

def __init__(self):
    self.menuBar.addcascademenu('Plugin','Ribosome')

    def doitall():
        first = load_frank_ribosome()
        make_pretty("")
        make_chains("all","","")
        make_chains("lsh","sticks","red")
        make_chains("ssh","sticks","blue")
        chain_color("modified")
        cmd.bg_color("black")
    self.menuBar.addmenuitem('Ribosome','command','Load from PDB',
                             label = 'Load from PDB',
                             command = lambda s=self : fetch_then_chains(s))
    self.menuBar.addmenuitem('Ribosome','command','Helices',label='Helices',command = lambda: helices())
    self.menuBar.addmenuitem('Ribosome', 'command', 'Load another session',
                             label = "Load another session",
                             command = lambda:  load_session(""))
    self.menuBar.addcascademenu('Ribosome','Split chains')
    self.menuBar.addmenuitem('Split chains', 'command', 'Frank chains',
                             label = 'Frank chains',
                             command = lambda: make_chains("frank","",""))
    self.menuBar.addmenuitem('Split chains', 'command', 'Random chains',
                             label = 'Random chains',
                             command = lambda: random_chains(""))
    self.menuBar.addmenuitem('Split chains', 'command', 'Saccharomyces helices',
                             label = 'Saccharomyces helices',
                             command = lambda: make_chains("saccharomyces","","red"))
    self.menuBar.addcascademenu('Ribosome','Color positions')
    self.menuBar.addmenuitem('Color positions', 'command', 'Modified Bases',
                             label = 'Modified Bases',
                             command = lambda: chain_color("modified"))
    self.menuBar.addmenuitem('Color positions', 'command', 'Custom file',
                             label = 'Custom file',
                             command = lambda: chain_color("custom"))
    self.menuBar.addcascademenu('Ribosome','Delete chains')
    self.menuBar.addmenuitem('Delete chains', 'command', 'Remove frank',
                             label = 'Remove frank',
                             command = lambda: delete_frank())
    self.menuBar.addmenuitem('Delete chains', 'command', 'Remove LSU',
                             label = 'Remove LSU',
                             command = lambda: delete_lsu())
    self.menuBar.addmenuitem('Delete chains', 'command', 'Remove SSU',
                             label = 'Remove SSU',
                             command = lambda: delete_ssu())
    self.menuBar.addmenuitem('Delete chains', 'command', 'Remove LSUp',
                             label = 'Remove LSUp',
                             command = lambda: delete_lsup())
    self.menuBar.addmenuitem('Delete chains', 'command', 'Remove SSUp',
                             label = 'Remove SSUp',
                             command = lambda: delete_ssup())
    self.menuBar.addmenuitem('Delete chains', 'command', 'Remove LSU helices',
                             label = 'Remove LSU helices',
                             command = lambda: delete_lsuh())
    self.menuBar.addmenuitem('Delete chains', 'command', 'Remove SSU helices',
                             label = 'Remove SSU helices',
                             command = lambda: delete_ssuh())
