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
    self.menuBar.addcascademenu('Ribosome','Delete chains')


