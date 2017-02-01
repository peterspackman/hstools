import tkinter as tk
import logging
from tkinter.filedialog import askopenfile
from .decompose import describe_surface
from .search import make_invariants

LOG = logging.getLogger('hs-csd-search')

class SearchApplication(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()
        self.cif = None

    def create_widgets(self):
        self.cif = tk.Button(self, text='Open CIF', command=self.open_cif)
        self.cif.pack(side="left")
        self.search_button = tk.Button(self, text='Search', command=self.search)
        self.search_button.pack(side="top")

        self.quit = tk.Button(self, text="Exit",
                              command=self.master.destroy)
        self.quit.pack(side="bottom")

    def search(self):
        LOG.debug("Searching")

    def open_cif(self):
        filename = askopenfile()
        if filename:
            self.cif = filename
        LOG.debug("cif: %s", self.cif)


def main():
    logging.basicConfig(level='DEBUG')
    root = tk.Tk()
    app = SearchApplication(master=root)
    app.mainloop()
